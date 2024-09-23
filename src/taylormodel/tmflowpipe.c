#include "tmflowpipe.h"

TaylorModel *computeTaylorPolynomial(ODEList *system, unsigned int order,
                                     unsigned int k) {
  assert(system != NULL);
  assert(order > 0);
  assert(k > 0);
  /* All the polynomial terms exceeding the truncation order
    would just get truncated anyways, so impose an explicit restriction. */
  assert(order <= k);

  /* The functions to seed each Lie derivation with. */
  TaylorModel *lieDerivativeSeed = initTaylorModel(system);
  /* The Taylor polynomials that will be built up in-place. */
  TaylorModel *polynomials = initTaylorModel(system);

  /* Start from i=1; case i=0 would be an order 0 Lie derivative. */
  for (unsigned int index = 1; index <= order; ++index) {
    TaylorModel *lieDeriv = lieDerivativeK(system, lieDerivativeSeed, index);

    TaylorModel *poly = polynomials;
    TaylorModel *deriv = lieDeriv;
    while (poly != NULL || deriv != NULL) {
      /* If one becomes NULL while the other does not,
        then there is a length mismatch. */
      assert((poly != NULL) == (deriv != NULL));

      /* Ensure each variable's derivative is added to that same variable's
        running Taylor polynomial. If this assertion fails, the variable
        ordering in the lists was messed up somehow. */
      assert(strcmp(poly->fun, deriv->fun) == 0);

      /* Make a copy to ensure the tree is decoupled from the derivative object.
       */
      ExpTree *derivExp = cpyExpTree(deriv->exp);

      /* fac(i) = gamma(i+1) */
      double factorial = tgamma(index + 1);
      char factorialStr[50];
      snprintf(factorialStr, sizeof(factorialStr), "%.15g", factorial);
      char indexStr[50];
      snprintf(indexStr, sizeof(indexStr), "%u", index);

      /* 1/i! */
      ExpTree *fac = newExpOp(EXP_DIV_OP, newExpLeaf(EXP_NUM, "1"),
                              newExpLeaf(EXP_NUM, factorialStr));
      /* t^i */
      ExpTree *tPow = newExpOp(EXP_EXP_OP, newExpLeaf(EXP_VAR, VAR_TIME),
                               newExpLeaf(EXP_NUM, indexStr));
      /* 1/i! * L^i(g) * t^i   where L^i(g) is the order i Lie derivative of
       * function g. */
      ExpTree *polyElement =
          newExpOp(EXP_MUL_OP, fac, newExpOp(EXP_MUL_OP, derivExp, tPow));

      /* Extend the polynomial in-place!. */
      poly->exp = newExpOp(EXP_ADD_OP, poly->exp, polyElement);

      poly = poly->next;
      deriv = deriv->next;
    }

    /* Cleanup */
    delTaylorModel(lieDeriv);
  }

  /* Cleanup */
  delTaylorModel(lieDerivativeSeed);

  return polynomials;
}

TaylorModel *lieDerivativeK(ODEList *system, TaylorModel *functions,
                            unsigned int order) {
  assert(system != NULL);
  assert(functions != NULL);

  /* Copy the input functions, so that the returned list is distinct from
    the input list. Then the lifetime of the returned instance TM is
    consistently the caller's responsibility. */
  if (order == 0)
    return cpyTaylorModel(functions);

  /* After pass i of the loop, the polynomial parts of the list of Taylor
    models represent a vector of order i Lie derivatives. So each pass
    raises the order of the Lie derivatives by one, until order k. */
  TaylorModel *ithLieDerivative = functions;
  for (unsigned int index = 0; index < order; ++index) {
    TaylorModel *old = ithLieDerivative;
    ithLieDerivative = lieDerivativeTaylorModel(system, ithLieDerivative);

    /* Cleanup: idx = 0 is the input, other intermediates may be deleted. */
    if (index > 0)
      delTaylorModel(old);
  }

  return ithLieDerivative;
}

TaylorModel *lieDerivativeTaylorModel(ODEList *system, TaylorModel *functions) {
  assert(system != NULL);
  assert(functions != NULL);

  ODEList *ode = system;
  TaylorModel *function = functions;
  TaylorModel *derived = NULL;

  /* Derive each of the functions individually w.r.t. the same ODe system. */
  while (ode != NULL || function != NULL) {
    /* If only one (XOR) is NULL but not the other, then there is a list length
      mismatch. Theoretically, they could have different lengths, but in our
      case their lengths are always the same. This assertion safeguards that
      loose assumption. */
    assert((ode != NULL) == (function != NULL));

    char *fun = strdup(function->fun);
    ExpTree *lieDeriv = lieDerivative(system, function->exp);
    Interval remainder = function->remainder;

    /* Post-process the result. */
    ExpTree *simplified = simplify(lieDeriv);
    delExpTree(lieDeriv);

    /* The initial tail should be NULL, since the list is extended head-first.
     */
    derived = newTMElem(derived, fun, simplified, remainder);

    ode = ode->next;
    function = function->next;
  }

  /* Reverse to ensure the output functions are
    ordered the same as the input functions. */
  return reverseTaylorModel(derived);
}

ExpTree *lieDerivative(ODEList *vectorField, ExpTree *function) {
  assert(vectorField != NULL);
  assert(function != NULL);

  /* Lf(g) = summ( d(g)/d(xi) * fi ) + d(g)/dt */
  ExpTree *lieDeriv = NULL;
  ODEList *ode = vectorField;
  while (ode != NULL) {
    /* d(g)/d(xi) * fi */
    ExpTree *dgdxi = derivative(function, ode->fun);
    ExpTree *fi = cpyExpTree(ode->exp);
    ExpTree *summElem = newExpOp(EXP_MUL_OP, dgdxi, fi);

    if (lieDeriv == NULL)
      lieDeriv = summElem;
    else
      lieDeriv = newExpOp(EXP_ADD_OP, lieDeriv, summElem);

    ode = ode->next;
  }
  /* Simply assume that a variable "t" exists. */
  /* d(g)/dt */
  lieDeriv = newExpOp(EXP_ADD_OP, lieDeriv, derivative(function, VAR_TIME));
  return lieDeriv;
}

TaylorModel *picardOperator(ODEList *vectorField, TaylorModel *functions) {
  assert(vectorField != NULL);
  assert(functions != NULL);

  /* Pf(g) = x0 + integral_0^t ( f(g(s), s) ds )
    Initialize the result to x0. */
  TaylorModel *picard = initTaylorModel(vectorField);
  TaylorModel *substitutedField = substituteTaylorModel(vectorField, functions);
  ExpTree *zero = newExpLeaf(EXP_NUM, "0");
  ExpTree *t = newExpLeaf(EXP_VAR, VAR_TIME);

  TaylorModel *function = picard;
  TaylorModel *substituted = substitutedField;

  /* Compute the complete Picard operator equation for each function. */
  while (function != NULL || substituted != NULL) {
    /* Enforce equal length lists. */
    assert((function != NULL) == (substituted != NULL));
    /* Ensure the function ordering was not messed up. */
    assert(strcmp(function->fun, substituted->fun) == 0);

    /* x0 + integral_0^t (f(g(t), t) dt) */
    ExpTree *x0 = function->exp;
    ExpTree *fg = substituted->exp;
    ExpTree *add =
        newExpOp(EXP_ADD_OP, x0, definiteIntegral(fg, VAR_TIME, zero, t));
    function->exp = add;

    function = function->next;
    substituted = substituted->next;
  }

  /* Clean */
  delTaylorModel(substitutedField);
  delExpTree(zero);
  delExpTree(t);
  return picard;
}

TaylorModel *picardOperatorTM(ODEList *vectorField, TaylorModel *functions) {
  assert(vectorField != NULL);
  assert(functions != NULL);

  // TODO: UNTESTED CODE: NEITHER MANUAL NOR AUTOMATED TESTS DONE

  /* Pf(g) = x0 + integral_0^t ( f(g(s), s) ds )
    Initialize the result to x0. */
  TaylorModel *picard = initTaylorModel(vectorField);
  TaylorModel *substitutedField = substituteTaylorModel(vectorField, functions);
  ExpTree *zero = newExpLeaf(EXP_NUM, "0");
  ExpTree *t = newExpLeaf(EXP_VAR, VAR_TIME);

  TaylorModel *function = picard;
  TaylorModel *substituted = substitutedField;

  /* Compute the complete Picard operator equation for each function. */
  while (function != NULL || substituted != NULL) {
    /* Enforce equal length lists. */
    assert((function != NULL) == (substituted != NULL));
    /* Ensure the function ordering was not messed up. */
    assert(strcmp(function->fun, substituted->fun) == 0);

    /* x0 + integral_0^t (f(g(t), t) dt) */
    ExpTree *x0 = function->exp;
    ExpTree *fg = cpyExpTree(substituted->exp);
    ExpTree *add =
        newExpOp(EXP_ADD_OP, x0, definiteIntegral(fg, VAR_TIME, zero, t));
    function->exp = add;

    function = function->next;
    substituted = substituted->next;
  }

  /* Clean */
  delTaylorModel(substitutedField);
  delExpTree(zero);
  delExpTree(t);
  return picard;
}

TaylorModel *substituteTaylorModel(ODEList *system, TaylorModel *functions) {
  assert(system != NULL);
  assert(functions != NULL);

  ExpTree *substituted = NULL;
  TaylorModel *function = functions;
  /* Substitute all functions into the current system component. */
  while (function != NULL) {
    ExpTree *source = substituted ? substituted : system->exp;
    substituted = substitute(source, function->fun, function->exp);

    /* Delete the old/intermediate tree, which always exists after the first
     * iteration. */
    if (source != system->exp)
      delExpTree(source);

    function = function->next;
  }

  char *fun = strdup(system->fun);
  Interval remainder = newInterval(0, 0);
  /* Base case: The tail/next of the last element is NULL. */
  if (system->next == NULL)
    return newTMElem(NULL, fun, substituted, remainder);
  /* Recursive case: The tail of the new element is everything built until now.
   */
  else
    return newTMElem(substituteTaylorModel(system->next, functions), fun,
                     substituted, remainder);
}

TaylorModel *initTaylorModel(ODEList *system) {
  /* Base case: The tail/next of the last element is NULL. */
  if (system == NULL)
    return NULL;

  /* Recursive case: The tail of the new element is everything built until now.
   */
  char *fun = strdup(system->fun);
  ExpTree *exp = newExpLeaf(EXP_VAR, system->fun);
  Interval remainder = newInterval(0, 0);
  return newTMElem(initTaylorModel(system->next), fun, exp, remainder);
}

Interval computeSafeRemainder(TaylorModel *polynomials) {
  return newInterval(0, 0);
}
