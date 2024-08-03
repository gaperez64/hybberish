#include "tmflowpipe.h"

TPList *computeTaylorPolynomial(ODEList *system, unsigned int order,
                                unsigned int k) {
  /* Enforce preconditions */
  assert(system != NULL);
  assert(order > 0);
  assert(k > 0);
  /* All the polynomial terms exceeding the truncation order
    would just get truncated anyways, so impose an explicit limit. */
  assert(order <= k);

  /* The functions to seed each Lie derivation with. */
  TPList *lieDerivativeSeed = initTPList(system);
  /* The Taylor polynomials that will be built up in-place. */
  TPList *polynomials = initTPList(system);

  /* Start from i=1; case i=0 would be an order 0 Lie derivative. */
  for (unsigned int index = 1; index <= order; ++index) {
    TPList *lieDeriv = lieDerivativeK(system, lieDerivativeSeed, index);

    TPList *poly = polynomials;
    TPList *deriv = lieDeriv;
    while (poly != NULL) {
      /* If one becomes NULL while the other does
        not, then there is a length mismatch. */
      assert(poly != NULL && deriv != NULL);

      /* Ensure each variable's derivative is added to that same variable's
        running Taylor polynomial. If this assertion fails, the variable
        ordering in the lists was messed up somehow. */
      assert(strcmp(poly->fun, deriv->fun) == 0);

      /* Make a copy to ensure the tree is decoupled from the derivative list
       * object. */
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
      ExpTree *tPow =
          newExpOp(EXP_EXP_OP, newExpLeaf(EXP_VAR, strdup(VAR_TIME)),
                   newExpLeaf(EXP_NUM, indexStr));
      /* 1/i! * L^i(g) * t^i   where L^i(g) is the order i Lie derivative of
       * function g. */
      ExpTree *polyElement =
          newExpOp(EXP_MUL_OP, fac, newExpOp(EXP_MUL_OP, derivExp, tPow));

      /* Extend the polynomial by adding the next element. */
      poly->exp = newExpOp(EXP_ADD_OP, poly->exp, polyElement);

      poly = poly->next;
      deriv = deriv->next;
    }

    /* Cleanup */
    delTPList(lieDeriv);
  }

  /* Cleanup */
  delTPList(lieDerivativeSeed);

  return polynomials;
}

TPList *lieDerivativeK(ODEList *system, TPList *functions, unsigned int order) {
  /* Enforce preconditions */
  assert(system != NULL);
  assert(functions != NULL);

  /* Copy so that the returned list's lifetime is always the caller's
   * responsibility. */
  if (order == 0)
    return cpyTPList(functions);

  TPList *ithLieDerivative = functions;
  for (unsigned int index = 0; index < order; ++index) {
    TPList *old = ithLieDerivative;
    ithLieDerivative = lieDerivativeTPList(system, ithLieDerivative);

    /* Cleanup: free intermediate derivative results;
      do not free at idx=0 since that is the input,
      and do not free at idx=order-1 since that is
      the output. */
    if (index > 0 && index < (order - 1)) {
      delTPList(old);
    }
  }

  return ithLieDerivative;
}

TPList *lieDerivativeTPList(ODEList *system, TPList *functions) {
  /* Enforce preconditions */
  assert(system != NULL);
  assert(functions != NULL);

  ODEList *ode = system;
  TPList *function = functions;
  TPList *derived = NULL;

  /* Derive each of the functions individually w.r.t. the same ODe system. */
  while (!(ode == NULL || function == NULL)) {
    /* If only one (XOR) is NULL but not the other, then there is a list length
      mismatch. Theoretically, they could have different lengths, but in our
      case their lengths are always the same. This assertion safeguards that
      loose assumption. */
    assert((ode != NULL) == (function != NULL));

    char *fun = strdup(function->fun);
    ExpTree *lieDeriv = lieDerivative(system, function->exp);

    /* Post-process the true Lie derivative result. */
    ExpTree *simplified = simplify(lieDeriv);
    delExpTree(lieDeriv);

    /* The initial tail should be NULL, since the list is extended head-first.
     */
    derived = newTPElem(derived, fun, simplified);

    ode = ode->next;
    function = function->next;
  }

  /* Reverse to ensure the output functions are
    ordered the same as the input functions. */
  return reverseTPList(derived);
}

ExpTree *lieDerivative(ODEList *vectorField, ExpTree *function) {
  /* Enforce preconditions */
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

TPList *initTPList(ODEList *system) {
  /* Base case: The tail/next of the last element is NULL. */
  if (system == NULL)
    return NULL;

  /* Recursive case: The tail of the new element is everything built until now.
   */
  char *fun = strdup(system->fun);
  ExpTree *exp = newExpLeaf(EXP_VAR, strdup(system->fun));
  return newTPElem(initTPList(system->next), fun, exp);
}

Interval computeSafeRemainder(TPList *polynomials) { return (Interval){0, 0}; }
