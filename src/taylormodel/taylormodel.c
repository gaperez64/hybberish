#include "taylormodel.h"

TaylorModel *newTaylorModel(char *const fun, ExpTree *const exp,
                            Interval remainder) {
  assert(fun != NULL);
  assert(exp != NULL);

  TaylorModel *list = (TaylorModel *)malloc(sizeof(TaylorModel));
  list->fun = fun;
  list->exp = exp;
  list->remainder = remainder;
  list->next = NULL;
  return list;
}

TaylorModel *appTMElem(TaylorModel *const tail, TaylorModel *head) {
  assert(head != NULL);
  assert(head->next == NULL);

  head->next = tail;
  return head;
}

TaylorModel *newTMElem(TaylorModel *const tail, char *const fun,
                       ExpTree *const exp, const Interval remainder) {
  TaylorModel *head = newTaylorModel(fun, exp, remainder);
  return appTMElem(tail, head);
}

void delTaylorModel(TaylorModel *list) {
  if (list->next != NULL)
    delTaylorModel(list->next);
  assert(list->fun != NULL);
  free(list->fun);
  assert(list->exp != NULL);
  delExpTree(list->exp);
  free(list);
}

void printTaylorModel(const TaylorModel *const list, FILE *where) {
  assert(list->fun != NULL);
  fprintf(where, "(p(%s) = ", list->fun);
  assert(list->exp != NULL);

  printExpTree(list->exp, where);
  fprintf(where, ", ");
  printInterval(&list->remainder, where);
  fprintf(where, ")");

  fprintf(where, "; ");
  if (list->next != NULL)
    printTaylorModel(list->next, where);
}

TaylorModel *cpyTaylorModel(const TaylorModel *const list) {
  /* Base case: The tail/next of the last element is NULL. */
  if (list == NULL)
    return NULL;

  /* Recursive case: The tail of the new element is everything built until now.
   */
  return appTMElem(cpyTaylorModel(list->next), cpyTaylorModelHead(list));
}

TaylorModel *cpyTaylorModelHead(const TaylorModel *const list) {
  char *fun = strdup(list->fun);
  ExpTree *exp = cpyExpTree(list->exp);
  return newTaylorModel(fun, exp, list->remainder);
}

TaylorModel *reverseTaylorModel(TaylorModel *const list) {
  /* Enforce preconditions */
  assert(list != NULL);

  /* Base case: push up the last element. */
  if (list->next == NULL)
    return list;

  /* Recursive case: Reverse tail and retrieve new head. */
  TaylorModel *lastElem = reverseTaylorModel(list->next);

  /* Only update upon backtracking, else would recurse on wrong or NULL ptrs. */
  assert(list->next != NULL);
  list->next->next = list; // reverse link: A->B  =>  A<-B
  list->next = NULL;       // sever own next: A->B  =>  A->NULL

  return lastElem;
}

Interval evaluateExpTree(const ExpTree *const tree,
                         const Domain *const domains) {
  assert(domains != NULL);
  assert(tree != NULL);

  switch (tree->type) {
  case EXP_NUM: {
    assert(tree->left == NULL);
    assert(tree->right == NULL);
    assert(tree->data != NULL);

    double value = atof(tree->data);
    return newInterval(value, value);
  }

  case EXP_VAR: {
    assert(tree->left == NULL);
    assert(tree->right == NULL);
    assert(tree->data != NULL);

    const Domain *dom = domains;
    while (dom != NULL) {
      assert(dom->var != NULL);

      if (strcmp(tree->data, dom->var) == 0)
        return dom->domain;

      dom = dom->next;
    }

    /* The expression tree contains a variable whose valuation is unknown. */
    assert(false);
    break;
  }

  case EXP_ADD_OP: {
    assert(tree->left != NULL);
    assert(tree->right != NULL);

    Interval left = evaluateExpTree(tree->left, domains);
    Interval right = evaluateExpTree(tree->right, domains);
    return addInterval(&left, &right);
  }

  case EXP_SUB_OP: {
    assert(tree->left != NULL);
    assert(tree->right != NULL);

    Interval left = evaluateExpTree(tree->left, domains);
    Interval right = evaluateExpTree(tree->right, domains);
    return subInterval(&left, &right);
  }

  case EXP_MUL_OP: {
    assert(tree->left != NULL);
    assert(tree->right != NULL);

    Interval left = evaluateExpTree(tree->left, domains);
    Interval right = evaluateExpTree(tree->right, domains);
    return mulInterval(&left, &right);
  }

  case EXP_DIV_OP: {
    assert(tree->left != NULL);
    assert(tree->right != NULL);

    Interval left = evaluateExpTree(tree->left, domains);
    Interval right = evaluateExpTree(tree->right, domains);
    return divInterval(&left, &right);
  }

  case EXP_NEG: {
    assert(tree->left != NULL);
    assert(tree->right == NULL);

    Interval left = evaluateExpTree(tree->left, domains);
    return negInterval(&left);
  }

  case EXP_EXP_OP: {
    assert(tree->left != NULL);
    assert(tree->right != NULL);

    /* Assume the exponent is always a natural number. */
    assert(tree->right->type == EXP_NUM);
    unsigned int exponent = (unsigned int)round(atof(tree->right->data));

    Interval left = evaluateExpTree(tree->left, domains);
    return pow2Interval(&left, exponent);
  }

  case EXP_FUN: {
    assert(tree->left != NULL);
    assert(tree->right == NULL);
    assert(tree->data != NULL);

    Interval left = evaluateExpTree(tree->left, domains);

    if (strcmp(tree->data, "sqrt") == 0)
      return sqrtInterval(&left);

    /* Unknown function, abort */
    assert(false);
    break;
  }

  /* Unknown operator or leaf to evaluate. */
  default:
    assert(false);
    break;
  }
}

double evaluateExpTreeReal(const ExpTree *const tree,
                           const Valuation *const values) {
  assert(values != NULL);
  assert(tree != NULL);

  switch (tree->type) {
  case EXP_NUM: {
    assert(tree->left == NULL);
    assert(tree->right == NULL);
    assert(tree->data != NULL);

    return atof(tree->data);
  }

  case EXP_VAR: {
    assert(tree->left == NULL);
    assert(tree->right == NULL);
    assert(tree->data != NULL);

    const Valuation *val = values;
    while (val != NULL) {
      assert(val->var != NULL);

      if (strcmp(tree->data, val->var) == 0)
        return val->val;

      val = val->next;
    }

    /* The expression tree contains a variable whose valuation is unknown. */
    assert(false);
    break;
  }

  case EXP_ADD_OP: {
    assert(tree->left != NULL);
    assert(tree->right != NULL);

    double left = evaluateExpTreeReal(tree->left, values);
    double right = evaluateExpTreeReal(tree->right, values);
    return left + right;
  }

  case EXP_SUB_OP: {
    assert(tree->left != NULL);
    assert(tree->right != NULL);

    double left = evaluateExpTreeReal(tree->left, values);
    double right = evaluateExpTreeReal(tree->right, values);
    return left - right;
  }

  case EXP_MUL_OP: {
    assert(tree->left != NULL);
    assert(tree->right != NULL);

    double left = evaluateExpTreeReal(tree->left, values);
    double right = evaluateExpTreeReal(tree->right, values);
    return left * right;
  }

  case EXP_DIV_OP: {
    assert(tree->left != NULL);
    assert(tree->right != NULL);

    double left = evaluateExpTreeReal(tree->left, values);
    double right = evaluateExpTreeReal(tree->right, values);
    return left / right;
  }

  case EXP_NEG: {
    assert(tree->left != NULL);
    assert(tree->right == NULL);

    double left = evaluateExpTreeReal(tree->left, values);
    return -left;
  }

  case EXP_EXP_OP: {
    assert(tree->left != NULL);
    assert(tree->right != NULL);

    /* Assume the exponent is always a natural number. */
    assert(tree->right->type == EXP_NUM);
    unsigned int exponent = (unsigned int)round(atof(tree->right->data));

    double left = evaluateExpTreeReal(tree->left, values);
    return pow(left, exponent);
  }

  case EXP_FUN: {
    assert(tree->left != NULL);
    assert(tree->right == NULL);
    assert(tree->data != NULL);

    double left = evaluateExpTreeReal(tree->left, values);

    if (strcmp(tree->data, "sqrt") == 0)
      return sqrt(left);

    /* Unknown function, abort */
    assert(false);
    break;
  }

  /* Unknown operator or leaf to evaluate. */
  default:
    assert(false);
    break;
  }
}

TaylorModel *evaluateExpTreeTM(const ExpTree *const tree,
                               const TaylorModel *const list,
                               const char *const fun,
                               const Domain *const variables,
                               const unsigned int k) {
  assert(list != NULL);
  assert(tree != NULL);
  assert(fun != NULL);

  switch (tree->type) {
  case EXP_NUM: {
    assert(tree->left == NULL);
    assert(tree->right == NULL);
    assert(tree->data != NULL);

    ExpTree *exp = cpyExpTree(tree);
    Interval remainder = newInterval(0, 0);
    TaylorModel *num = newTaylorModel(strdup(fun), exp, remainder);
    return num;
  }

  case EXP_VAR: {
    assert(tree->left == NULL);
    assert(tree->right == NULL);
    assert(tree->data != NULL);

    const TaylorModel *tm = list;
    while (tm != NULL) {
      assert(tm->fun != NULL);

      if (strcmp(tree->data, tm->fun) == 0) {
        /* Reorient the copied TM to correspond to the target fun/var. */
        TaylorModel *copied = cpyTaylorModelHead(tm);
        free(copied->fun);
        copied->fun = strdup(fun);
        return copied;
      }

      tm = tm->next;
    }

    /* The expression tree contains a variable without corresponding TM. */
    assert(false);
    break;
  }

  case EXP_ADD_OP: {
    assert(tree->left != NULL);
    assert(tree->right != NULL);

    TaylorModel *left = evaluateExpTreeTM(tree->left, list, fun, variables, k);
    TaylorModel *right =
        evaluateExpTreeTM(tree->right, list, fun, variables, k);
    TaylorModel *binop = addTM(left, right, variables, k);
    delTaylorModel(left);
    delTaylorModel(right);
    return binop;
  }

  case EXP_SUB_OP: {
    assert(tree->left != NULL);
    assert(tree->right != NULL);

    TaylorModel *left = evaluateExpTreeTM(tree->left, list, fun, variables, k);
    TaylorModel *right =
        evaluateExpTreeTM(tree->right, list, fun, variables, k);
    TaylorModel *binop = subTM(left, right, variables, k);
    delTaylorModel(left);
    delTaylorModel(right);
    return binop;
  }

  case EXP_MUL_OP: {
    assert(tree->left != NULL);
    assert(tree->right != NULL);

    TaylorModel *left = evaluateExpTreeTM(tree->left, list, fun, variables, k);
    TaylorModel *right =
        evaluateExpTreeTM(tree->right, list, fun, variables, k);
    TaylorModel *binop = mulTM(left, right, variables, k);
    delTaylorModel(left);
    delTaylorModel(right);
    return binop;
  }

  case EXP_DIV_OP: {
    assert(tree->left != NULL);
    assert(tree->right != NULL);

    TaylorModel *left = evaluateExpTreeTM(tree->left, list, fun, variables, k);
    TaylorModel *right =
        evaluateExpTreeTM(tree->right, list, fun, variables, k);
    TaylorModel *binop = divTM(left, right, variables, k);
    delTaylorModel(left);
    delTaylorModel(right);
    return binop;
  }

  case EXP_NEG: {
    assert(tree->left != NULL);
    assert(tree->right == NULL);

    TaylorModel *left = evaluateExpTreeTM(tree->left, list, fun, variables, k);
    TaylorModel *unop = negTM(left, variables, k);
    delTaylorModel(left);
    return unop;
  }

  case EXP_EXP_OP: {
    assert(tree->left != NULL);
    assert(tree->right != NULL);
    assert(tree->right->data != NULL);

    const unsigned int exponent = atou(tree->right->data);
    TaylorModel *left = evaluateExpTreeTM(tree->left, list, fun, variables, k);
    TaylorModel *binop = powTM(left, exponent, variables, k);
    delTaylorModel(left);

    return binop;
  }

  case EXP_FUN: {
    assert(tree->left != NULL);
    assert(tree->right == NULL);
    assert(tree->data != NULL);

    /* Unknown function, abort */
    assert(false);
    return NULL;
  }

  /* Unknown operator or leaf to evaluate. */
  default:
    assert(false);
    return NULL;
  }

  assert(false);
  return NULL;
}

TaylorModel *addTM(const TaylorModel *const left,
                   const TaylorModel *const right,
                   const Domain *const variables, const unsigned int k) {
  /* Require equal length lists: if only one is NULL
    then there is a list length mismatch. */
  assert((left == NULL) == (right == NULL));

  /* Base case: The tail/next of the last element is NULL. */
  if (left == NULL && right == NULL)
    return NULL;

  /* Only compose TMs that correspond to the same variable. */
  assert(left->fun != NULL && right->fun != NULL);
  assert(strcmp(left->fun, right->fun) == 0);

  /* Recursive case: The tail of the new element is everything built until now.
    (p1, I1) + (p2, I2) = (p1 + p2, I1 + I2) */
  char *fun = strdup(left->fun);
  ExpTree *exp =
      newExpOp(EXP_ADD_OP, cpyExpTree(left->exp), cpyExpTree(right->exp));
  Interval remainder = addInterval(&left->remainder, &right->remainder);
  TaylorModel *binaryOp = newTaylorModel(fun, exp, remainder);
  TaylorModel *truncated = truncateTM(binaryOp, variables, k);

  /* Clean */
  delTaylorModel(binaryOp);

  return appTMElem(addTM(left->next, right->next, variables, k), truncated);
}

TaylorModel *subTM(const TaylorModel *const left,
                   const TaylorModel *const right,
                   const Domain *const variables, const unsigned int k) {
  /* Require equal length lists: if only one is NULL
    then there is a list length mismatch. */
  assert((left == NULL) == (right == NULL));

  /* Base case: The tail/next of the last element is NULL. */
  if (left == NULL && right == NULL)
    return NULL;

  /* Only compose TMs that correspond to the same variable. */
  assert(left->fun != NULL && right->fun != NULL);
  assert(strcmp(left->fun, right->fun) == 0);

  /* Recursive case: The tail of the new element is everything built until now.
    (p1, I1) - (p2, I2) = (p1 - p2, I1 - I2) */
  char *fun = strdup(left->fun);
  ExpTree *exp =
      newExpOp(EXP_SUB_OP, cpyExpTree(left->exp), cpyExpTree(right->exp));
  Interval remainder = subInterval(&left->remainder, &right->remainder);
  TaylorModel *binaryOp = newTaylorModel(fun, exp, remainder);
  TaylorModel *truncated = truncateTM(binaryOp, variables, k);

  /* Clean */
  delTaylorModel(binaryOp);

  return appTMElem(subTM(left->next, right->next, variables, k), truncated);
}

TaylorModel *mulTM(const TaylorModel *const left,
                   const TaylorModel *const right,
                   const Domain *const variables, const unsigned int k) {
  /* Require equal length lists: if only one is NULL
    then there is a list length mismatch. */
  assert((left == NULL) == (right == NULL));

  /* Base case: The tail/next of the last element is NULL. */
  if (left == NULL && right == NULL)
    return NULL;

  /* Only compose TMs that correspond to the same variable. */
  assert(left->fun != NULL && right->fun != NULL);
  assert(strcmp(left->fun, right->fun) == 0);

  /* Recursive case: The tail of the new element is everything built until now.
    (p1, I1) * (p2, I2)
  = (p1 * p2 - pe, Int(pe) + Int(p1)*I2 + Int(p2)*I1 + I1*I2) */
  char *fun = strdup(left->fun);
  ExpTree *exp =
      newExpOp(EXP_MUL_OP, cpyExpTree(left->exp), cpyExpTree(right->exp));
  ExpTree *sumOfProds = toSumOfProducts(exp);

  Interval Intp1 = evaluateExpTree(left->exp, variables);
  Interval Intp2 = evaluateExpTree(right->exp, variables);
  Interval p1I2 = mulInterval(&Intp1, &right->remainder);
  Interval p2I1 = mulInterval(&Intp2, &left->remainder);
  Interval I1I2 = mulInterval(&left->remainder, &right->remainder);
  Interval remainder;
  remainder = addInterval(&p1I2, &p2I1);
  remainder = addInterval(&remainder, &I1I2);

  TaylorModel *binaryOp = newTaylorModel(fun, sumOfProds, remainder);
  TaylorModel *truncated = truncateTM(binaryOp, variables, k);

  /* Clean */
  delExpTree(exp);
  delTaylorModel(binaryOp);

  return appTMElem(mulTM(left->next, right->next, variables, k), truncated);
}

TaylorModel *divTM(const TaylorModel *const left,
                   const TaylorModel *const right,
                   const Domain *const variables, const unsigned int k) {
  /* Require equal length lists: if only one is NULL
    then there is a list length mismatch. */
  assert((left == NULL) == (right == NULL));

  /* Base case: The tail/next of the last element is NULL. */
  if (left == NULL && right == NULL)
    return NULL;

  /* Only compose TMs that correspond to the same variable. */
  assert(left->fun != NULL && right->fun != NULL);
  assert(strcmp(left->fun, right->fun) == 0);

  /* Int((p2, I2)) = Int(p2) + I2 */
  Interval enclosure = evaluateExpTree(right->exp, variables);
  enclosure = addInterval(&enclosure, &right->remainder);
  assert(!elemInterval(0, &enclosure));

  /* c = Mid(Int((p2, I2))) */
  double c = intervalMidpoint(&enclosure);
  char cStr[50];
  dtoa(cStr, sizeof(cStr), c);

  /* b = 1 / c  where c != 0 since 0 not in Int((p2, I2)) */
  double b = 1 / c;
  char bStr[50];
  dtoa(bStr, sizeof(bStr), b);

  /* pk(x) = 1/c * (1 - ((x - c) / c)^1 + ... + (-1)^k * ((x - c) / c)^k)
           =   b * (1 - ((x - c) * b)^1 + ... + (-1)^k * ((x - c) * b)^k) */
  ExpTree *inverseExp = newExpLeaf(EXP_NUM, "1");
  char *fun = left->fun;
  for (unsigned int it = 1; it <= k; ++it) {
    char expStr[50];
    utoa(expStr, sizeof(expStr), it);
    ExpTree *exponent = newExpLeaf(EXP_NUM, expStr);

    /* ((x - c) / c)^i
      BUT should not use DIV nodes, to avoid infinite recursion in evaluation.
    => ((x - c) * b)^i  should be used instead, where b = 1/c is evaluated. */
    ExpTree *cLeaf = newExpLeaf(EXP_NUM, cStr);
    ExpTree *bLeaf = newExpLeaf(EXP_NUM, bStr);
    ExpTree *xLeaf = newExpLeaf(EXP_VAR, fun);
    ExpTree *term = newExpOp(EXP_SUB_OP, xLeaf, cLeaf);
    term = newExpOp(EXP_MUL_OP, term, bLeaf);
    term = newExpOp(EXP_EXP_OP, term, exponent);

    /* Aggregate term into expression */
    ExpType opType = ((it % 2) == 0) ? EXP_ADD_OP : EXP_SUB_OP;
    inverseExp = newExpOp(opType, inverseExp, term);
  }
  ExpTree *bLeaf = newExpLeaf(EXP_NUM, bStr);
  inverseExp = newExpOp(EXP_MUL_OP, bLeaf, inverseExp);

  /* Setup leaves for remainder expression. */
  char expk1Str[12];
  char expk2Str[12];
  utoa(expk1Str, sizeof(expk1Str), k + 1);
  utoa(expk2Str, sizeof(expk2Str), k + 2);
  ExpTree *expk1 = newExpLeaf(EXP_NUM, expk1Str);
  ExpTree *expk2 = newExpLeaf(EXP_NUM, expk2Str);
  ExpTree *one = newExpLeaf(EXP_NUM, "1");
  ExpTree *cLeaf = newExpLeaf(EXP_NUM, cStr);
  ExpTree *xLeaf = newExpLeaf(EXP_VAR, fun);

  /* Compose TM (p2 - c, I2) */
  TaylorModel *rightMod = cpyTaylorModelHead(right);
  rightMod->exp = newExpOp(EXP_SUB_OP, rightMod->exp, cLeaf);

  /* Compose remainder expression. */
  /* 1 / x^(k+2) */
  ExpTree *factor1 = newExpOp(EXP_EXP_OP, xLeaf, expk2);
  factor1 = newExpOp(EXP_DIV_OP, one, factor1);
  ExpTree *sub = newExpOp(EXP_SUB_OP, cpyExpTree(xLeaf), cpyExpTree(cLeaf));
  /* (x - c)^(k+1) */
  ExpTree *factor2 = newExpOp(EXP_EXP_OP, sub, expk1);
  ExpTree *remExp = newExpOp(EXP_MUL_OP, factor1, factor2);
  /* (-1)^(k+1) * ... */
  if (((k + 1) % 2) == 1)
    remExp = newExpOp(EXP_NEG, remExp, NULL);

  /* Compute (p3, I3) by substituting (p2 - c, I2) for x in pk(x) */
  Interval inverseRemainder = evaluateExpTree(remExp, variables);
  printf("$$$$$$$$$$$$$$$$$$\\\n");
  fflush(stdout);
  TaylorModel *inverseTM =
      evaluateExpTreeTM(inverseExp, rightMod, right->fun, variables, k);
  printf("$$$$$$$$$$$$$$$$$$/\n");
  fflush(stdout);

  // TODO: /\/\/\ evaluateExpTM recursively calls divTM ???? Because DIV nodes are used above /\/\/\ 

  inverseTM->remainder = addInterval(&inverseTM->remainder, &inverseRemainder);

  /* Compose TM (p1, I1) / (p2, I2) = (p1, I1) * (p3, I3) */
  TaylorModel *binaryOp = mulTM(left, inverseTM, variables, k);
  TaylorModel *truncated = truncateTM(binaryOp, variables, k);

  /* Clean */
  delExpTree(remExp);
  delExpTree(inverseExp);
  delTaylorModel(rightMod);
  delTaylorModel(inverseTM);
  delTaylorModel(binaryOp);

  return appTMElem(divTM(left->next, right->next, variables, k), truncated);
}

TaylorModel *negTM(const TaylorModel *const list, const Domain *const variables,
                   const unsigned int k) {
  /* Base case: The tail/next of the last element is NULL. */
  if (list == NULL)
    return NULL;

  assert(list->fun != NULL);

  /* Recursive case: The tail of the new element is everything built until now.
    -(p, I) = (-p, -I) */
  char *fun = strdup(list->fun);
  ExpTree *exp = newExpOp(EXP_NEG, cpyExpTree(list->exp), NULL);
  Interval remainder = negInterval(&list->remainder);
  TaylorModel *unaryOp = newTaylorModel(fun, exp, remainder);
  TaylorModel *truncated = truncateTM(unaryOp, variables, k);

  /* Clean */
  delTaylorModel(unaryOp);

  return appTMElem(negTM(list->next, variables, k), truncated);
}

TaylorModel *powTM(const TaylorModel *const left, const unsigned int right,
                   const Domain *const variables, const unsigned int k) {
  /* Consistency with other TM arithmetic: NULL^k = NULL. */
  if (left == NULL)
    return NULL;

  /* For simplicity, disallow 0 exponent. */
  assert(right > 0);

  /* Unroll the integer exponent into successive multiplications.
    (p, I)^n = (p, I) * ... * (p, I) */
  TaylorModel *binaryOp = (TaylorModel *)left;
  for (unsigned int it = 1; it < right; ++it) {
    TaylorModel *intermediate = mulTM(left, binaryOp, variables, k);

    /* Clean: delete, except the first (input) iteration. */
    if (1 < it)
      delTaylorModel(binaryOp);

    binaryOp = intermediate;
  }
  TaylorModel *truncated = truncateTM(binaryOp, variables, k);

  /* Clean: except for exponent = 1, then this is the input TM. */
  if (right > 1)
    delTaylorModel(binaryOp);

  return truncated;
}

TaylorModel *intTM(const TaylorModel *const list,
                   const Interval *const intDomain, const char *const intVar,
                   const Domain *const variables, const unsigned int k) {
  assert(intDomain != NULL);
  assert(intVar != NULL);
  assert(variables != NULL);

  /* Base case: The tail/next of the last element is NULL. */
  if (list == NULL)
    return NULL;

  /* Definite integral bounds should be trees. */
  char lowerBoundStr[50];
  char upperBoundStr[50];
  snprintf(lowerBoundStr, sizeof(lowerBoundStr), "%.15g", intDomain->left);
  snprintf(upperBoundStr, sizeof(upperBoundStr), "%.15g", intDomain->right);
  ExpTree *lowerBound = newExpLeaf(EXP_NUM, lowerBoundStr);
  ExpTree *upperBound = newExpLeaf(EXP_NUM, upperBoundStr);

  /* Integration raises the degree of every single term by exactly one,
    so truncation before integration of terms of degree gte k is more efficient.
  */
  ExpTree *truncatedTerms = NULL;
  ExpTree *truncated = truncate2(list->exp, k - 1, &truncatedTerms);
  ExpTree *exp = definiteIntegral(truncated, intVar, lowerBound, upperBound);

  /* Il = (Int(pe) + I) * [ai, bi] */
  Interval enclosure; // Interval enclosure Int(pe)
  Interval remainder; // Updated remainder interval Il
  enclosure = evaluateExpTree(truncatedTerms, variables);
  remainder = addInterval(&enclosure, &list->remainder);
  remainder = mulInterval(&remainder, intDomain);

  /* Recursive case: The tail of the new element is everything built until now.
   */
  return newTMElem(intTM(list->next, intDomain, intVar, variables, k),
                   strdup(intVar), exp, remainder);
}

TaylorModel *truncateTM(const TaylorModel *const list,
                        const Domain *const variables, const unsigned int k) {
  /* Base case: The tail/next of the last element is NULL. */
  if (list == NULL)
    return NULL;

  assert(list->fun != NULL);
  assert(variables != NULL);

  /* Recursive case: The tail of the new element is everything built until now.
    trunc((p, I) = (p - pe, I + Int(pe))) where pe are the truncated terms and
    Int(pe) is their interval enclosure. */
  char *fun = strdup(list->fun);
  ExpTree *truncatedTerms = NULL;
  ExpTree *truncated = truncate2(list->exp, k, &truncatedTerms);
  Interval enclosure = (truncatedTerms != NULL)
                           ? evaluateExpTree(truncatedTerms, variables)
                           : newInterval(0, 0);
  Interval remainder = addInterval(&list->remainder, &enclosure);

  /* Clean. */
  if (truncatedTerms != NULL)
    delExpTree(truncatedTerms);

  return newTMElem(truncateTM(list->next, variables, k), fun, truncated,
                   remainder);
}
