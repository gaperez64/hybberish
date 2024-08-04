#include "transformations.h"



ExpTree* simplify(ExpTree *source) {
  ExpTree *simplified = simplifyOperators(source);

  return simplified;
}


ExpTree* toSumOfProducts(ExpTree *source) {
  assert(source != NULL);

  /* Base case: retain leaves. */
  if (source->left == NULL && source->right == NULL)
    return cpyExpTree(source);


  switch (source->type) {

  case EXP_MUL_OP: {
    assert(source->left != NULL);
    assert(source->right != NULL);

    ExpTree *leftSumOfProd = toSumOfProducts(source->left);
    ExpTree *rightSumOfProd = toSumOfProducts(source->right);

    bool leftIsDistributive = leftSumOfProd->type == EXP_ADD_OP || leftSumOfProd->type == EXP_SUB_OP;
    bool rightIsDistributive = rightSumOfProd->type == EXP_ADD_OP || rightSumOfProd->type == EXP_SUB_OP;

    /* Neither modified subtree is distributive, so simply retain them. */
    if (!leftIsDistributive && !rightIsDistributive)
      return newExpOp(source->type, leftSumOfProd, rightSumOfProd);

    /* Both are distributive, so distribute all elements of one across the other. */
    if (leftIsDistributive && rightIsDistributive) {
      ExpTree *res = distributeLeftDistributive(leftSumOfProd, rightSumOfProd);
      delExpTree(leftSumOfProd);
      delExpTree(rightSumOfProd);
      return res;
    }

    /* Only the left is distributive, so distribute the right across the left. */
    if (leftIsDistributive) {
      ExpTree *res = distributeLeft(rightSumOfProd, leftSumOfProd);
      delExpTree(leftSumOfProd);
      delExpTree(rightSumOfProd);
      return res;
    }

    /* Only the right is distributive, so distribute the left across the right. */
    ExpTree *res = distributeLeft(leftSumOfProd, rightSumOfProd);
    delExpTree(leftSumOfProd);
    delExpTree(rightSumOfProd);
    return res;
  }


  case EXP_DIV_OP:
  case EXP_ADD_OP:
  case EXP_SUB_OP:
  case EXP_EXP_OP:
  case EXP_FUN:
  case EXP_NEG: {
    ExpTree *leftSumOfProd = (source->left != NULL) ? toSumOfProducts(source->left) : NULL;
    ExpTree *rightSumOfProd = (source->right != NULL) ? toSumOfProducts(source->right) : NULL;
    char *dataCpy = (source->data != NULL) ? strdup(source->data) : NULL;
    return newExpTree(source->type, dataCpy, leftSumOfProd, rightSumOfProd);
  }


  default:
    assert(false);
    break;
  }

  return NULL;
}



/*
    Simplification Helper methods.
*/

ExpTree *simplifyOperators(ExpTree *source) {
  assert(source != NULL);

  /* Base case: retain leaves. */
  if (source->left == NULL && source->right == NULL)
    return cpyExpTree(source);

  /* Recursive case: simplify the subtrees/branches. */
  ExpTree *leftSimplified = source->left ?
    simplifyOperators(source->left) : NULL;
  ExpTree *rightSimplified = source->right ?
    simplifyOperators(source->right) : NULL;

  /* The left and/or right branch not existing (NULL ptr)
    is considered equivalent with being neutral. */
  bool leftIsNeutral  = true;
  bool rightIsNeutral = true;

  /* The left and/or right branch not existing (NULL ptr)
    is NOT considered equivalent with being absorbing. */
  bool leftIsAbsorbing  = false;
  bool rightIsAbsorbing = false;


  switch (source->type) {
  case EXP_ADD_OP:
  case EXP_SUB_OP:
    assert(source->left != NULL);
    assert(source->right != NULL);

    leftIsNeutral = isZeroExpTree(leftSimplified);
    rightIsNeutral = isZeroExpTree(rightSimplified);

    /* 0 +/- 0 = 0 */
    if (leftIsNeutral && rightIsNeutral) {
      delExpTree(rightSimplified);
      return leftSimplified;
    }

    /* a +/- b = a +/- b */
    if (!leftIsNeutral && !rightIsNeutral)
      return newExpOp(source->type, leftSimplified, rightSimplified);

    /* 0 +/- b = +/- b */
    if (leftIsNeutral) {
      delExpTree(leftSimplified);
      /* Edge case: 0 - b = -b */
      return (source->type == EXP_SUB_OP) ? newExpOp(EXP_NEG, rightSimplified, NULL) : rightSimplified;
    }

    /* a +/- 0 = a */
    delExpTree(rightSimplified);
    return leftSimplified;


  case EXP_MUL_OP:
    assert(source->left != NULL);
    assert(source->right != NULL);
    leftIsNeutral = isOneExpTree(leftSimplified);
    rightIsNeutral = isOneExpTree(rightSimplified);
    leftIsAbsorbing = isZeroExpTree(leftSimplified);
    rightIsAbsorbing = isZeroExpTree(rightSimplified);

    /* 0 * b = 0 */
    if (leftIsAbsorbing) {
      delExpTree(rightSimplified);
      return leftSimplified;
    }

    /* a * 0 = 0 */
    if (rightIsAbsorbing) {
      delExpTree(leftSimplified);
      return rightSimplified;
    }

    /* 0 + 0 = 0 */
    if (leftIsNeutral && rightIsNeutral) {
      delExpTree(rightSimplified);
      return leftSimplified;
    }
    
    /* a + b = a + b */
    if (!leftIsNeutral && !rightIsNeutral)
      return newExpOp(source->type, leftSimplified, rightSimplified);

    /* 0 + b = b */
    if (leftIsNeutral) {
      delExpTree(leftSimplified);
      return rightSimplified;
    }

    /* a + 0 = a */
    delExpTree(rightSimplified);
    return leftSimplified;


  case EXP_DIV_OP:
    assert(source->left != NULL);
    assert(source->right != NULL);
    leftIsNeutral = isOneExpTree(leftSimplified);
    rightIsNeutral = isOneExpTree(rightSimplified);
    leftIsAbsorbing = isZeroExpTree(leftSimplified);
    rightIsAbsorbing = isZeroExpTree(rightSimplified);

    /* a/0 is indeterminate */
    assert(!rightIsAbsorbing);

    /* 0/b = 0 */
    if (leftIsAbsorbing) {
      delExpTree(rightSimplified);
      return leftSimplified;
    }

    /* 1/1 = 1 */
    if (leftIsNeutral && rightIsNeutral) {
      delExpTree(rightSimplified);
      return leftSimplified;
    }

    /* a/b = a/b   and   1/b = 1/b */
    if (!rightIsNeutral)
      return newExpOp(source->type, leftSimplified, rightSimplified);

    /* a/1 = a */
    delExpTree(rightSimplified);
    return leftSimplified;


  case EXP_EXP_OP:
    assert(source->left != NULL);
    assert(source->right != NULL);

    bool leftIsZero = isZeroExpTree(leftSimplified);
    bool rightIsZero = isZeroExpTree(rightSimplified);
    bool leftIsOne = isOneExpTree(leftSimplified);
    bool rightIsOne = isOneExpTree(rightSimplified);


    /* Edge case: 0^0 is indeterminate. */
    assert(!(leftIsZero && rightIsZero));

    /* 1^0 = a^0 = 1 */
    if (rightIsZero) {
      delExpTree(leftSimplified);
      delExpTree(rightSimplified);
      return newOneExpTree();
    }

    /* 0^1 = 0^b = 0 */
    if (leftIsZero) {
      delExpTree(rightSimplified);
      return leftSimplified;
    }

    /* 1^1 = 1^b = 1 */
    if (leftIsOne) {
      delExpTree(rightSimplified);
      return leftSimplified;
    }

    /* a^1 = a */
    if (rightIsOne) {
      delExpTree(rightSimplified);
      return leftSimplified;
    }

    /* a^b */
    return newExpOp(EXP_EXP_OP, leftSimplified, rightSimplified);


  case EXP_NEG:
    assert(source->left != NULL);
    assert(source->right == NULL);

    /* -0 = 0 */
    if (isZeroExpTree(leftSimplified)) {
      return leftSimplified;
    }

    /* -a */
    return newExpOp(source->type, leftSimplified, NULL);


  case EXP_FUN:
    assert(source->left != NULL);
    assert(source->right == NULL);

    /* f(a) */
    return newExpTree(source->type, strdup(source->data), leftSimplified, NULL);


  default:
    assert(false);
    return NULL;
  }
}


bool isZeroExpTree(ExpTree *source) {
  return source != NULL && source->type == EXP_NUM && atof(source->data) == 0.0;
}

ExpTree * newZeroExpTree(void) {
  return newExpLeaf(EXP_NUM, "0");
}


bool isOneExpTree(ExpTree *source) {
  return source != NULL && source->type == EXP_NUM && atof(source->data) == 1.0;
}


ExpTree *newOneExpTree(void) {
  return newExpLeaf(EXP_NUM, "1");
}



/*
    Sum of Products Helper methods.
*/

ExpTree *distributeLeft(ExpTree *left, ExpTree *right) {
  /* Enforce preconditions */
  assert(left != NULL);
  assert(right != NULL);
  /* Right MUST be distributive. */
  assert(right->type == EXP_ADD_OP || right->type == EXP_SUB_OP);
  /* Left must NOT be distributive. */
  assert(!(left->type == EXP_ADD_OP || left->type == EXP_SUB_OP));


  ExpTree *leftSubDistributed;
  ExpTree *rightSubDistributed;

  /*   left * ((x + y) + ...) = left * (x + y) + ...   */
  if (right->left->type == EXP_ADD_OP || right->left->type == EXP_SUB_OP)
    leftSubDistributed = distributeLeft(left, right->left);
  else
    leftSubDistributed = newExpOp(EXP_MUL_OP, cpyExpTree(left), cpyExpTree(right->left));

  /*   left * (... + (y + z)) = left * ... + left * (y + z)   */
  if (right->right->type == EXP_ADD_OP || right->right->type == EXP_SUB_OP)
    rightSubDistributed = distributeLeft(left, right->right);
  else
    rightSubDistributed = newExpOp(EXP_MUL_OP, cpyExpTree(left), cpyExpTree(right->right));

  return newExpOp(right->type, leftSubDistributed, rightSubDistributed);
}


ExpTree *distributeLeftDistributive(ExpTree *left, ExpTree *right) {
  /* Enforce preconditions */
  assert(left != NULL);
  assert(right != NULL);
  /* Right MUST be distributive. */
  assert(right->type == EXP_ADD_OP || right->type == EXP_SUB_OP);
  /* Left MUST be distributive. */
  assert(left->type == EXP_ADD_OP || left->type == EXP_SUB_OP);


  ExpTree *leftDistributed;
  ExpTree *rightDistributed;

  /* The left->left subtree contains more terms to distribute over right. */
  if (left->left->type == EXP_ADD_OP || left->left->type == EXP_SUB_OP)
    leftDistributed = distributeLeftDistributive(left->left, right);
  /* Base case: distribute a non-distributive term across right. */
  else
    leftDistributed = distributeLeft(left->left, right);

  /* The left->right subtree contains more terms to distribute over right. */
  if (left->right->type == EXP_ADD_OP || left->right->type == EXP_SUB_OP)
    rightDistributed = distributeLeftDistributive(left->right, right);
  /* Base case: distribute a non-distributive term across right. */
  else
    rightDistributed = distributeLeft(left->right, right);

  /* Compose the distributions of all non-distributive subtrees. */
  return newExpOp(left->type, leftDistributed, rightDistributed);
}
