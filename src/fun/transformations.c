#include "transformations.h"



ExpTree* simplify(ExpTree *source) {
  ExpTree *simplified;
  ExpTree *temp = simplifyNeutralTerms(source);
  simplified = temp;

  temp = simplifyAbsorbingTerms(simplified);
  delExpTree(simplified);
  simplified = temp;

  return simplified;
}


ExpTree *simplifyNeutralTerms(ExpTree *source) {
  assert(source != NULL);

  /* Base case: retain leaves. */
  if (source->left == NULL && source->right == NULL)
    return cpyExpTree(source);

  /* Recursive case: determine the neutrality of the subtrees/branches. */
  ExpTree *leftSimplified = source->left ?
    simplifyNeutralTerms(source->left) : NULL;
  ExpTree *rightSimplified = source->right ?
    simplifyNeutralTerms(source->right) : NULL;

  /* The left and/or right branch not existing (NULL ptr)
    is considered equivalent with being neutral. */
  bool leftIsNeutral  = true;
  bool rightIsNeutral = true;


  switch (source->type) {
  /* binary operators: where zero (0) is neutral */
  case EXP_ADD_OP:
  case EXP_SUB_OP:
    assert(source->left != NULL);
    assert(source->right != NULL);

    leftIsNeutral = isZeroExpTree(leftSimplified);
    rightIsNeutral = isZeroExpTree(rightSimplified);

    if (leftIsNeutral && rightIsNeutral) {
      delExpTree(rightSimplified);
      return leftSimplified;
    }
    
    if (!leftIsNeutral && !rightIsNeutral)
      return newExpOp(source->type, leftSimplified, rightSimplified);

    if (leftIsNeutral) {
      delExpTree(leftSimplified);
      return (source->type == EXP_SUB_OP) ? newExpOp(EXP_NEG, rightSimplified, NULL) : rightSimplified;
    }

    /* Edge case: 0 - b = -b */
    delExpTree(rightSimplified);
    return leftSimplified;


  /* binary operators: where one (1) is neutral */
  case EXP_MUL_OP:
    assert(source->left != NULL);
    assert(source->right != NULL);
    leftIsNeutral = isOneExpTree(leftSimplified);
    rightIsNeutral = isOneExpTree(rightSimplified);

    if (leftIsNeutral && rightIsNeutral) {
      delExpTree(rightSimplified);
      return leftSimplified;
    }
    
    if (!leftIsNeutral && !rightIsNeutral)
      return newExpOp(source->type, leftSimplified, rightSimplified);

    if (leftIsNeutral) {
      delExpTree(leftSimplified);
      return rightSimplified;
    }

    delExpTree(rightSimplified);
    return leftSimplified;


  case EXP_DIV_OP:
    assert(source->left != NULL);
    assert(source->right != NULL);
    leftIsNeutral = isOneExpTree(leftSimplified);
    rightIsNeutral = isOneExpTree(rightSimplified);

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
    bool rightIsOne = isOneExpTree(rightSimplified);


    /* Edge case: 0^0 is indeterminate. */
    assert(!(leftIsZero && rightIsZero));

    /* a^1 = a */
    if (rightIsOne) {
      delExpTree(rightSimplified);
      return leftSimplified;
    }

    /* a^b is retained */
    return newExpOp(EXP_EXP_OP, leftSimplified, rightSimplified);


  /* remaining: never neutral; these operations should always be retained. */
  case EXP_NEG:
    assert(source->left != NULL);
    assert(source->right == NULL);
    return newExpOp(source->type, leftSimplified, NULL);
  case EXP_FUN:
    assert(source->left != NULL);
    assert(source->right == NULL);
    return newExpTree(source->type, strdup(source->data), leftSimplified, NULL);


  default:
    assert(false);
    return NULL;
  }

}


ExpTree *simplifyAbsorbingTerms(ExpTree *source) {
  assert(source != NULL);

  /* Base case: retain leaves. */
  if (source->left == NULL && source->right == NULL)
    return cpyExpTree(source);

  /* Recursive case: do simplification within subtrees/branches. */
  ExpTree *leftSimplified = source->left ?
    simplifyAbsorbingTerms(source->left) : NULL;
  ExpTree *rightSimplified = source->right ?
    simplifyAbsorbingTerms(source->right) : NULL;

  /* The left and/or right branch not existing (NULL ptr)
    is NOT considered equivalent with being absorbing. */
  bool leftIsAbsorbing  = false;
  bool rightIsAbsorbing = false;


  switch (source->type) {
  /* unary operators: where zero (0) is absorbing */
  case EXP_NEG:
    assert(source->left != NULL);
    assert(source->right == NULL);

    if (isZeroExpTree(leftSimplified)) {
      return leftSimplified;
    }
    return newExpOp(source->type, leftSimplified, NULL);


  /* binary operators: where zero (0) is absorbing */
  case EXP_MUL_OP:
    assert(source->left != NULL);
    assert(source->right != NULL);
    leftIsAbsorbing = isZeroExpTree(leftSimplified);
    rightIsAbsorbing = isZeroExpTree(rightSimplified);

    if (leftIsAbsorbing) {
      delExpTree(rightSimplified);
      return leftSimplified;
    }

    if (rightIsAbsorbing) {
      delExpTree(leftSimplified);
      return rightSimplified;
    }

    return newExpOp(source->type, leftSimplified, rightSimplified);


  case EXP_DIV_OP:
    assert(source->left != NULL);
    assert(source->right != NULL);
    leftIsAbsorbing = isZeroExpTree(leftSimplified);
    rightIsAbsorbing = isZeroExpTree(rightSimplified);

    /* a/0 is indeterminate */
    assert(!rightIsAbsorbing);

    if (leftIsAbsorbing) {
      delExpTree(rightSimplified);
      return leftSimplified;
    }

    return newExpOp(source->type, leftSimplified, rightSimplified);


  case EXP_EXP_OP:
    assert(source->left != NULL);
    assert(source->right != NULL);

    bool leftIsZero = isZeroExpTree(leftSimplified);
    bool rightIsZero = isZeroExpTree(rightSimplified);
    bool leftIsOne = isOneExpTree(leftSimplified);


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

    /* a^b is retained */
    return newExpOp(EXP_EXP_OP, leftSimplified, rightSimplified);


  /* remaining: no absorbing elements exist; subtrees should always be retained */
  case EXP_ADD_OP:
  case EXP_SUB_OP:
    assert(source->left != NULL);
    assert(source->right != NULL);
    return newExpOp(source->type, leftSimplified, rightSimplified);

  case EXP_FUN:
    assert(source->left != NULL);
    assert(source->right == NULL);
    return newExpTree(source->type, strdup(source->data), leftSimplified, NULL);


  default:
    assert(false);
    return NULL;
  }

}


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

ExpTree * newZeroExpTree() {
  return newExpLeaf(EXP_NUM, "0");
}


bool isOneExpTree(ExpTree *source) {
  return source != NULL && source->type == EXP_NUM && atof(source->data) == 1.0;
}


ExpTree *newOneExpTree() {
  return newExpLeaf(EXP_NUM, "1");
}
