#include "funexp.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
  /* to avoid silly warnings about unused parameters */
  (void)argc;
  (void)argv;

  /* Build leaves */
  ExpTree *x = newExpLeaf(EXP_VAR, strdup("x"));
  ExpTree *y = newExpLeaf(EXP_VAR, strdup("y"));
  ExpTree *z = newExpLeaf(EXP_NUM, strdup("z"));
  ExpTree *n1 = newExpLeaf(EXP_VAR, strdup("1"));
  ExpTree *n2 = newExpLeaf(EXP_NUM, strdup("2"));

  ExpTree *b1 = newExpLeaf(EXP_VAR, strdup("b"));
  ExpTree *b2 = newExpLeaf(EXP_VAR, strdup("b"));
  ExpTree *n4 = newExpLeaf(EXP_NUM, strdup("4"));
  ExpTree *a1 = newExpLeaf(EXP_VAR, strdup("a"));
  ExpTree *c1 = newExpLeaf(EXP_VAR, strdup("c"));
  ExpTree *m2 = newExpLeaf(EXP_NUM, strdup("2"));
  ExpTree *a2 = newExpLeaf(EXP_VAR, strdup("a"));

  /* (x + (y + (z + 1))) */
  ExpTree *sumLeft =
      newExpOp(EXP_ADD_OP, cpyExpTree(x),
               newExpOp(EXP_ADD_OP, cpyExpTree(y),
                        newExpOp(EXP_ADD_OP, cpyExpTree(z), cpyExpTree(n1))));
  /* ((x + y) + (z + 1)) */
  ExpTree *sumBalanced =
      newExpOp(EXP_ADD_OP, newExpOp(EXP_ADD_OP, cpyExpTree(x), cpyExpTree(y)),
               newExpOp(EXP_ADD_OP, cpyExpTree(z), cpyExpTree(n1)));

  /* (-b + sqrt((b^2) - ((4a)c))) / (2a) */
  ExpTree *exp = newExpOp(EXP_EXP_OP, b2, n2);
  ExpTree *foura = newExpOp(EXP_MUL_OP, n4, a1);
  ExpTree *fourac = newExpOp(EXP_MUL_OP, foura, c1);
  ExpTree *min = newExpOp(EXP_SUB_OP, exp, fourac);
  ExpTree *sqrt = newExpTree(EXP_FUN, strdup("sqrt"), min, NULL);
  ExpTree *neg = newExpOp(EXP_NEG, b1, NULL);
  ExpTree *sum = newExpOp(EXP_ADD_OP, neg, sqrt);
  ExpTree *twoa = newExpOp(EXP_MUL_OP, m2, a2);
  ExpTree *tree = newExpOp(EXP_DIV_OP, sum, twoa);

  /* Test exact tree equality. */
  {
    /* Equality */
    assert(isEqual(NULL, NULL));
    assert(isEqual(x, x));
    assert(isEqual(n2, n2));
    assert(isEqual(sumLeft, sumLeft));
    assert(isEqual(sumBalanced, sumBalanced));
    assert(isEqual(tree, tree));

    /* Inequality */
    assert(!isEqual(x, n2));
    assert(!isEqual(x, y));
    // Semantically equivalent but structurally different trees are not exactly
    // equal!
    assert(!isEqual(sumLeft, sumBalanced));
    assert(!isEqual(sumLeft, tree));
    assert(!isEqual(sumBalanced, tree));
  }
}
