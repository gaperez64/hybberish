#include "transformations.h"



void testSimplified(ExpTree *from, ExpTree *simplified, ExpTree *expected) {
  bool compare = isEqual(simplified, expected);

  printf("from:   ");
  if (from == NULL) printf("NULL");
  else printExpTree(from, stdout);
  printf("\nto:     ");
  if (simplified == NULL) printf("NULL");
  else printExpTree(simplified, stdout);
  printf("\nexpect: ");
  if (expected == NULL) printf("NULL");
  else printExpTree(expected, stdout);
  printf("\nequal:  %i", compare);
  printf("\n\n");
  fflush(stdout);

  assert(compare);
}


int main(int argc, char *argv[]) {
  /* to avoid silly warnings about unused parameters */
  (void)argc;
  (void)argv;
  

  ExpTree *exp;
  ExpTree *simpl;

  ExpTree *a = newExpLeaf(EXP_VAR, "a");
  ExpTree *b = newExpLeaf(EXP_VAR, "b");


  /* Test binary ADD (+) simplification */
  {
    /* Build expected results */
    ExpTree *zero = newExpLeaf(EXP_NUM, "0");
    ExpTree *one = newExpLeaf(EXP_NUM, "1");
    ExpTree *plus = newExpOp(EXP_ADD_OP, cpyExpTree(a), cpyExpTree(b));
    ExpTree *plusOne = newExpOp(EXP_ADD_OP, cpyExpTree(one), cpyExpTree(one));


    printf("=== ADD (+) ===\n");
    fflush(stdout);

    /* (a + b)  =>  (a + b) */
    exp = newExpOp(EXP_ADD_OP, cpyExpTree(a), cpyExpTree(b));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, plus);

    /* (0 + b)  =>  b */
    exp = newExpOp(EXP_ADD_OP, newZeroExpTree(), cpyExpTree(b));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, b);

    /* (a + 0)  =>  a */
    exp = newExpOp(EXP_ADD_OP, cpyExpTree(a), newZeroExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, a);

    /* (0 + 0)  =>  0 */
    exp = newExpOp(EXP_ADD_OP, newZeroExpTree(), newZeroExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, zero);

    /* (1 + 0)  =>  1 */
    exp = newExpOp(EXP_ADD_OP, newOneExpTree(), newZeroExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, one);

    /* (0 + 1)  =>  1 */
    exp = newExpOp(EXP_ADD_OP, newZeroExpTree(), newOneExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, one);

    /* (1 + 1)  =>  (1 + 1) */
    exp = newExpOp(EXP_ADD_OP, newOneExpTree(), newOneExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, plusOne);
  }

  /* Test binary SUB (-) simplification */
  {
    /* Compose expected results */
    ExpTree *zero = newExpLeaf(EXP_NUM, "0");
    ExpTree *one = newExpLeaf(EXP_NUM, "1");
    ExpTree *sub = newExpOp(EXP_SUB_OP, cpyExpTree(a), cpyExpTree(b));
    ExpTree *subOne = newExpOp(EXP_SUB_OP, cpyExpTree(one), cpyExpTree(one));
    ExpTree *neg = newExpOp(EXP_NEG, cpyExpTree(b), NULL);
    ExpTree *negOne = newExpOp(EXP_NEG, cpyExpTree(one), NULL);


    printf("=== SUB (-) ===\n");
    fflush(stdout);
    /* (a - b)  =>  (a - b) */
    exp = newExpOp(EXP_SUB_OP, cpyExpTree(a), cpyExpTree(b));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, sub);

    /* (0 - b)  =>  -b */
    exp = newExpOp(EXP_SUB_OP, newZeroExpTree(), cpyExpTree(b));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, neg);

    /* (a - 0)  =>  a */
    exp = newExpOp(EXP_SUB_OP, cpyExpTree(a), newZeroExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, a);

    /* (0 - 0)  =>  0 */
    exp = newExpOp(EXP_SUB_OP, newZeroExpTree(), newZeroExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, zero);

    /* (1 - 0)  =>  1 */
    exp = newExpOp(EXP_SUB_OP, newOneExpTree(), newZeroExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, one);

    /* (0 - 1)  =>  -1 */
    exp = newExpOp(EXP_SUB_OP, newZeroExpTree(), newOneExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, negOne);

    /* (1 - 1)  =>  (1 - 1) */
    exp = newExpOp(EXP_SUB_OP, newOneExpTree(), newOneExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, subOne);
  }

  /* Test binary MUL (*) simplification */
  {
    /* Compose expected results */
    ExpTree *zero = newExpLeaf(EXP_NUM, "0");
    ExpTree *one = newExpLeaf(EXP_NUM, "1");
    ExpTree *mul = newExpOp(EXP_MUL_OP, cpyExpTree(a), cpyExpTree(b));


    printf("=== MUL (*) ===\n");
    fflush(stdout);

    /* (a * b)  =>  (a * b) */
    exp = newExpOp(EXP_MUL_OP, cpyExpTree(a), cpyExpTree(b));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, mul);

    /* (0 * b)  =>  0 */
    exp = newExpOp(EXP_MUL_OP, newZeroExpTree(), cpyExpTree(b));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, zero);

    /* (a * 0)  =>  0 */
    exp = newExpOp(EXP_MUL_OP, cpyExpTree(a), newZeroExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, zero);

    /* (1 * b)  =>  b */
    exp = newExpOp(EXP_MUL_OP, newOneExpTree(), cpyExpTree(b));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, b);

    /* (a * 1)  =>  a */
    exp = newExpOp(EXP_MUL_OP, cpyExpTree(a), newOneExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, a);

    /* (1 * 1)  =>  1 */
    exp = newExpOp(EXP_MUL_OP, newOneExpTree(), newOneExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, one);

    /* (0 * 0)  =>  0 */
    exp = newExpOp(EXP_MUL_OP, newZeroExpTree(), newZeroExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, zero);

    /* (1 * 0)  =>  0 */
    exp = newExpOp(EXP_MUL_OP, newOneExpTree(), newZeroExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, zero);

    /* (0 * 1)  =>  0 */
    exp = newExpOp(EXP_MUL_OP, newZeroExpTree(), newOneExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, zero);
  }

  /* Test binary DIV (/) simplification */
  {
    /* Compose expected results */
    ExpTree *zero = newExpLeaf(EXP_NUM, "0");
    ExpTree *one = newExpLeaf(EXP_NUM, "1");
    ExpTree *div = newExpOp(EXP_DIV_OP, cpyExpTree(a), cpyExpTree(b));
    ExpTree *divOneL = newExpOp(EXP_DIV_OP, cpyExpTree(one), cpyExpTree(b));


    printf("=== DIV (/) ===\n");
    fflush(stdout);

    /* (a / b)  =>  (a / b) */
    exp = newExpOp(EXP_DIV_OP, cpyExpTree(a), cpyExpTree(b));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, div);

    /* (0 / b)  =>  0 */
    exp = newExpOp(EXP_DIV_OP, newZeroExpTree(), cpyExpTree(b));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, zero);

    /* (a / 0)  =>  indeterminate! */
    // exp = newExpOp(EXP_DIV_OP, cpyExpTree(a), newZeroExpTree());
    // simpl = simplifyOperators(exp);
    // testSimplified(exp, simpl, divZeroR);

    /* (1 / b)  =>  1 / b */
    exp = newExpOp(EXP_DIV_OP, newOneExpTree(), cpyExpTree(b));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, divOneL);

    /* (a / 1)  =>  a */
    exp = newExpOp(EXP_DIV_OP, cpyExpTree(a), newOneExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, a);

    /* (1 / 1)  =>  1 */
    exp = newExpOp(EXP_DIV_OP, newOneExpTree(), newOneExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, one);

    /* (0 / 0)  =>  indeterminate! */
    // exp = newExpOp(EXP_DIV_OP, newZeroExpTree(), newZeroExpTree());
    // simpl = simplifyOperators(exp);
    // testSimplified(exp, simpl, divZeroLR);

    /* (1 / 0)  =>  indeterminate! */
    // exp = newExpOp(EXP_DIV_OP, newOneExpTree(), newZeroExpTree());
    // simpl = simplifyOperators(exp);
    // testSimplified(exp, simpl, divOneZero);

    /* (0 / 1)  =>  0 */
    exp = newExpOp(EXP_DIV_OP, newZeroExpTree(), newOneExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, zero);
  }

  /* Test binary EXP (^, exponent or power) simplification */
  {
    /* Compose expected results */
    ExpTree *zero = newExpLeaf(EXP_NUM, "0");
    ExpTree *one = newExpLeaf(EXP_NUM, "1");
    ExpTree *pow = newExpOp(EXP_EXP_OP, cpyExpTree(a), cpyExpTree(b));


    printf("=== EXP (^) ===\n");
    fflush(stdout);

    /* (a ^ b)  =>  (a ^ b) */
    exp = newExpOp(EXP_EXP_OP, cpyExpTree(a), cpyExpTree(b));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, pow);

    /* (0 ^ b)  => 0 */
    exp = newExpOp(EXP_EXP_OP, newZeroExpTree(), cpyExpTree(b));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, zero);

    /* (0 ^ 0)  =>  indeterminate! */
    // exp = newExpOp(EXP_EXP_OP, newZeroExpTree(), newZeroExpTree());
    // simpl = simplifyAbsorbingTerms(exp);
    // testSimplified(exp, simpl);

    /* (a ^ 0)  =>  1 */
    exp = newExpOp(EXP_EXP_OP, cpyExpTree(a), newZeroExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, one);

    /* (1 ^ b)  =>  1 */
    exp = newExpOp(EXP_EXP_OP, newOneExpTree(), cpyExpTree(b));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, one);

    /* (a ^ 1)  =>  a */
    exp = newExpOp(EXP_EXP_OP, cpyExpTree(a), newOneExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, a);

    /* (1 ^ 1)  =>  1 */
    exp = newExpOp(EXP_EXP_OP, newOneExpTree(), newOneExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, one);

    /* (1 ^ 0)  =>  1 */
    exp = newExpOp(EXP_EXP_OP, newOneExpTree(), newZeroExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, one);

    /* (0 ^ 1)  =>  0 */
    exp = newExpOp(EXP_EXP_OP, newZeroExpTree(), newOneExpTree());
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, zero);
  }

  /* Test unary NEG (- or negation) simplification */
  {
    /* Compose expected results */
    ExpTree *zero = newExpLeaf(EXP_NUM, "0");
    ExpTree *nega = newExpOp(EXP_NEG, cpyExpTree(a), NULL);
    ExpTree *negOne = newExpOp(EXP_NEG, newExpLeaf(EXP_NUM, "1"), NULL);


    printf("=== NEG (-) ===\n");
    fflush(stdout);

    /* -a  =>  -a */
    exp = newExpOp(EXP_NEG, cpyExpTree(a), NULL);
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, nega);

    /* -1  =>  -1 */
    exp = newExpOp(EXP_NEG, newOneExpTree(), NULL);
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, negOne);

    /* -0  =>  0 */
    exp = newExpOp(EXP_NEG, newZeroExpTree(), NULL);
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, zero);
  }

  /* Test composite simplifications */
  {
    /* Compose expected results */
    ExpTree *zero = newExpLeaf(EXP_NUM, "0");
    ExpTree *one = newExpLeaf(EXP_NUM, "1");
    ExpTree *plusOne = newExpOp(EXP_ADD_OP, cpyExpTree(one), cpyExpTree(one));
    ExpTree *pow = newExpOp(EXP_EXP_OP, cpyExpTree(a), newExpOp(EXP_ADD_OP, newExpOp(EXP_SUB_OP, cpyExpTree(b), newOneExpTree()), cpyExpTree(b)));
    ExpTree *pow2 = newExpOp(EXP_EXP_OP, cpyExpTree(a), cpyExpTree(b));


    printf("=== COMPOSITE ===\n");
    fflush(stdout);
    /* ((1 * 1) + (1 * 1))  =>  (1 + 1) */
    ExpTree *mul = newExpOp(EXP_MUL_OP, newOneExpTree(), newOneExpTree());
    exp = newExpOp(EXP_ADD_OP, cpyExpTree(mul), cpyExpTree(mul));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, plusOne);

    /* ((0 + 0) + (0 + 0))  =>  0 */
    exp = newExpOp(EXP_ADD_OP, newExpOp(EXP_ADD_OP, newZeroExpTree(), newZeroExpTree()), newExpOp(EXP_ADD_OP, newZeroExpTree(), newZeroExpTree()));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, zero);

    /* ((1 * 1) * (1 * 1))  =>  1 */
    exp = newExpOp(EXP_MUL_OP, newExpOp(EXP_MUL_OP, newOneExpTree(), newOneExpTree()), newExpOp(EXP_MUL_OP, newOneExpTree(), newOneExpTree()));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, one);

    /* (((0 + 1) + 0) * ((1 * 1) + (1 + (0 + 0))))  =>  ((1 + 0) * (1 + 1))  =>  (1 * (1 + 1))  => (1 + 1) */
    exp = newExpOp(EXP_MUL_OP, newExpOp(EXP_ADD_OP, newExpOp(EXP_ADD_OP, newZeroExpTree(), newOneExpTree()), newZeroExpTree()), newExpOp(EXP_ADD_OP, cpyExpTree(mul), newExpOp(EXP_ADD_OP, newOneExpTree(), newExpOp(EXP_ADD_OP, newZeroExpTree(), newZeroExpTree()))));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, plusOne);

    /*  (1 * ((b / 1) - (b^0))) */
    ExpTree *exponentLeft = newExpOp(EXP_MUL_OP, newOneExpTree(), newExpOp(EXP_SUB_OP, newExpOp(EXP_DIV_OP, cpyExpTree(b), newOneExpTree()), newExpOp(EXP_EXP_OP, cpyExpTree(b), newZeroExpTree())));
    /* (-(1 * (0 + 0)) + (b / 1)) */
    ExpTree *exponentRight = newExpOp(EXP_ADD_OP, newExpOp(EXP_NEG, newExpOp(EXP_MUL_OP, newOneExpTree(), newExpOp(EXP_ADD_OP, newZeroExpTree(), newZeroExpTree())), NULL), newExpOp(EXP_DIV_OP, cpyExpTree(b), newOneExpTree()));
    /* a^(exponentLeft + exponentRight) = (a^((b - 1) + b)) */
    exp = newExpOp(EXP_EXP_OP, cpyExpTree(a), newExpOp(EXP_ADD_OP, exponentLeft, exponentRight));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, pow);

    /* (a^(b^(a^0)))  =>  (a^b) */
    exp = newExpOp(EXP_EXP_OP, cpyExpTree(a), newExpOp(EXP_EXP_OP, cpyExpTree(b), newExpOp(EXP_EXP_OP, cpyExpTree(a), newZeroExpTree())));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, pow2);

    /* ---(1 * -(0 + 0))  =>  ----0 */
    exp = newExpOp(EXP_NEG, newExpOp(EXP_NEG, newExpOp(EXP_NEG, newExpOp(EXP_MUL_OP, newOneExpTree(), newExpOp(EXP_NEG, newExpOp(EXP_ADD_OP, newZeroExpTree(), newZeroExpTree()), NULL)), NULL), NULL), NULL);
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, zero);
  }
}
