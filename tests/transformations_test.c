#include "transformations.h"

void testSimplified(ExpTree *from, ExpTree *simplified, ExpTree *expected) {
  bool compare = isEqual(simplified, expected);

  printf("from:   ");
  if (from == NULL)
    printf("NULL");
  else
    printExpTree(from, stdout);
  printf("\nto:     ");
  if (simplified == NULL)
    printf("NULL");
  else
    printExpTree(simplified, stdout);
  printf("\nexpect: ");
  if (expected == NULL)
    printf("NULL");
  else
    printExpTree(expected, stdout);
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
  ExpTree *c = newExpLeaf(EXP_VAR, "c");
  ExpTree *x = newExpLeaf(EXP_VAR, "x");
  ExpTree *y = newExpLeaf(EXP_VAR, "y");
  ExpTree *z = newExpLeaf(EXP_VAR, "z");

  /*
    Test Operator Simplification
  */
  printf("### Operator Simplification ###\n");
  fflush(stdout);

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
    printf("=== COMPOSITE ===\n");
    fflush(stdout);

    /* Compose expected results */
    ExpTree *zero = newExpLeaf(EXP_NUM, "0");
    ExpTree *one = newExpLeaf(EXP_NUM, "1");
    ExpTree *plusOne = newExpOp(EXP_ADD_OP, cpyExpTree(one), cpyExpTree(one));
    ExpTree *pow =
        newExpOp(EXP_EXP_OP, cpyExpTree(a),
                 newExpOp(EXP_ADD_OP,
                          newExpOp(EXP_SUB_OP, cpyExpTree(b), newOneExpTree()),
                          cpyExpTree(b)));
    ExpTree *pow2 = newExpOp(EXP_EXP_OP, cpyExpTree(a), cpyExpTree(b));

    /* ((1 * 1) + (1 * 1))  =>  (1 + 1) */
    ExpTree *mul = newExpOp(EXP_MUL_OP, newOneExpTree(), newOneExpTree());
    exp = newExpOp(EXP_ADD_OP, cpyExpTree(mul), cpyExpTree(mul));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, plusOne);

    /* ((0 + 0) + (0 + 0))  =>  0 */
    exp = newExpOp(EXP_ADD_OP,
                   newExpOp(EXP_ADD_OP, newZeroExpTree(), newZeroExpTree()),
                   newExpOp(EXP_ADD_OP, newZeroExpTree(), newZeroExpTree()));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, zero);

    /* ((1 * 1) * (1 * 1))  =>  1 */
    exp = newExpOp(EXP_MUL_OP,
                   newExpOp(EXP_MUL_OP, newOneExpTree(), newOneExpTree()),
                   newExpOp(EXP_MUL_OP, newOneExpTree(), newOneExpTree()));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, one);

    /* (((0 + 1) + 0) * ((1 * 1) + (1 + (0 + 0))))  =>  ((1 + 0) * (1 + 1))  =>
     * (1 * (1 + 1))  => (1 + 1) */
    exp = newExpOp(
        EXP_MUL_OP,
        newExpOp(EXP_ADD_OP,
                 newExpOp(EXP_ADD_OP, newZeroExpTree(), newOneExpTree()),
                 newZeroExpTree()),
        newExpOp(EXP_ADD_OP, cpyExpTree(mul),
                 newExpOp(EXP_ADD_OP, newOneExpTree(),
                          newExpOp(EXP_ADD_OP, newZeroExpTree(),
                                   newZeroExpTree()))));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, plusOne);

    /*  (1 * ((b / 1) - (b^0))) */
    ExpTree *exponentLeft = newExpOp(
        EXP_MUL_OP, newOneExpTree(),
        newExpOp(EXP_SUB_OP,
                 newExpOp(EXP_DIV_OP, cpyExpTree(b), newOneExpTree()),
                 newExpOp(EXP_EXP_OP, cpyExpTree(b), newZeroExpTree())));
    /* (-(1 * (0 + 0)) + (b / 1)) */
    ExpTree *exponentRight =
        newExpOp(EXP_ADD_OP,
                 newExpOp(EXP_NEG,
                          newExpOp(EXP_MUL_OP, newOneExpTree(),
                                   newExpOp(EXP_ADD_OP, newZeroExpTree(),
                                            newZeroExpTree())),
                          NULL),
                 newExpOp(EXP_DIV_OP, cpyExpTree(b), newOneExpTree()));
    /* a^(exponentLeft + exponentRight) = (a^((b - 1) + b)) */
    exp = newExpOp(EXP_EXP_OP, cpyExpTree(a),
                   newExpOp(EXP_ADD_OP, exponentLeft, exponentRight));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, pow);

    /* (a^(b^(a^0)))  =>  (a^b) */
    exp = newExpOp(
        EXP_EXP_OP, cpyExpTree(a),
        newExpOp(EXP_EXP_OP, cpyExpTree(b),
                 newExpOp(EXP_EXP_OP, cpyExpTree(a), newZeroExpTree())));
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, pow2);

    /* ---(1 * -(0 + 0))  =>  ----0 */
    exp = newExpOp(
        EXP_NEG,
        newExpOp(
            EXP_NEG,
            newExpOp(EXP_NEG,
                     newExpOp(EXP_MUL_OP, newOneExpTree(),
                              newExpOp(EXP_NEG,
                                       newExpOp(EXP_ADD_OP, newZeroExpTree(),
                                                newZeroExpTree()),
                                       NULL)),
                     NULL),
            NULL),
        NULL);
    simpl = simplifyOperators(exp);
    testSimplified(exp, simpl, zero);
  }

  /*
    Test Sum of Products transformation
  */
  printf("### To Sum of Products ###\n");
  fflush(stdout);

  /* Compose common subtrees. */
  ExpTree *xTy = newExpOp(EXP_MUL_OP, cpyExpTree(x), cpyExpTree(y));
  ExpTree *xTz = newExpOp(EXP_MUL_OP, cpyExpTree(x), cpyExpTree(z));
  ExpTree *yTz = newExpOp(EXP_MUL_OP, cpyExpTree(y), cpyExpTree(z));
  ExpTree *xTx = newExpOp(EXP_MUL_OP, cpyExpTree(x), cpyExpTree(x));
  ExpTree *yTy = newExpOp(EXP_MUL_OP, cpyExpTree(y), cpyExpTree(y));

  /* Test distributing the unary negative operator:
    -(x + y) = -x - y
  */
  {
    printf("=== Unary negative distribution ===\n\n");
    fflush(stdout);

    /* Compose expected results */
    ExpTree *zero = newExpLeaf(EXP_NUM, strdup("0"));
    ExpTree *powab = newExpOp(EXP_EXP_OP, cpyExpTree(a), cpyExpTree(b));
    ExpTree *sina = newExpTree(EXP_FUN, strdup("sin"), cpyExpTree(a), NULL);

    ExpTree *negDiv = newExpOp(EXP_DIV_OP,
                               newExpOp(EXP_SUB_OP, cpyExpTree(y),
                                        newExpOp(EXP_NEG, cpyExpTree(b), NULL)),
                               newExpOp(EXP_NEG, cpyExpTree(z), NULL));
    ExpTree *negAdd =
        newExpOp(EXP_NEG,
                 newExpOp(EXP_ADD_OP, cpyExpTree(sina),
                          newExpOp(EXP_NEG, cpyExpTree(negDiv), NULL)),
                 NULL);
    ExpTree *negBlocked =
        newExpOp(EXP_NEG, newExpOp(EXP_MUL_OP, cpyExpTree(a), negAdd), NULL);
    ExpTree *negSub =
        newExpOp(EXP_SUB_OP, newExpOp(EXP_NEG, cpyExpTree(sina), NULL),
                 newExpOp(EXP_NEG, cpyExpTree(negDiv), NULL));
    ExpTree *negBlockedTo =
        newExpOp(EXP_NEG, newExpOp(EXP_MUL_OP, cpyExpTree(a), negSub), NULL);

    ExpTree *aNegsina =
        newExpOp(EXP_NEG,
                 newExpOp(EXP_MUL_OP, cpyExpTree(a),
                          newExpOp(EXP_NEG, cpyExpTree(sina), NULL)),
                 NULL);
    ExpTree *aNegDiv = newExpOp(EXP_MUL_OP, cpyExpTree(a),
                                newExpOp(EXP_NEG, cpyExpTree(negDiv), NULL));
    ExpTree *negBlockedToDistr = newExpOp(EXP_ADD_OP, aNegsina, aNegDiv);

    ExpTree *apowNegab = newExpOp(
        EXP_EXP_OP, cpyExpTree(a),
        newExpOp(EXP_NEG, newExpOp(EXP_ADD_OP, cpyExpTree(a), cpyExpTree(b)),
                 NULL));
    ExpTree *nega = newExpOp(EXP_NEG, cpyExpTree(a), NULL);
    ExpTree *negaPlusb = newExpOp(
        EXP_NEG, newExpOp(EXP_ADD_OP, cpyExpTree(a), cpyExpTree(b)), NULL);
    ExpTree *add11 = newExpOp(EXP_ADD_OP, cpyExpTree(negaPlusb), cpyExpTree(b));
    ExpTree *sub11 = newExpOp(EXP_SUB_OP, nega, add11);
    ExpTree *add12 = newExpOp(EXP_ADD_OP, sub11, apowNegab);
    ExpTree *mulLeft = newExpOp(EXP_NEG, add12, NULL);
    ExpTree *negMul = newExpOp(
        EXP_NEG, newExpOp(EXP_MUL_OP, mulLeft, cpyExpTree(negaPlusb)), NULL);

    ExpTree *negaSub = newExpOp(
        EXP_SUB_OP, newExpOp(EXP_NEG, cpyExpTree(a), NULL), cpyExpTree(b));
    ExpTree *apowNegaSub =
        newExpOp(EXP_EXP_OP, cpyExpTree(a), cpyExpTree(negaSub));
    ExpTree *add21 = newExpOp(EXP_ADD_OP, cpyExpTree(negaSub), cpyExpTree(b));
    ExpTree *add22 = newExpOp(EXP_ADD_OP, cpyExpTree(a), add21);
    ExpTree *sub21 = newExpOp(EXP_SUB_OP, add22, apowNegaSub);
    ExpTree *negMulTo = newExpOp(
        EXP_NEG, newExpOp(EXP_MUL_OP, sub21, cpyExpTree(negaSub)), NULL);

    /* --x  =>  x */
    exp = newExpOp(EXP_NEG, newExpOp(EXP_NEG, cpyExpTree(x), NULL), NULL);
    simpl = distributeNeg(exp, false);
    testSimplified(exp, simpl, x);
    simpl = toSumOfProducts(exp);
    testSimplified(exp, simpl, x);

    /* -x  =>  -x */
    exp = newExpOp(EXP_NEG, cpyExpTree(x), NULL);
    simpl = distributeNeg(exp, false);
    testSimplified(exp, simpl, exp);
    simpl = toSumOfProducts(exp);
    testSimplified(exp, simpl, exp);

    /* -0  =>  -0 */
    exp = newExpOp(EXP_NEG, cpyExpTree(zero), NULL);
    simpl = distributeNeg(exp, false);
    testSimplified(exp, simpl, exp);
    simpl = toSumOfProducts(exp);
    testSimplified(exp, simpl, exp);

    /* -(a * b) = -(a * b) */
    exp = newExpOp(EXP_NEG, newExpOp(EXP_MUL_OP, cpyExpTree(a), cpyExpTree(b)),
                   NULL);
    simpl = distributeNeg(exp, false);
    testSimplified(exp, simpl, exp);
    simpl = toSumOfProducts(exp);
    testSimplified(exp, simpl, exp);

    /* -(a / b) = -(a / b) */
    exp = newExpOp(EXP_NEG, newExpOp(EXP_DIV_OP, cpyExpTree(a), cpyExpTree(b)),
                   NULL);
    simpl = distributeNeg(exp, false);
    testSimplified(exp, simpl, exp);
    simpl = toSumOfProducts(exp);
    testSimplified(exp, simpl, exp);

    /* -(a^b) != ((-a)^b) in many cases */
    exp = newExpOp(EXP_NEG, cpyExpTree(powab), NULL);
    simpl = distributeNeg(exp, false);
    testSimplified(exp, simpl, exp);
    simpl = toSumOfProducts(exp);
    testSimplified(exp, simpl, exp);

    /* -f(a) = -f(a) */
    exp = newExpOp(EXP_NEG, cpyExpTree(sina), NULL);
    simpl = distributeNeg(exp, false);
    testSimplified(exp, simpl, exp);
    simpl = toSumOfProducts(exp);
    testSimplified(exp, simpl, exp);

    /* ----(a^b) != ((-a)^b) in many cases */
    exp = newExpOp(
        EXP_NEG,
        newExpOp(
            EXP_NEG,
            newExpOp(EXP_NEG, newExpOp(EXP_NEG, cpyExpTree(powab), NULL), NULL),
            NULL),
        NULL);
    simpl = distributeNeg(exp, false);
    testSimplified(exp, simpl, powab);
    simpl = toSumOfProducts(exp);
    testSimplified(exp, simpl, powab);

    /* -(a * -(sin(a) + -((y - -b) / -z)))
    => -(a * (-sin(a) - -((y - -b) / -z))) */
    exp = negBlocked;
    simpl = distributeNeg(exp, false);
    testSimplified(exp, simpl, negBlockedTo);
    /* -(a * -(sin(a) + -((y - -b) / -z)))
    => (-(a * -sin(a)) + (a * -((y - -b) / -z))) */
    simpl = toSumOfProducts(exp);
    testSimplified(exp, simpl, negBlockedToDistr);

    /* -(-((-a - (-(a + b) + b)) + (a^-(a + b))) * -(a + b))
    => -(((a + ((-a - b) + b)) - (a^(-a - b))) * (-a - b)) */
    exp = negMul;
    simpl = distributeNeg(exp, false);
    testSimplified(exp, simpl, negMulTo);

    delExpTree(zero);
    delExpTree(powab);
    delExpTree(sina);
    delExpTree(negBlocked);
    delExpTree(negBlockedTo);
    delExpTree(negBlockedToDistr);
    delExpTree(negMul);
    delExpTree(negMulTo);
    delExpTree(negaSub);
    delExpTree(negaPlusb);
    delExpTree(negDiv);
  }

  /* Test single term distribution:
    x * (y1 + ... + yn)
    = x*y1 + ... + x*yn
  */
  {
    printf("=== Single term distribution ===\n\n");
    fflush(stdout);

    /* Compose expected results */
    ExpTree *zTzPz =
        newExpOp(EXP_MUL_OP, cpyExpTree(z),
                 newExpOp(EXP_ADD_OP, cpyExpTree(z), cpyExpTree(z)));
    ExpTree *xTzTz =
        newExpOp(EXP_MUL_OP, cpyExpTree(x),
                 newExpOp(EXP_MUL_OP, cpyExpTree(z), cpyExpTree(z)));
    ExpTree *rightBranch =
        newExpOp(EXP_ADD_OP, cpyExpTree(xTy), cpyExpTree(xTz));
    ExpTree *rightBranchLong = newExpOp(
        EXP_ADD_OP, cpyExpTree(xTy),
        newExpOp(EXP_SUB_OP, cpyExpTree(xTz),
                 newExpOp(EXP_ADD_OP, cpyExpTree(xTy),
                          newExpOp(EXP_ADD_OP, cpyExpTree(xTz),
                                   newExpOp(EXP_SUB_OP, cpyExpTree(xTy),
                                            cpyExpTree(xTz))))));
    ExpTree *rightNoRecurse =
        newExpOp(EXP_ADD_OP, cpyExpTree(xTy),
                 newExpOp(EXP_MUL_OP, cpyExpTree(x), cpyExpTree(zTzPz)));
    ExpTree *rightYesRecurse =
        newExpOp(EXP_ADD_OP, cpyExpTree(xTy),
                 newExpOp(EXP_ADD_OP, cpyExpTree(xTzTz), cpyExpTree(xTzTz)));
    ExpTree *balanced = newExpOp(
        EXP_ADD_OP,
        newExpOp(EXP_ADD_OP, cpyExpTree(xTx),
                 newExpOp(EXP_SUB_OP, cpyExpTree(xTy), cpyExpTree(xTy))),
        newExpOp(EXP_SUB_OP,
                 newExpOp(EXP_ADD_OP, cpyExpTree(xTy), cpyExpTree(xTy)),
                 cpyExpTree(xTz)));

    /* (x * (y + z))  =>  ((x * y) + (x * z)) */
    exp = newExpOp(EXP_MUL_OP, cpyExpTree(x),
                   newExpOp(EXP_ADD_OP, cpyExpTree(y), cpyExpTree(z)));
    simpl = distributeLeft(exp->left, exp->right);
    testSimplified(exp, simpl, rightBranch);
    simpl = toSumOfProducts(exp);
    testSimplified(exp, simpl, rightBranch);

    /* (x * (y + (z * (z + z))))  =>  ((x * y) + (x * (z * (z + z))))
      ==> The * operator in (z * (z + z)) blocks the singular distribution
          pass from distributing x into (z + z) */
    exp = newExpOp(EXP_MUL_OP, cpyExpTree(x),
                   newExpOp(EXP_ADD_OP, cpyExpTree(y), cpyExpTree(zTzPz)));
    simpl = distributeLeft(exp->left, exp->right);
    testSimplified(exp, simpl, rightNoRecurse);
    /* (x * (y + (z * (z + z))))  =>  ((x * y) + ((x * (z * z)) + (x * (z *
      z))))
      ==> The * operator in (z * (z + z)) does NOT block the recursive
      distribution passes from distributing x into (z + z) */
    simpl = toSumOfProducts(exp);
    testSimplified(exp, simpl, rightYesRecurse);

    /* (x * (y + (z - (y + (z + (y - z))))))  =>  ((x * y) + ((x * z) - ((x * y)
     * + ((x * z) + ((x * y) - (x * z)))))) */
    exp = newExpOp(
        EXP_MUL_OP, cpyExpTree(x),
        newExpOp(EXP_ADD_OP, cpyExpTree(y),
                 newExpOp(EXP_SUB_OP, cpyExpTree(z),
                          newExpOp(EXP_ADD_OP, cpyExpTree(y),
                                   newExpOp(EXP_ADD_OP, cpyExpTree(z),
                                            newExpOp(EXP_SUB_OP, cpyExpTree(y),
                                                     cpyExpTree(z)))))));
    simpl = distributeLeft(exp->left, exp->right);
    testSimplified(exp, simpl, rightBranchLong);
    simpl = toSumOfProducts(exp);
    testSimplified(exp, simpl, rightBranchLong);

    /* ((y + z) * x)  =>  ((x * y) + (x * z)) */
    exp =
        newExpOp(EXP_MUL_OP, newExpOp(EXP_ADD_OP, cpyExpTree(y), cpyExpTree(z)),
                 cpyExpTree(x));
    simpl = toSumOfProducts(exp);
    testSimplified(exp, simpl, rightBranch);

    /* (x * ((x + (y - y)) + ((y + y) - z)))  =>  (((x * x) + ((x * y) - (x *
     * y))) + (((x * y) + (x * y)) - (x * z))) */
    exp = newExpOp(
        EXP_MUL_OP, cpyExpTree(x),
        newExpOp(EXP_ADD_OP,
                 newExpOp(EXP_ADD_OP, cpyExpTree(x),
                          newExpOp(EXP_SUB_OP, cpyExpTree(y), cpyExpTree(y))),
                 newExpOp(EXP_SUB_OP,
                          newExpOp(EXP_ADD_OP, cpyExpTree(y), cpyExpTree(y)),
                          cpyExpTree(z))));
    simpl = distributeLeft(exp->left, exp->right);
    testSimplified(exp, simpl, balanced);
    simpl = toSumOfProducts(exp);
    testSimplified(exp, simpl, balanced);
  }

  /* Test recursive term distribution:
    (x1 + ... + xn) * (y1 + ... + ym)
    = x1*y1 + ... + xn*ym + ... + xn*y1 + ... + xn*ym
  */
  {
    printf("=== Distributive term distribution ===\n\n");
    fflush(stdout);

    /* Compose expected results */
    ExpTree *subRight = newExpOp(
        EXP_ADD_OP, newExpOp(EXP_SUB_OP, cpyExpTree(xTy), cpyExpTree(xTz)),
        newExpOp(EXP_SUB_OP, cpyExpTree(yTy), cpyExpTree(yTz)));
    ExpTree *subLeft = newExpOp(
        EXP_SUB_OP, newExpOp(EXP_ADD_OP, cpyExpTree(xTy), cpyExpTree(xTz)),
        newExpOp(EXP_ADD_OP, cpyExpTree(yTy), cpyExpTree(yTz)));

    ExpTree *xTz = newExpOp(EXP_MUL_OP, cpyExpTree(x), cpyExpTree(z));
    ExpTree *negx = newExpOp(EXP_NEG, cpyExpTree(x), NULL);
    ExpTree *negxTz = newExpOp(EXP_MUL_OP, cpyExpTree(negx), cpyExpTree(z));
    ExpTree *subNeg1 = newExpOp(
        EXP_SUB_OP,
        newExpOp(EXP_NEG, newExpOp(EXP_MUL_OP, cpyExpTree(x), cpyExpTree(y)),
                 NULL),
        cpyExpTree(xTz));
    ExpTree *subNeg2 =
        newExpOp(EXP_SUB_OP,
                 newExpOp(EXP_MUL_OP, cpyExpTree(x),
                          newExpOp(EXP_NEG, cpyExpTree(y), NULL)),
                 cpyExpTree(xTz));
    ExpTree *subNeg3 =
        newExpOp(EXP_SUB_OP,
                 newExpOp(EXP_MUL_OP, cpyExpTree(negx),
                          newExpOp(EXP_NEG, cpyExpTree(y), NULL)),
                 negxTz);

    /* ((x + y) * (y - z))  =>  (((x * y) - (x * z)) + ((y * y) - (y * z))) */
    exp =
        newExpOp(EXP_MUL_OP, newExpOp(EXP_ADD_OP, cpyExpTree(x), cpyExpTree(y)),
                 newExpOp(EXP_SUB_OP, cpyExpTree(y), cpyExpTree(z)));
    simpl = distributeLeftDistributive(exp->left, exp->right);
    testSimplified(exp, simpl, subRight);
    simpl = toSumOfProducts(exp);
    testSimplified(exp, simpl, subRight);

    /* ((x - y) * (y + z))  =>  (((x * y) + (x * z)) - ((y * y) + (y * z))) */
    exp =
        newExpOp(EXP_MUL_OP, newExpOp(EXP_SUB_OP, cpyExpTree(x), cpyExpTree(y)),
                 newExpOp(EXP_ADD_OP, cpyExpTree(y), cpyExpTree(z)));
    simpl = distributeLeftDistributive(exp->left, exp->right);
    testSimplified(exp, simpl, subLeft);
    simpl = toSumOfProducts(exp);
    testSimplified(exp, simpl, subLeft);

    /* -(x * (y + z))  =>  (-(x * y) - (x * z)) */
    exp = newExpOp(EXP_NEG,
                   newExpOp(EXP_MUL_OP, cpyExpTree(x),
                            newExpOp(EXP_ADD_OP, cpyExpTree(y), cpyExpTree(z))),
                   NULL);
    simpl = toSumOfProducts(exp);
    testSimplified(exp, simpl, subNeg1);

    /* (x * -(y + z))  =>  ((x * -y) - (x * z)) */
    exp = newExpOp(EXP_MUL_OP, cpyExpTree(x),
                   newExpOp(EXP_NEG,
                            newExpOp(EXP_ADD_OP, cpyExpTree(y), cpyExpTree(z)),
                            NULL));
    simpl = toSumOfProducts(exp);
    testSimplified(exp, simpl, subNeg2);

    /* (-x * -(y + z))  =>  ((-x * -y) - (-x * z)) */
    exp = newExpOp(EXP_MUL_OP, newExpOp(EXP_NEG, cpyExpTree(x), NULL),
                   newExpOp(EXP_NEG,
                            newExpOp(EXP_ADD_OP, cpyExpTree(y), cpyExpTree(z)),
                            NULL));
    simpl = toSumOfProducts(exp);
    testSimplified(exp, simpl, subNeg3);

    delExpTree(subRight);
    delExpTree(subLeft);
    delExpTree(xTz);
    delExpTree(negx);
    delExpTree(subNeg1);
    delExpTree(subNeg2);
    delExpTree(subNeg3);
  }

  {
    printf("=== Complex distribution ===\n\n");
    fflush(stdout);

    ExpTree *aTbTx =
        newExpOp(EXP_MUL_OP, cpyExpTree(a),
                 newExpOp(EXP_MUL_OP, cpyExpTree(b), cpyExpTree(x)));
    ExpTree *aTbTy =
        newExpOp(EXP_MUL_OP, cpyExpTree(a),
                 newExpOp(EXP_MUL_OP, cpyExpTree(b), cpyExpTree(y)));
    ExpTree *aTbTz =
        newExpOp(EXP_MUL_OP, cpyExpTree(a),
                 newExpOp(EXP_MUL_OP, cpyExpTree(b), cpyExpTree(z)));

    /* Compose expected results */
    ExpTree *twiceSingleTerms = newExpOp(
        EXP_ADD_OP,
        newExpOp(EXP_ADD_OP, cpyExpTree(aTbTx),
                 newExpOp(EXP_SUB_OP, cpyExpTree(aTbTy), cpyExpTree(aTbTy))),
        newExpOp(EXP_SUB_OP,
                 newExpOp(EXP_ADD_OP, cpyExpTree(aTbTy), cpyExpTree(aTbTy)),
                 cpyExpTree(aTbTz)));

    /* Watch out when deleting this mess of subtrees, some are copied and some
     * are not! */
    ExpTree *abx =
        newExpOp(EXP_MUL_OP, newExpOp(EXP_MUL_OP, cpyExpTree(a), cpyExpTree(b)),
                 cpyExpTree(x));
    ExpTree *aby =
        newExpOp(EXP_MUL_OP, newExpOp(EXP_MUL_OP, cpyExpTree(a), cpyExpTree(b)),
                 cpyExpTree(y));
    ExpTree *abya =
        newExpOp(EXP_MUL_OP, newExpOp(EXP_MUL_OP, cpyExpTree(a), cpyExpTree(b)),
                 newExpOp(EXP_MUL_OP, cpyExpTree(y), cpyExpTree(a)));
    ExpTree *abyb =
        newExpOp(EXP_MUL_OP, newExpOp(EXP_MUL_OP, cpyExpTree(a), cpyExpTree(b)),
                 newExpOp(EXP_MUL_OP, cpyExpTree(y), cpyExpTree(b)));
    ExpTree *abyc =
        newExpOp(EXP_MUL_OP, newExpOp(EXP_MUL_OP, cpyExpTree(a), cpyExpTree(b)),
                 newExpOp(EXP_MUL_OP, cpyExpTree(y), cpyExpTree(c)));
    ExpTree *abza =
        newExpOp(EXP_MUL_OP, newExpOp(EXP_MUL_OP, cpyExpTree(a), cpyExpTree(b)),
                 newExpOp(EXP_MUL_OP, cpyExpTree(z), cpyExpTree(a)));
    ExpTree *abzb =
        newExpOp(EXP_MUL_OP, newExpOp(EXP_MUL_OP, cpyExpTree(a), cpyExpTree(b)),
                 newExpOp(EXP_MUL_OP, cpyExpTree(z), cpyExpTree(b)));
    ExpTree *abzc =
        newExpOp(EXP_MUL_OP, newExpOp(EXP_MUL_OP, cpyExpTree(a), cpyExpTree(b)),
                 newExpOp(EXP_MUL_OP, cpyExpTree(z), cpyExpTree(c)));
    ExpTree *innerLSum1 =
        newExpOp(EXP_ADD_OP, newExpOp(EXP_ADD_OP, abya, cpyExpTree(abyb)),
                 newExpOp(EXP_ADD_OP, cpyExpTree(abyb), abyc));
    ExpTree *innerLSum2 =
        newExpOp(EXP_ADD_OP, newExpOp(EXP_ADD_OP, abza, cpyExpTree(abzb)),
                 newExpOp(EXP_ADD_OP, cpyExpTree(abzb), abzc));
    ExpTree *innerLSum = newExpOp(EXP_SUB_OP, innerLSum1, innerLSum2);
    ExpTree *leftLongTerm =
        newExpOp(EXP_SUB_OP, abx, newExpOp(EXP_ADD_OP, aby, innerLSum));
    ExpTree *acx =
        newExpOp(EXP_MUL_OP, newExpOp(EXP_MUL_OP, cpyExpTree(a), cpyExpTree(c)),
                 cpyExpTree(x));
    ExpTree *acy =
        newExpOp(EXP_MUL_OP, newExpOp(EXP_MUL_OP, cpyExpTree(a), cpyExpTree(c)),
                 cpyExpTree(y));
    ExpTree *acya =
        newExpOp(EXP_MUL_OP, newExpOp(EXP_MUL_OP, cpyExpTree(a), cpyExpTree(c)),
                 newExpOp(EXP_MUL_OP, cpyExpTree(y), cpyExpTree(a)));
    ExpTree *acyb =
        newExpOp(EXP_MUL_OP, newExpOp(EXP_MUL_OP, cpyExpTree(a), cpyExpTree(c)),
                 newExpOp(EXP_MUL_OP, cpyExpTree(y), cpyExpTree(b)));
    ExpTree *acyc =
        newExpOp(EXP_MUL_OP, newExpOp(EXP_MUL_OP, cpyExpTree(a), cpyExpTree(c)),
                 newExpOp(EXP_MUL_OP, cpyExpTree(y), cpyExpTree(c)));
    ExpTree *acza =
        newExpOp(EXP_MUL_OP, newExpOp(EXP_MUL_OP, cpyExpTree(a), cpyExpTree(c)),
                 newExpOp(EXP_MUL_OP, cpyExpTree(z), cpyExpTree(a)));
    ExpTree *aczb =
        newExpOp(EXP_MUL_OP, newExpOp(EXP_MUL_OP, cpyExpTree(a), cpyExpTree(c)),
                 newExpOp(EXP_MUL_OP, cpyExpTree(z), cpyExpTree(b)));
    ExpTree *aczc =
        newExpOp(EXP_MUL_OP, newExpOp(EXP_MUL_OP, cpyExpTree(a), cpyExpTree(c)),
                 newExpOp(EXP_MUL_OP, cpyExpTree(z), cpyExpTree(c)));
    ExpTree *innerRSum1 =
        newExpOp(EXP_ADD_OP, newExpOp(EXP_ADD_OP, acya, cpyExpTree(acyb)),
                 newExpOp(EXP_ADD_OP, cpyExpTree(acyb), acyc));
    ExpTree *innerRSum2 =
        newExpOp(EXP_ADD_OP, newExpOp(EXP_ADD_OP, acza, cpyExpTree(aczb)),
                 newExpOp(EXP_ADD_OP, cpyExpTree(aczb), aczc));
    ExpTree *innerRSum = newExpOp(EXP_SUB_OP, innerRSum1, innerRSum2);
    ExpTree *rightLongTerm =
        newExpOp(EXP_SUB_OP, acx, newExpOp(EXP_ADD_OP, acy, innerRSum));
    ExpTree *combinedLongTerm =
        newExpOp(EXP_ADD_OP, leftLongTerm, rightLongTerm);

    ExpTree *fSubexp =
        newExpOp(EXP_ADD_OP, newExpOp(EXP_MUL_OP, cpyExpTree(a), cpyExpTree(y)),
                 newExpOp(EXP_MUL_OP, cpyExpTree(a), cpyExpTree(z)));
    ExpTree *fDistributed = newExpOp(EXP_MUL_OP, cpyExpTree(x),
                                     newExpTree(EXP_FUN, "sin", fSubexp, NULL));

    /* (a * (EXP * b))
     = (a * (((x + (y - y)) + ((y + y) - z)) * b))
    => (((a * (b * x)) + ((a * (b * y)) - (a * (b * y)))) + (((a * (b * y)) + (a
    * (b * y))) - (a * (b * z)))) */
    exp = newExpOp(
        EXP_MUL_OP, cpyExpTree(a),
        newExpOp(
            EXP_MUL_OP,
            newExpOp(
                EXP_ADD_OP,
                newExpOp(EXP_ADD_OP, cpyExpTree(x),
                         newExpOp(EXP_SUB_OP, cpyExpTree(y), cpyExpTree(y))),
                newExpOp(EXP_SUB_OP,
                         newExpOp(EXP_ADD_OP, cpyExpTree(y), cpyExpTree(y)),
                         cpyExpTree(z))),
            cpyExpTree(b)));
    simpl = toSumOfProducts(exp);
    testSimplified(exp, simpl, twiceSingleTerms);

    /* OUTER = (a * (b + c))
       INNER = ((a + b) + (b + c))
       OUTER * (x - (y + ((y - z) * INNER)))
    */
    ExpTree *outerTerm =
        newExpOp(EXP_MUL_OP, cpyExpTree(a),
                 newExpOp(EXP_ADD_OP, cpyExpTree(b), cpyExpTree(c)));
    ExpTree *innerTerm =
        newExpOp(EXP_ADD_OP, newExpOp(EXP_ADD_OP, cpyExpTree(a), cpyExpTree(b)),
                 newExpOp(EXP_ADD_OP, cpyExpTree(b), cpyExpTree(c)));
    ExpTree *inner1 = newExpOp(EXP_SUB_OP, cpyExpTree(y), cpyExpTree(z));
    ExpTree *inner2 = newExpOp(EXP_MUL_OP, inner1, innerTerm);
    ExpTree *inner3 = newExpOp(EXP_ADD_OP, cpyExpTree(y), inner2);
    ExpTree *rightTerm = newExpOp(EXP_SUB_OP, cpyExpTree(x), inner3);
    exp = newExpOp(EXP_MUL_OP, outerTerm, rightTerm);
    simpl = toSumOfProducts(exp);
    testSimplified(exp, simpl, combinedLongTerm);

    /* (x * sin((a * (y + z))))  =>  (x * sin(((a * y) + (a * z)))) */
    ExpTree *distSubexp =
        newExpOp(EXP_MUL_OP, cpyExpTree(a),
                 newExpOp(EXP_ADD_OP, cpyExpTree(y), cpyExpTree(z)));
    ExpTree *f = newExpTree(EXP_FUN, "sin", distSubexp, NULL);
    exp = newExpOp(EXP_MUL_OP, cpyExpTree(x), f);
    simpl = toSumOfProducts(exp);
    testSimplified(exp, simpl, fDistributed);
  }

  /*
    Test variable substitution transformation
  */
  printf("### Variable substitution ###\n");
  fflush(stdout);

  {
    /* Compose expected results */
    ExpTree *addyz = newExpOp(EXP_ADD_OP, cpyExpTree(y), cpyExpTree(z));

    ExpTree *twoTz =
        newExpOp(EXP_MUL_OP, newExpLeaf(EXP_NUM, strdup("2")), cpyExpTree(z));
    ExpTree *pow = newExpOp(EXP_EXP_OP, cpyExpTree(twoTz), cpyExpTree(twoTz));
    ExpTree *sub = newExpOp(EXP_SUB_OP, pow, cpyExpTree(twoTz));
    ExpTree *mul = newExpOp(EXP_MUL_OP, cpyExpTree(a), sub);
    ExpTree *negsin =
        newExpOp(EXP_NEG,
                 newExpTree(EXP_FUN, strdup("sin"),
                            newExpOp(EXP_NEG, cpyExpTree(twoTz), NULL), NULL),
                 NULL);
    ExpTree *negDiv =
        newExpOp(EXP_NEG, newExpOp(EXP_DIV_OP, mul, negsin), NULL);

    /* The following notation is used for expressing the substitution of
      variable X by expression TARGET in expression SOURCE: (SOURCE)[x :=
      TARGET] */

    /* (x)[x := y]  =>  y */
    exp = cpyExpTree(x);
    simpl = substitute(exp, x->data, y);
    testSimplified(exp, simpl, y);
    delExpTree(simpl);
    delExpTree(exp);

    /* (3)[x := y]  =>  3 */
    exp = newExpLeaf(EXP_NUM, strdup("3"));
    simpl = substitute(exp, x->data, y);
    testSimplified(exp, simpl, exp);
    delExpTree(simpl);
    delExpTree(exp);

    /* ((x + z))[x := y]  =>  (y + z) */
    exp = newExpOp(EXP_ADD_OP, cpyExpTree(x), cpyExpTree(z));
    simpl = substitute(exp, x->data, y);
    testSimplified(exp, simpl, addyz);
    delExpTree(simpl);
    delExpTree(exp);

    /* (-((a * ((x^x) - x)) / -(sin(-x))))[x := (2 * z)]
    =>  -((a * (((2 * z)^(2 * z)) - (2 * z))) / -(sin(-(2 * z)))) */
    ExpTree *powX = newExpOp(EXP_EXP_OP, cpyExpTree(x), cpyExpTree(x));
    ExpTree *subX = newExpOp(EXP_SUB_OP, powX, cpyExpTree(x));
    ExpTree *mulX = newExpOp(EXP_MUL_OP, cpyExpTree(a), subX);
    ExpTree *negsinX =
        newExpOp(EXP_NEG,
                 newExpTree(EXP_FUN, strdup("sin"),
                            newExpOp(EXP_NEG, cpyExpTree(x), NULL), NULL),
                 NULL);
    exp = newExpOp(EXP_NEG, newExpOp(EXP_DIV_OP, mulX, negsinX), NULL);
    simpl = substitute(exp, x->data, twoTz);
    testSimplified(exp, simpl, negDiv);
    delExpTree(simpl);
    delExpTree(exp);

    delExpTree(addyz);
    delExpTree(twoTz);
    delExpTree(negDiv);
  }
}
