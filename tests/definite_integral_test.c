#include "funexp.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void test_integral(ExpTree *expr, char *var, ExpTree *expected,
                   ExpTree *lowerBound, ExpTree *upperBound) {
  /* evaluate integral */
  ExpTree *integral_expr = definiteIntegral(expr, var, lowerBound, upperBound);
  assert(integral_expr != NULL);
  bool compare = isEqual(expected, integral_expr);

  /* printing and comparison */
  printf("From:   ");
  printExpTree(expr, stdout);
  printf("\n");
  printf("To:     ");
  printExpTree(integral_expr, stdout);
  printf("\n");
  printf("Expect: ");
  printExpTree(expected, stdout);
  printf("\n");
  printf("Equal:  %i\n\n", compare);
  fflush(stdout); // Flush stdout
  assert(compare);
}

int main(int argc, char *argv[]) {
  /* to avoid silly warnings about unused parameters */
  (void)argc;
  (void)argv;

  /* Construct common leaves */
  ExpTree *x = newExpLeaf(EXP_VAR, strdup("x"));
  ExpTree *y = newExpLeaf(EXP_VAR, strdup("y"));
  ExpTree *a = newExpLeaf(EXP_VAR, strdup("a"));
  ExpTree *b = newExpLeaf(EXP_VAR, strdup("b"));
  ExpTree *half = newExpLeaf(EXP_NUM, strdup("0.5"));
  ExpTree *one = newExpLeaf(EXP_NUM, strdup("1"));
  ExpTree *two = newExpLeaf(EXP_NUM, strdup("2"));
  ExpTree *three = newExpLeaf(EXP_NUM, strdup("3"));
  ExpTree *four = newExpLeaf(EXP_NUM, strdup("4"));
  ExpTree *five = newExpLeaf(EXP_NUM, strdup("5"));

  /* integral of a constant */
  {
    ExpTree *fiveTa = newExpOp(EXP_MUL_OP, cpyExpTree(five), cpyExpTree(a));
    ExpTree *fiveTb = newExpOp(EXP_MUL_OP, cpyExpTree(five), cpyExpTree(b));
    ExpTree *res = newExpOp(EXP_SUB_OP, fiveTb, fiveTa);
    test_integral(five, "x", res, a, b);

    delExpTree(res);
  }

  /* integral of x w.r.t itself */
  {
    ExpTree *powa = newExpOp(EXP_EXP_OP, cpyExpTree(a), cpyExpTree(two));
    ExpTree *powb = newExpOp(EXP_EXP_OP, cpyExpTree(b), cpyExpTree(two));
    ExpTree *mula = newExpOp(EXP_MUL_OP, cpyExpTree(half), powa);
    ExpTree *mulb = newExpOp(EXP_MUL_OP, cpyExpTree(half), powb);
    ExpTree *res = newExpOp(EXP_SUB_OP, mulb, mula);

    test_integral(x, "x", res, a, b);
    delExpTree(res);
  }

  /* integral of x w.r.t another variable */
  {
    // TODO: "xNum"; The variables that are not integrated towards are converted
    // from EXP_VAR to EXP_NUM???
    ExpTree *xNum = newExpLeaf(EXP_NUM, strdup("x"));

    ExpTree *mula = newExpOp(EXP_MUL_OP, cpyExpTree(xNum), cpyExpTree(a));
    ExpTree *mulb = newExpOp(EXP_MUL_OP, cpyExpTree(xNum), cpyExpTree(b));
    ExpTree *res = newExpOp(EXP_SUB_OP, mulb, mula);
    test_integral(x, "y", res, a, b);

    delExpTree(xNum);
    delExpTree(res);
  }

  /* integral of a simple expression: x^2 */
  {
    ExpTree *pow2 = newExpOp(EXP_EXP_OP, cpyExpTree(x), cpyExpTree(two));
    ExpTree *numExp = newExpOp(EXP_ADD_OP, cpyExpTree(two), cpyExpTree(one));
    numExp = newExpOp(EXP_DIV_OP, cpyExpTree(one), numExp);
    ExpTree *powa3 = newExpOp(EXP_EXP_OP, cpyExpTree(a), cpyExpTree(three));
    ExpTree *powb3 = newExpOp(EXP_EXP_OP, cpyExpTree(b), cpyExpTree(three));
    ExpTree *mula = newExpOp(EXP_MUL_OP, cpyExpTree(numExp), powa3);
    ExpTree *mulb = newExpOp(EXP_MUL_OP, cpyExpTree(numExp), powb3);
    ExpTree *res = newExpOp(EXP_SUB_OP, mulb, mula);
    test_integral(pow2, "x", res, a, b);

    delExpTree(res);
    delExpTree(numExp);
    delExpTree(pow2);
  }

  /* integral of the polynomial: x^3 + 42x^2 + 10x - y */
  {
    // TODO: "yNum"; The variables that are not integrated towards are converted
    // from EXP_VAR to EXP_NUM???
    ExpTree *yNum = newExpLeaf(EXP_NUM, "y");
    ExpTree *num42 = newExpLeaf(EXP_NUM, "42");
    ExpTree *num10 = newExpLeaf(EXP_NUM, "10");
    ExpTree *polynomial = newExpOp(
        EXP_SUB_OP,
        newExpOp(
            EXP_ADD_OP,
            newExpOp(
                EXP_ADD_OP,
                newExpOp(EXP_EXP_OP, cpyExpTree(x), cpyExpTree(three)),
                newExpOp(EXP_MUL_OP, cpyExpTree(num42),
                         newExpOp(EXP_EXP_OP, cpyExpTree(x), cpyExpTree(two)))),
            newExpOp(EXP_MUL_OP, cpyExpTree(num10), cpyExpTree(x))),
        cpyExpTree(y));

    ExpTree *numExp1 = newExpOp(EXP_ADD_OP, cpyExpTree(three), cpyExpTree(one));
    numExp1 = newExpOp(EXP_DIV_OP, cpyExpTree(one), numExp1);
    ExpTree *numExp2 = newExpOp(EXP_ADD_OP, cpyExpTree(two), cpyExpTree(one));
    numExp2 = newExpOp(EXP_DIV_OP, cpyExpTree(one), numExp2);
    ExpTree *powa4 = newExpOp(EXP_EXP_OP, cpyExpTree(a), cpyExpTree(four));
    ExpTree *powa3 = newExpOp(EXP_EXP_OP, cpyExpTree(a), cpyExpTree(three));
    ExpTree *powa2 = newExpOp(EXP_EXP_OP, cpyExpTree(a), cpyExpTree(two));
    ExpTree *mula1 = newExpOp(EXP_MUL_OP, cpyExpTree(numExp1), powa4);
    ExpTree *mula2 = newExpOp(EXP_MUL_OP, cpyExpTree(num42), cpyExpTree(a));
    ExpTree *mula3 = newExpOp(EXP_MUL_OP, cpyExpTree(numExp2), powa3);
    ExpTree *mula4 = newExpOp(EXP_MUL_OP, mula2, mula3);
    ExpTree *adda1 = newExpOp(EXP_ADD_OP, mula1, mula4);
    ExpTree *mula5 = newExpOp(EXP_MUL_OP, cpyExpTree(num10), cpyExpTree(a));
    ExpTree *mula6 = newExpOp(EXP_MUL_OP, cpyExpTree(half), powa2);
    ExpTree *mula7 = newExpOp(EXP_MUL_OP, mula5, mula6);
    ExpTree *adda2 = newExpOp(EXP_ADD_OP, adda1, mula7);
    ExpTree *ya = newExpOp(EXP_MUL_OP, cpyExpTree(yNum), cpyExpTree(a));
    ExpTree *suba = newExpOp(EXP_SUB_OP, adda2, ya);
    ExpTree *powb4 = newExpOp(EXP_EXP_OP, cpyExpTree(b), cpyExpTree(four));
    ExpTree *powb3 = newExpOp(EXP_EXP_OP, cpyExpTree(b), cpyExpTree(three));
    ExpTree *powb2 = newExpOp(EXP_EXP_OP, cpyExpTree(b), cpyExpTree(two));
    ExpTree *mulb1 = newExpOp(EXP_MUL_OP, cpyExpTree(numExp1), powb4);
    ExpTree *mulb2 = newExpOp(EXP_MUL_OP, cpyExpTree(num42), cpyExpTree(b));
    ExpTree *mulb3 = newExpOp(EXP_MUL_OP, cpyExpTree(numExp2), powb3);
    ExpTree *mulb4 = newExpOp(EXP_MUL_OP, mulb2, mulb3);
    ExpTree *addb1 = newExpOp(EXP_ADD_OP, mulb1, mulb4);
    ExpTree *mulb5 = newExpOp(EXP_MUL_OP, cpyExpTree(num10), cpyExpTree(b));
    ExpTree *mulb6 = newExpOp(EXP_MUL_OP, cpyExpTree(half), powb2);
    ExpTree *mulb7 = newExpOp(EXP_MUL_OP, mulb5, mulb6);
    ExpTree *addb2 = newExpOp(EXP_ADD_OP, addb1, mulb7);
    ExpTree *yb = newExpOp(EXP_MUL_OP, cpyExpTree(yNum), cpyExpTree(b));
    ExpTree *subb = newExpOp(EXP_SUB_OP, addb2, yb);
    ExpTree *res = newExpOp(EXP_SUB_OP, subb, suba);
    test_integral(polynomial, "x", res, a, b);

    /* clean */
    delExpTree(polynomial);
    delExpTree(num10);
    delExpTree(num42);
    delExpTree(numExp1);
    delExpTree(numExp2);
    delExpTree(res);
    delExpTree(yNum);
  }

  delExpTree(x);
  delExpTree(y);
  delExpTree(a);
  delExpTree(b);
  delExpTree(half);
  delExpTree(one);
  delExpTree(two);
  delExpTree(three);
  delExpTree(four);
  delExpTree(five);

  return 0;
}
