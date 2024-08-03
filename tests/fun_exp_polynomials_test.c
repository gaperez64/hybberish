#include "funexp.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>



void testDegree(ExpTree *expression, unsigned int expectedDegree) {
  unsigned int degree = degreeMonomial(expression);

  printf("Expr:    ");
  printExpTree(expression, stdout);
  printf("\n");
  printf("Got deg: %u\n", degree);
  printf("Expect:  %u\n\n", expectedDegree);
  fflush(stdout);

  assert(degree == expectedDegree);
}



int main(int argc, char *argv[]) {
  /* to avoid silly warnings about unused parameters */
  (void)argc;
  (void)argv;

  /* Construct leaves */
  ExpTree *x = newExpLeaf(EXP_VAR, "x");
  ExpTree *y = newExpLeaf(EXP_VAR, "y");
  ExpTree *z = newExpLeaf(EXP_VAR, "z");
  ExpTree *num0 = newExpLeaf(EXP_NUM, "0");
  ExpTree *num2 = newExpLeaf(EXP_NUM, "2");
  ExpTree *num3 = newExpLeaf(EXP_NUM, "3");
  ExpTree *num1Dot = newExpLeaf(EXP_NUM, "1.");
  ExpTree *num1DotZero = newExpLeaf(EXP_NUM, "1.000");
  ExpTree *numDotZero = newExpLeaf(EXP_NUM, ".000");


  /* Constants */
  {
    testDegree(num0, 0);
    testDegree(num2, 0);
    testDegree(num1Dot, 0);
    testDegree(num1DotZero, 0);
    testDegree(numDotZero, 0);
  }

  /* Lone variables */
  {
    testDegree(x, 1);
    testDegree(y, 1);
    testDegree(z, 1);
  }

  /* x^0 */
  {
    ExpTree *expr = newExpOp(EXP_EXP_OP, cpyExpTree(x), cpyExpTree(num0));
    testDegree(expr, 0);

    delExpTree(expr);
  }

  /* x^.000 */
  {
    ExpTree *expr = newExpOp(EXP_EXP_OP, cpyExpTree(x), cpyExpTree(numDotZero));
    testDegree(expr, 0);

    delExpTree(expr);
  }

  /* z^1. */
  {
    ExpTree *expr = newExpOp(EXP_EXP_OP, cpyExpTree(z), cpyExpTree(num1Dot));
    testDegree(expr, 1);

    delExpTree(expr);
  }

  /* ((x^2) * (y^0)) */
  {
    ExpTree *expLeft = newExpOp(EXP_EXP_OP, cpyExpTree(x), cpyExpTree(num2));
    ExpTree *expRight = newExpOp(EXP_EXP_OP, cpyExpTree(y), cpyExpTree(num0));
    ExpTree *expr = newExpOp(EXP_MUL_OP, expLeft, expRight);
    testDegree(expr, 2);

    delExpTree(expr);
  }

  /* ((2 * x^3) * y^1.000) */
  {
    ExpTree *pow3 = newExpOp(EXP_EXP_OP, cpyExpTree(x), cpyExpTree(num3));
    ExpTree *pow1DotZero = newExpOp(EXP_EXP_OP, cpyExpTree(y), cpyExpTree(num1DotZero));
    ExpTree *multLeft = newExpOp(EXP_MUL_OP, cpyExpTree(num2), pow3);
    ExpTree *expr = newExpOp(EXP_MUL_OP, multLeft, pow1DotZero);
    testDegree(expr, 4);

    delExpTree(expr);
  }

  /* Negation */
  {
    ExpTree *negNum = newExpOp(EXP_NEG, cpyExpTree(num3), NULL);
    ExpTree *negVar = newExpOp(EXP_NEG, cpyExpTree(x), NULL);

    /* -(-(2 * x^3) * -(y^1.000)) */
    ExpTree *pow3 = newExpOp(EXP_EXP_OP, cpyExpTree(x), cpyExpTree(num3));
    ExpTree *pow1DotZero = newExpOp(EXP_NEG, newExpOp(EXP_EXP_OP, cpyExpTree(y), cpyExpTree(num1DotZero)), NULL);
    ExpTree *multLeft = newExpOp(EXP_NEG, newExpOp(EXP_MUL_OP, cpyExpTree(num2), pow3), NULL);
    ExpTree *negExp = newExpOp(EXP_NEG, newExpOp(EXP_MUL_OP, multLeft, pow1DotZero), NULL);

    testDegree(negNum, 0);
    testDegree(negVar, 1);
    testDegree(negExp, 4);

    delExpTree(negNum);
    delExpTree(negVar);
    delExpTree(negExp);
  }

  /* Clean */
  delExpTree(x);
  delExpTree(y);
  delExpTree(z);
  delExpTree(num2);
  delExpTree(num0);
  delExpTree(num3);
  delExpTree(num1Dot);
  delExpTree(num1DotZero);
}
