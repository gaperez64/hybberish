#include "funexp.h"
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void test_isLinear(ExpTree *expr, bool expected_result) {
  bool result = isLinear(expr);

  printf("Expression: ");
  printExpTree(expr, stdout);
  printf("\n");
  printf("Expected Result: %s\n", expected_result ? "true" : "false");
  printf("Result: %s\n", result ? "true" : "false");

  assert(result == expected_result);
}

int main(int argc, char *argv[]) {
  /* to avoid silly warnings about unused parameters */
  (void)argc;
  (void)argv;

  /* Test isLinear for sin(x) - should return true */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, strdup("x"));
    ExpTree *sin_x = newExpTree(EXP_FUN, strdup("sin"), cpyExpTree(x), NULL);
    test_isLinear(sin_x, true);
    delExpTree(x);
    delExpTree(sin_x);
  }

  /* Test isLinear for sin(x^2) - should return false */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, strdup("x"));
    ExpTree *x_squared =
        newExpOp(EXP_EXP_OP, cpyExpTree(x), newExpLeaf(EXP_NUM, "2"));
    ExpTree *sin_x_squared =
        newExpTree(EXP_FUN, strdup("sin"), cpyExpTree(x_squared), NULL);
    test_isLinear(sin_x_squared, false);
    delExpTree(x);
    delExpTree(x_squared);
    delExpTree(sin_x_squared);
  }

  /* Test isLinear for sin(2x^2) - should return false */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, strdup("x"));
    ExpTree *x_squared =
        newExpOp(EXP_EXP_OP, cpyExpTree(x), newExpLeaf(EXP_NUM, "2"));
    ExpTree *sin_2x_squared = newExpTree(
        EXP_FUN, strdup("sin"),
        newExpOp(EXP_MUL_OP, newExpLeaf(EXP_NUM, "2"), cpyExpTree(x_squared)),
        NULL);
    test_isLinear(sin_2x_squared, false);
    delExpTree(x);
    delExpTree(x_squared);
    delExpTree(sin_2x_squared);
  }

  return 0;
}
