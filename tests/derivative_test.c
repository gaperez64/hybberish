#include "funexp.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void test_derivative(ExpTree *expr, char *var, const char *expected_msg) {
  /* evaluate derivative */
  ExpTree *der = derivative(expr, var);
  assert(der != NULL);

  /* printing and comparison */
  char *buffer;
  size_t buflen;
  FILE *stream = open_memstream(&buffer, &buflen);
  assert(stream != NULL);
  printExpTree(der, stream);
  fclose(stream);
  printf("Expect: %s\n", expected_msg);
  printf("Got: %s\n", buffer);
  printf("!strcmp = %i\n", !strcmp(buffer, expected_msg));
  assert(!strcmp(buffer, expected_msg));
  fflush(stdout); // Flush stdout

  /* clean */
  free(buffer);
}

int main(int argc, char *argv[]) {
  /* to avoid silly warnings about unused parameters */
  (void)argc;
  (void)argv;

  /* derivative of a constant */
  test_derivative(newExpLeaf(EXP_NUM, strdup("5")), "x", "0");

  /* derivative of a variable w.r.t itself */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, strdup("x"));
    test_derivative(x, "x", "1");
    delExpTree(x);
  }

  /* derivative of a variable w.r.t another variable */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, strdup("x"));
    test_derivative(x, "y", "0");
    delExpTree(x);
  }

  /* derivative of a simple expression: x^2 */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, strdup("x"));
    test_derivative(newExpOp(EXP_EXP_OP, x, newExpLeaf(EXP_NUM, "2")), "x",
                    "((2 * 1) * (x^1))");
    delExpTree(x);
  }

  /* Construct the polynomial: x^3 + 42x^2 + 10x - y */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, strdup("x"));
    ExpTree *y = newExpLeaf(EXP_VAR, strdup("y"));

    /* Create copies of x and y */
    ExpTree *x_copy = cpyExpTree(x);
    ExpTree *y_copy = cpyExpTree(y);

    /* Create the polynomial expression using the copied variables */
    ExpTree *polynomial = newExpOp(
        EXP_SUB_OP,
        newExpOp(EXP_ADD_OP,
                 newExpOp(EXP_ADD_OP, newExpOp(EXP_EXP_OP, x_copy,
                                               newExpLeaf(EXP_NUM, "3")),
                          newExpOp(EXP_MUL_OP, newExpLeaf(EXP_NUM, "42"),
                                   newExpOp(EXP_EXP_OP, x_copy,
                                            newExpLeaf(EXP_NUM, "2")))),
                 newExpOp(EXP_MUL_OP, newExpLeaf(EXP_NUM, "10"), x_copy)),
        y_copy);

    /* Print the polynomial expression */
    printf("Polynomial expression: ");
    printExpTree(polynomial, stdout);
    printf("\n");

    /* Test derivative of the polynomial */
    test_derivative(polynomial, "x",
                    "(((((3 * 1) * (x^2)) + ((0 * (x^2)) + (42 "
                    "* ((2 * 1) * (x^1))))) + ((0 * x) + (10 * "
                    "1))) - 0)");

    /* Delete copied variables */
    delExpTree(x_copy);
    delExpTree(y_copy);
  }

  return 0;
}
