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
  delExpTree(der);
}

int main(int argc, char *argv[]) {
  /* to avoid silly warnings about unused parameters */
  (void)argc;
  (void)argv;

  /* derivative of a constant */
  {
    ExpTree *exp = newExpLeaf(EXP_NUM, "5");
    test_derivative(exp, "x", "0");

    delExpTree(exp);
  }

  /* derivative of a variable w.r.t itself */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, "x");
    test_derivative(x, "x", "1");
    delExpTree(x);
  }

  /* derivative of a variable w.r.t another variable */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, "x");
    test_derivative(x, "y", "0");
    delExpTree(x);
  }

  /* derivative of a simple expression: x^2 */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, "x");
    ExpTree *exp =
        newExpOp(EXP_EXP_OP, cpyExpTree(x), newExpLeaf(EXP_NUM, "2"));
    test_derivative(exp, "x", "((2 * 1) * (x^1))");

    delExpTree(x);
    delExpTree(exp);
  }

  /* Construct the polynomial: x^3 + 42x^2 + 10x - y */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, "x");
    ExpTree *y = newExpLeaf(EXP_VAR, "y");

    /* Create the polynomial expression */
    ExpTree *polynomial = newExpOp(
        EXP_SUB_OP,
        newExpOp(
            EXP_ADD_OP,
            newExpOp(
                EXP_ADD_OP,
                newExpOp(EXP_EXP_OP, cpyExpTree(x), newExpLeaf(EXP_NUM, "3")),
                newExpOp(EXP_MUL_OP, newExpLeaf(EXP_NUM, "42"),
                         newExpOp(EXP_EXP_OP, cpyExpTree(x),
                                  newExpLeaf(EXP_NUM, "2")))),
            newExpOp(EXP_MUL_OP, newExpLeaf(EXP_NUM, "10"), cpyExpTree(x))),
        cpyExpTree(y));

    /* Print the polynomial expression */
    printf("Polynomial expression: ");
    printExpTree(polynomial, stdout);
    printf("\n");

    /* Test derivative of the polynomial */
    test_derivative(polynomial, "x",
                    "(((((3 * 1) * (x^2)) + ((0 * (x^2)) + (42 "
                    "* ((2 * 1) * (x^1))))) + ((0 * x) + (10 * "
                    "1))) - 0)");

    /* Free memory */
    delExpTree(x);
    delExpTree(y);
    delExpTree(polynomial);
  }

  /* derivative of sin(x) */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, "x");
    ExpTree *sin_x = newExpTree(EXP_FUN, strdup("sin"), cpyExpTree(x), NULL);
    test_derivative(sin_x, "x", "(cos(x) * 1)");
    delExpTree(x);
    delExpTree(sin_x);
  }

  /* derivative of cos(x) */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, "x");
    ExpTree *cos_x = newExpTree(EXP_FUN, strdup("Cos"), cpyExpTree(x), NULL);
    test_derivative(cos_x, "x", "(-1 * (sin(x) * 1))");
    delExpTree(x);
    delExpTree(cos_x);
  }

  /* derivative of sine of a polynomial x^3 + 42x^2 + 10x - y */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, "x");
    ExpTree *y = newExpLeaf(EXP_VAR, "y");

    /* Construct the polynomial: x^3 + 42x^2 + 10x - y */
    ExpTree *polynomial = newExpOp(
        EXP_SUB_OP,
        newExpOp(
            EXP_ADD_OP,
            newExpOp(
                EXP_ADD_OP,
                newExpOp(EXP_EXP_OP, cpyExpTree(x), newExpLeaf(EXP_NUM, "3")),
                newExpOp(EXP_MUL_OP, newExpLeaf(EXP_NUM, "42"),
                         newExpOp(EXP_EXP_OP, cpyExpTree(x),
                                  newExpLeaf(EXP_NUM, "2")))),
            newExpOp(EXP_MUL_OP, newExpLeaf(EXP_NUM, "10"), cpyExpTree(x))),
        cpyExpTree(y));

    /* Print the polynomial expression */
    printf("Polynomial expression: ");
    printExpTree(polynomial, stdout);
    printf("\n");

    /* Test derivative of the sine of the polynomial */
    ExpTree *sine_polynomial =
        newExpTree(EXP_FUN, strdup("sin"), polynomial, NULL);
    test_derivative(sine_polynomial, "x",
                    "(cos(((((x^3) + (42 * (x^2))) + (10 * x)) - y)) * (((((3 "
                    "* 1) * (x^2)) + ((0 * (x^2)) + (42 * ((2 * 1) * (x^1))))) "
                    "+ ((0 * x) + (10 * 1))) - 0))");

    /* Free memory */
    delExpTree(x);
    delExpTree(y);
    delExpTree(sine_polynomial);
  }

  /* derivative of cosine of a polynomial x^3 + 42x^2 + 10x - y */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, "x");
    ExpTree *y = newExpLeaf(EXP_VAR, "y");

    /* Construct the polynomial: x^3 + 42x^2 + 10x - y */
    ExpTree *polynomial = newExpOp(
        EXP_SUB_OP,
        newExpOp(
            EXP_ADD_OP,
            newExpOp(
                EXP_ADD_OP,
                newExpOp(EXP_EXP_OP, cpyExpTree(x), newExpLeaf(EXP_NUM, "3")),
                newExpOp(EXP_MUL_OP, newExpLeaf(EXP_NUM, "42"),
                         newExpOp(EXP_EXP_OP, cpyExpTree(x),
                                  newExpLeaf(EXP_NUM, "2")))),
            newExpOp(EXP_MUL_OP, newExpLeaf(EXP_NUM, "10"), cpyExpTree(x))),
        cpyExpTree(y));

    /* Print the polynomial expression */
    printf("Polynomial expression: ");
    printExpTree(polynomial, stdout);
    printf("\n");

    /* Test derivative of the sine of the polynomial */
    ExpTree *cosine_polynomial =
        newExpTree(EXP_FUN, strdup("cos"), polynomial, NULL);
    test_derivative(cosine_polynomial, "x",
                    "(-1 * (sin(((((x^3) + (42 * (x^2))) + (10 * x)) - y)) * "
                    "(((((3 * 1) * (x^2)) + ((0 * (x^2)) + (42 * ((2 * 1) * "
                    "(x^1))))) + ((0 * x) + (10 * 1))) - 0)))");

    /* Free memory */
    delExpTree(x);
    delExpTree(y);
    delExpTree(cosine_polynomial);
  }

  /* Test derivative of sqrt(x) */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, "x");
    ExpTree *sqrt_x = newExpTree(EXP_FUN, strdup("sqrt"), cpyExpTree(x), NULL);
    test_derivative(sqrt_x, "x", "(0.5 * (1 / sqrt(x)))");
    delExpTree(x);
    delExpTree(sqrt_x);
  }

  /* Test derivative of sqrt(x^3) */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, "x");
    ExpTree *x_cubed =
        newExpOp(EXP_EXP_OP, cpyExpTree(x), newExpLeaf(EXP_NUM, "3"));
    ExpTree *sqrt_x_cubed =
        newExpTree(EXP_FUN, strdup("sqrt"), cpyExpTree(x_cubed), NULL);
    test_derivative(sqrt_x_cubed, "x",
                    "(0.5 * (((3 * 1) * (x^2)) / sqrt((x^3))))");
    delExpTree(x);
    delExpTree(x_cubed);
    delExpTree(sqrt_x_cubed);
  }

  return 0;
}
