#include "funexp.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void test_integral(ExpTree *expr, char *var, const char *expected_msg) {
  /* evaluate integral */
  ExpTree *integral_expr = integral(expr, var);
  assert(integral_expr != NULL);

  /* printing and comparison */
  char *buffer = NULL;
  size_t buflen = 0;
  FILE *stream = open_memstream(&buffer, &buflen);
  assert(stream != NULL);
  printExpTree(integral_expr, stream);
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

  /* integral of a constant */
  test_integral(newExpLeaf(EXP_NUM, strdup("5")), "x", "(5 * x)");

  /* integral of x w.r.t itself */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, strdup("x"));
    test_integral(x, "x", "(0.5 * (x^2))");
    delExpTree(x);
  }

  /* integral of x w.r.t another variable */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, strdup("x"));
    test_integral(x, "y", "(x * y)");
    delExpTree(x);
  }

  /* integral of a simple expression: x^2 */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, strdup("x"));
    test_integral(newExpOp(EXP_EXP_OP, x, newExpLeaf(EXP_NUM, "2")), "x",
                  "((1 / (2 + 1)) * (x^3))");
    delExpTree(x);
  }

  /* integral of the polynomial: x^3 + 42x^2 + 10x - y */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, strdup("x"));
    ExpTree *y = newExpLeaf(EXP_VAR, strdup("y"));

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

    /* Test integral of the polynomial */
    test_integral(polynomial, "x",
                  "(((((1 / (3 + 1)) * (x^4)) + ((42 * x) * ((1 / (2 + 1)) * "
                  "(x^3)))) + ((10 * x) * (0.5 * (x^2)))) - (y * x))");

    /* clean */
    delExpTree(x);
    delExpTree(y);
    delExpTree(polynomial);
  }

  /* Test integral of sin(x) */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, strdup("x"));
    ExpTree *sin_x = newExpTree(EXP_FUN, strdup("sin(x)"), cpyExpTree(x), NULL);
    test_integral(sin_x, "x", "(-1 * cos(x))");
    delExpTree(x);
    delExpTree(sin_x);
  }

  /* Test integral of sin(2x) */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, strdup("x"));
    ExpTree *sin_2x =
        newExpTree(EXP_FUN, strdup("sin(2x)"), cpyExpTree(x), NULL);
    test_integral(sin_2x, "x", "(-1/2.0)*cos(2.0)(x)");
    delExpTree(x);
    delExpTree(sin_2x);
  }

  /* Test integral of cos(x) */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, strdup("x"));
    ExpTree *cos_x = newExpTree(EXP_FUN, strdup("cos(x)"), cpyExpTree(x), NULL);
    test_integral(cos_x, "x", "(1 * sin(x))");
    delExpTree(x);
    delExpTree(cos_x);
  }

  /* Test integral of cos(nx) */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, strdup("x"));
    ExpTree *cos_2x =
        newExpTree(EXP_FUN, strdup("cos(2x)"), cpyExpTree(x), NULL);
    test_integral(cos_2x, "x", "(1/2.0)*sin(2.0)(x)");
    delExpTree(x);
    delExpTree(cos_2x);
  }

  /* Test integral of sqrt(x) */
  {
    ExpTree *x = newExpLeaf(EXP_VAR, strdup("x"));
    ExpTree *sqrt_x = newExpTree(EXP_FUN, strdup("sqrt"), cpyExpTree(x), NULL);
    test_integral(sqrt_x, "x", "((2/3) * (x^(3/2)))");
    delExpTree(x);
    delExpTree(sqrt_x);
  }

  return 0;
}
