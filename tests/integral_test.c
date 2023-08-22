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

  /* clean */
  free(buffer);
}

int main(int argc, char *argv[]) {
  /* to avoid silly warnings about unused parameters */
  (void)argc;
  (void)argv;

  /* variables */
  ExpTree *x = newExpLeaf(EXP_VAR, strdup("x"));
  ExpTree *y = newExpLeaf(EXP_VAR, strdup("y"));

  /* integral of a constant */
  test_integral(newExpLeaf(EXP_NUM, strdup("5")), "x", "(5 * x)");

  /* integral of x w.r.t itself */
  test_integral(x, "x", "(0.5 * (x^2))");

  /* integral of x w.r.t another variable */
  test_integral(x, "y", "(x * y)");

  /* integral of a simple expression: x^2 */
  test_integral(newExpOp(EXP_EXP_OP, x, newExpLeaf(EXP_NUM, "3")), "x",
                "((1 / (3 + 1)) * (x^(3 + 1)))");

  /* clean */
  delExpTree(x);
  delExpTree(y);

  return 0;
}
