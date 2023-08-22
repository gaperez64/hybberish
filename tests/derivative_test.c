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

  /* derivative of a constant */
  test_derivative(newExpLeaf(EXP_NUM, strdup("5")), "x", "0");

  /* derivative of a variable w.r.t itself */
  test_derivative(x, "x", "1");

  /* derivative of a variable w.r.t another variable */
  test_derivative(x, "y", "0");

  /* derivative of a simple expression: x^2 */
  test_derivative(newExpOp(EXP_EXP_OP, x, newExpLeaf(EXP_NUM, "2")), "x",
                  "((2 * 1) * (x^(2 - 1)))");

  /* clean */
  delExpTree(x);
  delExpTree(y);
  return 0;
}
