#include "funexp.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
  /* to avoid silly warnings about unused parameters */
  (void)argc;
  (void)argv;

  /* (-b + sqrt((b^2) - ((4a)c))) / (2a) */
  ExpTree *b1 = newExpLeaf(ODE_VAR, strdup("b"));
  ExpTree *b2 = newExpLeaf(ODE_VAR, strdup("b"));
  ExpTree *n2 = newExpLeaf(ODE_NUM, strdup("2"));
  ExpTree *n4 = newExpLeaf(ODE_NUM, strdup("4"));
  ExpTree *a1 = newExpLeaf(ODE_VAR, strdup("a"));
  ExpTree *c1 = newExpLeaf(ODE_VAR, strdup("c"));
  ExpTree *m2 = newExpLeaf(ODE_NUM, strdup("2"));
  ExpTree *a2 = newExpLeaf(ODE_VAR, strdup("a"));

  /* leaves ready, now build a tree */
  ExpTree *exp = newExpOp(ODE_EXP_OP, b2, n2);
  ExpTree *foura = newExpOp(ODE_MUL_OP, n4, a1);
  ExpTree *fourac = newExpOp(ODE_MUL_OP, foura, c1);
  ExpTree *min = newExpOp(ODE_SUB_OP, exp, fourac);
  ExpTree *sqrt = newExpTree(ODE_FUN, strdup("sqrt"), min, NULL);
  ExpTree *neg = newExpOp(ODE_NEG, b1, NULL);
  ExpTree *sum = newExpOp(ODE_ADD_OP, neg, sqrt);
  ExpTree *twoa = newExpOp(ODE_MUL_OP, m2, a2);
  ExpTree *tree = newExpOp(ODE_DIV_OP, sum, twoa);

  /* printing */
  char *buffer;
  size_t buflen;
  FILE *stream = open_memstream(&buffer, &buflen);
  assert(stream != NULL);
  const char msg[] = "((-b + sqrt(((b^2) - ((4 * a) * c)))) / (2 * a))";
  printExpTree(tree, stream);
  fclose(stream); /* close to flush and add a null byte*/
  printf("expect: |%s| = %lu\n", msg, strlen(msg));
  printf("got: |%s| = %lu\n", buffer, strlen(buffer));
  printf("!strcmp = %i\n", !strcmp(buffer, msg));
  assert(!strcmp(buffer, msg));
  fprintf(stderr, "done!\n");

  /* clean */
  delExpTree(tree);
  free(buffer);
  return 0;
}
