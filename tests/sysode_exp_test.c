#include "sysode.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
  /* to avoid silly warnings about unused parameters */
  (void)argc;
  (void)argv;

  /* (-b + sqrt((b^2) - ((4a)c))) / (2a) */
  ODEExpTree *b1 = newOdeExpLeaf(ODE_VAR, strdup("b"));
  ODEExpTree *b2 = newOdeExpLeaf(ODE_VAR, strdup("b"));
  ODEExpTree *n2 = newOdeExpLeaf(ODE_NUM, strdup("2"));
  ODEExpTree *n4 = newOdeExpLeaf(ODE_NUM, strdup("4"));
  ODEExpTree *a1 = newOdeExpLeaf(ODE_VAR, strdup("a"));
  ODEExpTree *c1 = newOdeExpLeaf(ODE_VAR, strdup("c"));
  ODEExpTree *m2 = newOdeExpLeaf(ODE_NUM, strdup("2"));
  ODEExpTree *a2 = newOdeExpLeaf(ODE_VAR, strdup("a"));

  /* leaves ready, now build a tree */
  ODEExpTree *exp = newOdeExpOp(ODE_EXP_OP, b2, n2);
  ODEExpTree *foura = newOdeExpOp(ODE_MUL_OP, n4, a1);
  ODEExpTree *fourac = newOdeExpOp(ODE_MUL_OP, foura, c1);
  ODEExpTree *min = newOdeExpOp(ODE_SUB_OP, exp, fourac);
  ODEExpTree *sqrt = newOdeExpTree(ODE_FUN, strdup("sqrt"), min, NULL);
  ODEExpTree *neg = newOdeExpOp(ODE_NEG, b1, NULL);
  ODEExpTree *sum = newOdeExpOp(ODE_ADD_OP, neg, sqrt);
  ODEExpTree *twoa = newOdeExpOp(ODE_MUL_OP, m2, a2);
  ODEExpTree *tree = newOdeExpOp(ODE_DIV_OP, sum, twoa);

  /* printing */
  char *buffer = (char *)malloc(100 * sizeof(char));
  FILE *stream = fmemopen(buffer, 100, "w");
  const char *msg = "((-b + sqrt(((b^2) - ((4 * a) * c)))) / (2 * a))";
  printOdeExpTree(tree, stderr);
  printOdeExpTree(tree, stream);
  fclose(stream); /* close to flush */
  printf("expect: |%s| = %lu\n", msg, strlen(msg));
  printf("got: |%s| = %lu\n", buffer, strlen(buffer));
  printf("!strcmp = %i\n", !strcmp(buffer, msg));
  assert(!strcmp(buffer, msg));

  /* clean */
  delOdeExpTree(tree);
  free(buffer);
  return 0;
}
