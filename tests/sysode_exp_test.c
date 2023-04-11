#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sysode.h"

int main(int argc, char *argv[]) {
  /* (-b + sqrt((b^2) - ((4a)c))) / (2a) */
  ODEExpTree *b1 = newOdeExpLeaf(ODE_VAR, strdup("b\0"));
  ODEExpTree *b2 = newOdeExpLeaf(ODE_VAR, strdup("b\0"));
  ODEExpTree *n2 = newOdeExpLeaf(ODE_NUM, strdup("2\0"));
  ODEExpTree *n4 = newOdeExpLeaf(ODE_NUM, strdup("4\0"));
  ODEExpTree *a1 = newOdeExpLeaf(ODE_VAR, strdup("a\0"));
  ODEExpTree *c1 = newOdeExpLeaf(ODE_VAR, strdup("c\0"));
  ODEExpTree *m2 = newOdeExpLeaf(ODE_NUM, strdup("2\0"));
  ODEExpTree *a2 = newOdeExpLeaf(ODE_VAR, strdup("a\0"));
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
  printOdeExpTree(tree, stream);
  fprintf(stream, "%c", '\0');
  printf("%s\n", buffer);
  assert(strcmp(buffer, "(-b + sqrt((b^2) - ((4 * a) * c))) / (2 * a)"));
  /* clean */
  delOdeExpTree(tree);
  free(buffer);
  return 0;
}
