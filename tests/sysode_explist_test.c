#include "sysode.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
  /* to avoid silly warnings about unused parameters */
  (void)argc;
  (void)argv;

  /* x' = -b; y' = sqrt(b^2 - 4ac); last' = 2at */
  ODEExpTree *b1 = newOdeExpLeaf(ODE_VAR, strdup("b"));
  ODEExpTree *b2 = newOdeExpLeaf(ODE_VAR, strdup("b"));
  ODEExpTree *n2 = newOdeExpLeaf(ODE_NUM, strdup("2"));
  ODEExpTree *n4 = newOdeExpLeaf(ODE_NUM, strdup("4"));
  ODEExpTree *a1 = newOdeExpLeaf(ODE_VAR, strdup("a"));
  ODEExpTree *c1 = newOdeExpLeaf(ODE_VAR, strdup("c"));
  ODEExpTree *m2 = newOdeExpLeaf(ODE_NUM, strdup("2"));
  ODEExpTree *a2 = newOdeExpLeaf(ODE_VAR, strdup("at"));

  /* leaves ready, now build the trees */
  ODEExpTree *t1 = newOdeExpOp(ODE_NEG, b1, NULL);

  ODEExpTree *exp = newOdeExpOp(ODE_EXP_OP, b2, n2);
  ODEExpTree *foura = newOdeExpOp(ODE_MUL_OP, n4, a1);
  ODEExpTree *fourac = newOdeExpOp(ODE_MUL_OP, foura, c1);
  ODEExpTree *min = newOdeExpOp(ODE_SUB_OP, exp, fourac);
  ODEExpTree *t2 = newOdeExpTree(ODE_FUN, strdup("sqrt"), min, NULL);

  ODEExpTree *t3 = newOdeExpOp(ODE_MUL_OP, m2, a2);

  /* tress ready, now creating the list */
  ODEList *list = newOdeList(strdup("x"), t1);
  list = newOdeElem(list, strdup("y"), t2);
  list = newOdeElem(list, strdup("last"), t3);

  /* printing */
  char *buffer = (char *)malloc(100 * sizeof(char));
  FILE *stream = fmemopen(buffer, 100, "w");
  const char *msg =
      "last' = (2 * at); y' = sqrt(((b^2) - ((4 * a) * c))); x' = -b";
  printOdeList(list, stream);
  fprintf(stream, "%c", '\0');
  printf("expect: %s\n", msg);
  printf("got: %s\n", buffer);
  assert(strcmp(buffer, msg));

  /* clean */
  delOdeList(list);
  free(buffer);
  return 0;
}
