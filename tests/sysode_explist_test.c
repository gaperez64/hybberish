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
  ExpTree *b1 = newExpLeaf(EXP_VAR, "b");
  ExpTree *b2 = newExpLeaf(EXP_VAR, "b");
  ExpTree *n2 = newExpLeaf(EXP_NUM, "2");
  ExpTree *n4 = newExpLeaf(EXP_NUM, "4");
  ExpTree *a1 = newExpLeaf(EXP_VAR, "a");
  ExpTree *c1 = newExpLeaf(EXP_VAR, "c");
  ExpTree *m2 = newExpLeaf(EXP_NUM, "2");
  ExpTree *a2 = newExpLeaf(EXP_VAR, "at");

  /* leaves ready, now build the trees */
  ExpTree *t1 = newExpOp(EXP_NEG, b1, NULL);
  ExpTree *exp = newExpOp(EXP_EXP_OP, b2, n2);
  ExpTree *foura = newExpOp(EXP_MUL_OP, n4, a1);
  ExpTree *fourac = newExpOp(EXP_MUL_OP, foura, c1);
  ExpTree *min = newExpOp(EXP_SUB_OP, exp, fourac);
  ExpTree *t2 = newExpTree(EXP_FUN, strdup("sqrt"), min, NULL);
  ExpTree *t3 = newExpOp(EXP_MUL_OP, m2, a2);

  /* tress ready, now creating the list */
  ODEList *list = newOdeList(strdup("x"), t1);
  list = newOdeElem(list, strdup("y"), t2);
  list = newOdeElem(list, strdup("last"), t3);

  /* printing */
  char buffer[100];
  FILE *stream = fmemopen(buffer, 100, "w");
  assert(stream != NULL);
  const char msg[] =
      "last' = (2 * at); y' = sqrt(((b^2) - ((4 * a) * c))); x' = -b; ";
  printOdeList(list, stream);
  fclose(stream); /* close to flush and write null byte */
  printf("expect: |%s| = %lu\n", msg, strlen(msg));
  printf("got: |%s| = %lu\n", buffer, strlen(buffer));
  printf("!strcmp = %i\n", !strcmp(buffer, msg));
  assert(!strcmp(buffer, msg));

  /* clean */
  delOdeList(list);
  return 0;
}
