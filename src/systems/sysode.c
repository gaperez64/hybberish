#include "sysode.h"
#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

ODEList *newOdeList(char *fun, ExpTree *exp) {
  ODEList *list = (ODEList *)malloc(sizeof(ODEList));
  list->fun = fun;
  list->exp = exp;
  list->next = NULL;
  return list;
}

ODEList *appOdeElem(ODEList *tail, ODEList *head) {
  assert(head != NULL);
  assert(head->next == NULL);
  head->next = tail;
  return head;
}

ODEList *newOdeElem(ODEList *tail, char *fun, ExpTree *exp) {
  ODEList *head = newOdeList(fun, exp);
  return appOdeElem(tail, head);
}

void delOdeList(ODEList *list) {
  if (list->next != NULL)
    delOdeList(list->next);
  assert(list->fun != NULL);
  free(list->fun);
  assert(list->exp != NULL);
  delExpTree(list->exp);
  free(list);
}

void printOdeList(ODEList *list, FILE *where) {
  assert(list->fun != NULL);
  fprintf(where, "%s' = ", list->fun);
  assert(list->exp != NULL);
  printExpTree(list->exp, where);
  fprintf(where, "; ");
  if (list->next != NULL)
    printOdeList(list->next, where);
}
