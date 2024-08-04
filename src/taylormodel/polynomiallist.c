#include "polynomiallist.h"
#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

TPList *newTPList(char *fun, ExpTree *exp) {
  TPList *list = (TPList *)malloc(sizeof(TPList));
  list->fun = fun;
  list->exp = exp;
  list->next = NULL;
  return list;
}

TPList *appTPElem(TPList *tail, TPList *head) {
  assert(head != NULL);
  assert(head->next == NULL);
  head->next = tail;
  return head;
}

TPList *newTPElem(TPList *tail, char *fun, ExpTree *exp) {
  TPList *head = newTPList(fun, exp);
  return appTPElem(tail, head);
}

void delTPList(TPList *list) {
  if (list->next != NULL)
    delTPList(list->next);
  assert(list->fun != NULL);
  free(list->fun);
  assert(list->exp != NULL);
  delExpTree(list->exp);
  free(list);
}

void printTPList(TPList *list, FILE *where) {
  assert(list->fun != NULL);
  fprintf(where, "p(%s) = ", list->fun);
  assert(list->exp != NULL);

  printExpTree(list->exp, where);

  fprintf(where, "; ");
  if (list->next != NULL)
    printTPList(list->next, where);
}

TPList *cpyTPList(TPList *list) {
  /* Base case: The tail/next of the last element is NULL. */
  if (list == NULL)
    return NULL;

  /* Recursive case: The tail of the new element is everything built until now.
   */
  char *fun = strdup(list->fun);
  ExpTree *exp = cpyExpTree(list->exp);
  return newTPElem(cpyTPList(list->next), fun, exp);
}

TPList *reverseTPList(TPList *list) {
  /* Enforce preconditions */
  assert(list != NULL);

  /* Base case: push up the last element. */
  if (list->next == NULL)
    return list;

  /* Recursive case: Retrieve new head. */
  TPList *lastElem = reverseTPList(list->next);

  /* While backtracking, reverse the old next link and sever your own next link.
   */
  assert(list->next != NULL);
  list->next->next = list; // reverse
  list->next = NULL;       // sever

  return lastElem;
}
