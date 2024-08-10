#include "taylormodel.h"
#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

TaylorModel *newTaylorModel(char *fun, ExpTree *exp, Interval remainder) {
  assert(fun != NULL);
  assert(exp != NULL);

  TaylorModel *list = (TaylorModel *)malloc(sizeof(TaylorModel));
  list->fun = fun;
  list->exp = exp;
  list->remainder = remainder;
  list->next = NULL;
  return list;
}

TaylorModel *appTMElem(TaylorModel *tail, TaylorModel *head) {
  assert(head != NULL);
  assert(head->next == NULL);

  head->next = tail;
  return head;
}

TaylorModel *newTMElem(TaylorModel *tail, char *fun, ExpTree *exp,
                       Interval remainder) {
  TaylorModel *head = newTaylorModel(fun, exp, remainder);
  return appTMElem(tail, head);
}

void delTaylorModel(TaylorModel *list) {
  if (list->next != NULL)
    delTaylorModel(list->next);
  assert(list->fun != NULL);
  free(list->fun);
  assert(list->exp != NULL);
  delExpTree(list->exp);
  free(list);
}

void printTaylorModel(TaylorModel *list, FILE *where) {
  assert(list->fun != NULL);
  fprintf(where, "(p(%s) = ", list->fun);
  assert(list->exp != NULL);

  printExpTree(list->exp, where);
  fprintf(where, ", ");
  printInterval(&list->remainder, where);
  fprintf(where, ")");

  fprintf(where, "; ");
  if (list->next != NULL)
    printTaylorModel(list->next, where);
}

TaylorModel *cpyTaylorModel(TaylorModel *list) {
  /* Base case: The tail/next of the last element is NULL. */
  if (list == NULL)
    return NULL;

  /* Recursive case: The tail of the new element is everything built until now.
   */
  char *fun = strdup(list->fun);
  ExpTree *exp = cpyExpTree(list->exp);
  return newTMElem(cpyTaylorModel(list->next), fun, exp, list->remainder);
}

TaylorModel *reverseTaylorModel(TaylorModel *list) {
  /* Enforce preconditions */
  assert(list != NULL);

  /* Base case: push up the last element. */
  if (list->next == NULL)
    return list;

  /* Recursive case: Reverse tail and retrieve new head. */
  TaylorModel *lastElem = reverseTaylorModel(list->next);

  /* Only update upon backtracking, else would recurse on wrong or NULL ptrs. */
  assert(list->next != NULL);
  list->next->next = list; // reverse link: A->B  =>  A<-B
  list->next = NULL;       // sever own next: A->B  =>  A->NULL

  return lastElem;
}
