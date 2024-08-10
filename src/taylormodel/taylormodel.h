#ifndef FUNEXP_LIST_H
#define FUNEXP_LIST_H

#include "funexp.h"
#include "interval.h"
#include "sysode.h"
#include <stdio.h>
#include <string.h>

/* A vector of Taylor Models, as a linked list.

  The models get derived from a system of ODEs. So, each model corresponds to a
  single variable solved for within the system.
*/
typedef struct TaylorModel {
  char *fun;
  ExpTree *exp;
  Interval remainder;
  struct TaylorModel *next;
} TaylorModel;

/* Create a new, single element list.

  Transfers ownership of all allocated arguments to the newly created instance.
*/
TaylorModel *newTaylorModel(char *fun, ExpTree *exp, Interval remainder);
/* Attach the second element as the head of the first list. */
TaylorModel *appTMElem(TaylorModel *tail, TaylorModel *head);
/* Allocate a new element and assign the passed member values. */
TaylorModel *newTMElem(TaylorModel *tail, char *fun, ExpTree *exp,
                       Interval remainder);

/* Create a copy of the given list. */
TaylorModel *cpyTaylorModel(TaylorModel *list);

/* Reverse the linked list in-place and return the new head. */
TaylorModel *reverseTaylorModel(TaylorModel *list);

/* Deallocate the given list. */
void delTaylorModel(TaylorModel *list);
/* Print a representation of the given list to the specified stream. */
void printTaylorModel(TaylorModel *, FILE *list);

#endif
