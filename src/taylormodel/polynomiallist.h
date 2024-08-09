#ifndef FUNEXP_LIST_H
#define FUNEXP_LIST_H

#include "funexp.h"
#include "sysode.h"
#include <stdio.h>
#include <string.h>

/* Lists of Taylor polynomials derived from a system of ODEs */
typedef struct TPList {
  char *fun;
  ExpTree *exp;
  struct TPList *next;
} TPList;

/* Create a new, single element list.

  Transfers ownership of allocated arguments to the newly created instance.
*/
TPList *newTPList(char *fun, ExpTree *exp);
/* Attach the second element as the head of the first list. */
TPList *appTPElem(TPList *tail, TPList *head);
/* Allocate a new element and assign the passed member values. */
TPList *newTPElem(TPList *tail, char *fun, ExpTree *exp);

/* Create a copy of the given list. */
TPList *cpyTPList(TPList *list);

/* Reverse the given list in-place and return the new head. */
TPList *reverseTPList(TPList *list);

/* Deallocate the given list. */
void delTPList(TPList *list);
/* Print a representation of the given list to the specified stream. */
void printTPList(TPList *, FILE *list);

#endif
