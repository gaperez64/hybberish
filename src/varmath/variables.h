#ifndef VARIABLES_H
#define VARIABLES_H

#include "interval.h"
#include <stdlib.h>

/* A vector of variables, each restricted to their own domain, as a linked list. */
typedef struct Domain {
  char *var;
  Interval domain;
  struct Domain *next;
} Domain;

/* Create a new, single element list.

  Transfers ownership of all allocated arguments to the newly created instance.
*/
Domain *newDomain(char *var, Interval domain);
/* Attach the second element as the head of the first list. */
Domain *appDomainElem(Domain *tail, Domain *head);
/* Allocate a new element and assign the passed member values. */
Domain *newDomainElem(Domain *tail, char *var, Interval domain);

/* Deallocate the given list. */
void delDomain(Domain *list);
/* Print a representation of the given list to the specified stream. */
void printDomain(Domain *list, FILE *where);

#endif
