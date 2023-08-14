#ifndef SYSODE_H
#define SYSODE_H

#include <stdio.h>

#include "funexp.h"

/* Lists of ODEs */
typedef struct ODEList {
  char *fun;
  ExpTree *exp;
  struct ODEList *next;
} ODEList;

ODEList *newOdeList(char *, ExpTree *);
ODEList *newOdeElem(ODEList *, char *, ExpTree *);
/* Append the second one to the first list,
 * It assumes the second one is a single node!
 */
ODEList *appOdeElem(ODEList *, ODEList *);
void delOdeList(ODEList *);
void printOdeList(ODEList *, FILE *);

#endif
