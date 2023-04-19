#ifndef SYSODE_H
#define SYSODE_H

#include <stdio.h>

/* Ordinary differential expression trees */
typedef enum ODEExpType {
  ODE_NUM,
  ODE_ADD_OP,
  ODE_SUB_OP,
  ODE_MUL_OP,
  ODE_DIV_OP,
  ODE_EXP_OP,
  ODE_NEG,
  ODE_VAR,
  ODE_FUN,
} ODEExpType;

typedef struct ODEExpTree {
  char *data;
  ODEExpType type;
  struct ODEExpTree *left;
  struct ODEExpTree *right;
} ODEExpTree;

ODEExpTree *newOdeExpLeaf(ODEExpType, char *);
ODEExpTree *newOdeExpOp(ODEExpType, ODEExpTree *, ODEExpTree *);
ODEExpTree *newOdeExpTree(ODEExpType, char *, ODEExpTree *, ODEExpTree *);
void delOdeExpTree(ODEExpTree *);
void printOdeExpTree(ODEExpTree *, FILE *);

/* Lists of ODEs */
typedef struct ODEList {
  char *fun;
  ODEExpTree *exp;
  struct ODEList *next;
} ODEList;

ODEList *newOdeList(char *, ODEExpTree *);
ODEList *newOdeElem(ODEList *, char *, ODEExpTree *);
/* Append the second one to the first list,
 * It assumes the first one is a single node!
 */
ODEList *appOdeElem(ODEList *, ODEList*);
void delOdeList(ODEList *);
void printOdeList(ODEList *, FILE *);

#endif
