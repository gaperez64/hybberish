#include <stdio.h>

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

typedef struct ODEList {
  char *var;
  ODEExpTree *exp;
  struct ODEList *next;
} ODEList;

ODEList *newOdeList(char *, ODEExpTree *);
