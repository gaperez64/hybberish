#ifndef FUNEXP_H
#define FUNEXP_H

#include <math.h>
#include <stdio.h>

/* Function expression trees */
typedef enum ExpType {
  EXP_NUM,
  EXP_ADD_OP,
  EXP_SUB_OP,
  EXP_MUL_OP,
  EXP_DIV_OP,
  EXP_EXP_OP,
  EXP_NEG,
  EXP_VAR,
  EXP_FUN,
} ExpType;

typedef struct ExpTree {
  char *data;
  ExpType type;
  struct ExpTree *left;
  struct ExpTree *right;
} ExpTree;

ExpTree *newExpLeaf(ExpType, char *);
ExpTree *newExpOp(ExpType, ExpTree *, ExpTree *);
ExpTree *newExpTree(ExpType, char *, ExpTree *, ExpTree *);
void delExpTree(ExpTree *);
void printExpTree(ExpTree *, FILE *);
ExpTree *derivative(ExpTree *, char *);
ExpTree *integral(ExpTree *expr, char *var);
ExpTree *cpyExpTree(ExpTree *src);

#endif
