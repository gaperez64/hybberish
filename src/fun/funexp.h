#ifndef FUNEXP_H
#define FUNEXP_H

#include <math.h>
#include <stdio.h>
#include <regex.h>
#include <stdbool.h>

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
bool isLinear(ExpTree *expr);

/* Verify the exact equality of the given trees.
  Both structure and content must match exactly. */
bool isEqual(ExpTree *expr1, ExpTree *expr2);

/* Compute the degree of the given monomial.

  A monomial is a polynomial of only a single term.
  The expression can only contain multiplication (*)
  unary negative (-) and exponentiation (^) with
  non-negative powers.
*/
unsigned int degreeMonomial(ExpTree *expr);



#endif
