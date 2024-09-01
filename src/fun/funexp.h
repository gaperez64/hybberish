#ifndef FUNEXP_H
#define FUNEXP_H

#include <math.h>
#include <regex.h>
#include <stdbool.h>
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

/* Construct a new leaf node.

  Internally duplicate the given name/data.
  Does NOT transfer ownership of any arguments to the newly created instance.
  Leaves have functional use for their data. So moving the string duplication
  to the 'constructor' internals simplifies most calls to this function.
*/
ExpTree *newExpLeaf(const ExpType, const char *const);

/* Construct an internal, operator node.

  All pointer arguments must be heap allocated or NULL.
  Transfers ownership of all allocated arguments to the newly created instance.
*/
ExpTree *newExpOp(ExpType, ExpTree *, ExpTree *);

/* Allocate and compose a new expression tree.

  All pointer arguments must be heap allocated or NULL.
  Transfers ownership of all allocated arguments to the newly created instance.
*/
ExpTree *newExpTree(ExpType, char *, ExpTree *, ExpTree *);
void delExpTree(ExpTree *);
void printExpTree(ExpTree *, FILE *);
ExpTree *derivative(ExpTree *, char *);
ExpTree *integral(ExpTree *expr, char *var);
ExpTree *definiteIntegral(ExpTree *expr, char *var, ExpTree *lowerBound,
                          ExpTree *upperBound);
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

/* A forward declaration of the substitute tree transformation. */
ExpTree *substitute(ExpTree *source, char *var, ExpTree *target);

#endif
