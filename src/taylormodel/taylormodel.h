#ifndef TAYLOR_MODEL_H
#define TAYLOR_MODEL_H

#include "funexp.h"
#include "interval.h"
#include "transformations.h"
#include "variables.h"
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
TaylorModel *newTaylorModel(char *const fun, ExpTree *const exp,
                            Interval remainder);
/* Attach the second element as the head of the first list. */
TaylorModel *appTMElem(TaylorModel *const tail, TaylorModel *head);
/* Allocate a new element and assign the passed member values. */
TaylorModel *newTMElem(TaylorModel *const tail, char *const fun,
                       ExpTree *const exp, const Interval remainder);

/* Create a copy of the given list. */
TaylorModel *cpyTaylorModel(const TaylorModel *const list);
/* Create a copy of only the head of the given list. Ignore the tail. */
TaylorModel *cpyTaylorModelHead(const TaylorModel *const list);

/* Reverse the linked list in-place and return the new head. */
TaylorModel *reverseTaylorModel(TaylorModel *const list);

/* Deallocate the given list. */
void delTaylorModel(TaylorModel *list);
/* Print a representation of the given list to the specified stream. */
void printTaylorModel(const TaylorModel *const list, FILE *where);

/* Perform interval-valued expression evaluation via interval arithmetic. */
Interval evaluateExpTree(const ExpTree *const tree,
                         const Domain *const domains);

/* Perform real-valued expression evaluation. */
double evaluateExpTreeReal(const ExpTree *const tree,
                           const Valuation *const values);

/* Perform expression evaluation via Taylor model arithmetic.

  tree: The tree to evaluate.
  list: The TMs to substitute the variables by
  fun: The variable w.r.t. which all TMs are constructed
  variables: The domain of each variable, for interval enclosure computation
  k: The truncation order
*/
TaylorModel *evaluateExpTreeTM(const ExpTree *const tree,
                               const TaylorModel *const list,
                               const char *const fun,
                               const Domain *const variables,
                               const unsigned int k);

/*
  Order k Taylor Model arithmetic
*/
TaylorModel *addTM(const TaylorModel *const left,
                   const TaylorModel *const right,
                   const Domain *const variables, const unsigned int k);
TaylorModel *subTM(const TaylorModel *const left,
                   const TaylorModel *const right,
                   const Domain *const variables, const unsigned int k);
TaylorModel *mulTM(const TaylorModel *const left,
                   const TaylorModel *const right,
                   const Domain *const variables, const unsigned int k);
TaylorModel *divTM(const TaylorModel *const left,
                   const TaylorModel *const right,
                   const Domain *const variables, const unsigned int k);
TaylorModel *negTM(const TaylorModel *const list, const Domain *const variables,
                   const unsigned int k);
TaylorModel *powTM(const TaylorModel *const left, const unsigned int right,
                   const Domain *const variables, const unsigned int k);
/* Definite integral of a Taylor Model.

  list: The Taylor model list to integrate
  variables: The list of variable valuations of all variables in each of the
  Taylor model polynomials intDomain: The interval w.r.t. which to integrate all
  Taylor models intVar: The variable w.r.t. which to integrate all Taylor models
  k: The order of the resulting Taylor model
*/
TaylorModel *intTM(const TaylorModel *const list,
                   const Interval *const intDomain, const char *const intVar,
                   const Domain *const variables, const unsigned int k);

/* Truncate the given Taylor model to order k.

  Truncates all terms of degree > k.
  Adds the interval enclosure of the truncated terms, w.r.t. the
  specified variables valuation, to the model's remainder interval. */
TaylorModel *truncateTM(const TaylorModel *const list,
                        const Domain *const variables, const unsigned int k);

#endif
