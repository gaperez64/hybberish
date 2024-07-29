/* Functions to transform an expression to a desirable form. */

#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H


#include "funexp.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>



/* Simplify the expression through algebraic manipulations. */
ExpTree* simplify(ExpTree *source);

/* Convert the expression to a sum of products. */
ExpTree* toSumOfProducts(ExpTree *source);

/* Convert the expression to Horner form. */
ExpTree* toHornerForm(ExpTree *source);

/* Truncate all terms of order strictly greater than k. */
ExpTree* truncate(ExpTree *source, unsigned int k, char *variable);



/*
    Helper methods.
*/

/* Simplify a given expression by applying any found absorbing and neutral elements to their operators. */
ExpTree *simplifyOperators(ExpTree *source);

/* Check if the given expression is a number leaf with data equivalent to '0'. */
bool isZeroExpTree(ExpTree *source);

/* Create a number leaf with data equivalent to '0'. */
ExpTree * newZeroExpTree();

/* Check if the given expression is a number leaf with data equivalent to '1'. */
bool isOneExpTree(ExpTree* source);

/* Create a number leaf with data equivalent to '1'. */
ExpTree *newOneExpTree();


#endif
