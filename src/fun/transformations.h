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
    Simplification Helper methods.
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



/*
    Sum of Products Helper methods.
*/

/* Apply the left-distributive property of multiplication (*) w.r.t.
  addition (+) and subtraction (-), recursively distributing the
  non-distributive left tree across the top-most, n-ary addition
  or subtraction that roots right.

  Note, right MUST be rooted by an addition or subtraction.
  Note, left must NOT be rooted by an addition or subtraction.

  Also, note that copies of left, not the original left object,
  are inserted into the result.
*/
ExpTree *distributeLeft(ExpTree *left, ExpTree *right);

/* Apply the left-distributive property of multiplication (*) w.r.t.
  addition (+) and subtraction (-). distributing each non-distributive
  subtree of the distributive left tree across the top-most, n-ary addition
  or subtraction that roots right, recursively.

  Note, right MUST be rooted by an addition or subtraction.
  Note, left MUST be rooted by an addition or subtraction.

  In other words, recurse down the subtrees of left to find
  all non-distributive subtrees, and distribute those subtrees
  across right.
*/
ExpTree *distributeLeftDistributive(ExpTree *left, ExpTree *right);



#endif
