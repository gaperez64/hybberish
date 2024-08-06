/* Functions to transform an expression to a desirable form. */

#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H

#include "funexp.h"
#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

/* Simplify the expression through algebraic manipulations. */
ExpTree *simplify(ExpTree *source);

/* Convert the expression to a sum of products. */
ExpTree *toSumOfProducts(ExpTree *source);

/* Convert the expression to Horner form. */
ExpTree *toHornerForm(ExpTree *source);

/* Truncate all terms of degree i, where i > k with k given. */
ExpTree *truncate(ExpTree *source, unsigned int k);

/* Substitute all variables with the given name inside a copy of the source tree
  by the target tree.

  Note, the lifetimes of all input instances are still the responsibility of the
  caller; copies are always used and or returned.
*/
ExpTree *substitute(ExpTree *source, char *var, ExpTree *target);

/*
    Simplification Helper methods.
*/

/* Simplify a given expression by applying any found absorbing and neutral
 * elements to their operators. */
ExpTree *simplifyOperators(ExpTree *source);

/* Check if the given expression is a number leaf with data equivalent to '0'.
 */
bool isZeroExpTree(ExpTree *source);

/* Create a number leaf with data equivalent to '0'. */
ExpTree *newZeroExpTree(void);

/* Check if the given expression is a number leaf with data equivalent to '1'.
 */
bool isOneExpTree(ExpTree *source);

/* Create a number leaf with data equivalent to '1'. */
ExpTree *newOneExpTree(void);

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

/* Distribute each encountered unary negative operator in the given expression
  tree.

  Suppose the given tree is "--x". Then passing true for the second argument,
  "adds" an implicit, new unary negative: input "--x" is interpreted as "---x"
  and returns "-x". Passing false simply distributes all existing unary negative
  operators: input "--x" is interpreted as "--x" and returns "x".

  unevenNegsFound: If true, then assume an uneven number of unary negative
  operators were encountered until now. In other words, distribute the unary
  negative. Else, assume an even number of unary negatives were encountered and
  consequently do NOT apply the operator. source: The tree to distribute over.
*/
ExpTree *distributeNeg(ExpTree *source, bool unevenNegsFound);

#endif
