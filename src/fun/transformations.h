/**
 * @file transformations.h
 * @author Thomas Gueutal (thomas.gueutal@student.uantwerpen.be)
 * @brief Functions to transform an expression tree to a desirable form.
 * @version 0.1
 * @date 2024-09-18
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H

#include "funexp.h"
#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief Simplify the expression through algebraic manipulations.
 * @details This is currently an alias for @ref simplifyOperators.
 *
 * @param[in] source The expression to simplify.
 * @return ExpTree* A newly heap-allocated, simplified expression tree.
 */
ExpTree *simplify(const ExpTree *source);

/**
 * @brief Convert the expression to a sum of products.
 * @details Distributes all multiplications in all levels of the tree.
 * The result is fully distributed, it is a sum of monomials.
 *
 * e.g. the following expressions and their distributions. \f{eqnarray*}{
 *          x * (y + z) & = & xy + xz \\
 *    (x + y) * (y - z) & = & xy - xz + yy - yz \\
 *    x * (y * (a + b)) & = & xya + xyb
 * \f}
 * @pre \p source may **not** be NULL.
 *
 * @param[in] source The expression to convert.
 * @return ExpTree* A newly heap-allocated expression tree; the input
 * expression as a sum of products.
 */
ExpTree *toSumOfProducts(const ExpTree *source);

/**
 * @brief Convert the expression to Horner form.
 * @details One expression may have multiple, algebraically equivalent Horner
 * forms. No guarantees are made about which one is selected internally.
 * @pre \p source may **not** be NULL.
 *
 * @param[in] source The expression to convert.
 * @return ExpTree* A newly heap-allocated expression tree; the input
 * expression as a Horner form.
 */
ExpTree *toHornerForm(const ExpTree *source);

/**
 * @brief Truncate all terms of degree i, where i > k.
 * @details e.g. the following symbolic truncation results. \f{eqnarray*}{
 *                 E & = & 1 + x + xy \\
 *    truncate(E, 1) & = & 1 + x
 * \f}
 * @pre \p source may **not** be NULL.
 *
 * @param[in] source The expression to truncate.
 * @param[in] k      The truncation degree k.
 * @return ExpTree* A newly heap-allocated expression tree; the truncated
 * input expression.
 */
ExpTree *truncate(const ExpTree *source, const unsigned int k);

/**
 * @brief Truncate all terms of degree i, where i > k.
 * @details Collects all truncated terms as a sum expression into
 * \p collectedTerms, in addition to truncating the input expression.
 *
 * e.g. the following symbolic truncation results. \f{eqnarray*}{
 *                    E  & = &  1 + x + xy - z^3 \\
 *    truncate(E, 1, C)  & = &  1 + x            \\
 *                    C  & = &  xy - z^3
 * \f}
 * @pre \p source may **not** be NULL.
 * @pre *\p collectedTerms (note the dereferencing) must be NULL,
 * but \p collectedTerms may **not** be NULL.
 * @post \p collectedTerms is NULL if no terms were truncated, else
 * it is not NULL.
 *
 * @param[in]  source         The expression to truncate.
 * @param[in]  k              The truncation degree k.
 * @param[out] collectedTerms The root of the sum expression that contains
 *                            all truncated terms.
 * @return ExpTree* A newly heap-allocated expression tree; the truncated
 * input expression.
 */
ExpTree *truncate2(const ExpTree *source, const unsigned int k,
                   ExpTree **collectedTerms);

/**
 * @brief Truncate all terms of degree i, where i > k.
 * @details This is an internal function and should normally not be called by
 * third-parties. Refer to @ref truncate and @ref truncate2 instead.
 *
 * Collects all truncated terms as a sum expression into
 * \p collectedTerms, in addition to truncating the input expression.
 * @pre \p source may **not** be NULL.
 * @pre Iff. \p collect is true, then *\p collectedTerms (note the
 * dereferencing) must be NULL, but \p collectedTerms may **not** be NULL.
 * Else there are no pre-conditions for \p collectedTerms.
 * @post \p collectedTerms is NULL if no terms were truncated or if \p collect
 * is false, else it is not NULL.
 *
 * @param[in] source          The expression to truncate.
 * @param[in] k               The truncation degree k.
 * @param[out] collectedTerms The root of the sum expression that contains all
 *                            truncated terms iff. collect is true.
 * @param[in] collect         If true, then collect the truncated terms, else
 *                            do not.
 * @return ExpTree* A newly heap-allocated expression tree; the truncated
 * input expression.
 */
ExpTree *truncateTerms(const ExpTree *source, const unsigned int k,
                       ExpTree **collectedTerms, const bool collect);

/**
 * @brief Substitute all variables with the given name in the source tree by
 * the target tree.
 * @details e.g. the following substitution. \f{eqnarray*}{
 *                      S  & = &  1 + 2x \\
 *                      T  & = &  1 + y  \\
 *    substitute(E, x, T)  & = &  1 + 2 * (1 + y)
 * \f}
 * @pre \p source, \p var and \p target may **not** be NULL.
 *
 * @param[in] source The expression to apply substitution to.
 * @param[in] var    The variable name to substitute in the source tree.
 * @param[in] target The expression to replace each occurrence of the
 *                   variable by.
 * @return ExpTree* A newly heap-allocated expression tree; the result of
 * applying substitution to the input expression.
 */
ExpTree *substitute(const ExpTree *source, const char *var,
                    const ExpTree *target);

/**
 * @brief Simplify a given expression by applying any found absorbing and
 * neutral elements to their operators.
 * @details This process is applied bottom-up. So new absorbing or neutral
 * elements that are created after applying deeper ones are pushed up to
 * the root of the expression. After this operation, no more absorbing or
 * neutral elements exist w.r.t. each operator node of the tree.
 *
 * e.g. the following applications. \f{eqnarray*}{
 *          0 + x  & = &  x     \\
 *          1 * x  & = &  x     \\
 *    0 * (x + y)  & = &  0     \\
 *    1 * (x + y)  & = &  x + y \\
 *    (0 * y) + x  & = &  x     \\
 * \f}
 *
 * @param[in] source The expression to simplify.
 * @return ExpTree* A newly heap-allocated expression tree; the result of
 * applying all neutral or absorbing elements in the input expression.
 */
ExpTree *simplifyOperators(const ExpTree *source);

/**
 * @brief Check if the given expression is a number leaf with data
 * equivalent to '0'.
 *
 * @return true  If the tree is a leaf matching the pattern.
 * @return false Else false.
 */
bool isZeroExpTree(const ExpTree *source);

/**
 * @brief Create a number leaf with data equivalent to '0'.
 *
 * @return ExpTree* A newly heap-allocated expression tree.
 */
ExpTree *newZeroExpTree(void);

/**
 * @brief Check if the given expression is a number leaf with
 * data equivalent to '1'.
 *
 * @return true  If the tree is a leaf matching the pattern.
 * @return false Else false.
 */
bool isOneExpTree(const ExpTree *source);

/**
 * @brief Create a number leaf with data equivalent to '1'.
 *
 * @return ExpTree* A newly heap-allocated expression tree.
 */
ExpTree *newOneExpTree(void);

/**
 * @brief Apply the left-distributive property of multiplication w.r.t.
 * addition and subtraction by distributing \p left across \p right.
 * @details Recursively distribute the non-distributive left operand across
 * the top-most, n-ary addition/subtraction that roots right.
 *
 * e.g. given left = x and right = (y + z), then distributing left over right
 * gives x * (y + z) = xy + xz. But left = (x + y) would not be allowed,
 * since left must be non-distributive.
 *
 * @pre Neither \p left nor \p right may be NULL
 * @pre \p right must be rooted by an addition or subtraction.
 * @pre \p left must **not** be rooted by an addition or subtraction, i.e.
 * left must be "non-distributive".
 *
 * @param[in] left  The expression to distribute.
 * @param[in] right The expression subject to distribution.
 * @return ExpTree* A newly heap-allocated expression tree; the result of
 * applying the distribution.
 */
ExpTree *distributeLeft(const ExpTree *left, const ExpTree *right);

/**
 * @brief Apply the left-distributive property of multiplication w.r.t.
 * addition and subtraction by distributing \p left across \p right.
 * @details Here left must be distributive, meaning it must be rooted by
 * an addition or by a subtraction. For each non-distributive subtree of
 * left, recursively distribute the subtree across the top-most, n-ary
 * addition/subtraction that roots right.
 *
 * In other words, recurse down the subtrees of left to find all
 * non-distributive subtrees, and distribute those subtrees across right.
 *
 * e.g. given left = (x + y) and right = (a - b), then distributing left over
 * right gives (x + y) * (a - b) = xa - xb + ya - yb. But left = x would not
 * be allowed, since left must be distributive.
 * @pre Neither \p left nor \p right may be NULL
 * @pre Both \p left and \p right must be rooted by an addition or subtraction.
 *
 * @param[in] left  The expression to distribute the subtrees of.
 * @param[in] right The expression subject to distribution.
 * @return ExpTree* A newly heap-allocated expression tree; the result of
 * applying the distribution.
 */
ExpTree *distributeLeftDistributive(const ExpTree *left, const ExpTree *right);

/**
 * @brief Distribute each encountered additive inverse (unary negation)
 * operator across its subtree, as deeply as possible.
 * @details The second, boolean parameter can be used to "add" an implicit
 * unary negation to the passed expression.
 *
 * e.g. suppose the given expression is "- - x". Then passing true for the
 * second argument, means input "- - x" is interpreted as "- - - x" and
 * returns "- x".
 * Passing false simply distributes all existing unary negative operators:
 * input "- - x" is interpreted as "- - x" and returns "x".
 * @pre \p source may **not** be NULL
 *
 * @param[in] source          The expression to distribute the negation
 *                            operators of.
 * @param[in] unevenNegsFound If true, then assume an uneven number of negation
 *                            operators have been seen until now.
 *                            Else assume an even number of negation operators
 *                            have been seen until now.
 * @return ExpTree* A newly heap-allocated expression tree; the result of
 * applying the distribution.
 */
ExpTree *distributeNeg(const ExpTree *source, const bool unevenNegsFound);

#endif
