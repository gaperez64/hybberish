/**
 * @file tmflowpipe.h
 * @author Thomas Gueutal (thomas.gueutal@student.uantwerpen.be)
 * @brief An implementation of Taylor model integration, for generating
 * flowpipe over-approximations given a system of ODEs.
 * @version 0.1
 * @date 2024-09-17
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef TM_FLOWPIPE_H
#define TM_FLOWPIPE_H

#include "funexp.h"
#include "interval.h"
#include "sysode.h"
#include "taylormodel.h"
#include "transformations.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

/// @brief The time variable; each ODE variable x is defined w.r.t. time as
///        \f$ \dot{x} \f$ = \f$ \frac{dx}{dt} \f$.
#define VAR_TIME "t"

/**
 * @brief Step 1 of TM integration: compute the vector of Taylor polynomials.
 * @details The Taylor polynomial \f$ p_l(\vec{x}_l, t) \f$ approximates
 * the true flow \f$ x(t) = \phi(\vec{x}_l, t) \f$ of the given ODEs in
 * the current TM integration iteration.
 * @pre The ODEs \p system may **not** be NULL.
 * @pre It must hold that \p order <= \p k. Else, the truncation would
 * simply strip away all polynomial terms of order > \p k.
 *
 * @param[in] system The system of ODEs whose true flows to over-approximate.
 * @param[in] order  The order of the generated Taylor polynomials with which
 *                   to approximate the true flow.
 * @param[in] k      The truncation order applied during TM arithmetic.
 * @return TaylorModel* A newly heap-allocated vector of Taylor models with as
 * polynomial parts the Taylor polynomials approximating the true flows of each
 * of the ODEs system's components.
 */
TaylorModel *computeTaylorPolynomial(ODEList *system, unsigned int order,
                                     unsigned int k);

/**
 * @brief Compute a **vector** of order k Lie derivatives.
 * @details By definition, higher order Lie derivatives can be computed
 * recursively:
 * \f$ L_f^{m+1}(g) = L_f(L_f^m(g)) \f$ for \f$ m = 1, 2, ..., k-1 \f$.
 *
 * Refer to @ref lieDerivativeTaylorModel for the semantics and implementation
 * of the first-order Lie derivative, and for the reasoning behind representing
 * the \p functions as a vector of Taylor models.
 * @pre The ODEs \p system and \p functions may **not** be NULL.
 *
 * @param[in] system    The system of ODEs to use in Lie derivation.
 * @param[in] functions The vector of functions to derive.
 * @param[in] order     k, the order of Lie derivative to compute.
 * @return TaylorModel* A newly heap-allocated vector of Taylor models with as
 * polynomial parts the order k Lie derivatives,
 * \f$ TM_i' = (L_f^k(p_i), I_i') \f$.
 */
TaylorModel *lieDerivativeK(ODEList *system, TaylorModel *functions,
                            unsigned int order);

/**
 * @brief Compute a **vector** of first-order Lie derivatives.
 * @details Here \p functions is a vector of Taylor models
 * \f$ TM_i = (p_i, I_i) \f$. The polynomial parts \f$ p_i \f$ of these models
 * represent the functions that the Lie derivative receives as input.
 * For each individual function \f$ p_i \f$ we compute the first-order Lie
 * derivative of that function, \f$ L_f(p_i) \f$, where \f$ f \f$ refers to
 * the entire ODEs \p system.
 * @pre The ODEs \p system and \p functions may **not** be NULL.
 *
 * @param[in] system    The system of ODEs to use in Lie derivation.
 * @param[in] functions The vector of functions to derive, expressed
 *                      as a vector of Taylor models.
 * @return TaylorModel* A newly heap-allocated vector of Taylor models with as
 * polynomial parts the first-order Lie derivatives,
 * \f$ TM_i' = (L_f(p_i), I_i') \f$.
 */
TaylorModel *lieDerivativeTaylorModel(ODEList *system, TaylorModel *functions);

/**
 * @brief Compute a single, first-order Lie derivative expression of a single
 * function expression w.r.t. an m-dimensional vector field.
 * @pre The ODEs \p vectorField and \p function may **not** be NULL.
 *
 * @param[in] vectorField The m-dimensional vector field f.
 * @param[in] function    The function expression g.
 * @return ExpTree* A newly heap-allocated, first-order Lie derivative
 * expression, \f$ L_f(g) \f$.
 */
ExpTree *lieDerivative(ODEList *vectorField, ExpTree *function);

/**
 * @brief Compute the picard operator for a **vector** of functions.
 * @details Here \p functions is a vector of Taylor models
 * \f$ TM_i = (p_i, I_i) \f$. The polynomial parts \f$ p_i \f$ of these models
 * represent the functions that the picard operator receives as input.
 * For each individual function \f$ p_i \f$ we compute the picard operator
 * result of that function, \f$ P_f(p_i) \f$, where \f$ f \f$ refers to
 * the entire ODEs \p vectorField. It computes
 * \f$ g_{j+1}(t) = P_f(g_j) = \vec{x} + \int_0^t f(g_j(t), t)dt \f$ for
 * \f$ j \ge 0 \f$ where \f$ g_j \f$ is the j-th order picard operator result.
 * Here j does *not* refer to the j-th component of any vector.
 * This approach strictly integrates w.r.t. dt and over interval [0, t].
 * @pre The ODEs \p vectorField and \p functions may **not** be NULL.
 * @post The remainder interval part of the output is indeterminate; this
 * implementation of the picard operator focuses purely on the polynomial
 * part, while the remainder part can take any value without consistency.
 *
 * @param[in] vectorField The m-dimensional vector field f.
 * @param[in] functions   The functions to apply the picard operator to.
 * @return TaylorModel* A newly heap-allocated vector of Taylor models with as
 * polynomial parts the picard operator results,
 * \f$ TM_i' = (P_f(p_i), I_i') \f$.
 */
TaylorModel *picardOperator(ODEList *vectorField, TaylorModel *functions);

/**
 * @brief Compute the TM extension of the picard operator for a **vector**
 * of functions.
 * @details Here \p functions is a vector of Taylor models
 * \f$ TM_i = (p_i, I_i) \f$.
 * Refer to @ref picardOperator for the semantics of computing the
 * polynomial part of the picard operator via algebraic operations.
 * The TM extension implies that the algebraic operations of the basic
 * picard operator are replaced by the Taylor model arithmetic equivalents.
 * So the TM extension of the picard operator must be evaluated via order k
 * TM arithmetic.
 * @pre The ODEs \p vectorField and \p functions may **not** be NULL.
 *
 * @param[in] vectorField The m-dimensional vector field f.
 * @param[in] functions   The Taylor models to apply the TM extension of the
 *                    picard operator to.
 * @return TaylorModel* A newly heap-allocated vector of Taylor models with as
 * polynomial parts the picard operator results and remainder parts computed
 * during order k TM arithmetic, \f$ TM_i' = (P_f(p_i), I_i') \f$.
 */
TaylorModel *picardOperatorTM(ODEList *vectorField, TaylorModel *functions);

/**
 * @brief Substitute all of the variables in each ODE by the corresponding
 * functions. This operation is applied elementwise to the ODEs.
 * @details Here \p functions is a vector of Taylor models
 * \f$ TM_i = (p_i, I_i) \f$. We denote the vector field \p system as follows
 * \f$ \dot{x} = f(\vec{x}, t) \f$.
 * Each of the functions \f$ g_j \f$ corresponds to a distinct variable
 * \f$ \vec{x}_j \f$ in the vector field f. So each substitution is
 * deterministic in function of the vector of functions \f$ g \f$.
 * @post The remainder interval part of the output is indeterminate; the
 * substitution procedure focuses purely on the polynomial part, while the
 * remainder part can take any value without consistency.
 * @pre The ODEs \p system and \p functions may **not** be NULL.
 *
 * @param[in] system    The m-dimensional vector field f.
 * @param[in] functions The functions g to subtitute into the vector field.
 * @return TaylorModel* A newly heap-allocated vector of Taylor models with as
 * polynomial parts the result of substituting the functions into the ODEs,
 * \f$ TM_i' = (f(g(\vec{x}, t), t), I_i') \f$.
 */
TaylorModel *substituteTaylorModel(ODEList *system, TaylorModel *functions);

/**
 * @brief Construct the identity polynomial list.
 * @details The identity polynomial matches each of the ODEs' variables'
 * derivatives \f$ \dot{x}_i \f$ to the variable itself, \f$ x_i \f$.
 *
 * i.e. given ODEs \f$
 *    x' =  ..., \quad
 *    y' =  ..., \quad
 *    z' =  ..., \quad ...
 * \f$
 * generate the corresponding polynomials \f$
 *    p_x(\vec{x}, t) = x, \quad
 *    p_y(\vec{x}, t) = y, \quad
 *    p_z(\vec{x}, t) = z, \quad ...
 * \f$.
 * @pre The ODEs \p system may **not** be NULL.
 *
 * @param[in] system The m-dimensional vector field f to extract the variables.
 * @return TaylorModel* A newly heap-allocated vector of Taylor models with as
 * polynomial parts identity polynomials, and as remainder the degenerate zero
 * interval, \f$ TM_i' = (x_i, [0, 0]) \f$.
 */
TaylorModel *initTaylorModel(ODEList *system);

/**
 * @brief Step 2 of TM integration: Compute a safe remainder interval, that
 * contains the true solutions to the system of ODEs, for each Taylor
 * polynomial approximation of the true solutions.
 * @pre The \p polynomials may **not** be NULL.
 *
 * @param[in] polynomials The polynomials generated in step 1 of
 *                        TM integration.
 * @return Interval The safe remainder intervals for the given polynomials.
 */
Interval computeSafeRemainder(TaylorModel *polynomials);

#endif
