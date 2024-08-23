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

/* The time variable, which is assumed to always exist in the context of TM
 * flowpipes. */
#define VAR_TIME "t"

/* Compute the Taylor polynomials approximating the given system of ODEs.

  system: The given system of ODEs
  order: The Taylor polynomial order
  k: The truncation order
*/
TaylorModel *computeTaylorPolynomial(ODEList *system, unsigned int order,
                                     unsigned int k);

/* Compute the Lie derivatives of the given functions
  and specified order, w.r.t. the given system of ODEs.

  system: The ODE system to use the vector field of
  functions: The initial set of functions to derive
  order: The order of the to compute Lie derivatives
*/
TaylorModel *lieDerivativeK(ODEList *system, TaylorModel *functions,
                            unsigned int order);
/* Compute the Lie derivative for a set of functions w.r.t. a vector field/ODE
  system.
*/
TaylorModel *lieDerivativeTaylorModel(ODEList *system, TaylorModel *functions);
/* Compute the Lie derivative of a function w.r.t. a vector field.  */
ExpTree *lieDerivative(ODEList *vectorField, ExpTree *function);

/* Apply the picard operator to the given functions w.r.t. the given ODEs.

  This approach strictly integrates w.r.t. the variable t and over interval [0,
  t].
*/
TaylorModel *picardOperator(ODEList *vectorField, TaylorModel *functions);
/* Apply the TM extension of the picard operator to the given Taylor models
  w.r.t. the given ODEs.

  This approach strictly integrates w.r.t. the variable t and over interval [0,
  t].
*/
TaylorModel *picardOperatorTM(ODEList *vectorField, TaylorModel *functions);
/* Substitute each of the functions, w.r.t. their assigned variable, into all
 * of the ODEs. */
TaylorModel *substituteTaylorModel(ODEList *system, TaylorModel *functions);

/* Construct the identity polynomial list.

  The identity polynomial matches each variable
  derivative to the variable itself, i.e. given ODEs
      x' = ...; y' = ...; z' = ...; ...
  generate the corresponding polynomials
      p(x) = x; p(y) = y; p(z) = z; ...
*/
TaylorModel *initTaylorModel(ODEList *system);

/* Compute a safe remainder interval, that contains the true solution
to the system of ODEs, for each Taylor polynomial.

  polynomial: The Taylor polynomial
*/
Interval computeSafeRemainder(TaylorModel *polynomials);

#endif
