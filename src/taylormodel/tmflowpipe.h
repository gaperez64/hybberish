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
