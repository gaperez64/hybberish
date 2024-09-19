/**
 * @file sysode.h
 * @author Guillermo A. Perez (guillermo.perez@uantwerpen.be)
 * @brief An interface for constructing and manipulating a system of ODEs.
 * @version 0.1
 * @date 2024-09-19
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef SYSODE_H
#define SYSODE_H

#include <stdio.h>

#include "funexp.h"

/* Lists of ODEs */

/**
 * @brief A system of ODEs as a linked list.
 * @details An ODE takes the form \f$ \dot{x} = f(\vec{x}, t) \f$ where
 * \f$ \dot{x} \f$. is the derivative of all ODE variables w.r.t. time,
 * \f$ f \f$ is the system's vector field and \f$ \vec{x} \f$ is a vector
 * the variables of the ODE.
 *
 * This can be expressed differently, in function of the separate
 * components of the m-dimensional ODE:
 * \f{eqnarray*}{
 *    \dot{x}_1  & = &  f_1(\vec{x}, t) \\
 *    \dot{x}_2  & = &  f_2(\vec{x}, t) \\
 *    ...                          \\
 *    \dot{x}_m  & = &  f_m(\vec{x}, t)
 * \f}
 */
typedef struct ODEList {
  /// The variable name of this ODE component.
  char *fun;
  /// The vector field of this ODE component.
  ExpTree *exp;
  /// @brief The next component of the ODEList.
  struct ODEList *next;
} ODEList;

/**
 * @brief Create a new, single element list.
 * @pre Both \p fun and \p exp may **not** be NULL and must be heap-allocated.
 * @post Transfers ownership of all heap-allocated arguments to the
 * newly created ODEList instance.
 *
 * @param fun The variable name of the ODE component.
 * @param exp The vector field of the ODE component.
 * @return ODEList* A heap-allocated ODE instance.
 */
ODEList *newOdeList(char *fun, ExpTree *exp);

/* Append the second one to the first list,
 * It assumes the second one is a single node!
 */

/**
 * @brief Attach the second element as the head of the first list.
 * @pre The \p head argument must be a single element list (not NULL).
 * @pre The \p tail argument must be NULL or heap-allocated.
 * @post Transfers ownership of \p tail to \p head. Ownership of
 * \p head is still the caller's.
 *
 * @param[in] tail The list to prepend a new head to.
 * @param[in] head The head to prepend to \p tail.
 * @return ODEList* The \p head pointer, with the tail attached.
 */
ODEList *appOdeElem(ODEList *tail, ODEList *head);

/**
 * @brief Allocate a new ODEList element and assign the passed member values.
 * @details The newly created element will be prepended to the passed tail
 * list if it exists. Else, the new element will have an empty tail, making
 * the result a single element list.
 *
 * @param tail The tail to prepend the new element to.
 * @param fun  The ODE variable to which the new element corresponds.
 * @param exp  The vector field of the new element.
 * @return ODEList* A heap-allocated, new ODEList element with the given tail.
 */
ODEList *newOdeElem(ODEList *tail, char *fun, ExpTree *exp);

/**
 * @brief Deallocate the given list.
 * @pre The given list must not be NULL.
 */
void delOdeList(ODEList *list);

/**
 * @brief Print a representation of the given list to the specified stream.
 * @pre Neither the given \p list nor stream ( \p where ) may be NULL.
 *
 * @param[in] list  The list to print.
 * @param[in] where The stream (destination) to print to.
 */
void printOdeList(ODEList *list, FILE *where);

#endif
