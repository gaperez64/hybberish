/**
 * @file variables.h
 * @author Thomas Gueutal (thomas.gueutal@student.uantwerpen.be)
 * @brief An interface for constructing and manipulating the valuation
 * of a vector of variables.
 * @details In this context, valuation refers to two concepts:
 *  1. A real-valued variable Valuation that fixes the value of
 *     a variable to a single, real value.
 *  2. An interval-valued variable Domain, that describes the
 *     domain of values a variable can take on.
 * @version 0.1
 * @date 2024-09-17
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef VARIABLES_H
#define VARIABLES_H

#include "interval.h"
#include <stdlib.h>

/**
 * @brief A class representing a vector of variables, each restricted
 * to their own domain, as a linked list.
 * @details For example, given an ordered set of variables {x, y},
 * a valid interval domain vector is ([-1, 1], [-0.5, 0.5]) where x &isin; [-1,
 * 1] and y &isin; [-0.5, 0.5]. This corresponds to a Domain linked list of the
 * form ("x", [-1, 1])-->("y", [-0.5, 0.5]).
 *
 * @invariant The member @ref Domain.var may never be NULL. This is only
 * guaranteed at construction if the correct constructor method,
 * @ref newDomain, is used.
 */
typedef struct Domain {
  /// @brief The name of the variable whose domain this is.
  char *var;
  /// @brief The interval domain of the variable. e.g. "x" in [-1, 1].
  Interval domain;
  /// @brief The next component of the Domain vector.
  struct Domain *next;
} Domain;

/**
 * @brief Create a new, single element linked list.
 * @pre \p var must **not** be NULL and heap-allocated.
 * @post Transfers ownership of all heap-allocated arguments to the
 * newly created Domain instance.
 *
 * @param[in] var    The name of the variable to construct the domain of.
 * @param[in] domain The interval domain of the variable.
 * @return Domain* A heap-allocated domain instance.
 */
Domain *newDomain(char *var, const Interval domain);

/**
 * @brief Attach the second element as the head of the first list.
 * @pre The \p head argument must be a single element list (not NULL).
 * @pre The \p tail argument must be NULL or heap-allocated.
 * @post Transfers ownership of \p tail to \p head. Ownership of
 * \p head is still the caller's.
 *
 * @param[in] tail The list to prepend a new head to.
 * @param[in] head The head to prepend to \p tail.
 * @return Domain* The \p head pointer, with the tail attached.
 */
Domain *appDomainElem(Domain *tail, Domain *head);

/**
 * @brief Allocate a new Domain element and assign the passed member values.
 * @details The newly created element will be prepended to the passed tail
 * list if it exists. Else, the new element will have an empty tail, making
 * the result a single element list. e.g. given \p tail ("y", [-0.5, 0.5]),
 * \p var "x" and \p domain [-1, 1], the result is
 * ("x", [-1, 1])-->("y", [-0.5, 0.5]).
 * @pre All pointer arguments must be heap-allocated or NULL.
 * @post Transfers ownership of all heap-allocated arguments to the
 * newly created Domain instance.
 *
 * @param[in] tail   The tail to prepend the newly created domain element to.
 * @param[in] var    The variable name of the new domain element.
 * @param[in] domain The interval domain of the new domain element.
 * @return Domain* A heap-allocated, new domain element with the given tail
 * attached.
 */
Domain *newDomainElem(Domain *tail, char *var, const Interval domain);

/**
 * @brief Deallocate the given list.
 * @pre The given list must not be NULL.
 */
void delDomain(Domain *list);

/**
 * @brief Print a representation of the given list to the specified stream.
 * @pre Neither the given \p list nor stream ( \p where ) may both be NULL.
 *
 * @param[in] list  The list to print.
 * @param[in] where The stream (destination) to print to.
 */
void printDomain(const Domain *list, FILE *where);

/**
 * @brief A class representing a vector of variables,
 * each with an exact, real-valued valuation, as a linked list.
 * @details For example, given an ordered set of variables {x, y},
 * a valid real valuation vector is (1, 0.5) where x = 1 and y = 0.5.
 * This corresponds to a Valuation linked list of the form
 * ("x", 1)-->("y", 0.5).
 *
 * @invariant The member @ref Valuation.var may never be NULL. This is only
 * guaranteed at construction if the correct constructor method,
 * @ref newValuation, is used.
 */
typedef struct Valuation {
  /// @brief The name of the variable whose valuation this is.
  char *var;
  /// @brief The real-valued valuation of the variable. e.g. "x" = 1.
  double val;
  /// @brief The next component of the Valuation vector.
  struct Valuation *next;
} Valuation;

/**
 * @brief Create a new, single element linked list.
 * @pre \p var must **not** be NULL and heap-allocated.
 * @post Transfers ownership of all heap-allocated arguments to the
 * newly created Valuation instance.
 *
 * @param[in] var The name of the variable to construct the valuation of.
 * @param[in] val The real-valued valuation of the variable.
 * @return Valuation* A heap-allocated valuation instance.
 */
Valuation *newValuation(char *var, const double val);

/**
 * @brief Attach the second element as the head of the first list.
 * @pre The \p head argument must be a single element list (not NULL).
 * @pre The \p tail argument must be NULL or heap-allocated.
 * @post Transfers ownership of \p tail to \p head. Ownership of
 * \p head is still the caller's.
 *
 * @param[in] tail The list to prepend a new head to.
 * @param[in] head The head to prepend to \p tail.
 * @return Valuation* The \p head pointer, with the tail attached.
 */
Valuation *appValuationElem(Valuation *tail, Valuation *head);

/**
 * @brief Allocate a new Valuation element and assign the passed member values.
 * @details The newly created element will be prepended to the passed tail
 * list if it exists. Else, the new element will have an empty tail, making
 * the result a single element list. e.g. given \p tail ("y", 0.5), \p var "x"
 * and \p val 1, the result is ("x", 1)-->("y", 0.5).
 * @pre All pointer arguments must be heap-allocated or NULL.
 * @post Transfers ownership of all heap-allocated arguments to the
 * newly created Valuation instance.
 *
 * @param[in] tail The tail to prepend the newly created valuation element to.
 * @param[in] var  The variable name of the new valuation element.
 * @param[in] val  The value of the new valuation element.
 * @return Valuation* A heap-allocated, new valuation element with the given
 * tail attached.
 */
Valuation *newValuationElem(Valuation *tail, char *var, const double val);

/**
 * @brief Deallocate the given list.
 * @pre The given list must not be NULL.
 */
void delValuation(Valuation *list);

/**
 * @brief Print a representation of the given list to the specified stream.
 * @pre Neither the given \p list nor stream ( \p where ) may both be NULL.
 *
 * @param[in] list  The list to print.
 * @param[in] where The stream (destination) to print to.
 */
void printValuation(const Valuation *list, FILE *where);

#endif
