/**
 * @file odeparse.h
 * @author Guillermo A. Perez (guillermo.perez@uantwerpen.be)
 * @brief A bison+ninja parser to read a multi-dimensional system of ODEs.
 * @version 0.1
 * @date 2023-04-19
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef ODEPARSE_H
#define ODEPARSE_H

#include "sysode.h"

/**
 * @brief Parse a string representation of a system of ODEs into a
 * structured form.
 * @details Parse the first parameter, containing a string representation
 * of a system of ODEs, and produce the second parameter, a structured
 * representation of the input, as output.
 *
 * @return int A return code.
 */
int parseOdeString(const char *, ODEList **);
void setOdeInputString(const char *);
void endOdeScan(void);

#endif
