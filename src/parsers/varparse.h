/**
 * @file varparse.h
 * @author Thomas Gueutal (thomas.gueutal@student.uantwerpen.be)
 * @brief A bison+ninja parser to read the domain specification of a
 * multi-dimensional vector of variables.
 * @version 0.1
 * @date 2023-04-19
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef VARPARSE_H
#define VARPARSE_H

#include "variables.h"

/**
 * @brief Parse a string representation of the domains of a vector of variables
 * into a structured form.
 * @details Parse the first parameter, containing a string representation
 * of the domains that a vector of variables is restricted to, and produce
 * the second parameter, a structured representation of the input, as output.
 *
 * @return int A return code.
 */
int parseVarString(const char *, Domain **);
void setVarInputString(const char *);
void endVarScan(void);

#endif
