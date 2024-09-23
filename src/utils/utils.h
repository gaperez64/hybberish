/**
 * @file utils.h
 * @author Thomas Gueutal (thomas.gueutal@student.uantwerpen.be)
 * @brief Utility functions that do not particularly belong anywhere.
 * @version 0.1
 * @date 2024-09-19
 *
 * @copyright Copyright (c) 2024
 *
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Convert a double value to a string representation, e.g. 1.5 to "1.5".
 *
 * @param[out] destination The destination buffer.
 * @param[in]  size        The size of the destination buffer.
 * @param[in]  value       The value to convert.
 */
void dtoa(char *const destination, const size_t size, const double value);

/**
 * @brief Convert a double value to a string representation, e.g. 1 to "1".
 *
 * @param[out] destination The destination buffer.
 * @param[in]  size        The size of the destination buffer.
 * @param[in]  value       The value to convert.
 */
void utoa(char *const destination, const size_t size, const unsigned int value);

/**
 * @brief Convert a string representation of an unsigned int, e.g. "1" to 1.
 * @details A signed integer string, such as "-1", will fail to convert.
 *
 * @param[in] source The char array to convert.
 * @return unsigned int
 */
unsigned int atou(const char *const source);
