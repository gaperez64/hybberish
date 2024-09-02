#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Convert a double value to a string representation, e.g. 1.5 to "1.5".

  The destination buffer size must be specified.
*/
void dtoa(char *const destination, const size_t size, const double value);

/* Convert a double value to a string representation, e.g. 1 to "1".

  The destination buffer size must be specified.
*/
void utoa(char *const destination, const size_t size, const unsigned int value);

/* Convert a string representation of an unsigned int, e.g. "1" to 1.

  A signed integer string, such as "-1", will fail to convert.
*/
unsigned int atou(const char *const source);
