#include "utils.h"

void dtoa(char *const destination, const size_t size, const double value) {
  snprintf(destination, size, "%.15g", value);
}

void utoa(char *const destination, const size_t size,
          const unsigned int value) {
  snprintf(destination, size, "%u", value);
}

unsigned int atou(const char *const source) {
  const int integer = atoi(source);
  /* Casting from a sign bit of '1' (-) would result in a changed value. */
  assert(integer >= 0);
  return (unsigned int)integer;
}
