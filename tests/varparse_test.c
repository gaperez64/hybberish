#include "variables.h"
#include "varparse.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
  /* to avoid silly warnings about unused parameters */
  (void)argc;
  (void)argv;

  /* Test the parsing of a string representing the domain values
    that a series of variables is restricted to.
    The result should be a Domain list that matches the input. */
  const char str[] =
      "x in [1, 1]; y in [ -1.4 , +3.1 ]; last in [-3.444, -3.333];";
  Domain *list;
  int res = parseVarString(str, &list);
  printDomain(list, stdout);
  fflush(stdout);
  printf("\n");
  assert(res == 0);

  /* printing */
  char buffer[100];
  FILE *stream = fmemopen(buffer, 100, "w");
  assert(stream != NULL);
  const char msg[] = "last in [-3.444000, -3.333000]; y in [-1.400000, "
                     "3.100000]; x in [1.000000, 1.000000]; ";
  printDomain(list, stream);
  fclose(stream); /* close to flush and write null byte */
  printf("\n");
  printf("expect: |%s| = %lu\n", msg, strlen(msg));
  printf("actual:    |%s| = %lu\n", buffer, strlen(buffer));
  printf("strcmp = %i\n", strcmp(buffer, msg));
  fflush(stdout);
  assert(strcmp(buffer, msg) == 0);

  /* clean */
  delDomain(list);
  return 0;
}
