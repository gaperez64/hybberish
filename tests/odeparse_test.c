#include "odeparse.h"
#include "sysode.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
  /* to avoid silly warnings about unused parameters */
  (void)argc;
  (void)argv;

  {
    const char str[] =
        "x' = -b; y' = sqrt((b^2) - 4 * a * c); last' = (2 * at);";
    ODEList *list;
    int res = parseOdeString(str, &list);
    assert(res == 0);

    /* printing */
    char buffer[100];
    FILE *stream = fmemopen(buffer, 100, "w");
    assert(stream != NULL);
    const char msg[] =
        "last' = (2 * at); y' = sqrt(((b^2) - ((4 * a) * c))); x' = -b; ";
    printOdeList(list, stream);
    fclose(stream); /* close to flush and write null byte */
    printf("expect: |%s| = %lu\n", msg, strlen(msg));
    fflush(stdout);
    printf("got: |%s| = %lu\n", buffer, strlen(buffer));
    fflush(stdout);
    printf("!strcmp = %i\n", !strcmp(buffer, msg));
    assert(!strcmp(buffer, msg));

    /* clean */
    delOdeList(list);
  }
  {
    const char str[] =
        "x' = -b; y' = sqrt((b^2) - 4 * a * c); last' = (2 * at);";
    ODEList *list;
    int res = parseOdeString(str, &list);
    assert(res == 0);

    /* printing */
    char buffer[100];
    FILE *stream = fmemopen(buffer, 100, "w");
    assert(stream != NULL);
    const char msg[] =
        "last' = (2 * at); y' = sqrt(((b^2) - ((4 * a) * c))); x' = -b; ";
    printOdeList(list, stream);
    fclose(stream); /* close to flush and write null byte */
    printf("expect: |%s| = %lu\n", msg, strlen(msg));
    fflush(stdout);
    printf("got: |%s| = %lu\n", buffer, strlen(buffer));
    fflush(stdout);
    printf("!strcmp = %i\n", !strcmp(buffer, msg));
    assert(!strcmp(buffer, msg));

    /* clean */
    delOdeList(list);
  }
  {
    const char str[] =
        "x' = -b; y' = sqrt((b^2) - 4 * a * c); last' = (2 * at);";
    ODEList *list;
    int res = parseOdeString(str, &list);
    assert(res == 0);

    /* printing */
    char buffer[100];
    FILE *stream = fmemopen(buffer, 100, "w");
    assert(stream != NULL);
    const char msg[] =
        "last' = (2 * at); y' = sqrt(((b^2) - ((4 * a) * c))); x' = -b; ";
    printOdeList(list, stream);
    fclose(stream); /* close to flush and write null byte */
    printf("expect: |%s| = %lu\n", msg, strlen(msg));
    fflush(stdout);
    printf("got: |%s| = %lu\n", buffer, strlen(buffer));
    fflush(stdout);
    printf("!strcmp = %i\n", !strcmp(buffer, msg));
    assert(!strcmp(buffer, msg));

    /* clean */
    delOdeList(list);
  }
  {
    const char str[] =
        "x' = -b; y' = sqrt((b^222) - 14 * acc * c); last' = (2 * at);";
    ODEList *list;
    int res = parseOdeString(str, &list);
    assert(res == 0);

    /* printing */
    char buffer[100];
    FILE *stream = fmemopen(buffer, 100, "w");
    assert(stream != NULL);
    const char msg[] =
        "last' = (2 * at); y' = sqrt(((b^2) - ((4 * a) * c))); x' = -b; ";
    printOdeList(list, stream);
    fclose(stream); /* close to flush and write null byte */
    printf("expect: |%s| = %lu\n", msg, strlen(msg));
    fflush(stdout);
    printf("got: |%s| = %lu\n", buffer, strlen(buffer));
    fflush(stdout);
    printf("!strcmp = %i\n", !strcmp(buffer, msg));
    // assert(!strcmp(buffer, msg));

    /* clean */
    delOdeList(list);
  }
  return 0;
}
