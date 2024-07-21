#include "interval.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>


bool approx(float got, float expected, float epsilon) {
  return fabs(got - expected) < epsilon;
}


/* Test if got = [a, b] = [expectedLeft, expectedRight]
  up to precision epsilon. and print out the results. */
void testInterval(Interval *got, float expectedLeft, float expectedRight, float epsilon) {
  assert(got != NULL);

  /* printing and comparison */
  bool compare = approx(got->left, expectedLeft, epsilon) &&
                 approx(got->right, expectedRight, epsilon);

  printf("Expect: [%f, %f]\n", expectedLeft, expectedRight);
  printf("Got:    [%f, %f]\n", got->left, got->right);
  printf("Compare = %i\n", compare);
  fflush(stdout); // Flush stdout
  assert(compare);
}


/* Test if got = expected, and print out the results. */
void testBool(bool got, bool expected) {
  /* printing and comparison */
  printf("Expect: %i\n", expected);
  printf("Got:    %i\n", got);
  fflush(stdout); // Flush stdout
  assert(got == expected);
}


/* Test if got = expected up to precision epsilon, and print out the results. */
void testFloat(float got, float expected, float epsilon) {
  /* printing and comparison */
  printf("Expect: %f\n", expected);
  printf("Got:    %f\n", got);
  fflush(stdout); // Flush stdout
  assert(fabs(expected - got) < epsilon);
}


/* Test if got is printed to the expected string, and print out the results. */
void testString(Interval *got, const char *expectedMsg) {
  /* printing and comparison */
  char *buffer = NULL;
  size_t buflen = 0;
  FILE *stream = open_memstream(&buffer, &buflen);
  assert(stream != NULL);
  printInterval(got, stream);
  fclose(stream);

  printf("Expect: \"%s\"\n", expectedMsg);
  printf("Got:    \"%s\"\n", buffer);
  printf("!strcmp = %i\n", !strcmp(buffer, expectedMsg));
  assert(!strcmp(buffer, expectedMsg));
  fflush(stdout); // Flush stdout

  /* clean */
  free(buffer);
}



int main(int argc, char *argv[]) {
  /* to avoid silly warnings about unused parameters */
  (void)argc;
  (void)argv;

  /* Setup */
  float eps = 0.0001;
  Interval iNeg   = newInterval(-2, -1);
  Interval iOrig  = newInterval(-1, 1);
  Interval iPos   = newInterval(1, 2);
  Interval iDegen = newInterval(12, 12);

  /* Enforce relationships within the test data,
    which are exploited in the testing. */
  assert(iNeg.right <= iOrig.left);
  assert(iOrig.right <= iPos.left);
  assert(iDegen.left == iDegen.right);


  /* Test binary addition. */
  {
    Interval res = addInterval(&iNeg, &iNeg);
    testInterval(&res, -4, -2, eps);

    res = addInterval(&iOrig, &iOrig);
    testInterval(&res, -2, 2, eps);

    res = addInterval(&iPos, &iPos);
    testInterval(&res, 2, 4, eps);

    res = addInterval(&iDegen, &iOrig);
    testInterval(&res, 11, 13, eps);
  }

  /* Test binary subtraction. */
  {
    Interval res = subInterval(&iNeg, &iNeg);
    testInterval(&res, -1, 1, eps);

    res = subInterval(&iOrig, &iOrig);
    testInterval(&res, -2, 2, eps);

    res = subInterval(&iPos, &iPos);
    testInterval(&res, -1, 1, eps);

    res = subInterval(&iDegen, &iNeg);
    testInterval(&res, 13, 14, eps);
  }

  /* Test binary multiplication. */
  {
    Interval res = mulInterval(&iNeg, &iNeg);
    testInterval(&res, 1, 4, eps);

    res = mulInterval(&iOrig, &iOrig);
    testInterval(&res, -1, 1, eps);

    res = mulInterval(&iPos, &iPos);
    testInterval(&res, 1, 4, eps);

    res = mulInterval(&iDegen, &iPos);
    testInterval(&res, 12, 24, eps);
  }

  /* Test binary division. */
  {
    Interval res = divInterval(&iNeg, &iNeg);
    testInterval(&res, 0.5, 2, eps);

    res = divInterval(&iOrig, &iOrig);
    testInterval(&res, -1, 1, eps);

    res = divInterval(&iPos, &iPos);
    testInterval(&res, 0.5, 2, eps);

    res = divInterval(&iDegen, &iNeg);
    testInterval(&res, -12, -6, eps);
  }

  /* Test binary interval equality. */
  {
    testBool(eqInterval(&iNeg, &iNeg, eps),   true);
    testBool(eqInterval(&iOrig, &iOrig, eps), true);
    testBool(eqInterval(&iPos, &iPos, eps),   true);

    /* Test inequality in a transitive fashion,
      instead of enumerating all possible options. */
    testBool(eqInterval(&iNeg, &iOrig, eps), false);
    testBool(eqInterval(&iOrig, &iPos, eps), false);
    testBool(eqInterval(&iPos, &iNeg, eps),  false);
  }

  /* Test binary interval membership. */
  {
    /* All intervals are closed, so should be contained in themselves. */
    testBool(inInterval(&iNeg, &iNeg),   true);
    testBool(inInterval(&iOrig, &iOrig), true);
    testBool(inInterval(&iPos, &iPos),   true);

    /*
        iNeg  = [a, b]
        iOrig = [c, d]
        iPos  = [e, f]
        encompassing = [a, f]
        where a <= b,c,d,e <= f
    */
    Interval encompassing = newInterval(iNeg.left, iPos.right);
    testBool(inInterval(&iNeg, &encompassing),  true);
    testBool(inInterval(&iOrig, &encompassing), true);
    testBool(inInterval(&iPos, &encompassing),  true);
    testBool(inInterval(&encompassing, &iNeg),  false);
    testBool(inInterval(&encompassing, &iOrig), false);
    testBool(inInterval(&encompassing, &iPos),  false);

    /*
        iNeg  = [a, b]
        leftEncompassing = [d, c]
        rightEncompassing = [c, e]
        where c in [a, b] but c != a & c != b
        so that d < a < c < b < e
    */
    float offset = 1.;
    Interval leftEncompassing = newInterval(
        iNeg.left - offset,
        (iNeg.left + iNeg.right) / 2.);
    Interval rightEncompassing = newInterval(
        (iNeg.left + iNeg.right) / 2.,
        iNeg.right + offset);
    testBool(inInterval(&iNeg, &leftEncompassing),  false);
    testBool(inInterval(&iNeg, &rightEncompassing), false);
  }

  /* Test interval properties. */
  {
    /* Width(I) */
    testFloat(intervalWidth(&iNeg),   (iNeg.right - iNeg.left)    , eps);
    testFloat(intervalWidth(&iOrig),  (iOrig.right - iOrig.left)  , eps);
    testFloat(intervalWidth(&iPos),   (iPos.right - iPos.left)    , eps);
    testFloat(intervalWidth(&iDegen), (iDegen.right - iDegen.left), eps);
    testFloat(intervalWidth(&iDegen), 0.                          , eps);

    /* Mid(I) */
    testFloat(intervalMidpoint(&iNeg),   (iNeg.left + iNeg.right) / 2.    , eps);
    testFloat(intervalMidpoint(&iOrig),  (iOrig.left + iOrig.right) / 2.  , eps);
    testFloat(intervalMidpoint(&iPos) ,  (iPos.left + iPos.right) / 2.    , eps);
    testFloat(intervalMidpoint(&iDegen), (iDegen.left + iDegen.right) / 2., eps);
    testFloat(intervalMidpoint(&iDegen), iDegen.left                      , eps);

    /* Mag(I) */
    testFloat(intervalMagnitude(&iNeg),   fmax(fabs(iNeg.left), fabs(iNeg.right))    , eps);
    testFloat(intervalMagnitude(&iOrig),  fmax(fabs(iOrig.left), fabs(iOrig.right))  , eps);
    testFloat(intervalMagnitude(&iPos),   fmax(fabs(iPos.left), fabs(iPos.right))    , eps);
    testFloat(intervalMagnitude(&iDegen), fmax(fabs(iDegen.left), fabs(iDegen.right)), eps);
    testFloat(intervalMagnitude(&iDegen), fabs(iDegen.left)                          , eps);

    /* Test interval degeneracy: I = [a, a] */
    testBool(intervalIsDegenerate(&iNeg, eps),   false);
    testBool(intervalIsDegenerate(&iOrig, eps),  false);
    testBool(intervalIsDegenerate(&iPos, eps),   false);
    testBool(intervalIsDegenerate(&iDegen, eps), true);
  }


  /* Test printing. */
  {
    testString(&iNeg,   strdup("[-2.000000, -1.000000]"));
    testString(&iOrig,  strdup("[-1.000000, 1.000000]"));
    testString(&iPos,   strdup("[1.000000, 2.000000]"));
    testString(&iDegen, strdup("[12.000000, 12.000000]"));
  }
}
