#include "interval.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

bool approx(float actual, float expected, float epsilon) {
  assert(epsilon >= 0);
  return fabs(actual - expected) < epsilon;
}

/* Test if actual = [a, b] = [expectedLeft, expectedRight]
  up to precision epsilon. and print out the results. */
void testInterval(Interval *actual, float expectedLeft, float expectedRight,
                  float epsilon) {
  assert(actual != NULL);
  assert(epsilon >= 0);

  /* printing and comparison */
  bool compare = approx(actual->left, expectedLeft, epsilon) &&
                 approx(actual->right, expectedRight, epsilon);

  printf("Expect: [%f, %f]\n", expectedLeft, expectedRight);
  printf("Actual:    [%f, %f]\n", actual->left, actual->right);
  printf("Compare = %i\n", compare);
  fflush(stdout); // Flush stdout
  assert(compare);
}

/* Test if actual = expected, and print out the results. */
void testBool(bool actual, bool expected) {
  /* printing and comparison */
  printf("Expect: %i\n", expected);
  printf("Actual:    %i\n", actual);
  fflush(stdout); // Flush stdout
  assert(actual == expected);
}

/* Test if actual = expected up to precision epsilon, and print out the results.
 */
void testFloat(float actual, float expected, float epsilon) {
  /* printing and comparison */
  printf("Expect: %f\n", expected);
  printf("Actual:    %f\n", actual);
  fflush(stdout); // Flush stdout
  assert(fabs(expected - actual) < epsilon);
}

/* Test if actual is printed to the expected string, and print out the results.
 */
void testString(Interval *actual, const char *expectedMsg) {
  /* printing and comparison */
  char *buffer = NULL;
  size_t buflen = 0;
  FILE *stream = open_memstream(&buffer, &buflen);
  assert(stream != NULL);
  printInterval(actual, stream);
  fclose(stream);

  printf("Expect: \"%s\"\n", expectedMsg);
  printf("Actual:    \"%s\"\n", buffer);
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
  Interval iNeg = newInterval(-2, -1);
  Interval iOrig = newInterval(-1, 1);
  Interval iPos = newInterval(1, 2);
  Interval iDegen = newInterval(12, 12);

  /* Enforce relationships within the test data,
    which are exploited in the testing. */
  assert(iNeg.right <= iOrig.left);
  assert(iOrig.right <= iPos.left);
  assert(iDegen.left == iDegen.right);

  /* Test binary addition for various interval combinations.
    Focus on the effects of positive values, negative values,
    or zero being part of an interval. */
  {
    printf("\n=== binary ADD (+) ===\n");
    fflush(stdout);

    Interval res = addInterval(&iNeg, &iNeg);
    testInterval(&res, -4, -2, eps);

    res = addInterval(&iOrig, &iOrig);
    testInterval(&res, -2, 2, eps);

    res = addInterval(&iPos, &iPos);
    testInterval(&res, 2, 4, eps);

    res = addInterval(&iDegen, &iOrig);
    testInterval(&res, 11, 13, eps);
  }

  /* Test binary subtraction for various interval combinations.
    Focus on the effects of positive values, negative values,
    or zero being part of an interval. */
  {
    printf("\n=== binary SUB (-) ===\n");
    fflush(stdout);

    Interval res = subInterval(&iNeg, &iNeg);
    testInterval(&res, -1, 1, eps);

    res = subInterval(&iOrig, &iOrig);
    testInterval(&res, -2, 2, eps);

    res = subInterval(&iPos, &iPos);
    testInterval(&res, -1, 1, eps);

    res = subInterval(&iDegen, &iNeg);
    testInterval(&res, 13, 14, eps);
  }

  /* Test binary multiplication for various interval combinations.
    Focus on the effects of positive values, negative values,
    or zero being part of an interval. */
  {
    printf("\n=== binary MUL (*) ===\n");
    fflush(stdout);

    Interval res = mulInterval(&iNeg, &iNeg);
    testInterval(&res, 1, 4, eps);

    res = mulInterval(&iOrig, &iOrig);
    testInterval(&res, -1, 1, eps);

    res = mulInterval(&iPos, &iPos);
    testInterval(&res, 1, 4, eps);

    res = mulInterval(&iDegen, &iPos);
    testInterval(&res, 12, 24, eps);
  }

  /* Test binary division for various interval combinations.
    Focus on the effects of positive values, negative values,
    or zero being part of an interval. */
  {
    printf("\n=== binary DIV (/) ===\n");
    fflush(stdout);

    Interval res = divInterval(&iNeg, &iNeg);
    testInterval(&res, 0.5, 2, eps);

    /* The divisor may not contain: [a, b] / [c, d]  =>  0 not in [c, d] */
    // res = divInterval(&iOrig, &iOrig);
    // testInterval(&res, -1, 1, eps);

    res = divInterval(&iPos, &iPos);
    testInterval(&res, 0.5, 2, eps);

    res = divInterval(&iDegen, &iNeg);
    testInterval(&res, -12, -6, eps);
  }

  /* Test (unary) additive inverse for various interval combinations.
    Focus on the effects of positive values, negative values,
    or zero being part of an interval. */
  {
    printf("\n=== unary NEG (-) ===\n");
    fflush(stdout);

    Interval res = negInterval(&iNeg);
    testInterval(&res, -iNeg.right, -iNeg.left, eps);

    res = negInterval(&iOrig);
    testInterval(&res, -iOrig.right, -iOrig.left, eps);

    res = negInterval(&iPos);
    testInterval(&res, -iPos.right, -iPos.left, eps);

    res = negInterval(&iDegen);
    testInterval(&res, -iDegen.right, -iDegen.left, eps);
  }

  /* Test (unary) square root function. */
  {
    printf("\n=== fun sqrt ===\n");
    fflush(stdout);

    Interval res = sqrtInterval(&iPos);
    testInterval(&res, sqrtl(iPos.left), sqrtl(iPos.right), eps);

    res = sqrtInterval(&iDegen);
    testInterval(&res, sqrtl(iDegen.left), sqrtl(iDegen.right), eps);
  }

  /* Test (binary) exponentiation for various interval combinations.
    Focus on the effects of positive values, negative values,
    or zero being part of an interval.

    Test both the smart and the naive exponentiation algorithm.
  */
  {
    printf("\n=== binary EXP (^) ===\n");
    fflush(stdout);

    Interval zero = newInterval(0, 0);
    unsigned int exponent;

    /* exponent = 0: special case */
    exponent = 0;
    Interval res = powInterval(&iNeg, exponent);
    testInterval(&res, 1, 1, eps);
    res = pow2Interval(&iNeg, exponent);
    testInterval(&res, 1, 1, eps);

    res = powInterval(&iOrig, exponent);
    testInterval(&res, 1, 1, eps);
    res = pow2Interval(&iOrig, exponent);
    testInterval(&res, 1, 1, eps);

    res = powInterval(&iPos, exponent);
    testInterval(&res, 1, 1, eps);
    res = pow2Interval(&iPos, exponent);
    testInterval(&res, 1, 1, eps);

    res = powInterval(&iDegen, exponent);
    testInterval(&res, 1, 1, eps);
    res = pow2Interval(&iDegen, exponent);
    testInterval(&res, 1, 1, eps);

    res = powInterval(&zero, exponent);
    testInterval(&res, 1, 1, eps);
    res = pow2Interval(&zero, exponent);
    testInterval(&res, 1, 1, eps);

    /* exponent = 1: ODD case */
    exponent = 1;
    res = powInterval(&iNeg, exponent);
    testInterval(&res, -2, -1, eps);
    res = pow2Interval(&iNeg, exponent);
    testInterval(&res, -2, -1, eps);

    res = powInterval(&iOrig, exponent);
    testInterval(&res, -1, 1, eps);
    res = pow2Interval(&iOrig, exponent);
    testInterval(&res, -1, 1, eps);

    res = powInterval(&iPos, exponent);
    testInterval(&res, 1, 2, eps);
    res = pow2Interval(&iPos, exponent);
    testInterval(&res, 1, 2, eps);

    res = powInterval(&iDegen, exponent);
    testInterval(&res, 12, 12, eps);
    res = pow2Interval(&iDegen, exponent);
    testInterval(&res, 12, 12, eps);

    res = powInterval(&zero, exponent);
    testInterval(&res, 0, 0, eps);
    res = pow2Interval(&zero, exponent);
    testInterval(&res, 0, 0, eps);

    /* exponent = 3: ODD case */
    exponent = 3;
    res = powInterval(&iNeg, exponent);
    testInterval(&res, -8, -1, eps);
    res = pow2Interval(&iNeg, exponent);
    testInterval(&res, -8, -1, eps);

    res = powInterval(&iOrig, exponent);
    testInterval(&res, -1, 1, eps);
    res = pow2Interval(&iOrig, exponent);
    testInterval(&res, -1, 1, eps);

    res = powInterval(&iPos, exponent);
    testInterval(&res, 1, 8, eps);
    res = pow2Interval(&iPos, exponent);
    testInterval(&res, 1, 8, eps);

    res = powInterval(&iDegen, exponent);
    testInterval(&res, 1728, 1728, eps);
    res = pow2Interval(&iDegen, exponent);
    testInterval(&res, 1728, 1728, eps);

    res = powInterval(&zero, exponent);
    testInterval(&res, 0, 0, eps);
    res = pow2Interval(&zero, exponent);
    testInterval(&res, 0, 0, eps);

    /* exponent = 2: EVEN case */
    exponent = 2;
    res = powInterval(&iNeg, exponent);
    testInterval(&res, 1, 4, eps);
    res = pow2Interval(&iNeg, exponent);
    testInterval(&res, 1, 4, eps);

    res = powInterval(&iOrig, exponent);
    testInterval(&res, 0, 1, eps);
    res = pow2Interval(&iOrig, exponent);
    testInterval(&res, -1, 1, eps);

    res = powInterval(&iPos, exponent);
    testInterval(&res, 1, 4, eps);
    res = pow2Interval(&iPos, exponent);
    testInterval(&res, 1, 4, eps);

    res = powInterval(&iDegen, exponent);
    testInterval(&res, 144, 144, eps);
    res = pow2Interval(&iDegen, exponent);
    testInterval(&res, 144, 144, eps);

    res = powInterval(&zero, exponent);
    testInterval(&res, 0, 0, eps);
    res = pow2Interval(&zero, exponent);
    testInterval(&res, 0, 0, eps);
  }

  /* Test binary interval equality. */
  {
    printf("\n=== binary EQ (=) ===\n");
    fflush(stdout);

    testBool(eqInterval(&iNeg, &iNeg, eps), true);
    testBool(eqInterval(&iOrig, &iOrig, eps), true);
    testBool(eqInterval(&iPos, &iPos, eps), true);

    /* Test inequality in a transitive fashion,
      instead of enumerating all possible options. */
    testBool(eqInterval(&iNeg, &iOrig, eps), false);
    testBool(eqInterval(&iOrig, &iPos, eps), false);
    testBool(eqInterval(&iPos, &iNeg, eps), false);
  }

  /* Test binary interval subset equal relationship.
    i.e. is interval I1 contained in interval I2? */
  {
    printf("\n=== binary subset (subseq) ===\n");
    fflush(stdout);

    /* All intervals are closed, so should be contained in themselves. */
    testBool(subeqInterval(&iNeg, &iNeg), true);
    testBool(subeqInterval(&iOrig, &iOrig), true);
    testBool(subeqInterval(&iPos, &iPos), true);

    /*
        iNeg  = [a, b]
        iOrig = [c, d]
        iPos  = [e, f]
        encompassing = [a, f]
        where a <= b,c,d,e <= f
    */
    Interval encompassing = newInterval(iNeg.left, iPos.right);
    testBool(subeqInterval(&iNeg, &encompassing), true);
    testBool(subeqInterval(&iOrig, &encompassing), true);
    testBool(subeqInterval(&iPos, &encompassing), true);
    testBool(subeqInterval(&encompassing, &iNeg), false);
    testBool(subeqInterval(&encompassing, &iOrig), false);
    testBool(subeqInterval(&encompassing, &iPos), false);

    /*
        iNeg  = [a, b]
        leftEncompassing = [d, c]
        rightEncompassing = [c, e]
        where c in [a, b] but c != a & c != b
        so that d < a < c < b < e
    */
    float offset = 1.;
    Interval leftEncompassing =
        newInterval(iNeg.left - offset, (iNeg.left + iNeg.right) / 2.);
    Interval rightEncompassing =
        newInterval((iNeg.left + iNeg.right) / 2., iNeg.right + offset);
    testBool(subeqInterval(&iNeg, &leftEncompassing), false);
    testBool(subeqInterval(&iNeg, &rightEncompassing), false);
  }

  /* Test binary interval membership.
    i.e. is a number x contained in interval I? */
  {
    printf("\n=== binary membership (elem of) ===\n");
    fflush(stdout);

    Interval i = newInterval(-1, 1);

    testBool(elemInterval(-30.0, &i), false);
    testBool(elemInterval(-1.00001, &i), false);
    testBool(elemInterval(-1, &i), true);
    testBool(elemInterval(0, &i), true);
    testBool(elemInterval(1, &i), true);
    testBool(elemInterval(1.00001, &i), false);
    testBool(elemInterval(30.0, &i), false);
  }

  /* Test interval properties. */
  {
    /* Width(I) */
    testFloat(intervalWidth(&iNeg), (iNeg.right - iNeg.left), eps);
    testFloat(intervalWidth(&iOrig), (iOrig.right - iOrig.left), eps);
    testFloat(intervalWidth(&iPos), (iPos.right - iPos.left), eps);
    testFloat(intervalWidth(&iDegen), (iDegen.right - iDegen.left), eps);
    testFloat(intervalWidth(&iDegen), 0., eps);

    /* Mid(I) */
    testFloat(intervalMidpoint(&iNeg), (iNeg.left + iNeg.right) / 2., eps);
    testFloat(intervalMidpoint(&iOrig), (iOrig.left + iOrig.right) / 2., eps);
    testFloat(intervalMidpoint(&iPos), (iPos.left + iPos.right) / 2., eps);
    testFloat(intervalMidpoint(&iDegen), (iDegen.left + iDegen.right) / 2.,
              eps);
    testFloat(intervalMidpoint(&iDegen), iDegen.left, eps);

    /* Mag(I) */
    testFloat(intervalMagnitude(&iNeg), fmax(fabs(iNeg.left), fabs(iNeg.right)),
              eps);
    testFloat(intervalMagnitude(&iOrig),
              fmax(fabs(iOrig.left), fabs(iOrig.right)), eps);
    testFloat(intervalMagnitude(&iPos), fmax(fabs(iPos.left), fabs(iPos.right)),
              eps);
    testFloat(intervalMagnitude(&iDegen),
              fmax(fabs(iDegen.left), fabs(iDegen.right)), eps);
    testFloat(intervalMagnitude(&iDegen), fabs(iDegen.left), eps);

    /* Test interval degeneracy: I = [a, a] */
    testBool(intervalIsDegenerate(&iNeg, eps), false);
    testBool(intervalIsDegenerate(&iOrig, eps), false);
    testBool(intervalIsDegenerate(&iPos, eps), false);
    testBool(intervalIsDegenerate(&iDegen, eps), true);
  }

  /* Test printing. */
  {
    testString(&iNeg, "[-2.000000, -1.000000]");
    testString(&iOrig, "[-1.000000, 1.000000]");
    testString(&iPos, "[1.000000, 2.000000]");
    testString(&iDegen, "[12.000000, 12.000000]");
  }
}
