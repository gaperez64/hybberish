#ifndef INTERVAL_H
#define INTERVAL_H

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

/* A closed interval I = [a, b], so that a <= b. */
typedef struct Interval {
  double left;
  double right;
} Interval;

/* Create a new interval that satisfies its definition. */
Interval newInterval(const double left, const double right);

void printInterval(const Interval *const source, FILE *where);

/*
  Interval operations.
*/
Interval addInterval(const Interval *const left, const Interval *const right);
Interval subInterval(const Interval *const left, const Interval *const right);
Interval mulInterval(const Interval *const left, const Interval *const right);
Interval divInterval(const Interval *const left, const Interval *const right);
Interval negInterval(const Interval *const source);
Interval sqrtInterval(const Interval *const source);
Interval powInterval(const Interval *const source, const unsigned int exponent);
Interval pow2Interval(const Interval *const source,
                      const unsigned int exponent);
bool eqInterval(const Interval *const left, const Interval *const right,
                const double epsilon);
bool subeqInterval(const Interval *const left, const Interval *const right);
bool elemInterval(const double elem, const Interval *const source);

/*
  Properties of an Interval I = [a, b].
*/
/* W(I) = b - a */
double intervalWidth(const Interval *const source);
/* Mid(I) = (a + b) / 2 */
double intervalMidpoint(const Interval *const source);
/* Mag(I) = max{|a|, |b|} */
double intervalMagnitude(const Interval *const source);

/* Interval I is degenerate if a = b.
  Due to precision concerns, an epsilon must be specified. */
bool intervalIsDegenerate(const Interval *const source, const double epsilon);

#endif
