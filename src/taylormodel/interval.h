#ifndef INTERVAL_H
#define INTERVAL_H


#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdbool.h>



/* A closed interval I = [a, b], so that a <= b. */
typedef struct Interval {
  float left;
  float right;
} Interval;



/* Create a new interval that satisfies its definition. */
Interval newInterval(float left, float right);

void printInterval(Interval *source, FILE *where);


/*
  Interval operations.
*/
Interval addInterval(Interval *left, Interval *right);
Interval subInterval(Interval *left, Interval *right);
Interval mulInterval(Interval *left, Interval *right);
Interval divInterval(Interval *left, Interval *right);
bool eqInterval(Interval *left, Interval *right, float epsilon);
bool inInterval(Interval *left, Interval *right);



/*
  Properties of an Interval I = [a, b].
*/
/* W(I) = b - a */
float intervalWidth(Interval *source);
/* Mid(I) = (a + b) / 2 */
float intervalMidpoint(Interval *source);
/* Mag(I) = max{|a|, |b|} */
float intervalMagnitude(Interval *source);

/* Interval I is degenerate if a = b.
  Due to precision concerns, an epsilon must be specified. */
bool intervalIsDegenerate(Interval *source, float epsilon);


#endif
