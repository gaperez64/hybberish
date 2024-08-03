#include "interval.h"

Interval newInterval(float left, float right) {
  assert(left <= right);
  return (Interval){left, right};
}

void printInterval(Interval *source, FILE *where) {
  fprintf(where, "[");
  fprintf(where, "%f", source->left);
  fprintf(where, ", ");
  fprintf(where, "%f", source->right);
  fprintf(where, "]");
}

Interval addInterval(Interval *left, Interval *right) {
  /* [a, b] + [c, d] = [a + c, b + d] */
  return (Interval){left->left + right->left, left->right + right->right};
}

Interval subInterval(Interval *left, Interval *right) {
  /* [a, b] - [c, d] = [a - d, b - c] */
  return (Interval){left->left - right->right, left->right - right->left};
}

Interval mulInterval(Interval *left, Interval *right) {
  /* [a, b] * [c, d] = [ min{a*c, a*d, b*c, b*d},
                         max{a*c, a*d, b*c, b*d} ] */
  float ac, ad, bc, bd;
  ac = left->left * right->left;
  ad = left->left * right->right;
  bc = left->right * right->left;
  bd = left->right * right->right;
  return (Interval){fmin(ac, fmin(ad, fmin(bc, bd))),
                    fmax(ac, fmax(ad, fmax(bc, bd)))};
}

Interval divInterval(Interval *left, Interval *right) {
  /* [a, b] / [c, d] = [a, b] * [1/d, 1/c] */
  assert(right->left != 0);
  assert(right->right != 0);
  Interval invertedRight = {1. / right->right, 1. / right->left};

  return mulInterval(left, &invertedRight);
}

bool eqInterval(Interval *left, Interval *right, float epsilon) {
  /* [a, b] = [c, d] iff. a=c and b=d */
  return fabs(left->left - right->left) < epsilon &&
         fabs(left->right - right->right) < epsilon;
}

bool inInterval(Interval *left, Interval *right) {
  /* [a, b] in [c, d] iff. c <= a <= b <= d
    where a <= b is assumed by construction. */
  return right->left <= left->left && left->right <= right->right;
}

float intervalWidth(Interval *source) { return source->right - source->left; }

float intervalMidpoint(Interval *source) {
  return (source->left + source->right) / 2.;
}

float intervalMagnitude(Interval *source) {
  return fmax(fabs(source->left), fabs(source->right));
}

bool intervalIsDegenerate(Interval *source, float epsilon) {
  /* Interval width should always be positive, so eps must be as well. */
  assert(epsilon > 0);
  return intervalWidth(source) < epsilon;
}
