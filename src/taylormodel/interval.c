#include "interval.h"

Interval newInterval(const double left, const double right) {
  assert(left <= right);
  return (Interval){left, right};
}

void printInterval(const Interval * const source, FILE *where) {
  assert(source != NULL);
  assert(where != NULL);
  fprintf(where, "[");
  fprintf(where, "%f", source->left);
  fprintf(where, ", ");
  fprintf(where, "%f", source->right);
  fprintf(where, "]");
}

Interval addInterval(const Interval * const left, const Interval * const right) {
  /* [a, b] + [c, d] = [a + c, b + d] */
  assert(left != NULL);
  assert(right != NULL);
  return (Interval){left->left + right->left, left->right + right->right};
}

Interval subInterval(const Interval * const left, const Interval * const right) {
  /* [a, b] - [c, d] = [a - d, b - c] */
  assert(left != NULL);
  assert(right != NULL);
  return (Interval){left->left - right->right, left->right - right->left};
}

Interval mulInterval(const Interval * const left, const Interval * const right) {
  /* [a, b] * [c, d] = [ min{a*c, a*d, b*c, b*d},
                         max{a*c, a*d, b*c, b*d} ] */
  assert(left != NULL);
  assert(right != NULL);
  double ac, ad, bc, bd;
  ac = left->left * right->left;
  ad = left->left * right->right;
  bc = left->right * right->left;
  bd = left->right * right->right;
  return (Interval){fmin(ac, fmin(ad, fmin(bc, bd))),
                    fmax(ac, fmax(ad, fmax(bc, bd)))};
}

Interval divInterval(const Interval * const left, const Interval * const right) {
  /* [a, b] / [c, d] = [a, b] * [1/d, 1/c] */
  assert(left != NULL);
  assert(right != NULL);
  /* not (c <= 0 <= d), i.e. 0 not in [c, d] */
  assert(!(right->left <= 0 && 0 <= right->right));
  Interval invertedRight = {1. / right->right, 1. / right->left};

  return mulInterval(left, &invertedRight);
}

Interval negInterval(const Interval * const source) {
  /* -[a, b] = [-b, -a] */
  assert(source != NULL);
  return (Interval){ -source->right, -source->left };
}

Interval sqrtInterval(const Interval * const source) {
  /* sqrt([a, b]) = [sqrt(a), sqrt(b)] */
  assert(source != NULL);
  assert(0 <= source->left && source->left <= source->right);
  return (Interval){ sqrtl(source->left), sqrtl(source->right) };
}

bool eqInterval(const Interval * const left, const Interval * const right, const double epsilon) {
  /* [a, b] = [c, d] iff. a=c and b=d */
  assert(left != NULL);
  assert(right != NULL);
  assert(epsilon >= 0);
  return fabs(left->left - right->left) < epsilon &&
         fabs(left->right - right->right) < epsilon;
}

bool inInterval(const Interval * const left, const Interval * const right) {
  /* [a, b] in [c, d] iff. c <= a <= b <= d
    where a <= b is assumed by construction. */
  assert(left != NULL);
  assert(right != NULL);
  return right->left <= left->left && left->right <= right->right;
}

double intervalWidth(const Interval * const source) {
  assert(source != NULL);
  return source->right - source->left;
}

double intervalMidpoint(const Interval * const source) {
  assert(source != NULL);
  return (source->left + source->right) / 2.;
}

double intervalMagnitude(const Interval * const source) {
  assert(source != NULL);
  return fmax(fabs(source->left), fabs(source->right));
}

bool intervalIsDegenerate(const Interval * const source, const double epsilon) {
  /* Interval width should always be positive, so eps must be as well. */
  assert(epsilon >= 0);
  return intervalWidth(source) < epsilon;
}
