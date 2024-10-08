#include "interval.h"

Interval newInterval(const double left, const double right) {
  assert(left <= right);
  return (Interval){left, right};
}

void printInterval(const Interval *const source, FILE *where) {
  assert(source != NULL);
  assert(where != NULL);
  fprintf(where, "[");
  fprintf(where, "%f", source->left);
  fprintf(where, ", ");
  fprintf(where, "%f", source->right);
  fprintf(where, "]");
}

Interval addInterval(const Interval *const left, const Interval *const right) {
  /* [a, b] + [c, d] = [a + c, b + d] */
  assert(left != NULL);
  assert(right != NULL);
  return newInterval(left->left + right->left, left->right + right->right);
}

Interval subInterval(const Interval *const left, const Interval *const right) {
  /* [a, b] - [c, d] = [a - d, b - c] */
  assert(left != NULL);
  assert(right != NULL);
  return newInterval(left->left - right->right, left->right - right->left);
}

Interval mulInterval(const Interval *const left, const Interval *const right) {
  /* [a, b] * [c, d] = [ min{a*c, a*d, b*c, b*d},
                         max{a*c, a*d, b*c, b*d} ] */
  assert(left != NULL);
  assert(right != NULL);
  double ac, ad, bc, bd;
  ac = left->left * right->left;
  ad = left->left * right->right;
  bc = left->right * right->left;
  bd = left->right * right->right;
  return newInterval(fmin(ac, fmin(ad, fmin(bc, bd))),
                     fmax(ac, fmax(ad, fmax(bc, bd))));
}

Interval divInterval(const Interval *const left, const Interval *const right) {
  /* [a, b] / [c, d] = [a, b] * [1/d, 1/c] */
  assert(left != NULL);
  assert(right != NULL);
  /* not (c <= 0 <= d), i.e. 0 not in [c, d] */
  assert(!(right->left <= 0 && 0 <= right->right));
  Interval invertedRight = newInterval(1. / right->right, 1. / right->left);

  return mulInterval(left, &invertedRight);
}

Interval negInterval(const Interval *const source) {
  /* -[a, b] = [-b, -a] */
  assert(source != NULL);
  return newInterval(-source->right, -source->left);
}

Interval sqrtInterval(const Interval *const source) {
  /* sqrt([a, b]) = [sqrt(a), sqrt(b)] */
  assert(source != NULL);
  assert(0 <= source->left && source->left <= source->right);
  return newInterval(sqrtl(source->left), sqrtl(source->right));
}

Interval powInterval(const Interval *const source,
                     const unsigned int exponent) {
  assert(source != NULL);

  /* Edge case: [a, b]^0 = [1, 1] always. */
  if (exponent == 0)
    return newInterval(1, 1);

  /* Note: this algorithm may produce more accurate results (more strict
    bounds) than simply applying interval multiplication multiple times.
    Note: finding [a, b]^n means computing x^n for all x in [a, b]. */
  double an, bn;
  an = pow(source->left, exponent);
  bn = pow(source->right, exponent);

  /* CASES: [a, b]^n for n is ODD, so for x in [a, b], x^n retains x's sign */
  if ((exponent % 2) == 1)
    return newInterval(an, bn);

  /* CASES: [a, b]^n for n is EVEN, so for x in [a, b], x^n is positive */
  /* 0 <= a <= b, so a^n <= b^n */
  if (source->left >= 0)
    return newInterval(an, bn);
  /* a <= b < 0, so b^n <= a^n */
  if (source->right < 0)
    return newInterval(bn, an);
  /* a < 0 <= b, so [a, 0]^n = [0, a^n] and [0, b]^n = [0, b^n].
    Notice that the input negative sub-interval becomes positive,
    hence the 0 lower bound in the output. */
  return newInterval(0, fmax(an, bn));
}

Interval pow2Interval(const Interval *const source,
                      const unsigned int exponent) {
  assert(source != NULL);

  if (exponent == 0)
    return newInterval(1, 1);

  Interval result = *source;
  /* Multiplications only happen from exponent > 1 onwards. */
  for (unsigned int it = 1; it < exponent; ++it) {
    result = mulInterval(source, &result);
  }
  return result;
}

bool eqInterval(const Interval *const left, const Interval *const right,
                const double epsilon) {
  /* [a, b] = [c, d] iff. a=c and b=d */
  assert(left != NULL);
  assert(right != NULL);
  assert(epsilon >= 0);
  return fabs(left->left - right->left) < epsilon &&
         fabs(left->right - right->right) < epsilon;
}

bool subeqInterval(const Interval *const left, const Interval *const right) {
  /* [a, b] subseteq [c, d] iff. c <= a <= b <= d
    where a <= b is assumed by construction. */
  assert(left != NULL);
  assert(right != NULL);
  return right->left <= left->left && left->right <= right->right;
}

bool elemInterval(const double elem, const Interval *const source) {
  /* x in [a, b] iff. a <= x <= b
    where a <= b is assumed by construction. */
  assert(source != NULL);
  return source->left <= elem && elem <= source->right;
}

double intervalWidth(const Interval *const source) {
  assert(source != NULL);
  return source->right - source->left;
}

double intervalMidpoint(const Interval *const source) {
  assert(source != NULL);
  return (source->left + source->right) / 2.;
}

double intervalMagnitude(const Interval *const source) {
  assert(source != NULL);
  return fmax(fabs(source->left), fabs(source->right));
}

bool intervalIsDegenerate(const Interval *const source, const double epsilon) {
  /* Interval width should always be positive, so eps must be as well. */
  assert(epsilon >= 0);
  return intervalWidth(source) < epsilon;
}
