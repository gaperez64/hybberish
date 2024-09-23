#include "funexp.h"
#include "sysode.h"
#include "taylormodel.h"
#include "tmflowpipe.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Print the ODEList. */
void printODETest(ODEList *got) {
  /* printing and comparison */
  char *buffer = NULL;
  size_t buflen = 0;
  FILE *stream = open_memstream(&buffer, &buflen);
  assert(stream != NULL);
  printOdeList(got, stream);
  fclose(stream);

  printf("Got ODE:    \"%s\"\n", buffer);
  fflush(stdout); // Flush stdout

  /* clean */
  free(buffer);
}

/* Print the ODEList. */
void printTPTest(TaylorModel *got) {
  /* printing and comparison */
  char *buffer = NULL;
  size_t buflen = 0;
  FILE *stream = open_memstream(&buffer, &buflen);
  assert(stream != NULL);
  printTaylorModel(got, stream);
  fclose(stream);

  printf("Got TP:     \"%s\"\n", buffer);
  fflush(stdout); // Flush stdout

  /* clean */
  free(buffer);
}

int main(int argc, char *argv[]) {
  /* to avoid silly warnings about unused parameters */
  (void)argc;
  (void)argv;

  {
    /* Setup */
    ExpTree *x = newExpLeaf(EXP_VAR, "x");
    ExpTree *y = newExpLeaf(EXP_VAR, "y");
    ExpTree *z = newExpLeaf(EXP_VAR, "z");
    /* construct an ODE chain 'var3 <- var2 <- var1' */
    ODEList *sys = newOdeElem(NULL, strdup(z->data), z);
    sys = newOdeElem(sys, strdup(y->data), y);
    sys = newOdeElem(sys, strdup(x->data), x);
    printODETest(sys);

    TaylorModel *list = initTaylorModel(sys);
    printf("(List) ");
    printTPTest(list);
    fflush(stdout);
    TaylorModel *cpy = cpyTaylorModel(list);
    printf("(Copy) ");
    printTPTest(cpy);
    fflush(stdout);
    cpy = reverseTaylorModel(cpy);
    printf("(Rev)  ");
    printTPTest(cpy);
    fflush(stdout);

    TaylorModel *poly = computeTaylorPolynomial(sys, 3, 3);
    printTPTest(poly);

    /* Clean */
    delOdeList(sys);
    delTaylorModel(list);
    delTaylorModel(cpy);
    delTaylorModel(poly);
  }

  {
    /* Setup */
    ExpTree *x = newExpLeaf(EXP_VAR, "x");
    ExpTree *y = newExpLeaf(EXP_VAR, "y");
    ExpTree *t = newExpLeaf(EXP_VAR, "t");

    {
      ExpTree *g = newExpOp(EXP_SUB_OP,
                            newExpOp(EXP_ADD_OP, cpyExpTree(t), cpyExpTree(x)),
                            cpyExpTree(y));
      printf("\ng   = ");
      printExpTree(g, stdout);
      printf("\n");

      /* x' = (x + y) */
      ExpTree *expX = newExpOp(EXP_ADD_OP, cpyExpTree(x), cpyExpTree(y));
      /* y' = (1 - (y^2)) */
      ExpTree *expY = newExpOp(
          EXP_SUB_OP, newExpLeaf(EXP_NUM, "1"),
          newExpOp(EXP_EXP_OP, cpyExpTree(y), newExpLeaf(EXP_NUM, "2")));
      printf("ODE = | x' | = | ");
      printExpTree(expX, stdout);
      printf(" |\n");
      printf("      | y' |   | ");
      printExpTree(expY, stdout);
      printf(" |\n");

      /* Construct an ODE chain. */
      ODEList *sys = newOdeElem(NULL, strdup(y->data), expY);
      sys = newOdeElem(sys, strdup(x->data), expX);
      printODETest(sys);

      TaylorModel *poly = computeTaylorPolynomial(sys, 2, 2);
      printTPTest(poly);

      printf("\n== Manual Lie Der calls ==\n");

      // ExpTree *lieDer = simplify(lieDerivative(sys, g));
      ExpTree *lieDer = lieDerivative(sys, g);
      printf("L(g):    ");
      printExpTree(lieDer, stdout);
      printf("\n");
      fflush(stdout);

      ExpTree *lieDer2nd = lieDerivative(sys, lieDer);
      ExpTree *simplified = simplify(lieDer2nd);
      printf("L^2(g):  ");
      printExpTree(simplified, stdout);
      printf("\n\n");
      fflush(stdout);

      /* Clean */
      delExpTree(g);
      delExpTree(lieDer);
      delExpTree(lieDer2nd);
      delExpTree(simplified);
      delOdeList(sys);
      delTaylorModel(poly);
    }

    {
      /* x' = (1 + y) */
      ExpTree *expX = newExpOp(EXP_ADD_OP, newOneExpTree(), cpyExpTree(y));
      /* y' = x^2 */
      ExpTree *expY =
          newExpOp(EXP_EXP_OP, cpyExpTree(x), newExpLeaf(EXP_NUM, "2"));
      printf("ODE = | x' | = | ");
      printExpTree(expX, stdout);
      printf(" |\n");
      printf("      | y' |   | ");
      printExpTree(expY, stdout);
      printf(" |\n");

      /* Construct an ODE chain. */
      ODEList *sys = newOdeElem(NULL, strdup(y->data), expY);
      sys = newOdeElem(sys, strdup(x->data), expX);
      printODETest(sys);

      /* Should (algebraically) simplify to for order=3:
        p(x) = 2x + 2xy;
        p(y) = 2 + 4y + 2y^3 + 2x^3
      */
      TaylorModel *functions = initTaylorModel(sys);
      TaylorModel *derived = lieDerivativeK(sys, functions, 3);
      printTPTest(derived);

      TaylorModel *poly = computeTaylorPolynomial(sys, 3, 3);
      printTPTest(poly);

      /* Clean */
      delOdeList(sys);
      delTaylorModel(functions);
      delTaylorModel(derived);
      delTaylorModel(poly);
    }

    /*
      Test Picard operator
    */
    printf("\n### Picard Operator ###\n");
    fflush(stdout);

    /* Test substitutions */
    {
      printf("=== Substitute functions into ODEs ===\n");
      fflush(stdout);

      /* x' = (1 + y) */
      ExpTree *odeX = newExpOp(EXP_ADD_OP, newOneExpTree(), cpyExpTree(y));
      /* y' = x^2 */
      ExpTree *odeY = newExpOp(
          EXP_NEG,
          newExpOp(EXP_EXP_OP, cpyExpTree(x), newExpLeaf(EXP_NUM, "2")), NULL);
      /* Construct an ODE chain. */
      ODEList *sys = newOdeElem(NULL, strdup(y->data), odeY);
      sys = newOdeElem(sys, strdup(x->data), odeX);

      /* x + s */
      ExpTree *funcX =
          newExpOp(EXP_ADD_OP, cpyExpTree(x), newExpLeaf(EXP_VAR, "t"));
      /* y */
      ExpTree *funcY = cpyExpTree(y);
      /* Construct a TaylorModel chain. */
      TaylorModel *functions =
          newTMElem(NULL, strdup(y->data), funcY, newInterval(0, 0));
      functions =
          newTMElem(functions, strdup(x->data), funcX, newInterval(0, 0));

      TaylorModel *substituted = substituteTaylorModel(sys, functions);
      printODETest(sys);
      printTPTest(functions);
      printTPTest(substituted);

      /* Clean */
      delTaylorModel(substituted);
      delTaylorModel(functions);
      delOdeList(sys);
    }

    /* Test actual picard operator */
    {
      printf("=== Picard operator ===\n");
      fflush(stdout);

      /* x' = (1 + y) */
      ExpTree *odeX = newExpOp(EXP_ADD_OP, newOneExpTree(), cpyExpTree(y));
      /* y' = x^2 */
      ExpTree *odeY =
          newExpOp(EXP_EXP_OP, cpyExpTree(x), newExpLeaf(EXP_NUM, "2"));
      /* Construct an ODE chain. */
      ODEList *sys = newOdeElem(NULL, strdup(y->data), odeY);
      sys = newOdeElem(sys, strdup(x->data), odeX);

      /* x + s */
      ExpTree *funcX =
          newExpOp(EXP_ADD_OP, cpyExpTree(x), newExpLeaf(EXP_VAR, "t"));
      /* y */
      ExpTree *funcY = cpyExpTree(y);
      /* Construct a TaylorModel chain. */
      TaylorModel *functions =
          newTMElem(NULL, strdup(y->data), funcY, newInterval(0, 0));
      functions =
          newTMElem(functions, strdup(x->data), funcX, newInterval(0, 0));

      TaylorModel *picard = picardOperator(sys, functions);
      printODETest(sys);
      printTPTest(functions);
      printTPTest(picard);

      Interval domx = newInterval(-1, 1);
      Interval doms = newInterval(0, 0.02);
      Interval third = newInterval(1. / 3., 1. / 3.);
      Interval J0 = newInterval(-0.1, 0.1);

      /* -(x^2)s - x(s^2) - (1/3)s^3 */
      Interval x2 = pow2Interval(&domx, 2);
      Interval s2 = pow2Interval(&doms, 2);
      Interval s3 = pow2Interval(&doms, 3);
      Interval thirds3 = mulInterval(&third, &s3);
      Interval x2s = mulInterval(&x2, &doms);
      Interval xs2 = mulInterval(&domx, &s2);
      Interval intpe = negInterval(&x2s);
      printf("Int(pe) = ");
      printInterval(&intpe, stdout);
      printf("\n");
      intpe = subInterval(&intpe, &xs2);
      printf("Int(pe) = ");
      printInterval(&intpe, stdout);
      printf("\n");
      intpe = subInterval(&intpe, &thirds3);
      printf("Int(pe) = ");
      printInterval(&intpe, stdout);
      printf("\n");

      printf("x^2 = ");
      printInterval(&x2, stdout);
      printf("\n");

      printf("s^2 = ");
      printInterval(&s2, stdout);
      printf("\n");

      printf("s^3 = ");
      printInterval(&s3, stdout);
      printf("\n\n");

      fflush(stdout);

      printf("-(x^2)s - x(s^2) - (1/3)(s^3)\n");
      printf("-");
      printInterval(&x2, stdout);
      printf("*");
      printInterval(&doms, stdout);
      printf(" - ");
      printInterval(&domx, stdout);
      printf("*");
      printInterval(&s2, stdout);
      printf(" - ");
      printInterval(&third, stdout);
      printf("*");
      printInterval(&s3, stdout);
      printf("\n");

      printf("-");
      printInterval(&x2s, stdout);
      printf(" - ");
      printInterval(&xs2, stdout);
      printf(" - ");
      printInterval(&thirds3, stdout);
      printf("\n");

      printInterval(&intpe, stdout);
      printf("\n");
      fflush(stdout);

      Interval intpeJ0 = addInterval(&intpe, &J0);
      Interval res = mulInterval(&intpeJ0, &doms);

      printf("\n(Int(pe) + J0) * [ai, bi]\n");
      printf("(");
      printInterval(&intpe, stdout);
      printf(" + ");
      printInterval(&J0, stdout);
      printf(") * ");
      printInterval(&doms, stdout);
      printf("\n");
      printInterval(&intpeJ0, stdout);
      printf(" * ");
      printInterval(&doms, stdout);
      printf("\n");
      printInterval(&res, stdout);
      printf("\n");
      fflush(stdout);

      /* Clean */
      delTaylorModel(picard);
      delTaylorModel(functions);
      delOdeList(sys);
    }

    delExpTree(x);
    delExpTree(y);
    delExpTree(t);
  }
}
