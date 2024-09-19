#include "funexp.h"
#include "odeparse.h"
#include "sysode.h"
#include "taylormodel.h"
#include "tmflowpipe.h"
#include "varparse.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
  /* to avoid silly warnings about unused parameters */
  (void)argc;
  (void)argv;

  /* Attach parsing to the flowpipe pipeline */
  {
    int code;                     // Parsing output code
    const unsigned int order = 3; // TM order
    const unsigned int k = 3;     // truncation order

    /* Setup */
    printf("=== Setup ===\n");
    fflush(stdout);

    char *odeStr = "x' = 1 + y; y' = x^2;";
    ODEList *odes;
    code = parseOdeString(odeStr, &odes);
    printOdeList(odes, stdout);
    printf("\n");
    fflush(stdout);
    assert(code == 0);

    char *varStr = "x in [-1, 1]; y in [-0.5, 0.5];";
    Domain *domains;
    code = parseVarString(varStr, &domains);
    printDomain(domains, stdout);
    printf("\n\n");
    fflush(stdout);
    assert(code == 0);

    /* TM flowpipe overapprox pipeline */
    printf("=== Flowpipe overapprox ===\n");
    fflush(stdout);
    TaylorModel *tms = computeTaylorPolynomial(odes, order, k);
    printTaylorModel(tms, stdout);
    printf("\n");
    fflush(stdout);

    // TODO: Add safe remainder computation ...

    /* Clean */
    delTaylorModel(tms);
    delDomain(domains);
    delOdeList(odes);
  }
}
