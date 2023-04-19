#ifndef ODEPARSE_H
#define ODEPARSE_H

#include "sysode.h"

int parseOdeString(const char *, ODEList **);
void setOdeInputString(const char *);
void endOdeScan(void);

#endif
