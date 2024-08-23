#ifndef VARPARSE_H
#define VARPARSE_H

#include "variables.h"

int parseVarString(const char *, Domain **);
void setVarInputString(const char *);
void endVarScan(void);

#endif
