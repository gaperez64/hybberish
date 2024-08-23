/* Scanner for a valuation of a set of variables  */

%{
#include "varparse.h"
#include "vars.tab.h"

void yyerror(char const *);
%}

%option yylineno
%option noyywrap
%option noinput
%option nounput
%option prefix="vars"

%%

"in" { return ELEMENT; }
";"  { return SCOLON; }
","  { return COMMA; }
"["  { return LBRAC; }
"]"  { return RBRAC; }

[+-]?([0-9]\.[0-9]*)        { varslval.num = atof(yytext); return FLOAT; }
[+-]?([0-9]|([1-9][0-9]+))  { varslval.num = atof(yytext); return INTEGER; }
[_a-zA-Z][_a-zA-Z0-9]*      { varslval.str = strdup(yytext); return IDENT; }

[ \t\n] { /* ignore white spaces */ }
.       { fprintf(stderr, "Unexpected symbol: %c\n", *yytext); return UNKNOWN; }

%%

void setVarInputString(const char *in) {
  yy_scan_string(in);
}

void endVarScan(void) {
  yy_delete_buffer(YY_CURRENT_BUFFER);
}
