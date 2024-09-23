/* Scanner for a system of ODEs */

%{
#include "odeparse.h"
#include "odes.tab.h"

void yyerror(char const *);
%}

%option yylineno
%option noyywrap
%option noinput
%option nounput
%option prefix="odes"

%%

"+" { return ADD; }
"-" { return SUB; }
"*" { return MUL; }
"/" { return DIV; }
"^" { return EXP; }
"=" { return EQUAL; }
"'" { return PRIME; }
";" { return SCOLON; }
"(" { return LPAR; }
")" { return RPAR; }

[0-9]\.[0-9]*          { odeslval.str = strdup(yytext); return FLOAT; }
[0-9]|([1-9][0-9]+)    { odeslval.str = strdup(yytext); return INTEGER; }
[_a-zA-Z][_a-zA-Z0-9]* { odeslval.str = strdup(yytext); return IDENT; }

[ \t\n] { /* ignore white spaces */ }
.       { fprintf(stderr, "Unexpected symbol: %c\n", *yytext); return UNKNOWN; }

%%

void setOdeInputString(const char *in) {
  yy_scan_string(in);
}

void endOdeScan(void) {
  yy_delete_buffer(YY_CURRENT_BUFFER);
}
