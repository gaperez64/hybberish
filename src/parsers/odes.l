/* Scanner for a system of ODEs */

%{
#include "odes.tab.h"

void yyerror (char const *);
%}

%option yylineno
%option noyywrap
%option reentrant
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

[0-9](\.[0-9]*)?       { odeslval = strdup(yytext); return NUMBER; }
[_a-zA-Z][_a-zA-Z0-9]* { odeslval = strdup(yytext); return IDENT; }

[ \t\n] { /* ignore white spaces */ }
.       { fprintf(stderr, "Unexpected symbol: %c\n", *yytext); return UNKNOWN; }

%%
