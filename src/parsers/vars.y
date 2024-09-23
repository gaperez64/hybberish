/* Parser for a valuation of a set of variables */

/* To be compatible with Bison 3.0.4 */
%define parse.error verbose
%define api.prefix {vars}

%code top {
  #include "varparse.h"
  #include "vars.lex.h"
  #include "variables.h"
  #include <stddef.h>
  #include <stdio.h>
}
%code {
  void yyerror(const char *);
  Domain *finlist;
}

/* tokens that will be used */
%union {
  Domain *list;
  Interval dom;
  double num;
  char *str;
}
%token SCOLON COMMA LBRAC RBRAC UNKNOWN
%token ELEMENT
%token <str> IDENT
%token <num> FLOAT INTEGER
%type <num> number
%type <dom> interval
%type <list> varlist vardef

%%

spec: varlist { finlist = $1; }
    ;

varlist: vardef          { $$ = $1; }
       | varlist vardef  { $$ = appDomainElem($1, $2); }
       ;

vardef: IDENT ELEMENT interval SCOLON  { $$ = newDomain($1, $3); }
      ;

interval: LBRAC number COMMA number RBRAC { $$ = newInterval($2, $4); }
        ;

number: FLOAT    { $$ = $1; }
      | INTEGER  { $$ = $1; }
      ;

%%

void yyerror(const char *str) {
  fprintf(stderr, "[line %d] Error: %s\n", varslineno, str);
}

int parseVarString(const char *in, Domain **ptlist) {
  setVarInputString(in);
  int rv = yyparse();
  endVarScan();
  *ptlist = finlist;
  return rv;
}
