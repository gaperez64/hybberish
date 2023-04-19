/* Parser for a system of ODEs */

%define api.value.type {char *}
%define parse.error detailed
%define api.prefix {odes}
%param { yyscan_t scanner }

%code top {
  #include <stdio.h>
  #include "odeparse.h"
  #include "odes.lex.h"
}
%code requires {
  typedef void* yyscan_t;
}
%code {
  void yyerror(const char *);
}

/* tokens that will be used */
%token NUMBER IDENT SCOLON LPAR RPAR UNKNOWN
%token ADD SUB MUL DIV EXP EQUAL PRIME

%%

odelist: odedef
       | odelist odedef
       ;

odedef: IDENT PRIME EQUAL sumofprods SCOLON
      ;

sumofprods: prod
          | sumofprods addsub prod
          ;

addsub: ADD | SUB;

prod: factor
    | prod muldiv factor
    ;

muldiv: MUL | DIV;

factor: term
      | SUB term
      | term EXP NUMBER
      ;

term: NUMBER
    | IDENT
    | IDENT LPAR sumofprods RPAR
    | LPAR sumofprods RPAR
    ;

%%

void yyerror(const char *str) {
  fprintf(stderr, "[line %d] Error: %s\n", yylineno, str);
}

int parseOdeString(const char *in) {
  setOdeInputString(in);
  int rv = yyparse();
  endOdeScan();
  return rv;
}
