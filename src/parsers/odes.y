/* Parser for a system of ODEs */

%define api.value.type {char*}
%define parse.lac full
%define parse.error detailed
%define api.prefix {odes}
%locations
%define api.pure
%param { yyscan_t scanner }

%code top {
  #include <stdio.h>
  #include "odes.lex.h"
}
%code requires {
  typedef void* yyscan_t;
}
%code {
  /*int yylex(YYSTYPE* yylvalp, YYLTYPE* yyllocp, yyscan_t scanner);*/
  void yyerror(YYLTYPE* yyllocp, yyscan_t unused, const char* msg);
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

void yyerror(YYLTYPE* yyllocp, yyscan_t unused, const char* msg) {
  fprintf(stderr, "[%d:%d]: %s\n",
          yyllocp->first_line, yyllocp->first_column, msg);
}
