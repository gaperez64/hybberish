/* Parser for a system of ODEs */

%define api.value.type {char *}
%define parse.error detailed
%define api.prefix {odes}

%code top {
  #include <stddef.h>
  #include <stdio.h>
  #include "odeparse.h"
  #include "odes.lex.h"
}
%code {
  void yyerror(const char *);
  ODEList *finlist;
}

/* tokens that will be used */
%union {
  ODEList *list;
  ODEExpTree *tree;
}
%token NUMBER IDENT SCOLON LPAR RPAR UNKNOWN
%token ADD SUB MUL DIV EXP EQUAL PRIME
%type <list> odelist odedef
%type <tree> sumofprods prod factor term

%%

spec: odelist { finlist = $1; }
    ;

odelist: odedef          { $$ = $1; }
       | odelist odedef  { $$ = appOdeElem($1, $2); }
       ;

odedef: IDENT PRIME EQUAL sumofprods SCOLON  { $$ = newOdeList($1, $4); }
      ;

sumofprods: prod                 { $$ = $1; }
          | sumofprods ADD prod  { $$ = newOdeExpOp(ODE_ADD_OP, $1, $3); }
          | sumofprods SUB prod  { $$ = newOdeExpOp(ODE_SUB_OP, $1, $3); }
          ;

prod: factor           { $$ = $1; }
    | prod MUL factor  { $$ = newOdeExpOp(ODE_MUL_OP, $1, $3); }
    | prod DIV factor  { $$ = newOdeExpOp(ODE_DIV_OP, $1, $3); }
    ;

factor: term             { $$ = $1; }
      | SUB term         { $$ = newOdeExpOp(ODE_NEG, $2, NULL); }
      | term EXP NUMBER  { $$ = newOdeExpOp(ODE_EXP_OP, $1, $3); }
      ;

term: NUMBER                      { $$ = newOdeExpLeaf(ODE_NUM, $1); }
    | IDENT                       { $$ = newOdeExpLeaf(ODE_VAR, $1); }
    | IDENT LPAR sumofprods RPAR  { $$ = newOdeExpTree(ODE_FUN, $1, $3, NULL); }
    | LPAR sumofprods RPAR        { $$ = $2; }
    ;

%%

void yyerror(const char *str) {
  fprintf(stderr, "[line %d] Error: %s\n", odeslineno, str);
}

int parseOdeString(const char *in, ODEList **ptlist) {
  setOdeInputString(in);
  int rv = yyparse();
  endOdeScan();
  *ptlist = finlist;
  return rv;
}
