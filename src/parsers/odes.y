/* Parser for a system of ODEs */

%define api.value.type {char*}
%define parse.lac full
%define parse.error detailed
%define api.prefix {odes}
%locations
%param { yyscan_t scanner }

%code top {
  #include <stdio.h>
}
%code requires {
  typedef void* yyscan_t;
}
%code {
  int yylex(YYSTYPE* yylvalp, YYLTYPE* yyllocp, yyscan_t scanner);
  void yyerror(YYLTYPE* yyllocp, yyscan_t unused, const char* msg);
}

/* tokens that will be used */
%token NUMBER IDENT SCOLON LPAR RPAR UNKNOWN
%token ADD SUB MUL DIV EXP EQUAL PRIME

%%

odelist: /* epsilon */
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

int main (int argc, char **argv) {
    yyparse();
}

void yyerror (char const *s) {
  fprintf (stderr, "%s\n", s);
}
