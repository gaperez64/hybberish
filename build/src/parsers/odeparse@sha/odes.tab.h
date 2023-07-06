/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

#ifndef YY_ODES_SRC_PARSERS_ODEPARSE_SHA_ODES_TAB_H_INCLUDED
# define YY_ODES_SRC_PARSERS_ODEPARSE_SHA_ODES_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef ODESDEBUG
# if defined YYDEBUG
#if YYDEBUG
#   define ODESDEBUG 1
#  else
#   define ODESDEBUG 0
#  endif
# else /* ! defined YYDEBUG */
#  define ODESDEBUG 0
# endif /* ! defined YYDEBUG */
#endif  /* ! defined ODESDEBUG */
#if ODESDEBUG
extern int odesdebug;
#endif

/* Token type.  */
#ifndef ODESTOKENTYPE
# define ODESTOKENTYPE
  enum odestokentype
  {
    SCOLON = 258,
    LPAR = 259,
    RPAR = 260,
    UNKNOWN = 261,
    ADD = 262,
    SUB = 263,
    MUL = 264,
    DIV = 265,
    EXP = 266,
    EQUAL = 267,
    PRIME = 268,
    NUMBER = 269,
    IDENT = 270
  };
#endif

/* Value type.  */
#if ! defined ODESSTYPE && ! defined ODESSTYPE_IS_DECLARED

union ODESSTYPE
{
#line 20 "../src/parsers/odes.y" /* yacc.c:1909  */

  ODEList *list;
  ODEExpTree *tree;
  char *str;

#line 84 "src/parsers/odeparse@sha/odes.tab.h" /* yacc.c:1909  */
};

typedef union ODESSTYPE ODESSTYPE;
# define ODESSTYPE_IS_TRIVIAL 1
# define ODESSTYPE_IS_DECLARED 1
#endif


extern ODESSTYPE odeslval;

int odesparse (void);

#endif /* !YY_ODES_SRC_PARSERS_ODEPARSE_SHA_ODES_TAB_H_INCLUDED  */
