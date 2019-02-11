/*
**  EVAL.C - A simple mathematical expression evaluator in C
**
**  operators supported: (
**                       )
**                       +
**                       -
**                       *
**                       /
**                       ^
**
**  limitations: 1 - No precedence rules are implemented.
**               2 - Numbers can be negated (e.g. "-13"), but not
**                   expressions (e.g. "-(13)").
**
**  Original Copyright 1991 by Bob Stout as part of
**  the MicroFirm Function Library (MFL)
**
**  This subset* version is hereby donated to the public domain.
**
**  *(The MFL version adds 150 lines of code, 5 level precedence,
**    logarithmic and transcendental operators, pi as a constant,
**    named variables, and fully understands negation.)
*/

/* 
 * Adapté sous la forme d'un object cpp: la forme originale utilisait 
 * des tableaux et variables globaux, ce qui en faisait un évaluateur 
 * non thread safe. 
 */

#ifndef EVAL_H
#define EVAL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define NUL '\0'

enum STATUS {R_ERROR = -2 /* range */, ERROR /* syntax */, SUCCESS};

extern char delims[];   /* Tokens */

class eval
{
 public:
  eval(): op_sptr(0), arg_sptr(0), parens(0) {}
  int evaluate(char *, double *);

 private:
  void       strupr(char *s);
  int        do_op(void);
  int        do_paren(void);
  void       push_op(char);
  void       push_arg(double);
  STATUS     pop_arg(double *);
  STATUS     pop_op(int *);
  char      *getexp(char *);
  char      *getop(char *);
  void       pack(char *);

  char       op_stack[256];  /* Operator stack       */
  double     arg_stack[256]; /* Argument stack       */
  char       token[256];     /* Token buffer         */
  int        op_sptr;        /* op_stack pointer     */
  int        arg_sptr;       /* arg_stack pointer    */
  int        parens;         /* Nesting level        */
  int        state ;         /* 0 = Awaiting expression
                                1 = Awaiting operator
                              */
};

#endif /* EVAL_H */
