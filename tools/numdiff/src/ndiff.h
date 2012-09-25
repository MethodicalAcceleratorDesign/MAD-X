#ifndef NDIFF_H
#define NDIFF_H

/*
 o---------------------------------------------------------------------o
 |
 | Numdiff
 |
 | Copyright (c) 2012+ Laurent Deniau
 |
 o---------------------------------------------------------------------o
  
   Purpose:
     numerical diff of files
     provides the main numdiff loop
 
 o---------------------------------------------------------------------o
*/

#include <stdio.h>
#include "types.h"

// ----- types

struct utest;
struct ndiff;
struct context;
struct constraint;

// ----- interface

#define T struct ndiff

T*    ndiff_alloc  (FILE *lhs, FILE *rhs, int n_);
void  ndiff_clear (T*);
void  ndiff_free  (T*);

// high level API
void  ndiff_loop     (T*, struct context*, int blank, int check);

// low level API
int   ndiff_skipLine (T*);
int   ndiff_readLine (T*);
int   ndiff_gotoLine (T*, const char *tag);
int   ndiff_fillLine (T*, const char *lhs, const char *rhs);

void  ndiff_diffLine (T*, int blank);
int   ndiff_nextNum  (T*, int blank); // return 0 if no number is found
int   ndiff_testNum  (T*, const struct context*, const struct constraint*);

int   ndiff_maxDisp  (T*, int max);   // return the previous value
void  ndiff_getInfo  (const T*, int *row_, int *col_, int *cnt_);
int   ndiff_feof     (const T*, int both);
int   ndiff_isempty  (const T*);

#undef T

// ----- testsuite

#ifndef NTEST

void ndiff_utest (struct utest*);

#endif // NTEST
#endif
