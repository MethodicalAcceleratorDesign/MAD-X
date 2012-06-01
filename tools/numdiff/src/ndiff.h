#ifndef NDIFF_H
#define NDIFF_H

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
void  ndiff_loop     (T*, struct context*, int debug);

// low level API
void  ndiff_diffLine (T*);

int   ndiff_skipLine (T*);
int   ndiff_readLine (T*);
int   ndiff_fillLine (T*, const char *lhs, const char *rhs);

int   ndiff_nextNum  (T*); // return 0 if no number is found
int   ndiff_testNum  (T*, const struct context*, const struct constraint*);

void  ndiff_getInfo  (const T*, int *row_, int *col_, int *cnt_);
int   ndiff_feof     (const T*);

#undef T

// ----- testsuite

#ifndef NTEST

void ndiff_utest (struct utest*);

#endif // NTEST
#endif
