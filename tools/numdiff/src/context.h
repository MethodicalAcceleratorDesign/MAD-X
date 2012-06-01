#ifndef CONTEXT_H
#define CONTEXT_H

#include <stdio.h>
#include "types.h"

// ----- types

struct utest;
struct context;
struct constraint;

// ----- interface

#define T struct context
#define C struct constraint

T*       context_alloc  (int n_);
void     context_clear  (T*); // reset + erase constraints
void     context_free   (T*);

// populate with constraints, invoke grow on need
T*       context_add     (T*, const C*);

// return 0 if no constraint are found, getInc requires increasing (row,col)
const C* context_getAt   (T*, int row, int col);
const C* context_getInc  (T*, int row, int col);

// return the contraint at the index
const C* context_getIdx  (const T*, int idx);
// return the index of the contraint
int      context_findIdx (const T*, const C*); // -1 for invalid constraint

// input/output context of (registered) constraints
T*       context_scan (      T*, FILE *fp);
void     context_print(const T*, FILE *fp); // for debug

#undef T
#undef C

// ----- testsuite

#ifndef NTEST

void context_utest (struct utest*);

#endif // NTEST
#endif
