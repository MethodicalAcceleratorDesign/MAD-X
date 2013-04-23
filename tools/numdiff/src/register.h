#ifndef REGISTER_H
#define REGISTER_H

/*
 o---------------------------------------------------------------------o
 |
 | Numdiff
 |
 | Copyright (c) 2012+ laurent.deniau@cern.ch
 | Gnu General Public License
 |
 o---------------------------------------------------------------------o
  
   Purpose:
     manage access to registers
     handle basic operations (-, /)
 
 o---------------------------------------------------------------------o
*/

#include "error.h"

// ----- constants

#define REG_MAX 16000

// ----- interface

static inline bool
reg_isvalid(int n)
{
  return n > 0 && n < REG_MAX;
}

static inline int
reg_encode(int n, int s, int i)
{
  ensure(reg_isvalid(n), "invalid register number %d", n);
  if (i) { n += REG_MAX; }
  if (s) { n  = -n;      }
  return n;
}

static inline int
reg_decode(int n, int *s, int *i)
{
  if (n < 0      ) { n  = -n;      if (s) *s = 1; }
  if (n > REG_MAX) { n -= REG_MAX; if (i) *i = 1; }
  ensure(reg_isvalid(n), "invalid register number %d", n);
  return n;
}

static inline double
reg_getval(const double *reg, int reg_n, int n, double val)
{
  if (n) {
    int s = 0, i = 0;

    if (n < 0      ) { s = 1; n  = -n;      }
    if (n > REG_MAX) { i = 1; n -= REG_MAX; }

    ensure(n > 0 && n <= reg_n, "invalid register number %d", n);
    val = reg[n-1];

    if (i) val = 1/val;
    if (s) val =  -val;
  }
  return val;
}

static inline void
reg_setval(double *reg, int reg_n, int n, double val)
{
  ensure(n > 0 && n <= reg_n, "invalid register number %d", n);
  reg[n-1] = val;
}

static inline void
reg_copy(double *reg, int reg_n, int from, int to)
{
    ensure(to   > 0 && to   <= reg_n, "invalid register %d", to  );
    ensure(from > 0 && from <= reg_n, "invalid register %d", from);
    reg[to-1] = reg_getval(reg, reg_n, from, 0);
}

static inline void
reg_move(double *reg, int reg_n, int from, int to, int cnt)
{
  if (!cnt)
    reg_copy(reg, reg_n, from, to);

  else {
    ensure(to   > 0 && to  +cnt <= reg_n, "invalid register range %d-%d", to  , to  +cnt);
    ensure(from > 0 && from+cnt <= reg_n, "invalid register range %d-%d", from, from+cnt);

    if (to <= from) // take care of overlapping registers
      for (int i=0; i < cnt; i++)    
        reg[to+i-1] = reg_getval(reg, reg_n, from+i, 0);
    else
      for (int i=cnt-1; i >= 0; i--)    
        reg[to+i-1] = reg_getval(reg, reg_n, from+i, 0);
  }
}

#endif

