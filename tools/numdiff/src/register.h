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
     handle arithmetic operations (+, -, *, /) and range (~)
 
 o---------------------------------------------------------------------o
*/

#include "error.h"

// ----- constants

#define REG_MAX 16000

// ----- interface

static inline bool
reg_isvalid(int rn)
{
  return rn > 0 && rn < REG_MAX;
}

static inline int
reg_encode(int rn, int s, int i)
{
  ensure(reg_isvalid(rn), "invalid register number %d", rn);
  if (i) { rn += REG_MAX; }
  if (s) { rn  = -rn;     }
  return rn;
}

static inline int
reg_decode(int rn, int *s, int *i)
{
  if (rn < 0      ) { rn  = -rn;     if (s) *s = 1; }
  if (rn > REG_MAX) { rn -= REG_MAX; if (i) *i = 1; }
  ensure(reg_isvalid(rn), "invalid register number %d", rn);
  return rn;
}

static inline double
reg_getval(const double *reg, int reg_n, int rn, double val)
{
  if (rn) {
    int s = 0, i = 0;

    if (rn < 0      ) { s = 1; rn  = -rn;     }
    if (rn > REG_MAX) { i = 1; rn -= REG_MAX; }

    ensure(rn > 0 && rn <= reg_n, "invalid register number %d", rn);
    val = reg[rn-1];

    if (i) val = 1/val;
    if (s) val =  -val;
  }
  return val;
}

static inline void
reg_setval(double *reg, int reg_n, int rn, double val)
{
  ensure(rn > 0 && rn <= reg_n, "invalid register number %d", rn);
  reg[rn-1] = val;
}

static inline void
reg_eval(double *reg, int reg_n, int dst, int src, int src2, int op)
{
  ensure(dst > 0 && dst < reg_n, "invalid register %d", dst);

  if (!op) { // nop = copy
    reg[dst] = reg_getval(reg, reg_n, src, 0);
  }
  else {
    ensure(src  > 0 && src  < reg_n, "invalid register %d", src);
    ensure(src2 > 0 && src2 < reg_n, "invalid register %d", src2);
    switch(op) {
    case '+': reg[dst] = reg[src] + reg[src2]; break;
    case '-': reg[dst] = reg[src] - reg[src2]; break;
    case '*': reg[dst] = reg[src] * reg[src2]; break;
    case '/': reg[dst] = reg[src] / reg[src2]; break;
    case '~':
      ensure(src < src2, "invalid range of registers %d~%d", src, src2);
      if (dst <= src) // take care of overlapping registers
        for (int i=0; i <= src2-src; i++)    
          reg[dst+i] = reg[src+i];
      else
        for (int i=src2-src; i >= 0; i--)    
          reg[dst+i] = reg[src+i];
      break;
    default:
      error("invalid register operation '%c'", op);
    }
  }
}

#endif

