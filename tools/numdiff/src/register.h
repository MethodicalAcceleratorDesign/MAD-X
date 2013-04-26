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

#include <stdlib.h>
#include <math.h>
#include "error.h"

// ----- constants

#define REG_MAX (1 << 13)
#define REG_UNARY_OP  "-/\\^|[]"
#define REG_BINARY_OP "+-*/%^<>~"

// ----- interface

static inline bool
reg_isvalid(short rn)
{
  return rn > 0 && rn < REG_MAX;
}

static inline short
reg_encode(short rn, char op)
{
  ensure(reg_isvalid(rn), "invalid register number %d", rn);

  switch(op) {
  case  0  : return 0 * REG_MAX + rn;
  case '-' : return 1 * REG_MAX + rn;
  case '/' : return 2 * REG_MAX + rn;
  case '\\': return 3 * REG_MAX + rn;
  case '^' : return 4 * REG_MAX + rn;
  case '|' : return 5 * REG_MAX + rn;
  case '[' : return 6 * REG_MAX + rn;
  case ']' : return 7 * REG_MAX + rn;
  default: error("invalid register unary operation '%c'", op); exit(EXIT_FAILURE);
  }
}

static inline short
reg_decode(short rn, char *op)
{
  switch(rn / REG_MAX) {
  case 0: if (op) *op =  0  ; break;
  case 1: if (op) *op = '-' ; break;
  case 2: if (op) *op = '/' ; break;
  case 3: if (op) *op = '\\'; break;
  case 4: if (op) *op = '^' ; break;
  case 5: if (op) *op = '|' ; break;
  case 6: if (op) *op = '[' ; break;
  case 7: if (op) *op = ']' ; break;
  default: error("invalid register unary operation '%c'", rn/REG_MAX);
  }
  rn &= REG_MAX-1;
  ensure(reg_isvalid(rn), "invalid register number %d", rn);
  return rn;
}

#ifdef __GNUC__
__attribute__((always_inline))
#endif
static inline double
reg_getval(const double *reg, short reg_n, short rn)
{
  short r = rn & (REG_MAX-1);
  ensure(r > 0 && r <= reg_n, "invalid register number %d", r);

  switch(rn / REG_MAX) {
  case 0: return reg[r-1];
  case 1: return -reg[r-1];
  case 2: return 1/reg[r-1];
  case 3: return -1/reg[r-1];
  case 4: return exp(reg[r-1]);
  case 5: return fabs(reg[r-1]);
  case 6: return floor(reg[r-1]);
  case 7: return ceil(reg[r-1]);
  default: error("invalid register unary operation '%c'", rn/REG_MAX); exit(EXIT_FAILURE);
  }
}

static inline void
reg_setval(double *reg, short reg_n, short rn, double val)
{
  ensure(rn > 0 && rn <= reg_n, "invalid register number %d", rn);
  reg[rn-1] = val;
}

#ifdef __GNUC__
__attribute__((always_inline))
#endif
static inline void
reg_eval(double *reg, short reg_n, short dst, short src, short src2, char op)
{
  ensure(dst > 0 && dst <= reg_n, "invalid register %d", dst);

  if (!op) { // nop = getval
    reg[dst-1] = reg_getval(reg, reg_n, src);
  }
  else {
    ensure(src  > 0 && src  <= reg_n, "invalid register %d", src);
    ensure(src2 > 0 && src2 <= reg_n, "invalid register %d", src2);

    switch(op) {
    case '+': reg[dst-1] = reg[src-1] + reg[src2-1]; break;
    case '-': reg[dst-1] = reg[src-1] - reg[src2-1]; break;
    case '*': reg[dst-1] = reg[src-1] * reg[src2-1]; break;
    case '/': reg[dst-1] = reg[src-1] / reg[src2-1]; break;
    case '%': reg[dst-1] = fmod(reg[src-1], reg[src2-1]); break;
    case '^': reg[dst-1] = pow(reg[src-1], reg[src2-1]); break;
    case '<': reg[dst-1] = reg[src-1] < reg[src2-1] ? reg[src-1] : reg[src2-1]; break;
    case '>': reg[dst-1] = reg[src-1] > reg[src2-1] ? reg[src-1] : reg[src2-1]; break;
    case '~':
      ensure(src < src2, "invalid range of registers %d~%d", src, src2);

      if (dst <= src) // take care of overlapping registers
        for (short i=0; i <= src2-src; i++)    
          reg[dst+i-1] = reg[src+i-1];
      else
        for (short i=src2-src; i >= 0; i--)    
          reg[dst+i-1] = reg[src+i-1];
      break;

    default:
      error("invalid register operation '%c'", op);
    }
  }
}

#endif

