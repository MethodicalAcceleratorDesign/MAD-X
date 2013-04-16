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
     create constraints content
     print, scan constraints from file
 
 o---------------------------------------------------------------------o
*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>

#include "constraint.h"
#include "utils.h"
#include "error.h"
#include "args.h"

#define T struct constraint
#define S struct slice

// ----- private

static void
printSlc(const S *s, FILE *out)
{
  if (slice_isFull(s)) {
    putc('*', out);
    return;
  }

  fprintf(out, "%u", slice_first(s));

  if (slice_isUnit(s)) return;

  if (slice_isInfinite(s)) {
    putc('-', out);
    putc('$', out);
  } else
   fprintf(out, "-%u", slice_last(s));

  if (slice_stride(s) != 1)
    fprintf(out, "/%u", slice_stride(s));    
}

static int
readSlcOrRng(S *s, FILE *in)
{
  int c, r = 1;
  uint first=0, last=0, stride=1;

  // skip spaces
  while((c = getc(in)) != EOF && isblank(c)) ;
  if (c == EOF) return EOF;

  // ('*'|num)
  if (c == '*') { last = UINT_MAX; goto finish; }
  else {
    ungetc(c, in);
    if (fscanf(in, "%u", &first) != 1) return EOF;
  }

  // (':'|'-')?
  c = getc(in);
       if (c == ':') r = 0;  // slice
  else if (c == '-') ;       // range
  else { ungetc(c, in); last = first; goto finish; }

  // ('$'|num)
  c = getc(in);
  if (c == '$') last = UINT_MAX;
  else {
    ungetc(c, in);
    if (fscanf(in, "%u", &last) != 1) return EOF;
  }

  // ('/'num)? 
  c = getc(in);
  if (c != '/') { ungetc(c, in); stride = 1; goto finish; }
  else
    if (fscanf(in, "%u", &stride) != 1) return EOF;

finish:
  if (r)
    *s = slice_initLastStride(first, last, stride);
  else
    *s = slice_initSizeStride(first, last, stride);

  trace("<-readSlcOrRng %u%c%u/%u", first, r ? '-' : ':', last, stride);

  return 0;
}

static bool
is_reg(int n)
{
  return n > 0 && n < 100;
}

static int
readEps(struct eps *e, FILE *in, int row)
{
  int c = 0, n = 0, rn = 0, cmd = eps_invalid;
  char str[16], buf[16], *end;

  while (1) {
    // parse next constraint
    *str = 0;
    n = fscanf(in, "%*[ \t]%16[^= \t\n\r!#]", str);
    str[sizeof str-1]=0;

    if (n == EOF || *str == 0) break;

         if (strcmp(str, "skip") == 0) {
      cmd |= eps_skip;  trace("[%d] skip", row);
    }
    else if (strcmp(str, "ign") == 0) {
      cmd |= eps_ign;   trace("[%d] ign", row);
    }
    else if (strcmp(str, "istr") == 0) {
      cmd |= eps_istr;  trace("[%d] istr", row);
    }
    else if (strcmp(str, "equ") == 0) {
      cmd |= eps_equ;   trace("[%d] equ", row);
    }
    else if (strcmp(str, "any") == 0) {
      cmd |= eps_any;   trace("[%d] any", row);
    }
    else if (strcmp(str, "all") == 0) {
      cmd &= ~eps_any;  trace("[%d] all", row);
    }
    else if (strcmp(str, "large") == 0) {
      cmd |= eps_large; trace("[%d] large", row);
    }
    else if (strcmp(str, "small") == 0) {
      cmd &= ~eps_large; trace("[%d] small", row);
    }
    else if (strcmp(str, "trace") == 0) {
      cmd |= eps_trace; trace("[%d] trace", row);
    }

// numeric constraints
    else if (strcmp(str, "scl") == 0 && (n = fscanf(in, "=%lf", &e->scl)) == 1) {
                        trace("[%d] scl=%g", row, e->scl);
      ensure(e->scl != 0.0, "invalid zero error scale (%s:%d)", option.cfg_file, row);
    }
    else if (strcmp(str, "scl") == 0 && (n = fscanf(in, "reg%d", &rn)) == 1) {
                        trace("[%d] scl=reg%d", row, rn); e->scl_reg = rn;
      ensure(is_reg(rn), "invalid register number (%s:%d)", option.cfg_file, row);
    }

    else if (strcmp(str, "abs") == 0 && (n = fscanf(in, "=%lf", &e->abs)) == 1) {
      cmd |= eps_abs;   trace("[%d] abs=%g", row, e->abs);      e->_abs = -e->abs;
      ensure(e->abs >= 0.0 && (cmd & eps_large || e->abs <= 1.0),
             "invalid absolute constraint (%s:%d)", option.cfg_file, row);
    }
    else if (strcmp(str, "abs") == 0 && (n = fscanf(in, "reg%d", &rn)) == 1) {
      cmd |= eps_abs;   trace("[%d] abs=reg%d", row, rn);  e->abs_reg = rn; e->_abs_reg = -rn;
      ensure(is_reg(rn), "invalid register number (%s:%d)", option.cfg_file, row);
    }

    else if (strcmp(str, "rel") == 0 && (n = fscanf(in, "=%lf", &e->rel)) == 1) {
      cmd |= eps_rel;   trace("[%d] rel=%g", row, e->rel);      e->_rel = -e->rel;
      ensure(e->rel >= 0.0 && (cmd & eps_large || e->rel <= 1.0),
             "invalid relative constraint (%s:%d)", option.cfg_file, row);
    }
    else if (strcmp(str, "rel") == 0 && (n = fscanf(in, "reg%d", &rn)) == 1) {
      cmd |= eps_rel;   trace("[%d] rel=reg%d", row, rn);  e->rel_reg = rn; e->_rel_reg = -rn;
      ensure(is_reg(rn), "invalid register number (%s:%d)", option.cfg_file, row);
    }

    else if (strcmp(str, "dig") == 0 && (n = fscanf(in, "=%lf", &e->dig)) == 1) {
      cmd |= eps_dig;   trace("[%d] dig=%g", row, e->dig);      e->_dig = -e->dig;
      ensure(e->dig >= 1.0, "invalid digital error (%s:%d)", option.cfg_file, row);
    }
    else if (strcmp(str, "dig") == 0 && (n = fscanf(in, "reg%d", &rn)) == 1) {
      cmd |= eps_dig;   trace("[%d] dig=reg%d", row, rn);  e->dig_reg = rn; e->_dig_reg = -rn;
      ensure(is_reg(rn), "invalid register number (%s:%d)", option.cfg_file, row);
    }

    else if (strcmp(str, "-abs") == 0 && (n = fscanf(in, "=%lf", &e->_abs)) == 1) {
      cmd |= eps_abs;   trace("[%d] -abs=%g", row, e->_abs);
      ensure(e->_abs <= 0.0 && (cmd & eps_large || e->_abs >= -1.0),
             "invalid absolute constraint (%s:%d)", option.cfg_file, row);
    }
    else if (strcmp(str, "-abs") == 0 && (n = fscanf(in, "reg%d", &rn)) == 1) {
      cmd |= eps_abs;   trace("[%d] -abs=reg%d", row, rn);  e->_abs_reg = rn;
      ensure(is_reg(rn), "invalid register number (%s:%d)", option.cfg_file, row);
    }

    else if (strcmp(str, "-rel") == 0 && (n = fscanf(in, "=%lf", &e->_rel)) == 1) {
      cmd |= eps_rel;   trace("[%d] -rel=%g", row, e->_rel);
      ensure(e->_rel <= 0.0 && (cmd & eps_large || e->_rel >= -1.0),
             "invalid relative constraint (%s:%d)", option.cfg_file, row);
    }
    else if (strcmp(str, "-rel") == 0 && (n = fscanf(in, "reg%d", &rn)) == 1) {
      cmd |= eps_rel;   trace("[%d] -rel=reg%d", row, rn);  e->_rel_reg = rn;
      ensure(is_reg(rn), "invalid register number (%s:%d)", option.cfg_file, row);
    }

    else if (strcmp(str, "-dig") == 0 && (n = fscanf(in, "=%lf", &e->_dig)) == 1) {
      cmd |= eps_dig;   trace("[%d] -dig=%g", row, e->_dig);
      ensure(e->_dig <= -1.0, "invalid digital error (%s:%d)", option.cfg_file, row);
    }
    else if (strcmp(str, "-dig") == 0 && (n = fscanf(in, "reg%d", &rn)) == 1) {
      cmd |= eps_dig;   trace("[%d] -dig=reg%d", row, rn);  e->_dig_reg = rn;
      ensure(is_reg(rn), "invalid register number (%s:%d)", option.cfg_file, row);
    }

// actions
    else if (strcmp(str, "omit") == 0 && (n = fscanf(in, "=%*['\"]%64[^'\"]%*['\"]", e->tag)) == 1) {
      e->tag[sizeof e->tag-1] = 0;
      cmd |= eps_omit | eps_equ;
      trace("[%d] omit='%s'", row, e->tag);
      ensure(*e->tag, "invalid empty tag (%s:%d)", option.cfg_file, row);
      ensure(!(cmd & (eps_goto | eps_gonum)), "omit tag conflicting with goto (%s:%d)", option.cfg_file, row);
    }
    else if (strcmp(str, "goto") == 0 && (n = fscanf(in, "=%*['\"]%64[^'\"]%*['\"]", e->tag)) == 1) {
      e->tag[sizeof e->tag-1] = 0; rn = 0;
      int is_num = (strtod(e->tag, &end), !*end) || (sscanf(e->tag, "reg%d", &rn) == 1 && is_reg(rn));
      cmd |= is_num ? eps_goto : eps_gonum | eps_istr;
      if (is_num && is_reg(rn)) e->gto_reg = rn;
      trace("[%d] goto='%s'%s", row, e->tag, cmd & eps_gonum ? " (num)" : "");
      ensure(*e->tag, "invalid empty tag (%s:%d)", option.cfg_file, row);
      ensure(!(cmd & eps_omit), "goto tag conflicting with omit (%s:%d)", option.cfg_file, row);
    }

// indirections
    else if (strncmp(str, "reg", 3) == 0 && (n = fscanf(in, "=%16[^ \t\n\r!#]", buf)) == 1) {
      buf[sizeof buf-1] = 0;
      rn = strtol(str+3, &end, 10);
      trace("[%d] reg%d=%s", row, rn, buf);
      ensure(is_reg(rn) && !*end, "invalid register number (%s:%d)", option.cfg_file, row);

           if (strcmp (buf, "lhs"   ) == 0) { e->lhs_reg = rn; cmd |= eps_lhs; }
      else if (strcmp (buf, "rhs"   ) == 0) { e->rhs_reg = rn; cmd |= eps_rhs; }
      else if (strncmp(buf, "reg", 3) == 0) { e->dst_reg = rn; cmd |= eps_cpy;
        int rn_p = 0, rn_q = 0;
        int k = sscanf(buf, "reg%d-%d", &rn_p, &rn_q);
        ensure((k == 1 && is_reg(rn_p) && rn_q == 0   ) ||
               (k == 2 && is_reg(rn_p) && is_reg(rn_q) && rn_p <= rn_q),
               "invalid registers range (%s:%d)", option.cfg_file, row);
        e->src_reg = rn_p;
        e->cnt_reg = rn_q ? rn_q-rn_p+1 : 1;
      }
      else error("invalid register assigment (%s:%d)", option.cfg_file, row);
    }

// invalid command
    else {
      cmd = eps_invalid;
      trace("[%d] invalid '%s'", row, str);
      break;
    }

    // next char
    ungetc((c = getc(in)), in);
    if (c == EOF || (isspace(c) && !isblank(c)) || c == '#' || c == '!') break; 
  }

  // cleanup non-persistant flags (e.g. large)
  e->cmd = (enum eps_cmd)(cmd & eps_mask);  // cast needed because of icc spurious warnings

  trace("<-readEps cmd = %d, str = '%s', c = '%c'", cmd, str, c);

  return cmd == eps_invalid || n == EOF ? EOF : 0;
}

// ----- interface

void
constraint_print(const T* cst, FILE *out)
{
  if (!out) out = stdout;
  if (!cst) { fprintf(out, "(null)"); return; }

  printSlc(&cst->row, out);
  putc(' ', out);
  printSlc(&cst->col, out);
  putc(' ', out);

  if (cst->eps.cmd & eps_any)    fprintf(out, "any ");
  if (cst->eps.cmd & eps_equ)    fprintf(out, "equ ");
  if (cst->eps.cmd & eps_ign)    fprintf(out, "ign ");
  if (cst->eps.cmd & eps_istr)   fprintf(out, "istr ");
  if (cst->eps.cmd & eps_skip)   fprintf(out, "skip ");
  if (cst->eps.cmd & eps_trace)  fprintf(out, "trace ");

  if (cst->eps.cmd & eps_omit)   fprintf(out, "omit='%s' ", cst->eps.tag);
  if (cst->eps.cmd & eps_goto)   fprintf(out, "goto='%s' ", cst->eps.tag);
  if (cst->eps.cmd & eps_gonum)  fprintf(out, "goto='%s' (num) ", cst->eps.tag);

  if (cst->eps.scl != 1.0)       fprintf(out, "scl=%g ", cst->eps.scl);
  if (cst->eps.scl_reg)          fprintf(out, "scl=reg%d ", cst->eps.scl_reg);

  if (cst->eps.cmd & eps_abs) {
    if (cst->eps.abs_reg) fprintf(out, "abs=reg%d ", cst->eps.abs_reg);
    else fprintf(out, cst->eps.abs == DBL_MIN ? "abs=eps " : "%sabs=%g ",
                      cst->eps.abs > 1        ? "large " : "", cst->eps.abs);

    if (cst->eps._abs_reg && cst->eps._abs_reg != -cst->eps.abs_reg)
      fprintf(out, "-abs=%sreg%d ", cst->eps._abs_reg < 0 ? "-" : "",
                                    cst->eps._abs_reg < 0 ? -cst->eps._abs_reg : cst->eps._abs_reg); 
    else if (cst->eps._abs != -cst->eps.abs)
        fprintf(out, cst->eps._abs == -DBL_MIN ? "-abs=-eps " : "%s-abs=%g ",
                     cst->eps._abs < -1        ? "large " : "", cst->eps._abs);
  }

  if (cst->eps.cmd & eps_rel) {
    if (cst->eps.rel_reg) fprintf(out, "rel=reg%d ", cst->eps.rel_reg);
    else fprintf(out, "%srel=%g ", cst->eps.rel > 1 ? "large " : "", cst->eps.rel);

    if (cst->eps._rel_reg && cst->eps._rel_reg != -cst->eps.rel_reg)
      fprintf(out, "-rel=%sreg%d ", cst->eps._rel_reg < 0 ? "-" : "",
                                    cst->eps._rel_reg < 0 ? -cst->eps._rel_reg : cst->eps._rel_reg); 
    else if (cst->eps._rel != -cst->eps.rel)
      fprintf(out, "%s-rel=%g ", cst->eps._rel < -1 ? "large " : "", cst->eps._rel);
  }

  if (cst->eps.cmd & eps_dig) {
    if (cst->eps.dig_reg) fprintf(out, "dig=reg%d ", cst->eps.dig_reg);
    else fprintf(out, "dig=%g ", cst->eps.dig);

    if (cst->eps._dig_reg && cst->eps._dig_reg != -cst->eps.dig_reg)
      fprintf(out, "-dig=%sreg%d ", cst->eps._dig_reg < 0 ? "-" : "",
                                    cst->eps._dig_reg < 0 ? -cst->eps._dig_reg : cst->eps._dig_reg); 
    else if (cst->eps._dig != -cst->eps.dig)
      fprintf(out, "-dig=%g ", cst->eps._dig);
  }   

  if (cst->eps.cmd & eps_lhs)  fprintf(out, "reg%d=lhs ", cst->eps.lhs_reg);
  if (cst->eps.cmd & eps_rhs)  fprintf(out, "reg%d=rhs ", cst->eps.rhs_reg);
  if (cst->eps.cmd & eps_cpy) {
    if (cst->eps.cnt_reg > 1)
      fprintf(out, "reg%d=reg%d-%d ", cst->eps.dst_reg, cst->eps.src_reg, cst->eps.src_reg+cst->eps.cnt_reg);
    else
      fprintf(out, "reg%d=reg%d ", cst->eps.dst_reg, cst->eps.src_reg);
  }
}

void
constraint_scan(T* cst, FILE *in, int *row)
{
  int c;
  assert(cst && row);

  *cst = (T){ .eps = { .cmd = eps_invalid, .scl=1.0 } };

  if (!in) in = stdin;

retry:

  while((c = getc(in)) != EOF && isblank(c)) ;

  // end of file
  if (c == EOF) return;

  ungetc(c, in);

  // comment or empty line
  if (c == '\n' || c == '\r' || c == '#' || c == '!') {
    if (skipLine(in, 0) == '\n') ++*row;
    goto retry;
  }

  cst->idx  = -1;
  cst->line = *row;

  ensure(readSlcOrRng(&cst->row, in      ) != EOF, "invalid row range (%s:%d)"   , option.cfg_file, *row);
  ensure(readSlcOrRng(&cst->col, in      ) != EOF, "invalid column range (%s:%d)", option.cfg_file, *row);
  ensure(readEps     (&cst->eps, in, *row) != EOF, "invalid constraint or command (%s:%d)", option.cfg_file, *row);

  // expand to all columns
  if (cst->eps.cmd & eps_skip || cst->eps.cmd & eps_goto)
    cst->col = slice_initAll();

  // adjust row count
  if (skipLine(in, 0) == '\n') ++*row;
}

