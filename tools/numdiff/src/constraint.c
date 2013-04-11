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

// ----- private

static void
printSlc(const struct slice *s, FILE *out)
{
  if (slice_first(s) <= 1 && slice_isInfinite(s) && slice_isDense(s)) {
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
readSlcOrRng(struct slice *s, FILE *in)
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

static int
readEps(struct eps *e, FILE *in, int row)
{
  int c = 0, n = 0, cmd = eps_invalid;
  char str[16], *end;

  while (1) {
    // parse next constraint
    *str = 0;
    n = fscanf(in, "%*[ \t]%10[^= \t\n\r!#]", str);
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
    else if (strcmp(str, "abs") == 0 && (n = fscanf(in, "=%lf", &e->abs)) == 1) {
      cmd |= eps_abs;   trace("[%d] abs=%g", row, e->abs);      e->_abs = -e->abs;
      ensure(e->abs >= 0.0 && (cmd & eps_large || e->abs <= 1.0),
             "invalid absolute constraint (%s:%d)", option.cfg_file, row);
    }
    else if (strcmp(str, "rel") == 0 && (n = fscanf(in, "=%lf", &e->rel)) == 1) {
      cmd |= eps_rel;   trace("[%d] rel=%g", row, e->rel);      e->_rel = -e->rel;
      ensure(e->rel >= 0.0 && (cmd & eps_large || e->rel <= 1.0),
             "invalid relative constraint (%s:%d)", option.cfg_file, row);
    }
    else if (strcmp(str, "dig") == 0 && (n = fscanf(in, "=%lf", &e->dig)) == 1) {
      cmd |= eps_dig;   trace("[%d] dig=%g", row, e->dig);      e->_dig = -e->dig;
      ensure(e->dig >= 1.0, "invalid digital error (%s:%d)", option.cfg_file, row);
    }

    else if (strcmp(str, "-abs") == 0 && (n = fscanf(in, "=%lf", &e->_abs)) == 1) {
      cmd |= eps_abs;   trace("[%d] -abs=%g", row, e->_abs);
      ensure(e->_abs <= 0.0 && (cmd & eps_large || e->_abs >= -1.0),
             "invalid absolute constraint (%s:%d)", option.cfg_file, row);
    }
    else if (strcmp(str, "-rel") == 0 && (n = fscanf(in, "=%lf", &e->_rel)) == 1) {
      cmd |= eps_rel;   trace("[%d] -rel=%g", row, e->_rel);
      ensure(e->_rel <= 0.0 && (cmd & eps_large || e->_rel >= -1.0),
             "invalid relative constraint (%s:%d)", option.cfg_file, row);
    }
    else if (strcmp(str, "-dig") == 0 && (n = fscanf(in, "=%lf", &e->_dig)) == 1) {
      cmd |= eps_dig;   trace("[%d] -dig=%g", row, e->_dig);
      ensure(e->_dig <= -1.0, "invalid digital error (%s:%d)", option.cfg_file, row);
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
      e->tag[sizeof e->tag-1] = 0;
      strtod(e->tag, &end);
      cmd |= *end ? eps_goto : eps_gonum | eps_istr; 
      trace("[%d] goto='%s'%s", row, e->tag, cmd & eps_gonum ? " (num)" : "");
      ensure(*e->tag, "invalid empty tag (%s:%d)", option.cfg_file, row);
      ensure(!(cmd & eps_omit), "goto tag conflicting with omit (%s:%d)", option.cfg_file, row);
    }
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

  if (cst->eps.cmd & eps_abs) {
    fprintf(out, cst->eps.abs == DBL_MIN ? "abs=eps " : "%sabs=%g ",
                 cst->eps.abs > 1 ? "large " : "", cst->eps.abs);
    if (cst->eps._abs != -cst->eps.abs)
      fprintf(out, cst->eps._abs == -DBL_MIN ? "-abs=-eps " : "%s-abs=%g ",
                   cst->eps._abs < -1 ? "large " : "", cst->eps._abs);
  }
  if (cst->eps.cmd & eps_rel) {
    fprintf(out, "%srel=%g ", cst->eps.rel > 1 ? "large " : "", cst->eps.rel);
    if (cst->eps._rel != -cst->eps.rel)
      fprintf(out, "%s-rel=%g ", cst->eps._rel < -1 ? "large " : "", cst->eps._rel);
  }
  if (cst->eps.cmd & eps_dig) {
    fprintf(out, "dig=%g ", cst->eps.dig);
    if (cst->eps._dig != -cst->eps.dig)
      fprintf(out, "-dig=%g ", cst->eps._dig);
  }   
}

void
constraint_scan(T* cst, FILE *in, int *row)
{
  int c;
  assert(cst && row);

  *cst = (T){ .eps = { .cmd = eps_invalid } };

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

  cst->line = *row;
  ensure(readSlcOrRng(&cst->row, in      ) != EOF, "invalid row range (%s:%d)"   , option.cfg_file, *row);
  ensure(readSlcOrRng(&cst->col, in      ) != EOF, "invalid column range (%s:%d)", option.cfg_file, *row);
  ensure(readEps     (&cst->eps, in, *row) != EOF, "invalid constraint or command (%s:%d)", option.cfg_file, *row);
  if (skipLine(in, 0) == '\n') ++*row;
}

