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
     numerical diff of files
     provides the main numdiff loop
 
 o---------------------------------------------------------------------o
*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>

#include "args.h"
#include "ndiff.h"
#include "constraint.h"
#include "context.h"
#include "error.h"
#include "utils.h"

#define T struct ndiff
#define C struct constraint

// ----- types

struct ndiff {
  // files
  FILE *lhs_f, *rhs_f;
  int   row_i,  col_i; // line, num-column

  // context
  struct context* cxt;

  // registers
  char reg[100][64];

  // options
  int blank, check;

  // diff counter
  int   cnt_i, max_i;

  // numbers counter
  long  num_i;

  // buffers
  int   lhs_i,  rhs_i; // char-columns
  int   buf_s;         // capacity
  char *lhs_b, *rhs_b;
};

// ----- private (parser helpers)

static inline int
is_separator (int c)
{
  return !c || isblank(c) || (ispunct(c) && !strchr(option.chr, c));
}

static inline int
is_number (char *buf)
{
  int i = 0;

  // sign
  if (buf[i] == '-' || buf[i] == '+' || buf[i] == ' ') i++;

  // dot
  if (buf[i] == '.') ++i;

  // digits
  return isdigit(buf[i]);
}

static inline char*
// assume that buf has been validated by is_number
backtrack_number (char *buf, const char *beg)
{
  if (*buf == ' ') ++buf;

  else
  if (*buf == '.') {
    if (buf > beg && (buf[-1] == '-' || buf[-1] == '+')) --buf;
  }

  else
  if (isdigit(*buf)) {
    if (buf > beg &&  buf[-1] == '.') --buf;
    if (buf > beg && (buf[-1] == '-' || buf[-1] == '+')) --buf;
  }

  return buf;
}

static inline int
// assume that buf has been validated by is_number and backtracked
is_number_start(char *buf, const char *beg)
{
  // number is at the beginning or is preceded by a separator
  return *buf == '-' || *buf == '+' || buf == beg || (buf > beg && is_separator(buf[-1]));
}

static inline int
parse_number (char *buf, int *d_, int *n_, int *e_, int *f_)
{
  int i = 0, d = 0, e = 0, n = 0;
  char c;

  // sign
  if (buf[i] == '-' || buf[i] == '+') i++;

  // drop leading zeros
  while(buf[i] == '0') i++;

  // digits
  while(isdigit(buf[i])) n++, i++;

  // dot
  if (buf[i] == '.') d = ++i;

  // decimals
  if (d) {
    // drop leading zeros
    if (!n) while(buf[i] == '0') i++;

    // digits
    while(isdigit(buf[i])) n++, i++;
  }

  // ensure at least ±# or ±#. or ±.#
  if(!(i > 0 && (isdigit(buf[i-1]) || (i > 1 &&  isdigit(buf[i-2])))))
    return 0;

  // exponent
  if (buf[i] == 'e' || buf[i] == 'E' || buf[i] == 'd' || buf[i] == 'D')
    c = buf[i], buf[i] = 'e', e = ++i;

  if (e) {
    // sign
    if (buf[i] == '-' || buf[i] == '+') i++;

    // digits
    while(isdigit(buf[i])) i++;

    // ensure e# or e±# otherwise backtrack
    if (!isdigit(buf[i-1]))
      i = e-1, buf[i] = c, e = 0;
  }

  if (n_) *n_ = n;
  if (d_) *d_ = d-1;
  if (e_) *e_ = e-1;
  if (f_) *f_ = d > 0 || e > 0;

  return i;
}

static inline void
skip_identifier(char *restrict *lhs, char *restrict *rhs, int strict)
{
  if (strict) {
    assert(lhs && rhs);
    while (**lhs == **rhs && !is_separator(**lhs)) ++*lhs, ++*rhs;
  }
  else {
    assert(lhs || rhs);
    if (lhs) while (!is_separator(**lhs)) ++*lhs;
    if (rhs) while (!is_separator(**rhs)) ++*rhs;
  }
}

static inline int
is_valid_omit(const char *lhs_p, const char *rhs_p, const T *dif, const char *tag)
{
  const char *p = tag+strlen(tag);

  while (--p >= tag && --lhs_p >= dif->lhs_b && --rhs_p >= dif->rhs_b)
    if (*p != *lhs_p || *p != *rhs_p) return false;

  return true;
}

// ----- private (ctor & dtor helpers)

static inline void
ndiff_setup (T *dif, int n)
{
  enum { min_alloc = 65536 };

  if (n < min_alloc) n = min_alloc;

  dif->lhs_b = malloc(n * sizeof *dif->lhs_b);
  dif->rhs_b = malloc(n * sizeof *dif->rhs_b);
  ensure(dif->lhs_b && dif->rhs_b, "out of memory");

  *dif = (T) {
    .lhs_f = dif->lhs_f, .rhs_f = dif->rhs_f,
    .lhs_b = dif->lhs_b, .rhs_b = dif->rhs_b,
    .blank = dif->blank, .check = dif->check,
    .cxt = dif->cxt,      
    .buf_s = n,
    .max_i = 25
  };
}

static inline void
ndiff_teardown (T *dif)
{
  free(dif->lhs_b);
  free(dif->rhs_b);

  *dif = (T) {
    .lhs_f = dif->lhs_f, .rhs_f = dif->rhs_f,
    .blank = dif->blank, .check = dif->check,
    .cxt = dif->cxt
  };
}

static inline void
ndiff_grow (T *dif, int n)
{
  if (n > dif->buf_s) { // enlarge on need
    dif->lhs_b = realloc(dif->lhs_b, n * sizeof *dif->lhs_b);
    dif->rhs_b = realloc(dif->rhs_b, n * sizeof *dif->rhs_b);
    ensure(dif->lhs_b && dif->rhs_b, "out of memory");
    dif->buf_s = n;
  }
}

static inline void
ndiff_reset_buf (T *dif)
{
  dif->lhs_i = dif->rhs_i = 0;
  dif->lhs_b[0] = dif->rhs_b[0] = 0;
}

// ----- private (error helpers)

static void
ndiff_error(const struct context *cxt,
            const C *c, const C *c2,
            int row, int col)
{
  warning("dual constraints differ at %d:%d", row, col);
  warning("getIncr select [#%d]", context_findIdx(cxt, c ));
  warning("getAt   select [#%d]", context_findIdx(cxt, c2));
  warning("rules list:");
  context_print(cxt, stderr);
  error("please report to mad@cern.ch");
}

static void
ndiff_header(void)
{
  if (option.test)
    warning("(*) files '" CSTR_RED("%s") "'|'" CSTR_RED("%s") "' from '%s' differ",
            option.lhs_file, option.rhs_file, option.test);
  else
    warning("(*) files '" CSTR_RED("%s") "'|'" CSTR_RED("%s") "' differ",
            option.lhs_file, option.rhs_file);
}

// -----------------------------------------------------------------------------
// ----- interface
// -----------------------------------------------------------------------------

T*
ndiff_alloc (FILE *lhs_f, FILE *rhs_f, struct context *cxt, int n_)
{
  assert(lhs_f && rhs_f);

  T *dif = malloc(sizeof *dif);
  ensure(dif, "out of memory");

  *dif = (T) { .lhs_f = lhs_f, .rhs_f = rhs_f, .cxt = cxt };

  ndiff_setup(dif, n_);
  return dif;
}

void
ndiff_free (T *dif)
{
  assert(dif);
  ndiff_teardown(dif);
  free(dif);
}

void
ndiff_clear (T *dif)
{
  assert(dif);
  ndiff_teardown(dif);
  ndiff_setup(dif, 0);
}

int
ndiff_skipLine (T *dif)
{
  assert(dif);
  int s1 = 0, s2 = 0;
  int c1, c2;

  ndiff_reset_buf(dif);

  c1 = skipLine(dif->lhs_f, &s1);
  c2 = skipLine(dif->rhs_f, &s2);

  dif->col_i  = 0;
  dif->row_i += 1;

  return c1 == EOF || c2 == EOF ? EOF : !EOF;
}

int
ndiff_fillLine (T *dif, const char *lhs_b, const char *rhs_b)
{
  assert(dif);
  assert(lhs_b && rhs_b);

  ndiff_reset_buf(dif);

  int s1 = strlen(lhs_b)+1; 
  int s2 = strlen(rhs_b)+1; 
  ndiff_grow(dif, imax(s1,s2));
  memcpy(dif->lhs_b, lhs_b, s1);
  memcpy(dif->rhs_b, rhs_b, s2);

  dif->col_i  = 0;
  dif->row_i += 1;

  return 0; // never fails
}

int
ndiff_readLine (T *dif)
{
  assert(dif);
  int s1 = 0, s2 = 0;
  int c1, c2, n = 0;

  trace("->readLine line %d", dif->row_i);

  ndiff_reset_buf(dif);

  while (1) {
    c1 = readLine(dif->lhs_f, dif->lhs_b+s1, dif->buf_s-s1, &n); s1 += n;
    c2 = readLine(dif->rhs_f, dif->rhs_b+s2, dif->buf_s-s2, &n); s2 += n;
    if (c1 == '\n' || c2 == '\n' || c1 == EOF || c2 == EOF) break;
    ndiff_grow(dif, 2*dif->buf_s);
  }

  dif->col_i  = 0;
  dif->row_i += 1;

  trace("  buffers: '%.25s'|'%.25s'", dif->lhs_b, dif->rhs_b);
  trace("<-readLine line %d", dif->row_i);

  return c1 == EOF || c2 == EOF ? EOF : !EOF;
}

int
ndiff_gotoLine (T *dif, const C *c)
{
  assert(dif && c);

  int c1=0, c2=0, i1=0, i2=0;

  trace("->gotoLine line %d", dif->row_i);

  // --- lhs ---
  while (1) {
    int s = 0, n = 0;

    dif->lhs_i    = 0;
    dif->lhs_b[0] = 0;

    if (c1 == EOF) break;

    while (1) {
      c1 = readLine(dif->lhs_f, dif->lhs_b+s, dif->buf_s-s, &n); s += n;
      if (c1 == '\n' || c1 == EOF) break;
      ndiff_grow(dif, 2*dif->buf_s);
    }

    i1 += 1;
    trace("  lhs[%d]: '%s'", dif->row_i+i1, dif->lhs_b);

    // search for tag
    if (strstr(dif->lhs_b, c->eps.tag)) break;
  }

  // --- rhs ---
  while (1) {
    int s = 0, n = 0;

    dif->rhs_i    = 0;
    dif->rhs_b[0] = 0;

    if (c2 == EOF) break;

    while (1) {
      c2 = readLine(dif->rhs_f, dif->rhs_b+s, dif->buf_s-s, &n); s += n;
      if (c2 == '\n' || c2 == EOF) break;
      ndiff_grow(dif, 2*dif->buf_s);
    }

    i2 += 1;
    trace("  rhs[%d]: '%s'", dif->row_i+i2, dif->rhs_b);

    // search for tag
    if (strstr(dif->rhs_b, c->eps.tag)) break;
  }

  dif->col_i  = 0;
  dif->row_i += i1 < i2 ? i1 : i2;

  // return with last lhs and rhs lines loaded if tag was found

  trace("  buffers: '%.25s'|'%.25s'", dif->lhs_b, dif->rhs_b);
  trace("<-gotoLine line %d (%+d|%+d)", dif->row_i, i1, i2);

  return c1 == EOF || c2 == EOF ? EOF : !EOF;
}

int
ndiff_gotoNum (T *dif, const C *c)
{
  assert(dif && c);

  if ((c->eps.cmd & eps_equ) && slice_isFull(&c->col))
    return ndiff_gotoLine(dif, c);

  int c1=0, c2=0, i1=0, i2=0;

  trace("->gotoNum line %d", dif->row_i);

  // --- lhs ---
  strcpy(dif->rhs_b, c->eps.tag);

  while (1) {
    int s = 0, n = 0;

    dif->lhs_i    = 0;
    dif->lhs_b[0] = 0;

    if (c1 == EOF) break;

    while (1) {
      c1 = readLine(dif->lhs_f, dif->lhs_b+s, dif->buf_s-s, &n); s += n;
      if (c1 == '\n' || c1 == EOF) break;
      ndiff_grow(dif, 2*dif->buf_s);
    }

    i1 += 1;
    trace("  lhs[%d]: '%s'", dif->row_i+i1, dif->lhs_b);

    // search for number
    int col = 0;
    for (dif->rhs_i=0; (col = ndiff_nextNum(dif, c)); dif->rhs_i=0) {
      if (slice_isElem(&c->col, col)) {
        if (ndiff_testNum(dif, c) == 0) goto lhs_done;
      }
      else
        dif->lhs_i += parse_number(dif->lhs_b+dif->lhs_i, 0,0,0,0);
    }
  }
lhs_done: ;

  // --- rhs ---
  char tag[sizeof c->eps.tag];
  memcpy(tag, dif->lhs_b, sizeof tag);
  strcpy(dif->lhs_b, c->eps.tag);
  C _c = *c;
  _c.eps.scl = -_c.eps.scl;

  while (1) {
    int s = 0, n = 0;

    dif->rhs_i    = 0;
    dif->rhs_b[0] = 0;

    if (c2 == EOF) break;

    while (1) {
      c2 = readLine(dif->rhs_f, dif->rhs_b+s, dif->buf_s-s, &n); s += n;
      if (c2 == '\n' || c2 == EOF) break;
      ndiff_grow(dif, 2*dif->buf_s);
    }

    i2 += 1;
    trace("  rhs[%d]: '%s'", dif->row_i+i2, dif->rhs_b);

    // search for number
    int col = 0;
    for (dif->lhs_i=0; (col = ndiff_nextNum(dif, c)); dif->lhs_i=0) {
      if (slice_isElem(&c->col, col)) {
        if (ndiff_testNum(dif, &_c) == 0) goto rhs_done;
      }
      else
        dif->rhs_i += parse_number(dif->rhs_b+dif->rhs_i, 0,0,0,0);
    }
  }
rhs_done: ;
  memcpy(dif->lhs_b, tag, sizeof tag);

  dif->lhs_i  = 0;
  dif->rhs_i  = 0;
  dif->col_i  = 0;
  dif->row_i += i1 < i2 ? i1 : i2;

  // return with last lhs and rhs lines loaded if number was found

  trace("  buffers: '%.25s'|'%.25s'", dif->lhs_b, dif->rhs_b);
  trace("<-gotoNum line %d (%+d|%+d)", dif->row_i, i1, i2);

  return c1 == EOF || c2 == EOF ? EOF : !EOF;
}

int
ndiff_nextNum (T *dif, const C *c)
{
  assert(dif);

  char *restrict lhs_p = dif->lhs_b+dif->lhs_i;
  char *restrict rhs_p = dif->rhs_b+dif->rhs_i;

  trace("->nextNum  line %d, column %d, char-column %d|%d", dif->row_i, dif->col_i, dif->lhs_i, dif->rhs_i);
  trace("  strings: '%.25s'|'%.25s'", lhs_p, rhs_p);

  if (ndiff_isempty(dif)) goto quit;

retry:

  // search for digits
  if (c->eps.cmd & eps_istr) {
    while (*lhs_p && !isdigit(*lhs_p)) ++lhs_p;
    while (*rhs_p && !isdigit(*rhs_p)) ++rhs_p;
  }
  // search for difference or digits
  else {
    while (*lhs_p && *lhs_p == *rhs_p && !isdigit(*lhs_p))
      ++lhs_p, ++rhs_p;

    // skip whitespaces differences
    if (dif->blank && (isblank(*lhs_p) || isblank(*rhs_p))) {
      while (isblank(*lhs_p)) ++lhs_p;
      while (isblank(*rhs_p)) ++rhs_p;
      goto retry;
    }
  }

  // end-of-line
  if (!*lhs_p && !*rhs_p)
    goto quit;

  // difference in not-a-number
  if (*lhs_p != *rhs_p && (!is_number(lhs_p) || !is_number(rhs_p)))
    goto quit_diff;

  // backtrack numbers
  lhs_p = backtrack_number(lhs_p, dif->lhs_b);
  rhs_p = backtrack_number(rhs_p, dif->rhs_b);

  trace("  backtracking numbers '%.25s'|'%.25s'", lhs_p, rhs_p);

  // at the start of a number?
  if (!is_number_start(lhs_p, dif->lhs_b) || !is_number_start(rhs_p, dif->rhs_b)) {
    if (c->eps.cmd & eps_istr) {
      if (!is_number_start(lhs_p, dif->lhs_b)) skip_identifier(&lhs_p, 0, false);
      if (!is_number_start(rhs_p, dif->rhs_b)) skip_identifier(&rhs_p, 0, false);
    }
    else {
      int strict = true;
      if (c->eps.cmd & eps_omit)
        strict = !is_valid_omit(lhs_p, rhs_p, dif, c->eps.tag);
      int j = strict ? 0 : strlen(c->eps.tag);
      trace("  %s strings '%.25s'|'%.25s'", strict ? "skipping" : "omitting", lhs_p-j, rhs_p-j);
      skip_identifier(&lhs_p, &rhs_p, strict);
    }
    goto retry;
  }

  // numbers found
  dif->lhs_i = lhs_p-dif->lhs_b;
  dif->rhs_i = rhs_p-dif->rhs_b;
  trace("  strnums: '%.25s'|'%.25s'", lhs_p, rhs_p);
  trace("<-nextNum  line %d, column %d, char-column %d|%d", dif->row_i, dif->col_i, dif->lhs_i, dif->rhs_i);
  return ++dif->num_i, ++dif->col_i;

quit_diff:
  dif->lhs_i = lhs_p-dif->lhs_b+1;
  dif->rhs_i = rhs_p-dif->rhs_b+1;
  if (!(c->eps.cmd & eps_gonum) && ++dif->cnt_i <= dif->max_i) {
    if (dif->cnt_i == 1) ndiff_header();
    warning("(%d) files differ at line %d and char-columns %d|%d",
            dif->cnt_i, dif->row_i, dif->lhs_i, dif->rhs_i);
    warning("(%d) strings: '%.25s'|'%.25s'", dif->cnt_i, lhs_p, rhs_p);
  }

quit:
  dif->lhs_i = lhs_p-dif->lhs_b+1;
  dif->rhs_i = rhs_p-dif->rhs_b+1;
  trace("<-nextNum  line %d, column %d, char-column %d|%d", dif->row_i, dif->col_i, dif->lhs_i, dif->rhs_i);
  return dif->col_i = 0;
}

int
ndiff_testNum (T *dif, const C *c)
{
  assert(dif && c);

  char *restrict lhs_p = dif->lhs_b+dif->lhs_i;
  char *restrict rhs_p = dif->rhs_b+dif->rhs_i;
  char *end;

  double lhs_d, rhs_d, dif_a=0, abs_a=0, min_a=0, pow_a=0;

  trace("->testNum  line %d, column %d, char-column %d|%d", dif->row_i, dif->col_i, dif->lhs_i, dif->rhs_i);
  trace("  strnums: '%.25s'|'%.25s'", lhs_p, rhs_p);

  // parse numbers
  int d1=0, d2=0, n1=0, n2=0, e1=0, e2=0, f1=0, f2=0;
  int l1 = parse_number(lhs_p, &d1, &n1, &e1, &f1);
  int l2 = parse_number(rhs_p, &d2, &n2, &e2, &f2);
  int ret = 0;

  // missing numbers
  if (!l1 || !l2) {
    l1 = l2 = 20;
    ret |= eps_ign;
    goto quit_diff;
  }

  // ignore difference
  if (c->eps.cmd & eps_ign) {
    trace("  ignoring numbers '%.25s'|'%.25s'", lhs_p, rhs_p);
    goto quit;
  }

  // omit difference
  if (c->eps.cmd & eps_omit) {
    if (is_valid_omit(lhs_p, rhs_p, dif, c->eps.tag)) {
      trace("  omitting numbers '%.25s'|'%.25s'", lhs_p, rhs_p);
      goto quit;
    }
  }

  // indirections
  // TODO

  // strict comparison...
  if (l1 == l2 && memcmp(lhs_p, rhs_p, l1) == 0)
    goto quit;

  // ...required
  if ((c->eps.cmd & eps_equ) && !(c->eps.cmd & eps_dra)) {
    ret |= eps_equ;
    goto quit_diff;
  }

  // convert numbers
  lhs_d = strtod(lhs_p, &end); assert(end == lhs_p+l1);
  rhs_d = strtod(rhs_p, &end); assert(end == rhs_p+l2);
  dif_a = (lhs_d - rhs_d) * c->eps.scl;
  abs_a = fabs(dif_a);
  min_a = fmin(fabs(lhs_d),fabs(rhs_d));
  pow_a = pow10(-imax(n1, n2));

  // if one number is zero -> relative becomes absolute
  if (!(min_a > 0)) min_a = 1.0;

  trace("  numdiff: |abs|=%.2g, |rel|=%.2g, ndig=%d", abs_a, abs_a/min_a, imax(n1, n2));   

  // absolute comparison
  if (c->eps.cmd & eps_abs)
    if (dif_a > c->eps.abs || dif_a < c->eps._abs) ret |= eps_abs;

  // relative comparison 
  if (c->eps.cmd & eps_rel)
    if (dif_a > c->eps.rel * min_a || dif_a < c->eps._rel * min_a) ret |= eps_rel;

  // input-specific relative comparison (does not apply to integers)
  if ((c->eps.cmd & eps_dig) && (f1 || f2))
    if (dif_a > c->eps.dig * min_a * pow_a || dif_a < c->eps._dig * min_a * pow_a) ret |= eps_dig;

  if ((c->eps.cmd & eps_any) && (ret & eps_dra) != (c->eps.cmd & eps_dra)) ret = 0;
  if (!ret) goto quit;

quit_diff:
  if (!(c->eps.cmd & eps_gonum) && ++dif->cnt_i <= dif->max_i) {
    if (dif->cnt_i == 1) ndiff_header();
    warning("(%d) files differ at line %d column %d between char-columns %d|%d and %d|%d",
            dif->cnt_i, dif->row_i, dif->col_i, dif->lhs_i+1, dif->rhs_i+1, dif->lhs_i+1+l1, dif->rhs_i+1+l2);

    char str[128];
    sprintf(str, "(%%d) numbers: '%%.%ds'|'%%.%ds'", l1,l2);
    warning(str, dif->cnt_i, lhs_p, rhs_p);

    if (ret & eps_ign)
      warning("(%d) one number is missing", dif->cnt_i);

    if (ret & eps_equ)
      warning("(%d) numbers strict representation differ", dif->cnt_i);

    if (ret & eps_abs)
      warning("(%d) absolute error (rule #%d, line %d: %g<=abs<=%g) |abs|=%.2g, |rel|=%.2g, ndig=%d",
              dif->cnt_i, context_findIdx(dif->cxt, c), context_findLine(dif->cxt, c),
              c->eps._abs, c->eps.abs, abs_a, abs_a/min_a, imax(n1, n2));   

    if (ret & eps_rel)
      warning("(%d) relative error (rule #%d, line %d: %g<=rel<=%g) |abs|=%.2g, |rel|=%.2g, ndig=%d",
              dif->cnt_i, context_findIdx(dif->cxt, c), context_findLine(dif->cxt, c),
              c->eps._rel, c->eps.rel, abs_a, abs_a/min_a, imax(n1, n2));   

    if (ret & eps_dig)
      warning("(%d) numdigit error (rule #%d, line %d: %g<=rel<=%g) |abs|=%.2g, |rel|=%.2g, ndig=%d",
              dif->cnt_i, context_findIdx(dif->cxt, c), context_findLine(dif->cxt, c),
              c->eps._dig*pow_a, c->eps.dig*pow_a, abs_a, abs_a/min_a, imax(n1, n2));   
 
  }
  ret = 1;

quit:
  if (c->eps.cmd & eps_lhs) {
    int rn = c->eps.dst_reg;
    memcpy(dif->reg[rn], lhs_p, l1);
    dif->reg[rn][l1] = 0;
  }
  if (c->eps.cmd & eps_rhs) {
    int rn = c->eps.dst_reg;
    memcpy(dif->reg[rn], rhs_p, l2);
    dif->reg[rn][l2] = 0;
  }

  dif->lhs_i += l1;
  dif->rhs_i += l2;
  trace("<-testNum  line %d, column %d, char-column %d|%d", dif->row_i, dif->col_i, dif->lhs_i, dif->rhs_i);

  return ret;
}

void
ndiff_option  (T *dif, const int *keep_, const int *blank_, const int *check_)
{
  assert(dif);
  
  if (keep_ ) dif->max_i = *keep_;
  if (blank_) dif->blank = *blank_; 
  if (check_) dif->check = *check_; 

  ensure(dif->max_i > 0, "number of kept diff must be positive");
}

void
ndiff_getInfo (const T *dif, int *row_, int *col_, int *cnt_, long *num_)
{
  assert(dif);

  if (row_) *row_ = dif->row_i;
  if (col_) *col_ = dif->col_i;
  if (cnt_) *cnt_ = dif->cnt_i;
  if (num_) *num_ = dif->num_i;
}

int
ndiff_feof (const T *dif, int both)
{
  assert(dif);

  return both ? feof(dif->lhs_f) && feof(dif->rhs_f)
              : feof(dif->lhs_f) || feof(dif->rhs_f);
}

int
ndiff_isempty (const T *dif)
{
  assert(dif);

  return !dif->lhs_b[dif->lhs_i] && !dif->rhs_b[dif->rhs_i];
}

// --- main ndiff loop --------------------------------------------------------

void
ndiff_loop(T *dif)
{
  assert(dif);

  const C *c, *c2;
  int row=0, col;
  int saved_level = logmsg_config.level;

  while(!ndiff_feof(dif, 0)) {
    ++row, col=0;

    c = context_getInc(dif->cxt, row, col);
    ensure(c, "invalid context");
    if (dif->check && c != (c2 = context_getAt(dif->cxt, row, col)))
      ndiff_error(dif->cxt, c, c2, row, col);

    // trace rule
    if (c->eps.cmd & eps_trace && c->eps.cmd & eps_sgg) {
      logmsg_config.level = trace_level;
      trace("~>active:  rule #%d, line %d, cmd = %d",
            context_findIdx(dif->cxt,c), context_findLine(dif->cxt,c), c->eps.cmd);
      logmsg_config.level = saved_level;
    }

    // skip this line
    if (c->eps.cmd & eps_skip) {
      ndiff_skipLine(dif);
      continue;
    }

    // goto or read line(s)
    if (c->eps.cmd & eps_goto) {
      ndiff_gotoLine(dif, c);
      ndiff_getInfo(dif, &row, 0, 0, 0);
    } else
    if (c->eps.cmd & eps_gonum) {
      ndiff_gotoNum(dif, c);
      ndiff_getInfo(dif, &row, 0, 0, 0);
    } else {
      ndiff_readLine(dif);
      if (ndiff_isempty(dif)) continue;
    }

    // for each number column, diff-chars between numbers
    while((col = ndiff_nextNum(dif, c))) {
      c = context_getInc(dif->cxt, row, col);
      ensure(c, "invalid context");
      if (dif->check && c != (c2 = context_getAt(dif->cxt, row, col)))
        ndiff_error(dif->cxt, c, c2, row, col); 

      // trace rule
      if (c->eps.cmd & eps_trace) {
        logmsg_config.level = trace_level;
        trace("~>active:  rule #%d, line %d, cmd = %d",
              context_findIdx(dif->cxt,c), context_findLine(dif->cxt,c), c->eps.cmd);
      }

      // check numbers
      ndiff_testNum(dif, c);

      // restore logmsg
      logmsg_config.level = saved_level;
    }
  }

  if (dif->blank) {
    skipSpace(dif->lhs_f, 0);
    skipSpace(dif->rhs_f, 0);
  }
}

#undef T
#undef C

// -----------------------------------------------------------------------------
// ----- testsuite
// -----------------------------------------------------------------------------

#ifndef NTEST

#include "utest.h"

#define T struct ndiff

// ----- debug

// ----- teardown

static T*
ut_teardown(T *dif)
{
  ndiff_clear(dif);
  return dif;
}

// ----- test

static void 
ut_testPow10(struct utest *utest, T* dif)
{
  (void)dif;

  for (int k = -100; k < 100; k++)
    UTEST(pow10(k) == pow(10, k));
}

static void
ut_testNul(struct utest *utest, T* dif)
{
  UTEST(dif != 0);
}

// ----- unit tests

static struct spec {
  const char *name;
  T*        (*setup)   (T*);
  void      (*test )   (struct utest*, T*);
  T*        (*teardown)(T*);
} spec[] = {
  { "power of 10",                          0        , ut_testPow10, 0           },
  { "empty input",                          0        , ut_testNul  , ut_teardown },
};
enum { spec_n = sizeof spec/sizeof *spec };

// ----- interface

void
ndiff_utest(struct utest *ut)
{
  assert(ut);
  T *dif = ndiff_alloc(stdout, stdout, 0, 0);

  utest_title(ut, "File diff");

  for (int k = 0; k < spec_n; k++) {
    utest_init(ut, spec[k].name);
    if (spec[k].setup)    dif = spec[k].setup(dif);
    spec[k].test(ut, dif);
    if (spec[k].teardown) dif = spec[k].teardown(dif);
    utest_fini(ut);
  }

  ndiff_free(dif);
}

#endif
