#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>

#include "ndiff.h"
#include "constraint.h"
#include "context.h"
#include "error.h"
#include "utils.h"

#define T struct ndiff

// ----- types

struct ndiff {
  // files
  FILE *lhs_f, *rhs_f;
  int   row_i,  col_i; // line, num-column

  // diff counter
  int   cnt_i;

  // buffers
  int   lhs_i,  rhs_i; // char-columns
  int   buf_s;         // capacity
  char *lhs_b, *rhs_b;
};

// ----- private (parser helpers)

static inline int
is_separator (int c)
{
  return isblank(c) || (ispunct(c) && c != '.' && c != '_');
}

static inline int
is_number_start(char *buf, const char *beg)
{
  // number starts by a '-' or is at the beginning or is preceded by a separator
  return *buf == '-' || buf == beg || (buf > beg && is_separator(buf[-1]));
}

static inline int
is_number (char *buf)
{
  int i = 0;

  // sign
  if (buf[i] == '+' || buf[i] == '-') i++;

  // digits
  if(isdigit(buf[i])) return 1;

  // dot
  if (buf[i] == '.') ++i;

  // decimals
  if(isdigit(buf[i])) return 1;

  return 0;
}

static inline char*
backtrace_number (char *buf, const char *beg)
{
  if (isdigit(*buf)) {
    if (buf > beg && buf[-1] == '.') --buf;
    if (buf > beg && buf[-1] == '-') --buf;
  }

  return buf;
}

static inline int
parse_number (char *buf, int *d_, int *n_, int *e_, int *f_)
{
  int i = 0, d = 0, n = 0, e = 0;

  // sign
  if (buf[i] == '+' || buf[i] == '-') i++;

  // digits
  while(isdigit(buf[i])) i++, n++;

  // dot
  if (buf[i] == '.') d = ++i;

  // decimals
  if (d) { 
    if (!n) n = 1;
    while(isdigit(buf[i])) i++, n++;
  }

  // ensure at least ±# or ±#. or ±.#
  if(!(i > 0 && (isdigit(buf[i-1]) || (i > 1 &&  isdigit(buf[i-2])))))
    return 0;

  // exponent
  if (buf[i] == 'e' || buf[i] == 'E' || buf[i] == 'd' || buf[i] == 'D')
    buf[i] = 'e', e = ++i;

  if (e) {
    // sign
    if (buf[i] == '+' || buf[i] == '-') i++;

    // digits
    while(isdigit(buf[i])) i++;

    // ensure e# or e±#
    if (!isdigit(buf[i-1])) return 0;
  }

  if (n_) *n_ = n;
  if (d_) *d_ = d-1;
  if (e_) *e_ = e-1;
  if (f_) *f_ = d > 0 || (e > 0 && buf[e] == '-');

  return i;
}

static inline int
skip_identifier_digits(const char *lhs, const char *rhs)
{
  int i = 0;

  while (lhs[i] == rhs[i] && (isdigit(lhs[i]) || lhs[i] == '.'))
    i++;

  return i;
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
    .buf_s = n
  };
}

static inline void
ndiff_teardown (T *dif)
{
  free(dif->lhs_b);
  free(dif->rhs_b);

  *dif = (T) {
    .lhs_f = dif->lhs_f, .rhs_f = dif->rhs_f
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
ndiff_error(const struct context    *cxt,
            const struct constraint *c,
            const struct constraint *c2,
            int row, int col)
{
  warning("dual constraints differ at %d:%d", row, col);
  warning("getIncr select [#%d]", context_findIdx(cxt, c )+1);
  warning("getAt   select [#%d]", context_findIdx(cxt, c2)+1);
  warning("rules list:");
  context_print(cxt, stderr);
  error("please report to mad@cern.ch");
}

// -----------------------------------------------------------------------------
// ----- interface
// -----------------------------------------------------------------------------

T*
ndiff_alloc (FILE *lhs_f, FILE *rhs_f, int n)
{
  assert(lhs_f && rhs_f);

  T *dif = malloc(sizeof *dif);
  ensure(dif, "out of memory");

  *dif = (T) { .lhs_f = lhs_f, .rhs_f = rhs_f };

  ndiff_setup(dif, n);
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
ndiff_readLine (T *dif)
{
  assert(dif);
  int s1 = 0, s2 = 0;
  int c1, c2, n = 0;

  ndiff_reset_buf(dif);

  while (1) {
    c1 = readLine(dif->lhs_f, dif->lhs_b+s1, dif->buf_s-s1, &n); s1 += n;
    c2 = readLine(dif->rhs_f, dif->rhs_b+s2, dif->buf_s-s2, &n); s2 += n;
    if (c1 == '\n' || c2 == '\n' || c1 == EOF || c2 == EOF) break;
    ndiff_grow(dif, 2*dif->buf_s);
  }

  dif->col_i  = 0;
  dif->row_i += 1;

  trace("<-readLine line %d", dif->row_i);
  trace("  buffers: '%.30s'|'%.30s'", dif->lhs_b, dif->rhs_b);

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

void
ndiff_diffLine (T *dif, int blank)
{
  assert(dif);

  char *lhs_p = dif->lhs_b+dif->lhs_i;
  char *rhs_p = dif->rhs_b+dif->rhs_i;

  trace("->diffLine line %d char-column %d|%d", dif->row_i, dif->lhs_i, dif->rhs_i);
  trace("  buffers: '%.30s'|'%.30s'", lhs_p, rhs_p);

retry:

  // fast search
  if (!strcmp(lhs_p, rhs_p)) {
    int n = strlen(lhs_p);
    dif->lhs_i += n;
    dif->rhs_i += n;
    return;
  }

  // slow search (find index)
  int i;
  for(i = 0; lhs_p[i] == rhs_p[i]; i++) ;

  lhs_p += i; dif->lhs_i += i;
  rhs_p += i; dif->rhs_i += i;
  i = 0;

  if (blank && (isblank(*lhs_p) || isblank(*rhs_p))) {
    while (isblank(*lhs_p)) ++lhs_p, ++dif->lhs_i;
    while (isblank(*rhs_p)) ++rhs_p, ++dif->rhs_i;
    goto retry;
  }

  dif->lhs_i += 1;
  dif->rhs_i += 1;
  dif->cnt_i += 1;
  warning("(%d) files differ at line %d at char-column %d|%d",
          dif->cnt_i, dif->row_i, dif->lhs_i, dif->rhs_i);
  warning("(%d) strings: '%.20s'|'%.20s'", dif->cnt_i, lhs_p, rhs_p);
}

int
ndiff_nextNum (T *dif, int blank)
{
  assert(dif);

  char *restrict lhs_p = dif->lhs_b+dif->lhs_i;
  char *restrict rhs_p = dif->rhs_b+dif->rhs_i;

retry:

  // search for difference or digits
  { int i = 0;

    while (lhs_p[i] && lhs_p[i] == rhs_p[i] && !isdigit(lhs_p[i]))
      i++;

    lhs_p += i; rhs_p += i;
  }

  // skip whitespaces differences
  if (blank && (isblank(*lhs_p) || isblank(*rhs_p))) {
    while (isblank(*lhs_p)) ++lhs_p;
    while (isblank(*rhs_p)) ++rhs_p;
    goto retry;
  }

  // end-of-line
  if (!*lhs_p && !*rhs_p)
    goto quit;

  // difference in not-a-number
  if (*lhs_p != *rhs_p && (!is_number(lhs_p) || !is_number(rhs_p)))
    goto quit_diff;

  // backtrace number
  lhs_p = backtrace_number(lhs_p, dif->lhs_b);
  rhs_p = backtrace_number(rhs_p, dif->rhs_b);

  // at the start of a number?
  if (!is_number_start(lhs_p, dif->lhs_b) || !is_number_start(rhs_p, dif->rhs_b)) {
    int i = skip_identifier_digits(lhs_p, rhs_p);
    lhs_p += i; rhs_p += i;
    if (!isdigit(*lhs_p)) goto retry;
    goto quit_diff;
  }

  // numbers found
  dif->lhs_i = lhs_p-dif->lhs_b;
  dif->rhs_i = rhs_p-dif->rhs_b;
  trace("<-nextNum line %d char-column %d|%d", dif->row_i, dif->lhs_i, dif->rhs_i);
  trace("  strnums: '%.20s'|'%.20s'", lhs_p, rhs_p);
  return ++dif->col_i;

quit_diff:
  dif->lhs_i = lhs_p-dif->lhs_b+1;
  dif->rhs_i = rhs_p-dif->rhs_b+1;
  dif->cnt_i += 1;
  warning("(%d) files differ at line %d and char-columns %d|%d",
          dif->cnt_i, dif->row_i, dif->lhs_i, dif->rhs_i);
  warning("(%d) strings: '%.20s'|'%.20s'", dif->cnt_i, lhs_p, rhs_p);
  return dif->col_i = 0;

quit:
  dif->lhs_i = lhs_p-dif->lhs_b+1;
  dif->rhs_i = rhs_p-dif->rhs_b+1;
  return dif->col_i = 0;
}

int
ndiff_testNum (T *dif, const struct context *cxt, const struct constraint *c)
{
  assert(dif && c);

  char *lhs_p = dif->lhs_b+dif->lhs_i;
  char *rhs_p = dif->rhs_b+dif->rhs_i;
  char *end;

  double lhs_d, rhs_d, dif_a, min_a;

  // parse numbers
  int d1, d2, n1, n2, e1, e2, f1, f2;
  int l1 = parse_number(lhs_p, &d1, &n1, &e1, &f1);
  int l2 = parse_number(rhs_p, &d2, &n2, &e2, &f2);
  int ret = 0;

  // invalid numbers
  if (!l1 || !l2) {
    l1 = l2 = 20;
    warning("(%d) one number is missing", dif->cnt_i+1);
    goto quit_diff;
  }

  // ignore difference
  if (c->eps.cmd == eps_ign)
    goto quit;

  // strict comparison...
  if (l1 == l2 && memcmp(lhs_p, rhs_p, l1) == 0)
    goto quit;

  // ...required or integer values
  if (c->eps.cmd == eps_equ || (c->eps.cmd != eps_abs && !f1 && !f2)) {
    warning("(%d) numbers strict representation differ", dif->cnt_i+1);
    goto quit_diff;
  }

  lhs_d = strtod(lhs_p, &end); assert(end == lhs_p+l1);
  rhs_d = strtod(rhs_p, &end); assert(end == rhs_p+l2);
  dif_a = fabs(lhs_d - rhs_d);
  min_a = fmin(fabs(lhs_d),fabs(rhs_d));

  // if one number is zero -> relative becomes absolute
  if (!(min_a > 0)) min_a = 1.0;

  // input-specific comparison
  if (c->eps.cmd & eps_dig) {
    double rel = c->eps.dig * pow10(-imax(n1, n2));
    if (dif_a > rel * min_a) {
      warning("(%d) numdigit error [rule #%d] rel = %g [|abs_err|=%.2g, |rel_err|=%.2g, ndig=%d]",
              dif->cnt_i+1, context_findIdx(cxt, c)+1, rel, dif_a, dif_a/min_a, imax(n1, n2));   
      ret = 1; 
    }
  }

  // relative comparison 
  if (c->eps.cmd & eps_rel)
    if (dif_a > c->eps.rel * min_a) {
      warning("(%d) relative error [rule #%d] rel = %g [|abs_err|=%.2g, |rel_err|=%.2g]",
              dif->cnt_i+1, context_findIdx(cxt, c)+1, c->eps.rel, dif_a, dif_a/min_a);   
      ret = 1; 
    }

  // absolute comparison
  if (c->eps.cmd & eps_abs)
    if (dif_a > c->eps.abs) {
      warning("(%d) absolute error [rule #%d] abs = %g [|abs_err|=%.2g]",
              dif->cnt_i+1, context_findIdx(cxt, c)+1, c->eps.abs, dif_a);   
      ret = 1;
    }

  if (!ret) goto quit;

quit_diff:
  ret = 1;
  dif->cnt_i += 1;
  warning("(%d) files differ at line %d column %d between char-columns %d|%d and %d|%d",
          dif->cnt_i, dif->row_i, dif->col_i, dif->lhs_i+1, dif->rhs_i+1, dif->lhs_i+1+l1, dif->rhs_i+1+l2);

  char str[128];
  sprintf(str, "(%%d) numbers: '%%.%ds'|'%%.%ds'", l1,l2);
  warning(str, dif->cnt_i, lhs_p, rhs_p);

quit:
  trace("<-testNum line %d char-column %d|%d", dif->row_i, dif->lhs_i, dif->rhs_i);
  trace("  numbers: [%d|%d] '%.30s'|'%.30s'", l1, l2, lhs_p, rhs_p);

  dif->lhs_i += l1;
  dif->rhs_i += l2;

  return ret;
}

void
ndiff_getInfo (const T *dif, int *row_, int *col_, int *cnt_)
{
  assert(dif);

  if (row_) *row_ = dif->row_i;
  if (col_) *col_ = dif->col_i;
  if (cnt_) *cnt_ = dif->cnt_i;
}

int
ndiff_feof (const T *dif)
{
  assert(dif);

  return feof(dif->lhs_f) || feof(dif->rhs_f);
}

void
ndiff_loop(struct ndiff *dif, struct context *cxt, int blank, int check)
{
  const struct constraint *c, *c2;
  int row = 0, col, n;

  while(!ndiff_feof(dif)) {
    ++row, col = 0;

    c = context_getInc(cxt, row, col);
    if (check && c != (c2 = context_getAt(cxt, row, col)))
      ndiff_error(cxt, c, c2, row, col); 

    // no constraint, diff-lines
    if (!c) {
      ndiff_readLine(dif);
      ndiff_diffLine(dif, blank);
      continue;
    }

    // skip this line
    if (c->eps.cmd == eps_skip) {
      ndiff_skipLine(dif);
      continue;
    }

    // normal constraint(s), read line
    ndiff_readLine(dif);

    // for each number column, diff-chars between numbers
    while((n = ndiff_nextNum(dif, blank))) {
      ++col;
      assert(n == col);

      c = context_getInc(cxt, row, col);
      if (check && c != (c2 = context_getAt(cxt, row, col)))
        ndiff_error(cxt, c, c2, row, col); 

      // no constraint means equal
      if (!c) {
        static const struct constraint cst_equ = { .eps = { .cmd = eps_equ } };
        c = &cst_equ;
      }

      ndiff_testNum(dif, cxt, c);
    }
  }
}

#undef T

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

  for (int k = -100; k < 100; k++) {
    UTEST(pow10(k) == pow(10, k)); 
//    if (pow10(k) != pow(10, k))
//      fprintf(stderr, "pow10(%d)=%g != pow(10,%d)=%g, err=%g\n", k, pow10(k), k, pow(10,k), fabs(pow10(k)-pow(10,k)));
  }
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
  T *dif = ndiff_alloc(stdout, stdout, 0);

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
