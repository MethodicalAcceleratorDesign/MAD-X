#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "constraint.h"
#include "utils.h"
#include "error.h"

#define T struct constraint

// ----- constants

const char * const eps_cmd_cstr[] = {
  [eps_invalid]             = "invalid",
  [eps_dig]                 = "dig",
  [eps_rel]                 = "rel",
  [eps_rel|eps_dig]         = "rel&dig",
  [eps_abs]                 = "abs",
  [eps_abs|eps_dig]         = "abs&dig",
  [eps_abs|eps_rel]         = "abs&rel",
  [eps_abs|eps_rel|eps_dig] = "abs&rel&dig",
  [eps_equ]                 = "equ",
  [eps_ign]                 = "ign",
  [eps_skip]                = "skip",
  [eps_goto]                = "goto",
};

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
  int c = 0, r = 0;
  char str[16];

  e->cmd = eps_invalid;

  while (1) {
    // parse next constraint
    *str = 0;
    r = fscanf(in, "%*[ \t]%10[^= \t\n\r#!]", str);
    if (r == EOF || *str == 0) break;

         if (strcmp(str, "skip") == 0) {
      e->cmd |= eps_skip; trace("[%d] skip", row);
    }
    else if (strcmp(str, "ign") == 0) {
      e->cmd |= eps_ign;  trace("[%d] ign", row);
    }
    else if (strcmp(str, "equ") == 0) {
      e->cmd |= eps_equ;  trace("[%d] equ", row);
    }
    else if (strcmp(str, "dig") == 0 && (r = fscanf(in, "=%lf", &e->dig)) == 1) {
      e->cmd |= eps_dig;  trace("[%d] dig=%g", row, e->dig);
    }
    else if (strcmp(str, "rel") == 0 && (r = fscanf(in, "=%lf", &e->rel)) == 1) {
      e->cmd |= eps_rel;  trace("[%d] rel=%g", row, e->rel);
    }
    else if (strcmp(str, "abs") == 0 && (r = fscanf(in, "=%lf", &e->abs)) == 1) {
      e->cmd |= eps_abs;  trace("[%d] abs=%g", row, e->abs);
    }
    else if (strcmp(str, "goto") == 0 && (r = fscanf(in, "='%32[^']'", e->tag)) == 1) {
      e->cmd |= eps_goto; e->tag[sizeof e->tag-1] = 0;
                          trace("[%d] goto='%s'", row, e->tag);
    }
    else {
                          trace("[%d] invalid '%s'", row, str);
      e->cmd = eps_invalid;
      break;
    }

    // next char
    ungetc((c = getc(in)), in);
    if (c == EOF || (isspace(c) && !isblank(c)) || c == '#' || c == '!') break; 
  }

  trace("<-readEps cmd = %d, str = '%s', c = '%c'", e->cmd, str, c);

  return e->cmd == eps_invalid || r == EOF ? EOF : 0;
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

  if (cst->eps.cmd >= eps_dig && cst->eps.cmd < eps_equ) {
    if (cst->eps.cmd & eps_dig) fprintf(out, "dig=%g ", cst->eps.dig);    
    if (cst->eps.cmd & eps_rel) fprintf(out, "rel=%g ", cst->eps.rel);    
    if (cst->eps.cmd & eps_abs) fprintf(out, "abs=%g ", cst->eps.abs);    
  } else
  if (cst->eps.cmd == eps_goto) {
    fprintf(out, "goto='%s' ", cst->eps.tag);
  } else
    fprintf(out, "%s", eps_cstr(cst->eps.cmd));
}

void
constraint_scan(T* cst, FILE *in, int *row)
{
  int c;
  assert(cst && row);

  cst->eps.cmd = eps_invalid;

  if (!in) in = stdin;

retry:

  while((c = getc(in)) != EOF && isblank(c)) ;

  // end of file
  if (c == EOF) return;

  ungetc(c, in);

  // comment or empty line
  if (c == '\n' || c == '\r' || c == '!' || c == '#') {
    if (skipLine(in, 0) == '\n') ++*row;
    goto retry;
  }

  ensure(readSlcOrRng(&cst->row, in      ) != EOF, "invalid row constraint line %d"   , *row);
  ensure(readSlcOrRng(&cst->col, in      ) != EOF, "invalid column constraint line %d", *row);
  ensure(readEps     (&cst->eps, in, *row) != EOF, "invalid eps constraint line %d"   , *row);
  if (skipLine(in, 0) == '\n') ++*row;
}

