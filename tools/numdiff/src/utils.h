#ifndef UTILS_H
#define UTILS_H

/*
 o---------------------------------------------------------------------o
 |
 | Numdiff
 |
 | Copyright (c) 2012+ Laurent Deniau
 |
 o---------------------------------------------------------------------o
  
   Purpose:
     provides utilities
 
 o---------------------------------------------------------------------o
*/

#include <stdio.h>
#include <ctype.h>
#include <math.h>

// extern functions

FILE* open_indexedFile(const char* str, int idx, const char *ext, int optext, int required);
void  accum_summary(int total, int failed);

// inline functions

#if !__STDC__ || __STDC_VERSION__ < 199901L
static inline int
isblank(int c)
{
  return c == ' ' || c == '\t';
}

static inline double
fmin(double a, double b)
{
  return a < b ? a : b;
}

static inline double
fmax(double a, double b)
{
  return a > b ? a : b;
}
#endif

static inline int
imin (int a, int b)
{
  return a<b ? a : b;
}

static inline int
imax (int a, int b)
{
  return a>b ? a : b;
}

static inline double 
pow10(int i)
{
  extern const double *const pow10_table99;
  return -100 < i && i < 100 ? pow10_table99[i] : pow(10, i);
}

// ----- public (read helpers)

static inline int
skipSpace (FILE *fp, int *i_)
{
  int c = 0, i = 0;

  while ((c = getc(fp)) != EOF) {
    if (!isspace(c)) break;
    i++;
  }

  if (i_) *i_ = i;

  return c; 
}

static inline int
skipLine (FILE *fp, int *i_)
{
  int c = 0, i = 0;

  while ((c = getc(fp)) != EOF) {
    if (c == '\n') break;                 // \n   : Unix, Linux, MacOSX
    if (c == '\r') {
      if ((c = getc(fp)) != '\n')         // \r\n : Windows
        ungetc(c, fp);                    // \r   : Mac (old)
      c = '\n'; break;
    }
    i++;
  }

  if (i_) *i_ = i;

  return c;
}

static inline int
readLine (FILE *fp, char *buf, int n, int *i_)
{
  int c = 0, i = 0;

  while (i < n-1 && (c = getc(fp)) != EOF) {
    if (c == '\n') break;                 // \n   : Unix, Linux, MacOSX
    if (c == '\r') {
      if ((c = getc(fp)) != '\n')         // \r\n : Windows
        ungetc(c, fp);                    // \r   : Mac (old)
      c = '\n'; break;
    }
    buf[i++] = c;
  }
  buf[i] = 0;

  if (i_) *i_ = i;

  return c;
}

#endif
