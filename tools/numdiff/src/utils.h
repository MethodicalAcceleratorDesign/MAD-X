#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <math.h>

// extern functions

FILE* open_indexedFile(const char* str, int idx, const char *fmt, int strict);

// inline functions

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
