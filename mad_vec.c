#include <math.h>
#include "mad_wrap_f.h"

double
vdot(int* n, double* v1, double* v2)
  /* returns dot product of vectors v1 and v2 */
{
  int i;
  double dot = 0;
  for (i = 0; i < *n; i++)  dot += v1[i] * v2[i];
  return dot;
}

double
vmod(int* n, double* v)
{
  int i;
  double mod = 0;
  for (i = 0; i < *n; i++)  mod += v[i] * v[i];
  return sqrt(mod);
}

void
zero_double(double* a, int n)
  /* sets first n values in double array a to zero */
{
  int j;
  for (j = 0; j < n; j++)  a[j] = 0;
}


