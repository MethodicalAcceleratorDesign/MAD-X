#include "madx.h"

double
vdot_(int* n, double* v1, double* v2)
  /* returns dot product of vectors v1 and v2 */
{
  int i;
  double dot = 0;
  for (i = 0; i < *n; i++)  dot += v1[i] * v2[i];
  return dot;
}

double
vmod_(int* n, double* v)
{
  int i;
  double mod = 0;
  for (i = 0; i < *n; i++)  mod += v[i] * v[i];
  return sqrt(mod);
}

void
zero_double_(double* a, int n)
  /* sets first n values in double array a to zero */
{
  int j;
  for (j = 0; j < n; j++)  a[j] = 0;
}


