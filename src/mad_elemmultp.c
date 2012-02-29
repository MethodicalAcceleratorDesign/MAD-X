#include "madx.h"

double
mult_par(const char* par, struct element* el)
  /* returns multipole parameter for par = "k0l" or "k0sl" etc. */
{
  double val = zero;
  char* p;
  if (*par == 'k' && (p = strchr(par, 'l')) != NULL)
  {
    *p = '\0';  /* suppress trailing l */
    int skew = 0;
    if ((p = strchr(par, 's')) != NULL)
    {
      skew = 1; *p = '\0';
    }
    int k = 0;
    sscanf(&par[1], "%d", &k);
    double vect[FIELD_MAX];
    int l;
    if (skew) l = element_vector(el, "ksl", vect);
    else      l = element_vector(el, "knl", vect);
    if (k < l) val = vect[k];
  }
  return val;
}


