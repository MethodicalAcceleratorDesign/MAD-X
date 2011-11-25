#include "madx.h"

double
mult_par(char* par, struct element* el)
  /* returns multipole parameter for par = "k0l" or "k0sl" etc. */
{
  char tmp[12];
  char* p;
  double val = zero, vect[FIELD_MAX];
  int k = 0, l, skew = 0;
  strcpy(tmp, par);
  if (*tmp == 'k' && (p = strchr(tmp, 'l')) != NULL)
  {
    *p = '\0';  /* suppress trailing l */
    if ((p = strchr(tmp, 's')) != NULL)
    {
      skew = 1; *p = '\0';
    }
    sscanf(&tmp[1], "%d", &k);
    if (skew) l = element_vector(el, "ksl", vect);
    else      l = element_vector(el, "knl", vect);
    if (k < l) val = vect[k];
  }
  return val;
}


