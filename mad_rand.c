#include "madx.h"

// private functions

static void
irngen(void)
  /* creates random number for frndm() */
{
  int i, j;
  for (i = 0; i < NJ_RAND; i++)
  {
    if ((j = irn_rand[i] - irn_rand[i+NR_RAND-NJ_RAND]) < 0) j += MAX_RAND;
    irn_rand[i] = j;
  }
  for (i = NJ_RAND; i < NR_RAND; i++)
  {
    if ((j = irn_rand[i] - irn_rand[i-NJ_RAND]) < 0) j += MAX_RAND;
    irn_rand[i] = j;
  }
  next_rand = 0;
}

// interface

void
init55(int seed)
  /* initializes random number algorithm */
{
  int i, ii, k = 1, j = abs(seed)%MAX_RAND;
  irn_rand[NR_RAND-1] = j;
  for (i = 0; i < NR_RAND-1; i++)
  {
    ii = (ND_RAND*(i+1))%NR_RAND;
    irn_rand[ii-1] = k;
    if ((k = j - k) < 0) k += MAX_RAND;
    j = irn_rand[ii-1];
  }
  /* warm up */
  for (i = 0; i < 3; i++) irngen();
}


double
frndm(void)
  /* returns random number r with 0 <= r < 1 from flat distribution */
{
  const double one = 1;
  double scale = one / MAX_RAND;
  if (next_rand == NR_RAND) irngen();
  return scale*irn_rand[next_rand++];
}

double
grndm(void)
  /* returns random number x from normal distribution */
{
  double xi1 = 2*frndm()-one, xi2=2*frndm()-one, zzr;
  while ((zzr = xi1*xi1+xi2*xi2) > one)
  {
    xi1 = 2*frndm()-one; xi2=2*frndm()-one;
  }
  zzr = sqrt(-2*log(zzr)/zzr);
  return xi1*zzr;
}

double
tgrndm(double cut)
  /* returns random variable from normal distribution cat at 'cut' sigmas */
{
  double ret = zero;
  if (cut > zero)
  {
    ret = grndm();
    while (fabs(ret) > fabs(cut))  ret = grndm();
  }
  return ret;
}


