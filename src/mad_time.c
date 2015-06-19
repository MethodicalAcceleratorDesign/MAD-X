#include "madx.h"
#include <time.h>

#ifndef _WIN32
#include <sys/time.h>  /* for gettimeofday */
#else
#include <sys/timeb.h> /* for ftime */
#endif

void
time_stamp(const char* place)
{
  time_t now;
  int k, l;

  (void)place;
  time(&now);    /* get system time */
  k = (int)now - (int)start_time;
  l = (int)now - (int)last_time;
  last_time = now;
  fprintf(prt_file, "sec.s since start: %d   since last call: %d\n", k, l);
}

#if 0 // not used... (commented in micado)
static double
fextim(void)
{
   double mytime;

   #ifndef _WIN32 /* gettimeofday available */
     struct timeval tp;
     gettimeofday(&tp,0);
     mytime = (double)(tp.tv_sec%10000) + 1.e-6 * tp.tv_usec; /* seconds from epoch, modulo 10 000 */
   #else /* use old ftime */
     struct timeb tp;
     ftime(&tp);
     mytime = (double)(tp.time%10000) + 0.001*tp.millitm;
   #endif

   /* printf("Time now:  %-6.3f\n",mytime);    */
   return mytime;
}
#endif
