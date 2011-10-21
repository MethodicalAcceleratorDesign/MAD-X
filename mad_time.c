#include "madx.h"
#include <time.h>

#ifndef _WIN32
#include <sys/time.h> /* for gettimeofday */
#endif

void
time_stamp(char* place)
{
  time_t now;
  int k, l;
  time(&now);    /* get system time */
  k = (int)now - (int)start_time;
  l = (int)now - (int)last_time;
  last_time = now;
  fprintf(prt_file, "sec.s since start: %d   since last call: %d\n", k, l);
}

float
fextim(void)
{
   float mytime;
   struct timeval tp;
   struct timezone tzp;
   gettimeofday(&tp,&tzp);
   mytime = (float)(tp.tv_sec%10000) + 1.e-6 * tp.tv_usec; /* seconds from epoch, modulo 10 000 */
   /* printf("Time now:  %-6.3f\n",mytime);    */
   return(mytime);
}


