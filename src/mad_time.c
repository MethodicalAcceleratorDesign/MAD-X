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
