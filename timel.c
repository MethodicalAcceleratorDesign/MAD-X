/* timex, timest, etc Eric McIntosh 28/9/99 */
#include <sys/types.h>
#include <time.h>
#include <sys/times.h>
#include <sys/param.h>
#include <sys/time.h>
#include <sys/resource.h>

#ifndef HZ
#define HZ 60.;
#endif

struct tms tps;
static float timlim;
static time_t timstart, timlast;
static int tml_init = 1;
float deftim = 999.;

                   /*  local routine called by timst, and time_init */
static void time_st(timl)
float timl;
{
    times(&tps);
    timlim = timl;
    timstart =  tps.tms_utime+tps.tms_cutime+tps.tms_stime+tps.tms_cstime;
    timlast  = timstart;
    tml_init = 0;
    return;
}
                   /*  local routine to start by default  */
static void time_init()
{
	struct rlimit rlimit;
	float  maxtime;

	maxtime=deftim;

	if (getrlimit(RLIMIT_CPU, &rlimit)==0) {
		if ( rlimit.rlim_cur != RLIM_INFINITY )
		   maxtime = (float) rlimit.rlim_cur;
	}	

	time_st(maxtime);
	return;
}

void timest_(timl)
float *timl;
{
 float  maxtime;

 if (tml_init != 0) {

    maxtime = *timl;
    time_st(maxtime);
 }
 return;
}
void timex_(tx)
float *tx;
{
   time_t timnow;
   if (tml_init) {
       time_init();
       *tx = 0.;
   }
   else {
       times(&tps);
       timnow = tps.tms_utime+tps.tms_cutime+tps.tms_stime+tps.tms_cstime;
       *tx = (float) (timnow - timstart) / HZ;
   }
   return;
}

void timed_(td)
float *td;
{
   time_t timnow;
   if (tml_init) {
       time_init();
       *td = timlim;
   }
   else {
       times(&tps);
       timnow = tps.tms_utime+tps.tms_cutime+tps.tms_stime+tps.tms_cstime;
       *td = (float) (timnow - timlast) / HZ;
       timlast = timnow;
   }
   return;
}

void timel_(tl)
float *tl;
{
   time_t timnow;
   if (tml_init) {
       time_init();
       *tl = timlim;
   }
   else {
       times(&tps);
       timnow = tps.tms_utime+tps.tms_cutime+tps.tms_stime+tps.tms_cstime;
       *tl = timlim - (float) (timnow - timstart) / HZ;
   }
   return;
}

