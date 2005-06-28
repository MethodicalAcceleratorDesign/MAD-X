#include <stdio.h>
#include <time.h>
#define cbyt     cbyt_
#define mydtime  mydtime_
  void cbyt(int* source, int* s_pos, int* target, int* t_pos, int* n)
/* inserts n_bit byte from source at position s_pos in target at t_pos.
   Attention: least significant bit is #1, positions are those of the least
   significant bit in byte */
{
  int mask = 1;
  int so = *source, tg = *target;
  so >>= (*s_pos - 1);
  mask <<= *n;
  mask -= 1;
  so &= mask;
  mask <<= (*t_pos - 1);  so <<= (*t_pos - 1);
  tg &= ~mask; tg |= so;
  *target = tg;
}

  void mydtime(int* year, int* month, int* day, int* hour, int* minute,
               int* sec)
{
  time_t _time;
  struct tm* tm;
  time(&_time); /* initialize timing */
  tm = localtime(&_time); /* split system time */
  *year = tm->tm_year%100;
  *month = tm->tm_mon;
  *day = tm->tm_mday;
  *hour = tm->tm_hour;
  *minute = tm->tm_min;
  *sec = tm->tm_sec;
}
