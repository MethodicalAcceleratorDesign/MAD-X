/* timex, timest, etc Eric McIntosh 28/9/99 */
#include <sys/types.h>
#include <time.h>
#include <sys/times.h>
#include <sys/param.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdio.h>
#include <unistd.h>

int usleep_(int *timl)
{
  int blah;
  printf("sleep time: %d\n",*timl);
  blah = sleep(*timl);
  return blah;
}
