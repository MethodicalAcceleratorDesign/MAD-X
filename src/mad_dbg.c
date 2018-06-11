#include <stdio.h>
#include "mad_dbg.h"

/* special assert that loops to let the debugger to catch the process through
   its PID and perform a backtrace. */

static volatile int __assert_foo = 0;

void __assert_fail(const char *assertion, const char *file, int line,
                   const char *function)
{
  fprintf(stderr, "MY assertion failed: %s in %s, %s:%d\n",
          assertion, function, file, line);
  while (__assert_foo == 0); // LOOP!!!
}
