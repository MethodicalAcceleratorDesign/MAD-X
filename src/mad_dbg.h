#ifndef MAD_DBG_H
#define MAD_DBG_H

// always enable assertion
#undef  NDEBUG
#define NDEBUG 1
#include <assert.h>

#ifdef _ASSERT_DBG

/* special assert that loops to let the debugger to catch the process through
   its PID and perform a backtrace. */

#undef  assert
#define assert(c) ((void)( (c) || (__assert_fail(#c, __FILE__, __LINE__, __func__),1) ))

#endif // MAD_ASSERT_DBG

void __assert_fail(const char *assertion, const char *file, int line, const char *function);

#endif // MAD_DBG_H
