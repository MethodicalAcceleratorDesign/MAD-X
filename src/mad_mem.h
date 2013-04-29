#ifndef MAD_MEM_H
#define MAD_MEM_H

#ifdef _USEGC

#include <gc.h>

#define mymalloc(fn, sz)          mad_mem_check_ptr(fn, GC_MALLOC(sz)        )
#define mycalloc(fn, n, sz)       mad_mem_check_ptr(fn, GC_MALLOC((n)*(sz))  )
#define myrealloc(fn, p, sz)      mad_mem_check_ptr(fn, GC_REALLOC((p),(sz)) )
#define mymalloc_atomic(fn, sz)   mad_mem_check_ptr(fn, GC_MALLOC_ATOMIC(sz) )
#define myfree(fn, p)

#else

#include <stdlib.h>

#define mymalloc(fn, sz)          mad_mem_check_ptr(fn, malloc(sz)           )
#define mycalloc(fn, n, sz)       mad_mem_check_ptr(fn, calloc((n),(sz))     )
#define myrealloc(fn, p, sz)      mad_mem_check_ptr(fn, realloc((p),(sz))    )
#define mymalloc_atomic(fn, sz)   mad_mem_check_ptr(fn, malloc(sz)           )
#define myfree(fn, p)

#endif

// SIGSEGV handler
void  mad_mem_handler(int sig);

// --- inliners ---------------------------------------------------------------

static inline void*
mad_mem_check_ptr(const char *caller, void *ptr)
{
  if (!ptr)
    fatal_error("memory overflow, called from routine:", caller);

  return ptr;
}

#endif // MAD_MEM_H

