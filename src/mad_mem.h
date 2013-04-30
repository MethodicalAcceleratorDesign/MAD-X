#ifndef MAD_MEM_H
#define MAD_MEM_H

#ifdef _USEGC

// #define GC_DEBUG
#include <gc.h>
#include <string.h>

#define mymalloc(fn, sz)            myptrchk(fn, GC_MALLOC(sz)        )
#define mycalloc(fn, n, sz)         myptrchk(fn, GC_MALLOC((n)*(sz))  )
#define myrealloc(fn, p, sz)        myptrchk(fn, GC_REALLOC((p),(sz)) )
#define mymalloc_atomic(fn, sz)     myptrchk(fn, GC_MALLOC_ATOMIC(sz) )
#define mycalloc_atomic(fn, n, sz)  memset(myptrchk(fn, GC_MALLOC_ATOMIC((n)*(sz))),0,(n)*(sz))
#define myfree(fn, p)               (void)(GC_FREE(p), fn)
#define mymemcheck()                GC_gcollect()

#else

#include <stdlib.h>

#define mymalloc(fn, sz)            myptrchk(fn, malloc(sz)        )
#define mycalloc(fn, n, sz)         myptrchk(fn, calloc((n),(sz))  )
#define myrealloc(fn, p, sz)        myptrchk(fn, realloc((p),(sz)) )
#define mymalloc_atomic(fn, sz)     myptrchk(fn, malloc(sz)        )
#define mycalloc_atomic(fn, n, sz)  myptrchk(fn, calloc((n),(sz))  )
#define myfree(fn, p)               (void)(free(p), fn)
#define mymemcheck()              

#endif

// SIGSEGV handler
void  mad_mem_handler(int sig);

// --- inliners ---------------------------------------------------------------

static inline void*
myptrchk(const char *caller, void *ptr)
{
  if (!ptr)
    fatal_error("memory overflow, called from routine:", caller);

  return ptr;
}

#endif // MAD_MEM_H

