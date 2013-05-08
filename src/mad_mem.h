#ifndef MAD_MEM_H
#define MAD_MEM_H

#ifdef _USEGC

// #define GC_DEBUG

#include <gc.h>

#define mymalloc(fn, sz)            myptrchk(fn, GC_MALLOC_IGNORE_OFF_PAGE(sz))
#define mymalloc_atomic(fn, sz)     myptrchk(fn, GC_MALLOC_ATOMIC_IGNORE_OFF_PAGE(sz))
#define myrealloc(fn, p, sz)        myptrchk(fn, GC_REALLOC((p),(sz)))
#define myfree(fn, p)               ((void)((void)fn, (p)=0))

#define mycalloc(fn, n, sz)         memset(mymalloc(fn, (n)*(sz)), 0, (n)*(sz))
#define mycalloc_atomic(fn, n, sz)  memset(mymalloc_atomic(fn, (n)*(sz)), 0, (n)*(sz))
#define myrecalloc(fn, p, osz, sz)  ((void*)((char*)memset((char*)myptrchk(fn,GC_REALLOC((p),(sz)))+(osz),0,(sz)-(osz))-(osz)))

#else

#define mymalloc(fn, sz)            myptrchk(fn, malloc(sz))
#define mymalloc_atomic(fn, sz)     myptrchk(fn, malloc(sz))
#define myrealloc(fn, p, sz)        myptrchk(fn, realloc((p),(sz)))
#define myfree(fn, p)               ((void)(free(p), (void)fn, (p)=0))

#define mycalloc(fn, n, sz)         myptrchk(fn, calloc((n),(sz)))
#define mycalloc_atomic(fn, n, sz)  myptrchk(fn, calloc((n),(sz)))
#define myrecalloc(fn, p, osz, sz)  ((void*)((char*)memset((char*)myptrchk(fn,realloc((p),(sz)))+(osz),0,(sz)-(osz))-(osz)))

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
