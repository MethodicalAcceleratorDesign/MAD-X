/* src/config.h.  Generated from config.h.in by configure.  */
/* src/config.h.in.  Generated from configure.ac by autoheader.  */

/* Inline assembly available (only VC/x86_64) */
/* #undef AO_ASM_X64_AVAILABLE */

/* Assume Windows Server 2003, Vista or later target (only VC/x86) */
/* #undef AO_ASSUME_VISTA */

/* Assume target is not old AMD Opteron chip (only x86_64) */
/* #undef AO_CMPXCHG16B_AVAILABLE */

/* Define to avoid C11 atomic intrinsics even if available. */
/* #undef AO_DISABLE_GCC_ATOMICS */

/* Force test_and_set to use SWP instruction instead of LDREX/STREX (only arm
   v6+) */
/* #undef AO_FORCE_USE_SWP */

/* Force compare_and_swap definition via fetch_compare_and_swap */
/* #undef AO_GENERALIZE_ASM_BOOL_CAS */

/* No pthreads library available */
/* #undef AO_NO_PTHREADS */

/* Assume target is not sparc v9+ (only sparc) */
/* #undef AO_NO_SPARC_V9 */

/* Assume ancient MS VS Win32 headers (only VC/arm v6+, VC/x86) */
/* #undef AO_OLD_STYLE_INTERLOCKED_COMPARE_EXCHANGE */

/* Prefer C11 atomic intrinsics over assembly-based implementation even in
   case of inefficient implementation (do not use assembly for any atomic_ops
   primitive if C11/GCC atomic intrinsics available) */
/* #undef AO_PREFER_BUILTIN_ATOMICS */

/* Prefer generalized definitions to direct assembly-based ones */
/* #undef AO_PREFER_GENERALIZED */

/* Trace AO_malloc/free calls (for debug only) */
/* #undef AO_TRACE_MALLOC */

/* Assume single-core target (only arm v6+) */
/* #undef AO_UNIPROCESSOR */

/* Assume Win32 _Interlocked primitives available as intrinsics (only VC/arm)
   */
/* #undef AO_USE_INTERLOCKED_INTRINSICS */

/* Use nanosleep() instead of select() (only if atomic operations are
   emulated) */
/* #undef AO_USE_NANOSLEEP */

/* Do not block signals in compare_and_swap (only if atomic operations are
   emulated) */
/* #undef AO_USE_NO_SIGNALS */

/* Use Pentium 4 'mfence' instruction (only x86) */
/* #undef AO_USE_PENTIUM4_INSTRS */

/* Emulate atomic operations via slow and async-signal-unsafe pthread locking
   */
/* #undef AO_USE_PTHREAD_DEFS */

/* Prefer GCC built-in CAS intrinsics in favor of inline assembly (only
   gcc/x86, gcc/x86_64) */
/* #undef AO_USE_SYNC_CAS_BUILTIN */

/* Use Win32 Sleep() instead of select() (only if atomic operations are
   emulated) */
/* #undef AO_USE_WIN32_PTHREADS */

/* Emulate double-width CAS via pthread locking in case of no hardware support
   (only gcc/x86_64, the emulation is unsafe) */
/* #undef AO_WEAK_DOUBLE_CAS_EMULATION */

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the `getpagesize' function. */
#define HAVE_GETPAGESIZE 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have a working `mmap' system call. */
#define HAVE_MMAP 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/param.h> header file. */
#define HAVE_SYS_PARAM_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Define to disable assertion checking. */
#define NDEBUG 1

/* Name of package */
#define PACKAGE "libatomic_ops"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "https://github.com/ivmai/libatomic_ops/issues"

/* Define to the full name of this package. */
#define PACKAGE_NAME "libatomic_ops"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "libatomic_ops 7.6.2"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "libatomic_ops"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "7.6.2"

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Version number of package */
#define VERSION "7.6.2"

/* Indicates the use of pthreads (NetBSD). */
/* #undef _PTHREADS */

/* Required define if using POSIX threads. */
#define _REENTRANT 1
