/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file   : mdb.h
 * purpose: definitions for general use, for mdblib, and for mdbmth.
 *
 * Michael Borland, 1988
 $Log: not supported by cvs2svn $
 Revision 1.92  2006/02/08 23:02:24  soliday
 Added the random_oag, gauss_rn_oag and gauss_rn_lim_oag functions.

 Revision 1.91  2005/11/16 19:04:09  shang
 added wofz() -- complex error function for mdbmth

 Revision 1.90  2005/11/04 22:46:59  soliday
 Updated code to be compiled by a 64 bit processor.

 Revision 1.89  2005/04/07 19:33:21  borland
 Added WILDCARD_MATCH ability to match_string.  Used by sddsxref -wildMatch option.

 Revision 1.88  2005/03/03 17:22:19  soliday
 Updated to work on WIN32

 Revision 1.87  2005/02/02 16:07:27  soliday
 Moved a few routines from mdblib to mdbcommon

 Revision 1.86  2004/12/21 20:14:24  shang
 added routines for finding files between dates and sort_and_return_index() function.

 Revision 1.85  2004/12/17 20:34:51  soliday
 Updated declaration for rawread, tt_attach, tt_detach

 Revision 1.84  2004/12/03 17:42:28  soliday
 Put the sys/types.h and sys/stat.h includes back in because Linux was
 complaining without them.

 Revision 1.83  2004/12/02 23:02:01  soliday
 The section of code inside the _MATCH_STRING_ ifdef in both the match_string.h
 and mdb.h is now the same. So it will not matter which one is included first.

 Revision 1.82  2004/11/04 16:12:23  shang
 added functions for reading file links, gettting file state and checking if
 file is modified.

 Revision 1.81  2004/07/16 16:28:48  shang
 added countLimit argument to despikeData()

 Revision 1.80  2004/04/01 14:39:44  borland
 Added SIMPLEX_VERBOSE_LEVEL1 and SIMPLEX_VERBOSE_LEVEL2 flags for simplexMin()
 and simplexMinimization().

 Revision 1.79  2004/02/27 16:29:15  borland
 Added prototype for trapizoidIntegration1().

 Revision 1.78  2003/12/19 19:31:50  soliday
 Added strcmp_nh funciton.

 Revision 1.77  2003/11/07 16:47:45  borland
 Added prototype for solveQuadratic().

 Revision 1.76  2003/10/07 16:21:57  shang
 added modified bessel functions: dbeskv

 Revision 1.75  2003/08/28 22:38:38  soliday
 Added some definitions for vxWorks

 Revision 1.74  2003/08/28 15:50:14  soliday
 If MIN or MAX are already defined it will not redefine them.

 Revision 1.73  2003/07/23 16:21:20  soliday
 Added isinf for WIN32

 Revision 1.72  2003/07/22 21:09:45  soliday
 Removed reference to IEEE.

 Revision 1.71  2003/07/22 20:01:18  soliday
 Added support for Kylix.

 Revision 1.70  2003/07/09 19:16:51  soliday
 Fixed bessel function prototypes.

 Revision 1.69  2003/06/18 19:52:32  borland
 Removed the query() and query_e() macros, and replaced with queryn()
 and queryn_e() macros.  The former used gets() whereas the latter do
 not.

 Revision 1.68  2003/06/14 21:35:00  borland
 Added prototypes for double-precision bessel functions in mdbmth/dbessel.c

 Revision 1.67  2003/01/21 18:55:04  borland
 simplexMin() has another argument: a factor by which to multiply the
 parameter ranges at the end of a pass to get the initial step sizes for
 the next pass.  A value of 1.0 reproduces the old behavior.

 Revision 1.66  2003/01/16 20:07:05  soliday
 Added optimAbort

 Revision 1.65  2003/01/15 22:59:18  borland
 Added SIMPLEX_START_FROM_VERTEX1 flag for simplexMin().
 Added simplexDivisor argument for simplexMin(); default value is 3
 to reproduce old behavior.

 Revision 1.64  2003/01/08 22:42:26  borland
 Added prototype for binaryArraySearch().

 Revision 1.63  2003/01/08 19:33:18  borland
 Updated prototype for binaryIndexSearch().

 Revision 1.62  2002/10/28 17:10:47  shang
 added sorting function declarations

 Revision 1.61  2002/09/24 21:03:59  borland
 Added prototypes for randomSampleMin() and randomWalkMin().  Added
 target argument for grid_sample_opt() and grid_search_min().

 Revision 1.60  2002/09/09 19:33:30  soliday
 Added the missing TouchFile declaration.

 Revision 1.59  2002/08/15 16:52:08  soliday
 *** empty log message ***

 Revision 1.58  2002/08/14 15:40:15  soliday
 Added Open License

 Revision 1.57  2002/07/24 20:41:51  shang
 added search path related and file generation related functions

 Revision 1.56  2002/07/13 21:06:11  borland
 Added SIMPLEX_RANDOM_SIGNS macro.

 Revision 1.55  2002/06/26 16:28:05  soliday
 Fixed the round definition for negative numbers.

 Revision 1.54  2002/06/19 23:12:14  borland
 Fixed spelling of Savitzky-Golay routines.

 Revision 1.53  2002/06/18 13:43:40  borland
 Added prototype for trapazoidIntegration().

 Revision 1.52  2002/02/26 03:11:01  borland
 Added prototypes for NAFF functions (fftpackC.h) and
 1d parabolic optimization routine (mdb.h).  These are both due
 to removing the NAFF functions from sddsnaff.c

 Revision 1.51  2002/02/18 17:43:25  borland
 Added prototype for simplexMinAbort.

 Revision 1.50  2002/01/07 21:32:27  borland
 Added prototype for approximate_percentiles() function.

 Revision 1.49  2001/10/12 21:07:03  soliday
 Fixed function declaration for WIN32 and Linux.

 Revision 1.48  2001/09/28 19:51:45  shang
 add prototype of substituteTagValue()

 Revision 1.47  2001/07/31 20:47:25  borland
 Added prototype for OneDScanOptimize() and updated prototype for
 simplexMin().

 Revision 1.46  2001/07/13 14:51:19  soliday
 Added replaceString declaration.

 Revision 1.45  2001/05/31 03:18:27  borland
 Changes for simplexMin().

 Revision 1.44  2000/11/21 19:08:11  borland
 Updated prototype for powellMin().

 Revision 1.43  2000/11/06 17:59:15  borland
 Added prototype for powellMin()

 Revision 1.42  2000/11/04 17:49:03  borland
 Changed prototype for nextHaltonSequencePoint and added prototype for
 startHaltonSequence.

 Revision 1.41  2000/11/02 21:26:00  borland
 Revised prototypes for zeroIntHalve, zeroNewton, zeroInterp, and
 nextHaltonSequencePoint.

 Revision 1.40  2000/11/02 19:43:58  borland
 Added prototypes for nextHaltonSequencePoint and randomizeOrder.

 Revision 1.39  2000/10/31 21:25:27  soliday
 Fixed the declaration of computeMode so that it works with WIN32.

 Revision 1.38  2000/10/11 21:45:55  soliday
 Changed definition of isinf so that the sunmath library is no longer needed.

 Revision 1.37  2000/10/07 01:16:02  borland
 Added prototype for computeMedian function.

 Revision 1.36  2000/08/17 21:27:22  soliday
 computeCorrelations is now used by elegant so declaration was changed so that
 it can be called inside a DLL.

 Revision 1.35  2000/08/10 21:10:11  soliday
 Added definition for isinf on Solaris with gcc

 Revision 1.34  2000/08/09 21:59:58  borland
 Changed prototype for savitzky-golay smoother.

 Revision 1.33  2000/04/19 17:00:10  soliday
 Borland C no longer includes mdbtc.h

 Revision 1.32  2000/04/17 20:25:22  soliday
 Added option to define binaryInsert with Borland C.

 Revision 1.31  2000/04/17 19:27:11  soliday
 Removed binaryInsert prototype when compiling with Borland C.

 Revision 1.30  2000/04/13 16:10:27  soliday
 Changed WIN32 to _WIN32

 Revision 1.29  2000/04/11 16:19:24  soliday
 Modified prototypes to work with new mdbcommon library.

 Revision 1.28  2000/04/06 22:24:53  soliday
 Added support for Borland C.

 Revision 1.27  2000/03/27 20:25:59  borland
 Added prototype for random_4().

 Revision 1.26  2000/01/18 20:47:22  soliday
 Renamed compress to compressString to avoid a conflict with ZLIB.

 Revision 1.25  1999/09/14 18:04:58  soliday
 Added export commands for WIN32 dll files.

 Revision 1.24  1999/08/03 17:55:42  soliday
 Added keep_alloc_record declaration

 Revision 1.23  1999/07/29 21:23:23  borland
 Added prototype and macro for differential equation routine.

 Revision 1.22  1999/07/22 15:35:00  soliday
 Added macros for fopen modes.

 Revision 1.21  1999/07/09 14:24:42  soliday
 Borland added shiftedLinearCorrelationCoefficient

 Revision 1.20  1999/07/01 19:25:34  borland
 Added prototypes for Savitzky-Golay filters.

 Revision 1.19  1999/05/25 18:45:59  soliday
 Altered sleep macro for WIN32, also defined popen and pclose for WIN32

 Revision 1.18  1999/05/04 14:46:42  borland
 Added some WIN32-specific conditional compilation statements.

 Revision 1.17  1999/01/07 21:45:39  borland
 Modified prototypes for simplexMin and simplexMinimization.

 Revision 1.16  1998/08/26 14:49:41  borland
 Treatment of IEEE math function isinf is now uniform.  If on solaris
 and sunmath is missing, then modify mdb.h.

 Revision 1.15  1997/03/27 22:20:30  borland
 Added TimeToEpochText().

 Revision 1.14  1997/02/05 20:50:48  saunders
 Added 'extern "C" {}' and renamed some arguments in func prototypes.

 Revision 1.13  1997/02/03 21:21:09  borland
 Added prototype and definitions for renameRobust().

 Revision 1.12  1996/10/22 18:48:17  borland
 Added prototypes for poisson statistics significance level routines.

 Revision 1.11  1996/10/07 17:28:46  borland
 Changed prototype for despikeData to reflect long integer return value.

 Revision 1.10  1996/08/26 20:08:47  borland
 Added prototypes for new mdblib routines tokenIsInteger() and tokenIsNumber().

 Revision 1.9  1996/08/16 20:03:18  borland
 Added prototype for normSigLevel().

 * Revision 1.8  1996/03/28  04:58:43  borland
 * Changed time conversion routine prototypes (long integers are now short
 * integers).
 *
 * Revision 1.7  1996/03/19  23:59:51  borland
 * Added prototypes for new time conversion functions.
 *
 * Revision 1.6  1995/12/12  03:16:46  borland
 * Changed prototype for linearCorrelationCoefficient(); added prototype for
 * compute_percentiles().
 *
 * Revision 1.5  1995/12/02  02:15:15  borland
 * Added prototype for strslide().
 *
 * Revision 1.4  1995/11/13  16:08:35  borland
 * Added prototype for function mtime().
 *
 * Revision 1.3  1995/09/12  03:19:44  borland
 * Added prototypes for wild_match_ci, strchr_ci, strcmp_ci
 *
 * Revision 1.2  1995/09/05  21:15:14  saunders
 * First test release of the SDDS1.5 package.
 *
 */
#ifndef _MDB_
#define _MDB_ 1

#include <math.h>
#include <float.h>
#include <limits.h>
#include <time.h>
#if defined(_WIN32)
typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
#define PRId32 "ld"
#define SCNd32 "ld"
#define INT32_MAX (2147483647)
#else
#include <inttypes.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#if defined(SUNOS4)
#define SUN_SPARC 1
#endif

#if defined(linux)
#define LINUX 1
#endif

  /*
#if !(defined(IEEE_MATH) && (defined(SUNOS4) || defined(SOLARIS) || defined(LINUX)))
#define isinf(x) (0)
#endif
  */

#if defined(SOLARIS)
#include <ieeefp.h>
#define isinf(x) ((x==x) && !finite(x))
#endif
#if defined(_WIN32)
#define isnan(x) _isnan(x)
#define isinf(x) (0)
#endif
#if defined(vxWorks)
#include <private/mathP.h>
#define isinf(x) isInf(x)
#define isnan(x) isNan(x)
#endif

#include <string.h>
#include <stdio.h>

#define FOPEN_WRITE_MODE "wb"
#define FOPEN_READ_MODE  "rb"
#define FOPEN_READ_AND_WRITE_MODE "r+b"

#undef epicsShareFuncMDBLIB
#undef epicsShareFuncMDBMTH
#undef epicsShareFuncMDBCOMMON
#if (defined(_WIN32) && !defined(__CYGWIN32__)) || (defined(__BORLANDC__) && defined(__linux__))
#if defined(EXPORT_MDBLIB)
#define epicsShareFuncMDBLIB  __declspec(dllexport)
#else
#define epicsShareFuncMDBLIB
#endif
#if defined(EXPORT_MDBMTH)
#define epicsShareFuncMDBMTH  __declspec(dllexport)
#else
#define epicsShareFuncMDBMTH
#endif
#if defined(EXPORT_MDBCOMMON)
#define epicsShareFuncMDBCOMMON  __declspec(dllexport)
#else
#define epicsShareFuncMDBCOMMON
#endif
#else
#define epicsShareFuncMDBLIB
#define epicsShareFuncMDBMTH
#define epicsShareFuncMDBCOMMON
#endif
  
#if defined(SUNOS4) && defined(GNU_C)
/* prototypes for functions not defined in stdio: */
extern int printf(const char *format_spec, ...);
extern int fprintf(FILE *file_ptr, const char *format_spec, ...);
/* int sprintf(char *str, const char *format_spec, ...); */
extern int scanf(const char *format_spec, ...);
extern int fscanf(FILE *file_ptr, const char *format_spec, ...);
extern int sscanf(char *str, const char *format_spec, ...);
extern int fputs(const char *string, FILE *file_ptr);
extern int puts(const char *string);
extern int fputc(char c, FILE *file_ptr);
extern int fclose(FILE *file_ptr);
extern int close(int file_descriptor);
extern void perror(char *s);
extern int fseek(FILE *file_ptr, long offset, int direction);
extern int fread(void *data, int size, int number, FILE *fp);
extern int fwrite(void *data, int size, int number, FILE *fp);
extern int fflush(FILE *file_ptr);
/* prototypes for functions not fully prototyped in math.h: */
extern double   acos(double x);
extern double   asin(double x);
extern double   atan(double x);
extern double   atan2(double y, double x);
extern double   ceil(double x);
extern double   cos(double x);
extern double   cosh(double x);
extern double   exp(double x);
extern double   fabs(double x);
extern double   floor(double x);
extern double   fmod(double x, double y);
extern double   frexp(double value, int *expo);
extern double   ldexp(double value, int expo);
extern double   log(double x);
extern double   log10(double x);
extern double   modf(double value, double *iptr);
extern double   pow(double x, double y);
extern double   sin(double x);
extern double   sinh(double x);
extern double   sqrt(double x);
extern double   tan(double x);
extern double   tanh(double x);
#endif

/* double-precision Bessel functions */
#define EPS 1.0e-16
#define FPMIN 1.0e-30
#define MAXIT 10000
#define XMIN 2.0
#define PI 3.141592653589793
epicsShareFuncMDBMTH double dbesi0(double x);
epicsShareFuncMDBMTH double dbesi1(double x);
epicsShareFuncMDBMTH double dbesj0(double x);
epicsShareFuncMDBMTH double dbesj1(double x);
epicsShareFuncMDBMTH double dbesk0(double x);
epicsShareFuncMDBMTH double dbesk1(double x);
epicsShareFuncMDBMTH double dbesy0(double x);
epicsShareFuncMDBMTH double dbesy1(double x);
/*modified bessel function dbeskv */
epicsShareFuncMDBMTH double chebev(double a, double b, double c[], int m, double x);
epicsShareFuncMDBMTH void beschb(double x, double *gam1, double *gam2, double *gampl, double *gammi);
epicsShareFuncMDBMTH void dbeskv(double x, double xnu, double *ri, double *rk, double *rip, double *rkp);

#include <stdlib.h>

epicsShareFuncMDBLIB long PackSuffixType(char *filename, char **unpackedName, unsigned long mode);
epicsShareFuncMDBLIB void *array_1d(long size_of_elem, long lower_index, long upper_index);
epicsShareFuncMDBLIB void **array_2d(long size_of_elem, long lower1, long upper1, long lower2,
                long upper2);
epicsShareFuncMDBLIB int free_array_1d(void *array, long size_of_elem, long lower_index,
                long upper_index);
epicsShareFuncMDBLIB int free_array_2d(void **array, long size_of_elem, long lower1, long upper1,
                long lower2, long upper2);
epicsShareFuncMDBLIB void **zarray_2d(long size, long n1, long n2);
epicsShareFuncMDBLIB void **resize_zarray_2d(long size, long old_n1, long old_n2,
            void **array, long n1, long n2);
epicsShareFuncMDBLIB int free_zarray_2d(void **array, long n1, long n2);
epicsShareFuncMDBLIB void zero_memory(void *memory, long n_bytes);
epicsShareFuncMDBLIB int tfree(void *ptr);
epicsShareFuncMDBLIB void keep_alloc_record(char *filename);

epicsShareFuncMDBLIB void fill_int_array(int *array, long n, int value);
epicsShareFuncMDBLIB void fill_short_array(short *array, long n, short value);
epicsShareFuncMDBLIB void fill_long_array(long *array, long n, long value);
epicsShareFuncMDBLIB void fill_float_array(float *array, long n, float value);
epicsShareFuncMDBLIB void fill_double_array(double *array, long n, double value);

epicsShareFuncMDBLIB void *tmalloc(unsigned long size_of_block);
epicsShareFuncMDBLIB void *trealloc(void *ptr, unsigned long size_of_block);

/* String-related macro definitions: */
#define chop_nl(m_s) ( ((m_s)[strlen(m_s)-1]=='\n') ? (m_s)[strlen(m_s)-1]=0 : 0)

#define queryn(s, t, n) ((*(t)=0),fputs(s,stdout),fgets(t,n,stdin),chop_nl(t))
#define queryn_e(s, t, n) ((*(t)=0),fputs(s,stderr),fgets(t,n,stdin),chop_l(t))

#define is_yes(c) ((c)=='y' || (c)=='Y')
#define is_no(c) ((c)=='n' || (c)=='N')

/*   -- Data-scanning routines: */
epicsShareFuncMDBLIB extern long   query_long(char *prompt, long default_value);
epicsShareFuncMDBLIB  extern int    query_int(char *prompt, int default_value);
epicsShareFuncMDBLIB  extern short  query_short(char *prompt, short default_value);
epicsShareFuncMDBLIB extern double query_double(char *prompt, double default_value);
epicsShareFuncMDBLIB extern float  query_float(char *prompt, float default_value);
epicsShareFuncMDBLIB extern int   get_double(double *target, char *source);
epicsShareFuncMDBLIB extern int   get_long(long *target, char *source);
epicsShareFuncMDBLIB extern int   get_short(short *target, char *source);
epicsShareFuncMDBLIB extern int   get_int(int *target, char *source);
epicsShareFuncMDBLIB extern int   get_float(float *target, char *source);
epicsShareFuncMDBLIB extern char  *get_token(char *source);
epicsShareFuncMDBLIB extern char  *get_token_buf(char *source, char *buffer, long buffer_length);
epicsShareFuncMDBLIB extern char  *get_token_t(char *source, char *token_delimiters);
epicsShareFuncMDBLIB extern char  *get_token_tq(char *source, char *token_start,
                           char *token_end, char *quote_start, char *quote_end);
epicsShareFuncMDBLIB long tokenIsInteger(char *token);
epicsShareFuncMDBLIB long tokenIsNumber(char *token);

/*   -- String routines: */
epicsShareFuncMDBLIB extern char *trim_spaces(char *s);
epicsShareFuncMDBLIB  extern char *replace_chars(char *string, char *from, char *to);
epicsShareFuncMDBLIB  extern char *rcdelete(char *string, char c_lower, char c_upper);
epicsShareFuncMDBLIB extern char *compressString(char *string, char *chars_to_compress);
epicsShareFuncMDBLIB extern char *delete_chars(char *s, char *t);
epicsShareFuncMDBLIB extern char *delete_bounding(char *string, char *chars_to_delete);
epicsShareFuncMDBLIB extern char *str_toupper(char *string);
epicsShareFuncMDBLIB extern char *str_tolower(char *string);
epicsShareFuncMDBLIB extern long is_blank(char *string);
epicsShareFuncMDBLIB extern char *str_in(char *string, char *sub_string);
epicsShareFuncMDBLIB  extern char *str_inn(char *string, char *sub_string, long n_char_to_check);
epicsShareFuncMDBLIB  extern char *insert(char *place_to_insert, char *string_to_insert);
epicsShareFuncMDBLIB extern char *pad_with_spaces(char *s, int n_spaces);
epicsShareFuncMDBLIB extern char *cp_str(char **target, char *source);
epicsShareFuncMDBLIB  extern char *cpn_str(char **target, char *source, long n_characters);
epicsShareFuncMDBLIB extern long edit_string(char *text, char *edit);
epicsShareFuncMDBLIB extern void edit_strings(char **string, long strings, char *buffer, char *edit);
epicsShareFuncMDBLIB char *strslide(char *s, long distance);

/* ---search path routines-- */
epicsShareFuncMDBLIB extern void setSearchPath(char *input);
epicsShareFuncMDBLIB extern char *findFileInSearchPath(char *filename);

/* --file stat routines-- */
#include <sys/types.h>
#include <sys/stat.h>
epicsShareFuncMDBLIB extern char *dir_name (const char *path);
epicsShareFuncMDBLIB extern char *read_file_link(const char *filename);
epicsShareFuncMDBLIB extern char *read_file_lastlink(const char *filename);
epicsShareFuncMDBLIB extern char *read_last_link_to_file(const char *filename);

epicsShareFuncMDBLIB extern long get_file_stat(const char *filename, const char *lastlink, struct stat *filestat);
epicsShareFuncMDBLIB extern long file_is_modified(const char *inputfile, char **final_file, struct stat *input_stat);

/* -- find files routines ---*/
#include <ctype.h>
#if !defined(_WIN32)
#include <dirent.h>
#endif
epicsShareFuncMDBCOMMON extern short make_four_digit_year (short year);
epicsShareFuncMDBCOMMON extern long is_leap_year (short year);
epicsShareFuncMDBCOMMON extern char **find_files_between_dates(char *directory, char *rootname, char *suffix, 
                            short startYear, short startMonth, short startDay, short startJDay,
                            short endYear, short endMonth, short endDay, short endJDay,
                            char *filter, char **extensionList, long extensions,
                            long tailsOnly, long *files, long increaseOrder);
void sort_files_by_start_time(char *directory, long isTail, char **fileList, long files, long increaseOrder);

epicsShareFuncMDBCOMMON extern char **ls_dir (char *path, char *matchstr, long tailsOnly, long *files);

#if !defined(__BORLANDC__) || defined(DefineBinaryInsert)
epicsShareFuncMDBLIB long binaryInsert(void **array, long members, void *newMember, 
             int (*compare)(void *c1, void *c2), int32_t *duplicate);
#endif
epicsShareFuncMDBLIB long binaryIndexSearch(void **array, long members, void *key, 
                       int (*compare)(void *c1, void *c2), long bracket);
epicsShareFuncMDBLIB long binaryArraySearch(void *array, size_t elemSize, long members, void *key, 
                                            int (*compare)(void *c1, void *c2), long bracket);

/* sort routines (previously sort.h) */
#if !defined(_MDBSORT_INCLUDED_)
#define _MDBSORT_INCLUDED_ 1 

/*following structs and function are moved from sddsxref.c for quick sorting. */
typedef struct {
  char *stringKey;
  double doubleKey;
  long rowIndex;
} KEYED_INDEX;

typedef struct {
  KEYED_INDEX **equivalent;
  long equivalents, nextIndex;
} KEYED_EQUIVALENT;

epicsShareFuncMDBLIB extern int CompareStringKeyedIndex(const void *ki1, const void *ki2);
epicsShareFuncMDBLIB extern int CompareDoubleKeyedIndex(const void *ki1, const void *ki2);
epicsShareFuncMDBLIB extern int CompareStringKeyedGroup(void *kg1, void *kg2);
epicsShareFuncMDBLIB extern int CompareDoubleKeyedGroup(void *kg1, void *kg2);
epicsShareFuncMDBLIB extern KEYED_EQUIVALENT **MakeSortedKeyGroups(long *keyGroups, long keyType, void *data, long points);
epicsShareFuncMDBLIB extern long FindMatchingKeyGroup(KEYED_EQUIVALENT **keyGroup, long keyGroups, long keyType,void *searchKeyData, long reuse);
epicsShareFuncMDBLIB extern long *sort_and_return_index(void *data, long type, long rows, long increaseOrder);

/* sort routines (previously sort.h) */
epicsShareFuncMDBLIB extern int double_cmpasc(const void *a, const void *b);
extern int double_cmpdes(const void *a, const void *b);
extern void double_copy(void *a, void *b);
extern int float_cmpasc(const void *a, const void *b);
extern int float_cmpdes(const void *a, const void *b);
extern void float_copy(void *a, void *b);
epicsShareFuncMDBLIB extern int long_cmpasc(const void *a, const void *b);
extern int long_cmpdes(const void *a, const void *b);
extern void long_copy(void *a, void *b);
epicsShareFuncMDBLIB extern int string_cmpasc(const void *a, const void *b);
extern int string_cmpdes(const void *a, const void *b);
epicsShareFuncMDBLIB extern void string_copy(void *a, void *b);
epicsShareFuncMDBLIB extern int row_compare(const void *a, const void *b);
extern void row_copy(void *a, void *b);
epicsShareFuncMDBLIB extern void set_up_row_sort(int sort_by_column, size_t n_columns,
    size_t element_size,
    int (*compare)(const void *a, const void *b));
epicsShareFuncMDBLIB extern int unique(void *base, size_t n_items, size_t size,
    int (*compare)(const void *a, const void *b),
    void (*copy)(void *a, void *b));
#endif

/* string array matching (previously match_string.h): */
#if !defined(_MATCH_STRING_)
epicsShareFuncMDBLIB extern long match_string(char *string, char **option_list, long n_options,
                        long match_mode_flags);

epicsShareFuncMDBLIB extern int strncmp_case_insensitive(char *s1, char *s2, long n);
epicsShareFuncMDBLIB extern int strcmp_case_insensitive(char *s1, char *s2);
#if defined(_WIN32)
#if defined(__BORLANDC__)
#define strcasecmp(s, t) stricmp(s, t)
#define strncasecmp(s, t, n) strnicmp(s, t, n)
#else
#define strcasecmp(s, t) _stricmp(s, t)
#define strncasecmp(s, t, n) _strnicmp(s, t, n)
#endif
#endif
#define _MATCH_STRING_ 1


#define DCL_STYLE_MATCH 0
#define UNIQUE_MATCH DCL_STYLE_MATCH
#define CASE_SENSITIVE 1
#define MATCH_WHOLE_STRING 2
#define RETURN_FIRST_MATCH 8
#define EXACT_MATCH (CASE_SENSITIVE|MATCH_WHOLE_STRING|RETURN_FIRST_MATCH)
#define WILDCARD_MATCH 16

#endif

epicsShareFuncMDBLIB extern char *clean_filename(char *filename);
epicsShareFuncMDBLIB extern long fexists(char *filename);
#define RENAME_OVERWRITE 0x0001UL
extern long renameRobust(char *oldName, char *newName, unsigned long flags);

extern char *exp_notation(double x, long n1, long n2);
extern void add_to_headers(char **header, long n_headers, char **item,
    long min_width, long format_index);
long replaceFile(char *file, char *replacement);
epicsShareFuncMDBLIB long replaceFileAndBackUp(char *file, char *replacement);
epicsShareFuncMDBLIB void add_to_standard_headers(char *name_header, char *unit_header,
    char *printf_string, char *new_name, char *new_unit, char *new_format,
    long min_width);
extern long format_length(char *format_specifier);
extern char *sbinary(char *s, int len, long n);
epicsShareFuncMDBLIB extern long bitsSet(unsigned long data);
epicsShareFuncMDBLIB extern void interpret_escapes(char *s);
epicsShareFuncMDBLIB extern int replace_string(char *target, char *source, char *orig, char *newOne);
epicsShareFuncMDBLIB int replace_stringn(char *t, char *s, char *orig, char *repl, long count_limit);
epicsShareFuncMDBLIB void interpret_escaped_quotes(char *s);
epicsShareFuncMDBLIB extern int replaceString(char *t, char *s, char *orig, 
			 char *repl, long count_limit, long here);

extern char **wild_list(int *n_names_ret, int **origin_ret, char **item_list,
        int num_items);
epicsShareFuncMDBLIB int wild_match(char *string, char *tmplate);
epicsShareFuncMDBLIB int wild_match_ci(char *string, char *tmplate);
char *strchr_ci(char *s, char c);
epicsShareFuncMDBLIB int strcmp_ci(const char *s, const char *t);
epicsShareFuncMDBLIB char *expand_ranges(char *tmplate);
epicsShareFuncMDBLIB int has_wildcards(char *tmplate);
epicsShareFuncMDBLIB char *unescape_wildcards(char *tmplate);
epicsShareFuncMDBLIB int strcmp_nh(const char *s, const char *t);

/*   -- Routines for flagging and aborting on errors: */
epicsShareFuncMDBLIB extern void bomb(char *error_message, char *usage_message);
epicsShareFuncMDBLIB extern long bombre(char *error_message, char *usage_message, long return_value);
extern long err_mess(long status, char *routine_name, char *message);
extern long err_mess_sys(long status, char *routine_name, char *message);
extern void fatal_err(long error_code, char *error_message);

/*   -- IO routines: */
extern char *ffgets(char *target, long target_length, FILE *file_pointer);
epicsShareFuncMDBLIB extern FILE *fopen_e(char *file_name, char *open_mode, long error_mode);
#define FOPEN_EXIT_ON_ERROR 0
#define FOPEN_RETURN_ON_ERROR 1
#define FOPEN_INFORM_OF_OPEN  2
#define FOPEN_SAVE_IF_EXISTS  4
epicsShareFuncMDBLIB extern void backspace(long n);

/*   -- Run-time statistics: */
epicsShareFuncMDBLIB extern void init_stats(void);
epicsShareFuncMDBLIB extern void report_stats(FILE *fp, char *label);
epicsShareFuncMDBLIB extern double delapsed_time();
epicsShareFuncMDBLIB extern long memory_count(void);

/* terminal IO routines: */
epicsShareFuncMDBLIB extern int tt_attach(void);
epicsShareFuncMDBLIB extern char rawread(void);
epicsShareFuncMDBLIB extern void tt_detach(void);

/*   -- Miscellaneous routines: */
extern long log_usage(char *log_file_name, char *program_name);
epicsShareFuncMDBLIB extern char *tmpname(char *target);
epicsShareFuncMDBLIB char *mtime(void);
short IsLeapYear(short year);
short JulianDayFromMonthDay(short month, short day, short year, short *julianDay);
short MonthDayFromJulianDay(short julianDay, short year, short *month, short *day);
epicsShareFuncMDBLIB short TimeEpochToBreakdown(short *year, short *jDay, short *month, short *day, double *hour, double epochTime);
epicsShareFuncMDBLIB short TimeEpochToText(char *text, double epochTime);
epicsShareFuncMDBLIB short TimeBreakdownToEpoch(short year, short jDay, short month, short day, double hour, double *epochTime);


/********************* end of mdblib routines *****************************/

/************************ mathematical stuff *******************************/
#ifndef _MDB_MTH_
#define _MDB_MTH_ 1

#define FABS(x) fabs(x)
#define SIGN(x) ((x)<0?-1:((x)>0?1:0))
#define IS_NEGATIVE(x) (SIGN(x)==-1)
#define IS_POSITIVE(x) (SIGN(x)==1)

/* mdbmth routines */
epicsShareFuncMDBMTH long gaussianQuadrature(double (*fn)(), double a, double b, long n, double err, double *result);
epicsShareFuncMDBCOMMON extern int fixcount(char *filename, long n_points);
epicsShareFuncMDBMTH extern long factorial(long n);
epicsShareFuncMDBMTH extern double dfactorial(long n);
epicsShareFuncMDBMTH extern double g_int(double (*function)(), double lower_limit,
                    double upper_limit, long n_panels, double error_limit);
epicsShareFuncMDBMTH extern double simpson(double (*function)(), double lower_limit,
                    double upper_limit, long n_panels);
epicsShareFuncMDBMTH long trapazoidIntegration(double *x, double *y, long n, double *integral);
epicsShareFuncMDBMTH long trapazoidIntegration1(double *x, double *y, long n, double *integral);
epicsShareFuncMDBMTH extern double ipow(double base, long power);
epicsShareFuncMDBMTH extern double amod(double x, double m);
epicsShareFuncMDBMTH extern double zeroInterp(double (*function)(), double value,
                                              double x_initial, double x_final,
                                              double x_step, double effective_zero);
epicsShareFuncMDBMTH extern double zeroIntHalve(double (*function)(), double value,
                                                double x_initial, double x_final,
                                                double x_step, double effective_zero);
epicsShareFuncMDBMTH extern double zeroNewton(double (*function)(), 
                                              double value, double x_initial, double dx_deriv,
                                              long n_passes, double effective_zero);
#if defined(__BORLANDC__)
epicsShareFuncMDBMTH extern double poly2(double *coef, long n_coefs, double x);
#else
epicsShareFuncMDBMTH extern double poly(double *coef, long n_coefs, double x);
#endif
epicsShareFuncMDBMTH extern double dpoly(double *coef, long n_coefs, double x);
epicsShareFuncMDBMTH extern double polyp(double *coef, long *power, long n_coefs, double x);
epicsShareFuncMDBMTH extern double dpolyp(double *coef, long *power, long n_coefs, double x);
epicsShareFuncMDBMTH extern int solveQuadratic(double a, double b, double c, double *solution);
epicsShareFuncMDBMTH extern double K_cei(double k);
epicsShareFuncMDBMTH extern double E_cei(double k);
epicsShareFuncMDBMTH extern double dK_cei(double k);
epicsShareFuncMDBMTH extern double dE_cei(double k);
epicsShareFuncMDBMTH extern double celliptic(double kc, double p, double a, double b);
epicsShareFuncMDBMTH extern float drand(long dummy);
epicsShareFuncMDBMTH extern double rdrand(double lower_limit, double upper_limit);
epicsShareFuncMDBMTH extern void r_theta_rand(double *r, double *theta, double r_min,
                           double r_max);
epicsShareFuncMDBMTH extern double random_1(long iseed);
epicsShareFuncMDBMTH extern double random_2(long iseed);
epicsShareFuncMDBMTH extern double random_3(long iseed);
epicsShareFuncMDBMTH extern double random_4(long iseed);
epicsShareFuncMDBMTH extern double urandom_gauss(long iseed);
epicsShareFuncMDBMTH extern double gauss_rn(long iseed, double (*urandom)(long iseed1));
epicsShareFuncMDBMTH extern double gauss_rn_lim(double mean, double sigma, double limit_in_sigmas, double (*urandom)(long iseed));

epicsShareFuncMDBMTH extern double random_oag(long iseed, long increment);
epicsShareFuncMDBMTH extern double gauss_rn_oag(long iseed, long increment, double (*urandom)(long iseed1, long increment));
epicsShareFuncMDBMTH extern double gauss_rn_lim_oag(double mean, double sigma, double limit_in_sigmas, long increment, double (*urandom)(long iseed, long increment));


epicsShareFuncMDBMTH extern long randomizeOrder(char *ptr, long size, long length, long iseed, double (*urandom)(long iseed1));
epicsShareFuncMDBMTH extern double nextHaltonSequencePoint(long ID);
epicsShareFuncMDBMTH extern int32_t startHaltonSequence(int32_t *radix, double value);
epicsShareFuncMDBMTH extern long convertSequenceToGaussianDistribution(double *data, long points, double limit);
epicsShareFuncMDBMTH extern double KS_Qfunction(double lambda);
extern double twoVariableKStest(double *d1, long n1, double *d2, long n2, double *MaxCDFerror);
epicsShareFuncMDBCOMMON extern double linearCorrelationCoefficient(double *data1, double *data2, 
                                           short *accept1, short *accept2, long rows, long *count);
epicsShareFuncMDBCOMMON extern double linearCorrelationSignificance(double r, long rows);
epicsShareFuncMDBCOMMON extern double shiftedLinearCorrelationCoefficient(double *data1, double *data2, 
                                                  short *accept1, short *accept2,
                                                  long rows, long *count, long shift);
epicsShareFuncMDBMTH extern double betaInc(double x, double a, double b);
epicsShareFuncMDBMTH extern double betaComp(double a, double b);
epicsShareFuncMDBMTH extern double gammaP(double a, double x);
epicsShareFuncMDBMTH extern double gammaQ(double a, double x);
epicsShareFuncMDBMTH extern double tTailSigLevel(double t0, long nu, long tails);
epicsShareFuncMDBMTH extern double FSigLevel(double var1, double var2, long nu1, long nu2);
epicsShareFuncMDBMTH extern double rSigLevel(double r0, long nu);
epicsShareFuncMDBMTH double ChiSqrSigLevel(double ChiSquared0, long nu);
epicsShareFuncMDBMTH double normSigLevel(double z0, long tails);
epicsShareFuncMDBMTH double poissonSigLevel(long n, double n0);

epicsShareFuncMDBMTH extern long is_prime(long number);
epicsShareFuncMDBMTH extern long smallest_factor(long number);
epicsShareFuncMDBMTH extern long next_prime_factor(long *number);
epicsShareFuncMDBMTH long largest_prime_factor(long number);
epicsShareFuncMDBMTH extern void copy_sp_array(float  **c, float  *o, long n);
epicsShareFuncMDBMTH extern void copy_dp_array(double **c, double *o, long n);
epicsShareFuncMDBMTH extern double bessel_Jn(double z, long n);
epicsShareFuncMDBMTH extern double bessel_Yn(double z, long n);
epicsShareFuncMDBMTH extern double Ai(double x), Bi(double x), Aip(double x),
        Bip(double x);
epicsShareFuncMDBMTH extern int wofz(double *xi, double *yi, double *u, double *v, long *flag);


epicsShareFuncMDBCOMMON long lsfn(double *xd, double *yd, double *sy,
    long nd, long nf, double *coef, double *s_coef,
    double *chi, double *diff);
epicsShareFuncMDBCOMMON long lsfp(double *xd, double *yd, double *sy,
    long n_pts, long n_terms, long *power, double *coef, double *s_coef,
    double *chi, double *diff);
epicsShareFuncMDBCOMMON long lsfg(double *xd, double *yd, double *sy,
    long n_pts, long n_terms, int32_t *order,
    double *coef, double *s_coef, double *chi, double *diff,
    double (*fn)(double x, long ord));

/*functions for generation file names, moved from SDDSepics.c, May 8, 2002 */
/*i.e. the declarations of functions in logfile_generation.c */
epicsShareFuncMDBCOMMON extern double computeYearStartTime(double StartTime);
epicsShareFuncMDBCOMMON extern void getTimeBreakdown(double *Time, double *Day, 
                                                     double *Hour, double *JulianDay,
                                                     double *Year, double *Month, 
                                                     char **TimeStamp);
epicsShareFuncMDBCOMMON extern void makeTimeBreakdown(double Time, double *ptrTime,
                                                      double *ptrDay, double *ptrHour,
                                                      double *ptrJulianDay, double *ptrYear,
                                                      double *ptrMonth, char **ptrTimeStamp);
epicsShareFuncMDBCOMMON extern char *makeTimeStamp(double Time);
epicsShareFuncMDBCOMMON extern double getTimeInSecs(void);
epicsShareFuncMDBCOMMON extern double getHourOfDay(void);
epicsShareFuncMDBCOMMON extern char *getHourMinuteSecond ();
epicsShareFuncMDBCOMMON extern void checkGenerationFileLocks(char *match_date);
epicsShareFuncMDBCOMMON extern char *MakeGenerationFilename(char *rootname, long digits, 
                                                            char *delimiter, char *lastFile);
epicsShareFuncMDBCOMMON extern char *MakeDailyGenerationFilename(char *rootname, long digits,
                                                                 char *delimiter, long timetag);
epicsShareFuncMDBCOMMON extern char *MakeSCRDailyTimeGenerationFilename(char *rootname);
#define DEFAULT_GENERATIONS_DIGITS 4
#define USE_TIMETAG   0x0010UL
epicsShareFuncMDBCOMMON extern void usleepSystemIndependent(long usec);


#define DIFFEQ_EXIT_COND_FAILED -4
#define DIFFEQ_ZERO_STEPSIZE -3
#define DIFFEQ_CANT_TAKE_STEP -2
#define DIFFEQ_OUTSIDE_INTERVAL -1
#define DIFFEQ_XI_GT_XF 0
#define DIFFEQ_SOLVED 1
#define DIFFEQ_SOLVED_ALREADY 1
#define DIFFEQ_ZERO_FOUND 2
#define DIFFEQ_END_OF_INTERVAL 3
epicsShareFuncMDBMTH char *diffeq_result_description(long return_code);

epicsShareFuncMDBMTH extern long rk_odeint(
    double *y_i, void (*derivs)(double *yp, double *y, double x),
    long n_eq, double *accuracy, long *accmode, double *tiny, long *misses,
    double *x0, double xf, double x_accuracy, double h_start, double h_max,
    double *h_rec, double (*exfn)(double *yp, double *y, double x),
    double exit_accuracy, long n_to_skip,
    void (*store_data)(double *yp, double *y, double x, double exval)  );
epicsShareFuncMDBMTH extern long rk_odeint1(
    double *y_i, void (*derivs)(double *yp, double *y, double x),
    long n_eq, double *accuracy, long *accmode,
    double *tiny, long *misses, double *x0, double xf, double x_accuracy,
    double h_start, double h_max, double *h_rec );
epicsShareFuncMDBMTH extern long rk_odeint2(
    double *y_i, void (*derivs)(double *qp, double *q, double t), long n_eq,
    double *accuracy, long *accmode,
    double *tiny, long *misses, double *x0, double xf, double x_accuracy,
    double h_start, double h_max, double *h_rec, double exit_value,
    long i_exit_value, double exit_accuracy, long n_to_skip );
epicsShareFuncMDBMTH extern long rk_odeint3(
    double *y_i, void (*derivs)(double *yp, double *y, double x), long n_eq,
    double *accuracy, long *accmode, double *tiny, long *misses,
    double *x0, double xf, double x_accuracy, double h_start, double h_max,
    double *h_rec, double (*exfn)(double *yp, double *y, double x),
    double exit_accuracy);
epicsShareFuncMDBMTH extern long rk_odeint4(
    double *y_i, void (*derivs)(double *qp, double *q, double t), long n_eq,
    double *accuracy, long *accmode,
    double *tiny, long *misses, double *x0, double xf, double x_accuracy,
    double h_start, double h_max, double *h_rec, double exit_value,
    long i_exit_value, double exit_accuracy, long n_to_skip,
    void (*store_data)(double *qp, double *q, double t, double exfn) );
epicsShareFuncMDBMTH extern long rk_odeint_na(
    double *y_i, void (*derivs)(double *yp, double *y, double x),
    long n_eq, double *null1, long *null2, double *null3, long *null4,
    double *x0, double xf, double dummy1,
    double h, double dummy2, double *dummy3 );
epicsShareFuncMDBMTH extern long rk_odeint3_na(
    double *y_i, void (*derivs)(double *yp, double *y, double x), long n_eq,
    double *accuracy, long *accmode, double *tiny, long *misses,
    double *x0, double xf, double x_accuracy, double h_start, double h_max,
    double *h_rec, double (*exfn)(double *yp, double *y, double x),
    double exit_accuracy);
epicsShareFuncMDBMTH extern long bs_odeint(
    double *y_i, void (*derivs)(double *yp, double *y, double x), long n_eq,
    double *accuracy, long *accmode, double *tiny, long *misses,
    double *x0, double xf, double x_accuracy, double h_start, double h_max,
    double *h_rec,
    double (*exfn)(double *yp, double *y, double x), double exit_accuracy, long n_to_skip,
    void (*store_data)(double *qp, double *q, double t, double exfn_value)  );
epicsShareFuncMDBMTH extern long bs_odeint1(
    double *y_i, void (*derivs)(double *yp, double *y, double x), long n_eq,
    double *accuracy, long *accmode, double *tiny, long *misses,
    double *x0, double xf, double x_accuracy, double h_start, double h_max,
    double *h_rec );
epicsShareFuncMDBMTH extern long bs_odeint2(
    double *y_i, void (*derivs)(double *qp, double *q, double t), long n_eq,
    double *accuracy, long *accmode, double *tiny, long *misses,
    double *x0, double xf, double x_accuracy, double h_start, double h_max,
    double *h_rec,
    double exit_value, long i_exit_value, double exit_accuracy, long n_to_skip );
epicsShareFuncMDBMTH extern long bs_odeint3(
    double *y_i, void (*derivs)(double *yp, double *y, double x), long n_eq,
    double *accuracy, long *accmode, double *tiny, long *misses,
    double *x0, double xf, double x_accuracy, double h_start, double h_max,
    double *h_rec, double (*exfn)(double *yp, double *y, double x),
    double exit_accuracy);
epicsShareFuncMDBMTH extern long bs_odeint4(
    double *y_i, void (*derivs)(double *qp, double *q, double t), long n_eq,
    double *accuracy, long *accmode,
    double *tiny, long *misses, double *x0, double xf, double x_accuracy,
    double h_start, double h_max, double *h_rec, double exit_value,
    long i_exit_value, double exit_accuracy, long n_to_skip,
    void (*store_data)(double *qp, double *q, double t, double exfn) );
epicsShareFuncMDBMTH extern long bss_odeint(
    float *y_i, void (*derivs)(float *yp, float *y, float x), long n_eq,
    float *accuracy, long *accmode, float *tiny, long *misses,
    float *x0, float xf, float x_accuracy, float h_start, float h_max,
    float *h_rec, float (*exfn)(float *yp, float *y, float x),
    float exit_accuracy, long n_to_skip,
    void (*store_data)(float *yp, float *y, float x, float exval)  );
epicsShareFuncMDBMTH extern long bss_odeint1(
    float *y_i, void (*derivs)(float *yp, float *y, float x), long n_eq,
    float *accuracy, long *accmode, float *tiny, long *misses, float *x0,
    float xf, float x_accuracy, float h_start, float h_max, float *h_rec );
epicsShareFuncMDBMTH extern long bss_odeint3(
    float *y_i, void (*derivs)(float *yp, float *y, float x), long n_eq,
    float *accuracy, long *accmode, float *tiny, long *misses,
    float *x0, float xf, float x_accuracy, float h_start, float h_max,
    float *h_rec, float (*exfn)(float *yp, float *y, float x),
    float exit_accuracy);
epicsShareFuncMDBMTH extern long rks_odeint(
    float *y_i, void (*derivs)(float *dydx, float *y, float x),
    long n_eq, float *accuracy, long *accmode, float *tiny, long *misses,
    float *x0, float xf, float x_accuracy, float h_start, float h_max,
    float *h_rec, float (*exit_func)(float *dydx, float *y, float x),
    float exit_accuracy, long n_to_skip,
    void (*store_data)(float *dydx, float *y, float x, float exval)  );
epicsShareFuncMDBMTH extern long rks_odeint1(
    float *y_i, void (*derivs)(float *yp, float *y, float x), long n_eq,
    float *accuracy, long *accmode, float *tiny, long *misses, float *x0,
    float xf, float x_accuracy, float h_start, float h_max, float *h_rec);
epicsShareFuncMDBMTH extern long rks_odeint2(
    float *y_i, void (*derivs)(float *yp, float *y, float x), long n_eq,
    float *accuracy, long *accmode, float *tiny, long *misses, float *x0,
    float xf, float x_accuracy, float h_start, float h_max, float *h_rec,
    float exit_value, long i_exit_value, float exit_accuracy, long n_to_skip);
epicsShareFuncMDBMTH extern long rks_odeint_na(
    float *y_i, void (*derivs)(float *yp, float *y, float x),
    long n_eq, float *null1, long *null2, float *null3, long *null4,
    float *x0, float xf, float dummy1,
    float h, float dummy2, float *dummy3 );
epicsShareFuncMDBMTH extern long stiff_odeint1(
    double *y_i, void (*derivs)(double *yp, double *y, double x),
    long n_eq, double *accuracy, long *accmode,
    double *tiny, long *misses, double *x0, double xf, double x_accuracy,
    double h_start, double h_max, double *h_rec );
epicsShareFuncMDBMTH void mmid(double *y, double *dydx, long nvar, double xs, double htot,
    long nstep, double *yout, void (*derivs)(double *dydxa, double *ya, double xa)
    );
epicsShareFuncMDBMTH void mmid2(double *y, double *dydx, long nvar, double xs, double htot,
    long nstep, double *yout, void (*derivs)(double *dydxa, double *ya, double xa)
    );
epicsShareFuncMDBMTH extern long mmid_odeint3_na(
    double *y_i, void (*derivs)(double *yp, double *y, double x), long n_eq,
    double *accuracy, long *accmode, double *tiny, long *misses,
    double *x0, double xf, double x_accuracy, double h_start, double h_max,
    double *h_rec, double (*exfn)(double *yp, double *y, double x),
    double exit_accuracy);

epicsShareFuncMDBMTH void smoothData(double *data, long rows, long smoothPoints, long smoothPasses);
epicsShareFuncMDBMTH long despikeData(double *data, long rows, long neighbors, long passes, long averageOf,
    double threshold, long countLimit);
epicsShareFuncMDBCOMMON void SavitzkyGolayCoefficients(double *coef, long maxCoefs,
                              long order, long nLeft, long nRight,
                              long derivativeOrder, long wrapAround);
epicsShareFuncMDBCOMMON long SavitzkyGolaySmooth(double *data, long rows,
                        long order, long nLeft, long nRight, long derivativeOrder);
epicsShareFuncMDBCOMMON void TouchFile(char *filename);

#define SavitzyGolaySmooth(data, rows, order, nLeft, nRight, derivativeOrder) \
SavitzkyGolaySmooth(data, rows, order, nLeft, nRight, derivativeOrder)

epicsShareFuncMDBMTH extern long optimAbort(long abort);
epicsShareFuncMDBMTH extern double minc(double (*fn)(double *param), double *x, double *dx,
    double *dx_lim, double *xlo, double *xhi, long np, long ns_max,
    long p_flag);

epicsShareFuncMDBMTH void set_argument_offset(double offset);
epicsShareFuncMDBMTH void set_argument_scale(double scale);
epicsShareFuncMDBMTH double dtcheby(double x, long n);
epicsShareFuncMDBMTH double tcheby(double x, long n);
epicsShareFuncMDBMTH double ipower(double x, long n);
epicsShareFuncMDBMTH double dipower(double x, long n);
epicsShareFuncMDBMTH double eval_sum(double (*fn)(double x, long ord), double *coef, int32_t *order, long n_coefs, double x0);
epicsShareFuncMDBMTH long powellMin(double *yReturn, double *xGuess, double *dxGuess, double *xLowerLimit,
                                    double *xUpperLimit, long dims, double target, double tolerance,
                                    double (*func)(double *x, long *invalid), 
                                    void (*report)(double ymin, double *xmin, long pass, long evals, long dims),
                                    long maxPasses, long maxEvaluations, long linMinIterations);
#define SIMPLEX_NO_1D_SCANS        0x0001UL
#define SIMPLEX_RANDOM_SIGNS       0x0002UL
#define SIMPLEX_START_FROM_VERTEX1 0x0004UL
#define SIMPLEX_VERBOSE_LEVEL1     0x0008UL
#define SIMPLEX_VERBOSE_LEVEL2     0x0010UL
epicsShareFuncMDBMTH long simplexMinAbort(long abort);
epicsShareFuncMDBMTH long simplexMin(double *yReturn, double *xGuess, double *dxGuess, double *xLowerLimit,
                double *xUpperLimit, short *disable,
                long dimensions, double target, double tolerance, double (*func)(double *x, long *invalid), 
                void (*report)(double ymin, double *xmin, long pass, long n_evals, long n_dim),
                long maxEvaluations, long maxPasses, long maxDivisions, double divisorFactor, 
                double passRangeFactor, unsigned long flags);
epicsShareFuncMDBMTH long simplexMinimization(double **simplexVector, double *fValue, double *coordLowerLimit,
                         double *coordUpperLimit, short *disable, long dimensions, long activeDimensions,
                         double target, double tolerance, long tolerance_mode,
                         double (*function)(double *x, long *invalid), long maxEvaluations, 
                         long *evaluations, unsigned long flags);
#define ONEDSCANOPTIMIZE_REFRESH 0x0001UL
epicsShareFuncMDBMTH long OneDScanOptimize (double *yReturn, double *xGuess,double *dxGuess,
                 double *xLowerLimit, double *xUpperLimit,short *disable,
                 long dimensions, 
                 double target,              /* will return if any value is <= this */
                 double tolerance,           /* <0 means fractional, >0 means absolute */
                 double (*func)(double *x, long *invalid), 
                 void (*report)(double ymin, double *xmin, long pass, long evals, long dims),
                 long maxSteps,
                 long maxDivsions, long maxRepeats,
                 unsigned long flags);                                            

epicsShareFuncMDBMTH long OneDParabolicOptimization
(double *yReturn, double *xGuess, double dx,
 double xLower, double xUpper, 
 double (*func)(double x, long *invalid),
 long maxCycles, double dxLimit, double tolerance,
 long maximize);

epicsShareFuncMDBMTH long grid_search_min(double *best_result, double *best_x, double *lower, double *upper,
    double *step, long n_dimen, double target, double (*func)(double *x, long *invalid));
epicsShareFuncMDBMTH long grid_sample_min(double *best_result, double *best_x, double *lower, double *upper,
    double *step, long n_dimen, double target, double (*func)(double *x, long *invalid), double sample_fraction);
epicsShareFuncMDBMTH long randomSampleMin(double *best_result, double *best_x,
                                          double *lower, double *upper, long n_dimen,
                                          double target,
                                          double (*func)(double *x, long *invalid), long nSamples);
epicsShareFuncMDBMTH long randomWalkMin(double *best_result, double *best_x,
                                          double *lower, double *upper, double *range, long n_dimen,
                                          double target,
                                          double (*func)(double *x, long *invalid), long nSamples);

epicsShareFuncMDBMTH long advance_values(double *value, long *value_index, double *initial, double *step, long n_values,
    long *counter, long *max_count, long n_indices);
epicsShareFuncMDBMTH long advance_counter(long *counter, long *max_count, long n_indices);

epicsShareFuncMDBMTH extern long compute_average(double *value, double *x, long n);
epicsShareFuncMDBMTH extern long compute_middle(double *value, double *x, long n);
epicsShareFuncMDBMTH extern long compute_median(double *value, double *x, long n);
epicsShareFuncMDBMTH extern long compute_percentile(double *value, double *x, long n, double percentile);
epicsShareFuncMDBMTH extern long compute_percentiles(double *value, double *percent, long values, double *x, long n);
epicsShareFuncMDBMTH extern long approximate_percentiles(double *value, double *percent, long values, double *x, long n, long bins);

epicsShareFuncMDBMTH extern long find_average(double *value, double *x, long n);
epicsShareFuncMDBMTH extern long find_middle(double *value, double *x, long n);
epicsShareFuncMDBMTH extern long find_median(double *value, double *x, long n);
epicsShareFuncMDBMTH extern long find_percentile(double *value, double *x, long n, double percentile);
epicsShareFuncMDBMTH extern long find_median_of_row(double *value, double **x, long index, long n);

epicsShareFuncMDBMTH extern long make_histogram(double *hist, long n_bins, double lo, double hi, double *data,
    long n_pts, long new_start);
epicsShareFuncMDBMTH extern long make_histogram_weighted(double *hist, long n_bins, double lo, double hi, double *data,
    long n_pts, long new_start, double *weight);
epicsShareFuncMDBMTH long computeMode(double *result, double *data, long pts, double binSize, long bins);

epicsShareFuncMDBMTH long findCrossingPoint(long start, double *data, long points, double level, long direction,
                       double *interpData, double *result);
epicsShareFuncMDBMTH long findTopBaseLevels(double *top, double *base, double *data, long points,
                       long bins, double sigmasRequired);

epicsShareFuncMDBMTH extern double standardDeviation(double *x, long n);
epicsShareFuncMDBMTH long unweightedLinearFit(double *xData, double *yData, long nData, double *slope, double *intercept, double *variance);
epicsShareFuncMDBMTH extern double rmsValue(double *y, long n);
epicsShareFuncMDBMTH extern double arithmeticAverage(double *y, long n);
epicsShareFuncMDBMTH extern double meanAbsoluteDeviation(double *y, long n);
epicsShareFuncMDBMTH extern long computeMoments(double *mean, double *rms, double *standardDev,
          double *meanAbsoluteDev, double *x, long n);
epicsShareFuncMDBMTH extern long computeCorrelations(double *C11, double *C12, double *C22, double *x, double *y, long n);
epicsShareFuncMDBMTH extern long computeWeightedMoments(double *mean, double *rms, double *standardDev,
          double *meanAbsoluteDev, double *x, double *w, long n);
extern long accumulateMoments(double *mean, double *rms, double *standardDev,
          double *x, long n, long reset);
extern long accumulateWeightedMoments(double *mean, double *rms, double *standardDev,
          double *x, double *w, long n, long reset);
epicsShareFuncMDBMTH extern double weightedAverage(double *y, double *w, long n);
epicsShareFuncMDBMTH extern double weightedRMS(double *y, double *w, long n);
epicsShareFuncMDBMTH extern double weightedMAD(double *y, double *w, long n);
epicsShareFuncMDBMTH extern double weightedStDev(double *y, double *w, long n);

epicsShareFuncMDBMTH extern int find_min_max(double *min, double *max, double *list, long n);
epicsShareFuncMDBMTH extern int index_min_max(long *imin, long *imax, double *list, long n);
extern int assign_min_max(double *min, double *max, double val);
extern int find_min_max_2d(double *min, double *max, double **value,
            long n1, long n2);
extern int find_min_max_2d_float(float *min, float *max, float **value,
    long n1, long n2);
extern int find_min(double *min, double *loc, double *c1, double *c2, long n);
extern int find_max(double *max, double *loc, double *c1, double *c2, long n);
epicsShareFuncMDBMTH extern double max_in_array(double *array, long n);
epicsShareFuncMDBMTH extern double min_in_array(double *array, long n);

/*interpolate functions from interp.c */
typedef struct {
  double value;
  unsigned long flags;
} OUTRANGE_CONTROL;
#define OUTRANGE_VALUE       0x00000001
#define OUTRANGE_SKIP        0x00000002
#define OUTRANGE_SATURATE    0x00000004
#define OUTRANGE_EXTRAPOLATE 0x00000008
#define OUTRANGE_ABORT       0x00000010
#define OUTRANGE_WARN        0x00000020
#define OUTRANGE_WRAP        0x00000040
epicsShareFuncMDBMTH extern double interpolate(double *f, double *x, long n, double xo, 
                                               OUTRANGE_CONTROL *belowRange, 
                                               OUTRANGE_CONTROL *aboveRange, 
                                               long order, unsigned long *returnCode, long M);
epicsShareFuncMDBMTH extern double interp(double *y, double *x, long n, double x0, long warn, long order, long *returnCode);
int interpolate_minimum(double *fmin, double *zmin, double *value, double z_lo,
    double z_hi, long n);
epicsShareFuncMDBMTH double LagrangeInterp(double *x, double *f, long order, double x0, long *returnCode);

epicsShareFuncMDBLIB extern void substituteTagValue(char *input, long buflen, 
                        char **macroTag, char **macroValue, long macros); 

#define iceil(x) ((int)ceil(x))
#define round(x) ( x < 0.0 ? ((int)((x)-.5)) : ((int)((x)+.5)) )

#ifndef MIN
#define MIN(x,y) ( ((x)>(y)) ? (y) : (x))
#endif
#ifndef MAX
#define MAX(x,y) ( ((x)<(y)) ? (y) : (x))
#endif

#define SWAP_LONG(x, y) {long tmp_swap_long; tmp_swap_long=(x); (x)=(y); (y)=tmp_swap_long; }
#define SWAP_INT(x, y) {int tmp_swap_int; tmp_swap_int=(x); (x)=(y); (y)=tmp_swap_int; }
#define SWAP_SHORT(x, y) {short tmp_swap_short; tmp_swap_short=(x); (x)=(y); (y)=tmp_swap_short; }
#define SWAP_DOUBLE(x, y) {double tmp_swap_double; tmp_swap_double=(x); (x)=(y); (y)=tmp_swap_double; }
#define SWAP_FLOAT(x, y) {float tmp_swap_float; tmp_swap_float=(x); (x)=(y); (y)=tmp_swap_float; }
#define SWAP_PTR(x, y) {void *tmp_swap_ptr; tmp_swap_ptr=(x); (x)=(y); (y)=tmp_swap_ptr; }

#define INTERPOLATE(y1, y2, x1, x2, x0) (((y2)-(y1))/((x2)-(x1))*((x0)-(x1)) + (y1))

#define sqr(x)  ipow(x, 2)
#define pow2(x) sqr(x)
#define pow3(x) ipow(x, 3)
#define pow4(x) ipow(x, 4)
#define pow5(x) ipow(x, 5)

/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file   : constants.h
 * purpose: mathematical and physical constants
 *
 * Michael Borland, 1991
 $Log: not supported by cvs2svn $
 Revision 1.7  2005/01/27 20:15:36  borland
 Restored previous version.  Changes resulted in too many annoying discrepancies.

 Revision 1.5  2002/08/14 15:40:14  soliday
 Added Open License

 Revision 1.4  1999/01/06 16:45:13  emery
 Added classical radius of electron.

 Revision 1.3  1997/08/25 19:24:05  borland
 Undefines macro PI before trying to define it.

 * Revision 1.2  1995/09/05  21:14:58  saunders
 * First test release of the SDDS1.5 package.
 *
 */

/* mathematical constants from Abramowitz and Stegun */
#undef PI
#define PI   3.141592653589793238462643
#define PIx2 6.283185307179586476925287
#define PIo2 1.570796326794896619231322
#define SQRT2 1.4142135623730950488

/* physical constants from 1988 Particle Properties Data Booklet */
    /* speed of light */
#define c_cgs	(2.99792458e10)
#define c_mks   (2.99792458e8)
    /* unit charge */
#define e_cgs   (4.8032068e-10)
#define e_mks   (1.60217733e-19)
    /* electron mass */
#define me_cgs (9.1093897e-28)
#define me_mks (9.1093897e-31)
#define me_mev (0.51099906)
    /* classical electron radius */
#define re_cgs (2.81794092E-13)
#define re_mks (2.81794092E-15)

    /* Planck's constant */
#define h_cgs (6.6260755e-27)
#define hbar_cgs (h_cgs/PIx2)
#define h_mks (6.6260755e-34)
#define hbar_mks (h_mks/PIx2)
    /* Boltzmann's constant */
#define k_boltzmann_cgs (1.380658e-16)
#define k_boltzmann_mks (1.380658e-23)
    /* mu and epsilon */
#define mu_o (4*PI*1e-7)
#define epsilon_o (1.e7/(4*PI*sqr(c_mks)))


#endif  /* _MDB_MTH_ */

#if defined(_WIN32)
#include <windows.h>
#define sleep(sec) Sleep(sec * 1000)
#define popen(x,y) _popen(x,y)
#define pclose(x) _pclose(x)
#endif

/* machine-specific include file: */
#ifdef VAX_VMS
#include "mdbvax.h"
#endif
#ifdef SUNOS4
#include "mdbsunos4.h"
#endif
#if defined(__TURBOC__) && !defined(__BORLANDC__)
#include "mdbtc.h"
#endif

#ifdef __cplusplus
}
#endif 

#endif /* _MDB_ */

/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/*
 $Log: not supported by cvs2svn $
 Revision 1.19  2005/02/22 17:30:21  shang
 added scanItemListLong() function to deal with unsigned long long type flags

 Revision 1.18  2004/03/16 23:26:35  borland
 Added SCANITEMLIST_IGNORE_VALUELESS macro.

 Revision 1.17  2003/07/22 20:01:18  soliday
 Added support for Kylix.

 Revision 1.16  2002/08/14 15:40:17  soliday
 Added Open License

 Revision 1.15  2002/03/22 22:53:18  soliday
 Replaced free_scanargs with free_scanargs2

 Revision 1.14  2002/03/21 23:10:47  soliday
 Added free_scanargs2

 Revision 1.13  2002/03/07 01:18:53  soliday
 Added parse_string

 Revision 1.12  2002/01/28 16:49:49  soliday
 Added free_scanargs.

 Revision 1.11  2000/07/19 16:10:10  soliday
 Added ability to call from C++

 Revision 1.10  2000/04/11 16:19:32  soliday
 Modified prototypes to work with new mdbcommon library.

 Revision 1.9  2000/01/18 19:59:55  soliday
 Added support for ZLIB.

 Revision 1.8  1999/09/14 18:06:21  soliday
 Added export commands for WIN32 dll files.

 Revision 1.7  1999/07/22 16:22:06  soliday
 Added contains_keyword_phrase

 Revision 1.6  1996/05/29 21:44:48  borland
 Added mode flags for scanItemLists().

 * Revision 1.5  1996/02/14  01:02:08  borland
 * Added prototype for scanItemList().
 *
 * Revision 1.4  1996/01/21  00:15:54  borland
 * Added bit flag definitions and prototypes for new versions of unpacking
 * routines.
 *
 * Revision 1.3  1996/01/19  00:18:11  borland
 * SDDS.h: Added popenUsed to SDDS_LAYOUT structure.
 * scan.h: Added prototypes for unpacking functions.
 *
 * Revision 1.2  1995/09/05  21:15:39  saunders
 * First test release of the SDDS1.5 package.
 *
*/

/* define structure for use with scanargs(), scanlist() */
#if !defined(SCAN_INCLUDED)
#define SCAN_INCLUDED 1

#undef epicsShareFuncMDBLIB
#undef epicsShareFuncMDBCOMMON
#if (defined(_WIN32) && !defined(__CYGWIN32__)) || (defined(__BORLANDC__) && defined(__linux__))
#if defined(EXPORT_MDBLIB)
#define epicsShareFuncMDBLIB  __declspec(dllexport)
#else
#define epicsShareFuncMDBLIB
#endif
#if defined(EXPORT_MDBCOMMON)
#define epicsShareFuncMDBCOMMON  __declspec(dllexport)
#else
#define epicsShareFuncMDBCOMMON
#endif
#else
#define epicsShareFuncMDBLIB
#define epicsShareFuncMDBCOMMON
#endif

#if defined(zLib)
#include "zlib.h"
#endif

#ifdef __cplusplus 
extern "C" {
#endif

typedef struct {
    long arg_type;	/* type of argument */
    long n_items;	/* number of items in list */
    char **list;	/* the list */
    } SCANNED_ARG;

/* possible values for arg_type */
#define OPTION 		1
#define A_LIST 		2

epicsShareFuncMDBCOMMON extern int scanargs(SCANNED_ARG **scanned, int argc, char **argv);
  /*
epicsShareFuncMDBCOMMON extern void free_scanargs(SCANNED_ARG *scanned, int argc);
  */
epicsShareFuncMDBCOMMON extern void free_scanargs(SCANNED_ARG **scanned, int argc);
epicsShareFuncMDBCOMMON extern int scanargsg(SCANNED_ARG **scanned, int argc, char **argv);
extern int parseList(char ***list, char *string);
extern int parse_string(char ***list, char *string);
epicsShareFuncMDBCOMMON long processPipeOption(char **item, long items, unsigned long *flags);
#define USE_STDIN      0x0001UL
#define USE_STDOUT     0x0002UL
#define DEFAULT_STDIN  0x0004UL
#define DEFAULT_STDOUT 0x0008UL

epicsShareFuncMDBCOMMON void processFilenames(char *programName, char **input, char **output, unsigned long pipeFlags, long noWarnings, long *tmpOutputUsed);


#include <stdio.h>
#define UNPACK_REQUIRE_SDDS 0x00000001UL
#define UNPACK_USE_PIPE     0x00000002UL
long PackSuffixType(char *filename, char **unpackedName,  unsigned long mode);
epicsShareFuncMDBLIB FILE *UnpackFopen(char *filename, unsigned long mode, short *popenUsed, char **tmpFileUsed);
#if defined(zLib)
epicsShareFuncMDBLIB gzFile *UnpackGZipOpen(char *filename);
#endif

/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: SDDS.h
 * purpose: SDDS data types
 *          format
 *
 * Michael Borland, 1993
 $Log: not supported by cvs2svn $
 Revision 1.6  2002/08/14 15:40:13  soliday
 Added Open License

 Revision 1.5  1999/02/19 22:52:26  borland
 Added SDDS_ANY_INTEGER_TYPE for use with SDDS_CheckXXX procedures.

 Revision 1.4  1997/12/19 16:55:46  borland
 Fixed SDDS_RowCount macro (more parentheses).  Added prototype for
 SDDS_Malloc.  Added new "type": SDDS_ANY_FLOATING_TYPE.

 * Revision 1.3  1995/09/06  14:12:01  saunders
 * First test release of SDDS1.5
 *
 */

#if !defined(_SDDSTYPES_)

#define _SDDSTYPES_ 1

#define SDDS_DOUBLE    1
#define SDDS_FLOAT     2
#define SDDS_LONG      3
#define SDDS_SHORT     4
#define SDDS_STRING    5
#define SDDS_CHARACTER 6
#define SDDS_NUM_TYPES 6
#define SDDS_INTEGER_TYPE(type) ((type)==SDDS_LONG || (type)==SDDS_SHORT)
#define SDDS_FLOATING_TYPE(type) ((type)==SDDS_DOUBLE || (type)==SDDS_FLOAT)
#define SDDS_NUMERIC_TYPE(type) (SDDS_INTEGER_TYPE(type) || SDDS_FLOATING_TYPE(type))
#define SDDS_VALID_TYPE(type) (type>=1 && type<=SDDS_NUM_TYPES)

/* used by SDDS_Check*() routines  */
#define SDDS_ANY_NUMERIC_TYPE (SDDS_NUM_TYPES+1)
#define SDDS_ANY_FLOATING_TYPE   (SDDS_NUM_TYPES+2)
#define SDDS_ANY_INTEGER_TYPE   (SDDS_NUM_TYPES+3)

#endif

epicsShareFuncMDBLIB extern long scanItemList(unsigned long *flags, char **item, long *items, unsigned long mode, ...);
epicsShareFuncMDBLIB extern long scanItemListLong(unsigned long long *flags, char **item, long *items, unsigned long mode, ...);
/* usage: scanItemList(&flags, item, &items, mode,
               <keyword>, <SDDS-type>, <pointer>, <number-required>, <set-flag>, etc.
               NULL)
 */
#define SCANITEMLIST_UNKNOWN_VALUE_OK    0x00000001UL
#define SCANITEMLIST_UNKNOWN_KEYVALUE_OK 0x00000002UL
#define SCANITEMLIST_REMOVE_USED_ITEMS   0x00000004UL
#define SCANITEMLIST_IGNORE_VALUELESS    0x00000008UL

epicsShareFuncMDBLIB extern long contains_keyword_phrase(char *string);
/* Obsolete: */
epicsShareFuncMDBLIB extern long scan_item_list(unsigned long *flags, char **item, long *items, ...);

/* usage: scan_item_list(&flags, item, &items, 
               <keyword>, <SDDS-type>, <pointer>, <number-required>, <set-flag>, etc.
               NULL)
 */

#ifdef __cplusplus
}
#endif

#endif

/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file   : table.h
 * purpose: provide definitions for use with routines get_table()
 *	    and put_table(), which read and write data in dpl format
 * definition of dpl format:
 *   -Files are ordinary text format; fortran carraige control is not
 *     recommended.
 *   -Lines in file:
 *     1:  label for x-axis (independent variable)
 *     2:  label for y-axis (dependent variable)
 *     3:  label for plot title
 *     4:  label for top of plot
 *     5:  N: integer number of data points that follow
 *     6:  x[0]       y[0]    {sigma_y[0]  |  {sigma_x[0]  sigma_y[0] } }
 *                      .
 *                      .
 *                      .
 *     N+5:  x[N-1]     y[N-1]  {sigma_y[N-1]  |  {sigma_x[N-1]  sigma_y[N-1] } }
 *     [EOF]
 *   -The data points are in free format, with no restriction except that
 *    non-data text should not contain the characters ., +, -, or 0-9.
 *   -Any line beginning with '!' will be ignored.
 *   -Lines beyond N+5 will be ignored.
 *
 * Michael Borland, 1988
 $Log: not supported by cvs2svn $
 Revision 1.7  2005/11/04 22:47:00  soliday
 Updated code to be compiled by a 64 bit processor.

 Revision 1.6  2003/07/22 20:01:18  soliday
 Added support for Kylix.

 Revision 1.5  2002/08/14 15:40:18  soliday
 Added Open License

 Revision 1.4  2000/04/11 16:19:45  soliday
 Modified prototypes to work with new mdbcommon library.

 Revision 1.3  1999/09/14 18:06:29  soliday
 Added export commands for WIN32 dll files.

 Revision 1.2  1995/09/05 21:15:43  saunders
 First test release of the SDDS1.5 package.

 */
#include <stdio.h>

#ifndef _TABLE_INCLUDED_
#define _TABLE_INCLUDED_ 1

#undef epicsShareFuncMDBCOMMON
#undef epicsShareFuncSDDS
#if (defined(_WIN32) && !defined(__CYGWIN32__)) || (defined(__BORLANDC__) && defined(__linux__))
#if defined(EXPORT_MDBCOMMON)
#define epicsShareFuncMDBCOMMON  __declspec(dllexport)
#else
#define epicsShareFuncMDBCOMMON
#endif
#if defined(EXPORT_SDDS)
#define epicsShareFuncSDDS  __declspec(dllexport)
#else
#define epicsShareFuncSDDS
#endif
#else
#define epicsShareFuncMDBCOMMON
#define epicsShareFuncSDDS
#endif

/* control bit-flags for get_table() */
#define SWAP 1
#define REVERSE 2
#define REORDER_ASCENDING 4
#define REORDER_DESCENDING 8
#define SAVE_SIGMA_ARRAYS 16
#define READ_LABELS_ONLY 32
#define SDDS_NOCOMPRESS_NAMES 64

typedef struct {
	double *c1, *c2;           /* arrays of data in cols 1 & 2 */
        double *s1, *s2;           /* sigmas of data in cols 1 & 2 */
    	char *xlab, *ylab;         /* axis labels */
        char  *topline, *title;    /* other plot labels */
        long flags;                 /* data description bit-flags */
#define SIGMA_X_PRESENT 1
#define SIGMA_Y_PRESENT 2
	long n_data;                /* number of data points */
	} TABLE;

typedef struct {
	float *c1, *c2;           /* arrays of data in cols 1 & 2 */
        float *s1, *s2;           /* sigmas of data in cols 1 & 2 */
    	char *xlab, *ylab;         /* axis labels */
        char  *topline, *title;    /* other plot labels */
        long flags;                 /* data description bit-flags */
#define SIGMA_X_PRESENT 1
#define SIGMA_Y_PRESENT 2
	long n_data;                /* number of data points */
	} TABLE_FLOAT;

epicsShareFuncMDBCOMMON extern long get_table(TABLE *tab, char *file, long sample_interval, long flags);
extern void put_table(char *file, TABLE *tab, char *format);
extern long get_table_float(TABLE_FLOAT *tab, char *file, 
                           long sample_interval, long flags);
extern void put_table_float(char *file, TABLE_FLOAT *tab, char *format);
extern double *double_array_from_float(float *f_array, long n_elements);
extern float  *float_array_from_double(double *d_array, long n_elements);
extern char *fgets_skip(char *s, long slen, FILE *fp, char skip_char, long skip_lines);
extern int fixcount(char *filename, long n_points);

epicsShareFuncSDDS int32_t SDDS_ReadIntoMplTable(TABLE *mpl_data, char *file, int32_t sample_interval, int32_t mpl_flags, char *SDDS_tags);
epicsShareFuncSDDS int32_t SDDS_WriteMplTable(TABLE *mpl_data, char *file);

#endif   /* _TABLE_INCLUDED_ */






/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: SDDS.h
 * purpose: define structures and prototypes for use with SDDS
 *          format
 *
 * Michael Borland, 1993
 $Log: not supported by cvs2svn $
 Revision 1.71  2006/03/16 18:50:02  soliday
 Added the SDDS_BreakIntoLockedFile function.

 Revision 1.70  2005/11/04 22:46:59  soliday
 Updated code to be compiled by a 64 bit processor.

 Revision 1.69  2005/01/13 16:52:51  shang
 added is_string argument to is_memory() function to get the memory data type.

 Revision 1.68  2004/11/15 16:59:19  soliday
 Moved and renamed the procedures to write non native binary data from
 sddsendian to SDDSlib.

 Revision 1.67  2004/05/13 18:16:08  soliday
 Added SDDS_AppendToArrayVararg

 Revision 1.66  2004/02/27 21:57:57  soliday
 Modified SDDS_SetRowsOfInterest

 Revision 1.65  2004/02/27 21:14:08  soliday
 Modified SDDS_RegisterProgramName

 Revision 1.64  2004/02/05 22:56:15  soliday
 Added missing declaration for SDDS_CheckTabularData

 Revision 1.63  2003/10/17 16:11:29  soliday
 Added SDDS_VerifyColumnExists, SDDS_VerifyParameterExists,
 and SDDS_VerifyArrayExists

 Revision 1.62  2003/09/23 21:49:23  soliday
 Added the SDDS_MATCH_EXCLUDE_STRING definition.

 Revision 1.61  2003/09/11 22:31:35  soliday
 Added SDDS_GetParameterAsLong and SDDS_ConvertToLong

 Revision 1.60  2003/07/22 20:01:16  soliday
 Added support for Kylix.

 Revision 1.59  2003/02/15 22:54:16  borland
 Added prototype for SDDS_DoFSync().

 Revision 1.58  2002/08/14 15:40:12  soliday
 Added Open License

 Revision 1.57  2002/07/25 00:10:37  borland
 Added definition of SDDS_PRINT_BUFLEN.

 Revision 1.56  2002/07/24 20:23:13  shang
 added SDDS_InitializeInputFromSearchPath()

 Revision 1.55  2002/06/15 02:28:00  borland
 Added SDDS_ForceActive().

 Revision 1.54  2002/06/05 16:20:41  shang
 added SDDS_GetParameterAsString()

 Revision 1.53  2002/05/26 02:07:12  borland
 Added autoRecovered variable to SDDS_DATASET to allow recording fact that
 auto-read-recovery has occurred.

 Revision 1.52  2002/02/18 19:27:47  soliday
 Added the ability to enable or disable fsync.

 Revision 1.51  2002/01/15 22:46:05  soliday
 Added SetRowCountMode, UpdateRowCount, SetAutoReadRecovery

 Revision 1.50  2002/01/11 17:08:19  soliday
 Added SDDS_SyncDataSet

 Revision 1.49  2002/01/07 21:33:16  borland
 Added prototypes for SDDS_GetRowLimit() and SDDS_SetRowLimit().

 Revision 1.48  2002/01/07 19:31:23  soliday
 Fixed declaration for SDDS_SetNoRowCounts.

 Revision 1.47  2001/12/22 21:21:32  soliday
 Fixed RW_ASSOCIATES typo.

 Revision 1.46  2001/12/21 03:29:39  soliday
 Fixed some declarations.

 Revision 1.45  2001/11/30 15:35:13  borland
 Changes by H. Shang: addition of SDDS_GotoPage() and changes required to
 support this.

 Revision 1.44  2001/02/19 18:01:25  borland
 Added prototype for SDDS_GetParameterByIndex().

 Revision 1.43  2000/11/27 17:15:24  borland
 Fixed typo in previous change.

 Revision 1.42  2000/11/27 17:06:27  borland
 Added prototype for SDDS_FreeStringData().

 Revision 1.41  2000/10/31 21:24:53  soliday
 Fixed the declarations of some new functions so that they work with WIN32.

 Revision 1.40  2000/10/24 01:16:45  borland
 Added prototype for SDDS_GetArrayInDoubles().

 Revision 1.39  2000/09/19 20:14:11  borland
 Added function prototypes and changed some structures (SDDS_DATASET and
 SDDS_LAYOUT) to support reading of little/big endian files.

 Revision 1.38  2000/06/20 18:02:42  borland
 Added prototype for SDDS_SetDefaultIOBufferSize.

 Revision 1.37  2000/06/20 14:28:37  soliday
 Moved some definitions from SDDS_internal.h to SDDS.h so that sddsendian
 can access them.

 Revision 1.36  2000/04/14 16:36:12  soliday
 gzFile now is defined even when we are not using the zLib library.

 Revision 1.35  2000/04/13 17:04:58  soliday
 Fixed prototype for SDDS_NumberOfErrors and SDDS_ClearErrors

 Revision 1.34  2000/03/30 19:58:35  borland
 Added SDDS_DisconnectFile() and SDDS_ReconnectFile(), plus disconnected member
 to SDDS_LAYOUT structure.

 Revision 1.33  2000/03/30 17:26:57  soliday
 Added SDDS_InitializeAppendToPage

 Revision 1.32  2000/03/27 20:25:47  borland
 Added prototype for SDDS_ShortenTable().

 Revision 1.31  2000/02/25 23:15:12  evans
 Added SDDS_Free.  Memory allocated with SDDS_Malloc and others should
 be freed with SDDS_Free.  It is necessary for WIN32 when degugging
 SDDS applications using the release version of the SDDS library
 because the memory allocation models for debug and release are
 inconsistent.

 Revision 1.30  2000/01/18 20:06:46  soliday
 Added support for ZLIB.

 Revision 1.29  1999/11/03 22:42:47  soliday
 Added WriteBinaryString declaration

 Revision 1.28  1999/09/16 22:03:15  soliday
 Modified the way global variables are exported and imported to a WIN32 dll

 Revision 1.27  1999/09/14 18:03:28  soliday
 Added export commands for WIN32 dll files.

 Revision 1.26  1999/06/01 19:17:05  borland
 Added SDDS_Calloc.

 Revision 1.25  1999/04/14 13:58:47  borland
 Fixed some function types and returns.

 Revision 1.24  1999/02/04 21:13:19  soliday
 Added prototype for SDDS_GetColumnsInLong

 Revision 1.23  1998/12/16 21:20:35  borland
 Modified prototypes for SDDS_TransferAll* functions and added required
 macros.

 Revision 1.22  1998/08/25 01:45:47  borland
 Moved SDDS_FILEBUFFER_SIZE definition here from SDDS_binary.c

 Revision 1.21  1998/06/30 21:30:48  borland
 Added prototype for SDDS_FreeStringArray().

 Revision 1.20  1997/12/19 16:55:45  borland
 Fixed SDDS_RowCount macro (more parentheses).  Added prototype for
 SDDS_Malloc.  Added new "type": SDDS_ANY_FLOATING_TYPE.

 Revision 1.19  1997/07/08 14:41:17  borland
 Added SDDS_RowCount macro.

 * Revision 1.18  1996/07/05  16:36:02  borland
 * Modified prototypes for SDDS_PrintTypedValue and SDDS_SprintTypedValue
 * to include new mode argument.
 *
 * Revision 1.17  1996/03/27  17:51:08  borland
 * Added prototype and macros for SDDS_GetErrorMessages().
 *
 * Revision 1.16  1996/03/21  20:07:32  borland
 * Added mode argument to prototype for SDDS_ReadPageSparse().
 *
 * Revision 1.15  1996/03/21  16:39:50  borland
 * Added definitions for major revision of SDDS library.  Added file buffer
 * structure to SDDS_DATASET.  Added prototype for SDDS_ReadPageSparse().
 *
 * Revision 1.14  1996/02/12  17:18:55  borland
 * Added prototype for SDDS_SetAutoCheckMode() and mode flags to go with it.
 *
 * Revision 1.13  1996/01/20  05:38:03  borland
 * Added prototypes for SDDS_AppendLayout().
 *
 * Revision 1.12  1996/01/19  00:18:09  borland
 * SDDS.h: Added popenUsed to SDDS_LAYOUT structure.
 * scan.h: Added prototypes for unpacking functions.
 *
 * Revision 1.11  1996/01/18  04:31:56  borland
 * Added prototype for SDDS_GetInternalColumn()
 *
 * Revision 1.10  1996/01/14  01:56:27  borland
 * Added prototype and flag definitions for SDDS_FilterByNumScan().
 *
 * Revision 1.9  1996/01/10  17:23:57  borland
 * Added prototype for SDDS_LockFile().
 *
 * Revision 1.8  1996/01/09  23:22:58  borland
 * Added prototypes for SDDS_TransferAll*Definitions routines.
 *
 * Revision 1.7  1995/12/12  10:02:28  borland
 * Added prototype for SDDS_SetNoRowCounts().
 *
 * Revision 1.6  1995/12/12  04:30:38  borland
 * Added layout_written field to SDDS_LAYOUT struct; added  file_had_data
 * field to SDDS_DATASET struct; both in support of write/append to
 * files with no_row_counts=1.
 *
 * Revision 1.5  1995/11/13  16:08:13  borland
 * Added first_row_in_mem and writing_page entries to SDDS_DATASET structure
 * to support partially flushed tabular data.
 *
 * Revision 1.4  1995/11/08  04:25:53  borland
 * Added prototypes for new routines FreeArray, DefineParameterLikeArray,
 * and DefineColumnLikeArray.
 *
 * Revision 1.3  1995/09/12  03:18:41  borland
 * Added flag bit #define's for name validity control.
 * Added SDDS_CI_* #defines for SDDS_SetRowsOfInterest case-insensitive
 *   mode.
 * Added SDDS_NOCASE_COMPARE for SDDS_MatchRowsOfInterest case-insenstive
 *   mode.
 *
 * Revision 1.2  1995/09/05  21:14:46  saunders
 * First test release of the SDDS1.5 package.
 *
 */

#if !defined(_SDDS_)
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#if defined(zLib)
#include "zlib.h"
#endif

#ifdef __cplusplus 
extern "C" {
#endif


#define _SDDS_ 1

#define SDDS_VERSION 1

/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: SDDS.h
 * purpose: SDDS data types
 *          format
 *
 * Michael Borland, 1993
 $Log: not supported by cvs2svn $
 Revision 1.6  2002/08/14 15:40:13  soliday
 Added Open License

 Revision 1.5  1999/02/19 22:52:26  borland
 Added SDDS_ANY_INTEGER_TYPE for use with SDDS_CheckXXX procedures.

 Revision 1.4  1997/12/19 16:55:46  borland
 Fixed SDDS_RowCount macro (more parentheses).  Added prototype for
 SDDS_Malloc.  Added new "type": SDDS_ANY_FLOATING_TYPE.

 * Revision 1.3  1995/09/06  14:12:01  saunders
 * First test release of SDDS1.5
 *
 */

#if !defined(_SDDSTYPES_)

#define _SDDSTYPES_ 1

#define SDDS_DOUBLE    1
#define SDDS_FLOAT     2
#define SDDS_LONG      3
#define SDDS_SHORT     4
#define SDDS_STRING    5
#define SDDS_CHARACTER 6
#define SDDS_NUM_TYPES 6
#define SDDS_INTEGER_TYPE(type) ((type)==SDDS_LONG || (type)==SDDS_SHORT)
#define SDDS_FLOATING_TYPE(type) ((type)==SDDS_DOUBLE || (type)==SDDS_FLOAT)
#define SDDS_NUMERIC_TYPE(type) (SDDS_INTEGER_TYPE(type) || SDDS_FLOATING_TYPE(type))
#define SDDS_VALID_TYPE(type) (type>=1 && type<=SDDS_NUM_TYPES)

/* used by SDDS_Check*() routines  */
#define SDDS_ANY_NUMERIC_TYPE (SDDS_NUM_TYPES+1)
#define SDDS_ANY_FLOATING_TYPE   (SDDS_NUM_TYPES+2)
#define SDDS_ANY_INTEGER_TYPE   (SDDS_NUM_TYPES+3)

#endif
#if defined(_WIN32)
typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
#define PRId32 "ld"
#define SCNd32 "ld"
#define INT32_MAX (2147483647)
#else
#include <inttypes.h>
#endif

#undef epicsShareFuncSDDS
#undef epicsShareExtern
#if (defined(_WIN32) && !defined(__CYGWIN32__)) || (defined(__BORLANDC__) && defined(__linux__))
#if defined(EXPORT_SDDS)
#define epicsShareFuncSDDS  __declspec(dllexport)
#define epicsShareExtern extern __declspec(dllexport)
#else
#define epicsShareFuncSDDS
#define epicsShareExtern extern __declspec(dllimport)
#endif
#else
#define epicsShareFuncSDDS
#define epicsShareExtern extern
#endif


#define RW_ASSOCIATES 0

#define SDDS_BINARY 1
#define SDDS_ASCII  2
#define SDDS_NUM_DATA_MODES 2

  /*
extern char *SDDS_type_name[SDDS_NUM_TYPES];
extern int32_t SDDS_type_size[SDDS_NUM_TYPES];
extern char *SDDS_data_mode[SDDS_NUM_DATA_MODES];
  */

epicsShareExtern char *SDDS_type_name[SDDS_NUM_TYPES];
epicsShareExtern int32_t SDDS_type_size[SDDS_NUM_TYPES];
extern char *SDDS_data_mode[SDDS_NUM_DATA_MODES];

/* this shouldn't be changed without changing buffer sizes in namelist routines: */
#define SDDS_MAXLINE 1024

#define SDDS_PRINT_BUFLEN 16834

#define SDDS_NORMAL_DEFINITION 0
#define SDDS_WRITEONLY_DEFINITION 1

typedef struct {
    char *name, *symbol, *units, *description, *format_string;
    int32_t type;
    char *fixed_value;
    /* these are used internally and are not set by the user: */
    int32_t definition_mode, memory_number;
    } PARAMETER_DEFINITION;
#define SDDS_PARAMETER_FIELDS 7

typedef struct {
    char *name, *symbol, *units, *description, *format_string;
    int32_t type, field_length;
    /* these are used internally and are not set by the user: */
    int32_t definition_mode, memory_number, pointer_number;
    } COLUMN_DEFINITION;
#define SDDS_COLUMN_FIELDS 7

typedef struct {
    char *name, *symbol, *units, *description, *format_string, *group_name;
    int32_t type, field_length, dimensions;
    } ARRAY_DEFINITION;
#define SDDS_ARRAY_FIELDS 9

typedef struct {
    char *name, *filename, *path, *description, *contents;
    int32_t sdds;
    } ASSOCIATE_DEFINITION;
#define SDDS_ASSOCIATE_FIELDS 6

typedef struct {
  int32_t index;
  char *name;
} SORTED_INDEX;
int SDDS_CompareIndexedNames(void *s1, void *s2);
int SDDS_CompareIndexedNamesPtr(const void *s1, const void *s2);

typedef struct {
    int32_t mode, lines_per_row, no_row_counts, fixed_row_count, fsync_data;
    int32_t additional_header_lines;
    } DATA_MODE;

#if !defined(zLib)
  typedef void *voidp;
  typedef voidp gzFile;
#endif

typedef struct {
    int32_t n_columns, n_parameters, n_associates, n_arrays;

    char *description;
    char *contents;
    int32_t version;
    short layout_written;

    DATA_MODE data_mode;
    COLUMN_DEFINITION *column_definition;
    PARAMETER_DEFINITION *parameter_definition;
    ARRAY_DEFINITION *array_definition;
    ASSOCIATE_DEFINITION *associate_definition;
    SORTED_INDEX **column_index, **parameter_index, **array_index;

    char *filename;
    FILE *fp;
    gzFile *gzfp;
    short disconnected;
    short gzipFile;
    short popenUsed;
    uint32_t byteOrderDeclared;
    } SDDS_LAYOUT;

typedef struct {
    ARRAY_DEFINITION *definition;
    int32_t *dimension, elements;
    /* Address of array of data values, stored contiguously.
     * For STRING data, the "data values" are actually the addresses of the individual strings.
     */
    void *data;
    void *pointer;
    } SDDS_ARRAY;

typedef struct {
  char *data, *buffer;
  int32_t bytesLeft, bufferSize;
} SDDS_FILEBUFFER ;

#define SDDS_FILEBUFFER_SIZE  262144

typedef struct {
    SDDS_LAYOUT layout, original_layout;
    short swapByteOrder;
    SDDS_FILEBUFFER fBuffer;
    int32_t page_number;

    short mode; /*file mode*/
    short page_started; /*0 or 1 for page not started or already started*/
    int32_t pages_read; /*the number of pages read so far*/
    int32_t endOfFile_offset;/*the offset in the end of the file*/
    int32_t *pagecount_offset; /*the offset of each read page */ 
    int32_t rowcount_offset;  /* ftell() position of row count */
    int32_t n_rows_written;   /* number of tabular data rows written to disk */
    int32_t last_row_written; /* virtual index of last row written */
    int32_t first_row_in_mem; /* virtual index of first row in memory */
    short writing_page;    /* 1/0 for writing/not-writing page */
    int32_t n_rows_allocated; /* number of tabular data rows for which space is alloc'd */
    int32_t n_rows;           /* number of rows stored in memory */
    int32_t *row_flag;        /* accept/reject flags for rows stored in memory */
    short file_had_data;   /* indicates that file being appended to already had some data in it (i.e.,
                            * more than just a header.  Affects no_row_counts=1 output.
                            */
    short autoRecover;
    short autoRecovered;
    
    int32_t n_of_interest;
    int32_t *column_order;          /* column_order[i] = internal index of user's ith column */
    int32_t *column_flag;           /* column_flag[i] indicates whether internal ith column has been selected */

    /* array of SDDS_ARRAY structures for storing array data */
    SDDS_ARRAY *array;

    /* array for parameter data.  The address of the data for the ith parameter
     * is parameter[i].  So *(<type-name> *)parameter[i] gives the data itself.  For type
     * SDDS_STRING the "data itself" is actually the address of the string, and the type-name
     * is "char *".
     */
    void **parameter;

    /* array for tabular data.  The starting address for the data for the ith column
     * is data[i].  So ((<type-name> *)data[i])[j] gives the jth data value.  
     * For type SDDS_STRING the "data value" is actually the address of the string, and
     * the type-name is "char *".
     */
    void **data;
    } SDDS_DATASET;

typedef SDDS_DATASET SDDS_TABLE;

/* prototypes for routines to prepare and write SDDS files */
epicsShareFuncSDDS extern int32_t SDDS_InitializeOutput(SDDS_DATASET *SDDS_dataset, int32_t data_mode,
						     int32_t lines_per_row, char *description,
						     char *contents, char *filename);
epicsShareFuncSDDS extern int32_t SDDS_InitializeAppend(SDDS_DATASET *SDDS_dataset, char *filename);
epicsShareFuncSDDS extern int32_t SDDS_InitializeAppendToPage(SDDS_DATASET *SDDS_dataset, char *filename, 
                                                           int32_t updateInterval,
							   int32_t *rowsPresentReturn);
epicsShareFuncSDDS extern int32_t SDDS_DisconnectFile(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_ReconnectFile(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_FreeStringData(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_Terminate(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern void SDDS_SetTerminateMode(uint32_t mode);
#define TERMINATE_DONT_FREE_TABLE_STRINGS 0x0001
#define TERMINATE_DONT_FREE_ARRAY_STRINGS 0x0002
epicsShareFuncSDDS extern int32_t SDDS_SetRowCountMode(SDDS_DATASET *SDDS_dataset, uint32_t mode);
#define SDDS_VARIABLEROWCOUNT 0x0001UL
#define SDDS_FIXEDROWCOUNT 0x0002UL
#define SDDS_NOROWCOUNT 0x0004UL
epicsShareFuncSDDS extern int32_t SDDS_SetAutoReadRecovery(SDDS_DATASET *SDDS_dataset, uint32_t mode);
#define SDDS_NOAUTOREADRECOVER 0x0001UL
#define SDDS_AUTOREADRECOVER 0x0002UL
epicsShareFuncSDDS extern int32_t SDDS_UpdateRowCount(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern void SDDS_DisableFSync(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern void SDDS_EnableFSync(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_DoFSync(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_DefineParameter(SDDS_DATASET *SDDS_dataset, char *name, char *symbol, char *units, char *description, 
                           char *format_string, int32_t type, char *fixed_value);
epicsShareFuncSDDS extern int32_t SDDS_DefineParameter1(SDDS_DATASET *SDDS_dataset, char *name, char *symbol, char *units, char *description, 
                           char *format_string, int32_t type, void *fixed_value);
epicsShareFuncSDDS extern int32_t SDDS_DefineColumn(SDDS_DATASET *SDDS_dataset, char *name, char *symbol, char *units, char *description, 
                        char *format_string, int32_t type, int32_t field_length);
epicsShareFuncSDDS extern int32_t SDDS_DefineArray(SDDS_DATASET *SDDS_dataset, char *name, char *symbol, char *units, char *description, 
                        char *format_string, int32_t type, int32_t field_length, int32_t dimensions, char *group_name);
epicsShareFuncSDDS extern int32_t SDDS_DefineAssociate(SDDS_DATASET *SDDS_dataset, char *name,
                                 char *filename, char *path, char *description, char *contents, int32_t sdds);
epicsShareFuncSDDS extern int32_t SDDS_IsValidName(char *name, char *dataClass);
epicsShareFuncSDDS extern int32_t SDDS_SetNameValidityFlags(uint32_t flags);
#define SDDS_ALLOW_ANY_NAME 0x0001UL
#define SDDS_ALLOW_V15_NAME 0x0002UL
epicsShareFuncSDDS extern int32_t SDDS_DefineSimpleColumn(SDDS_DATASET *SDDS_dataset, char *name, char *unit, int32_t type);
epicsShareFuncSDDS extern int32_t SDDS_DefineSimpleParameter(SDDS_DATASET *SDDS_dataset, char *name, char *unit, int32_t type);
epicsShareFuncSDDS extern int32_t SDDS_DefineSimpleColumns(SDDS_DATASET *SDDS_dataset, int32_t number, char **name, char **unit, int32_t type);
epicsShareFuncSDDS extern int32_t SDDS_DefineSimpleParameters(SDDS_DATASET *SDDS_dataset, int32_t number, char **name, char **unit, int32_t type);
epicsShareFuncSDDS extern int32_t SDDS_SetNoRowCounts(SDDS_DATASET *SDDS_dataset, int32_t value);
epicsShareFuncSDDS extern int32_t SDDS_WriteLayout(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_EraseData(SDDS_DATASET *SDDS_dataset);

epicsShareFuncSDDS extern int32_t SDDS_ProcessColumnString(SDDS_DATASET *SDDS_dataset, char *string, int32_t mode);
epicsShareFuncSDDS extern int32_t SDDS_ProcessParameterString(SDDS_DATASET *SDDS_dataset, char *string, int32_t mode);
epicsShareFuncSDDS extern int32_t SDDS_ProcessArrayString(SDDS_DATASET *SDDS_dataset, char *string);
epicsShareFuncSDDS extern int32_t SDDS_ProcessAssociateString(SDDS_DATASET *SDDS_dataset, char *string);

epicsShareFuncSDDS extern int32_t SDDS_InitializeCopy(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source, char *filename, char *filemode);
epicsShareFuncSDDS extern int32_t SDDS_CopyLayout(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source);
epicsShareFuncSDDS extern int32_t SDDS_AppendLayout(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source, uint32_t mode);
epicsShareFuncSDDS extern int32_t SDDS_CopyPage(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source);
#define SDDS_CopyTable(a, b) SDDS_CopyPage(a, b)
epicsShareFuncSDDS extern int32_t SDDS_CopyParameters(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source);
epicsShareFuncSDDS extern int32_t SDDS_CopyArrays(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source);
epicsShareFuncSDDS extern int32_t SDDS_CopyColumns(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source);
epicsShareFuncSDDS extern int32_t SDDS_CopyRowsOfInterest(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source);
epicsShareFuncSDDS extern int32_t SDDS_CopyRow(SDDS_DATASET *SDDS_target, int32_t target_row, SDDS_DATASET *SDDS_source, int32_t source_srow);
epicsShareFuncSDDS extern int32_t SDDS_CopyRowDirect(SDDS_DATASET *SDDS_target, int32_t target_row, SDDS_DATASET *SDDS_source, int32_t source_row);
epicsShareFuncSDDS extern int32_t SDDS_CopyAdditionalRows(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source);

epicsShareFuncSDDS extern void SDDS_DeferSavingLayout(int32_t mode);
epicsShareFuncSDDS extern int32_t SDDS_SaveLayout(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_RestoreLayout(SDDS_DATASET *SDDS_dataset);

#define SDDS_BY_INDEX 1
#define SDDS_BY_NAME  2

epicsShareFuncSDDS extern int32_t SDDS_StartPage(SDDS_DATASET *SDDS_dataset, int32_t expected_n_rows);
#define SDDS_StartTable(a, b) SDDS_StartPage(a, b)
epicsShareFuncSDDS extern int32_t SDDS_ClearPage(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_LengthenTable(SDDS_DATASET *SDDS_dataset, int32_t n_additional_rows);
epicsShareFuncSDDS extern int32_t SDDS_ShortenTable(SDDS_DATASET *SDDS_dataset, int32_t rows);
#define SDDS_SET_BY_INDEX SDDS_BY_INDEX
#define SDDS_SET_BY_NAME SDDS_BY_NAME
#define SDDS_PASS_BY_VALUE 4
#define SDDS_PASS_BY_REFERENCE 8
#define SDDS_PASS_BY_STRING 16
epicsShareFuncSDDS extern int32_t SDDS_SetParameters(SDDS_DATASET *SDDS_dataset, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_SetRowValues(SDDS_DATASET *SDDS_dataset, int32_t mode, int32_t row, ...);
epicsShareFuncSDDS extern int32_t SDDS_WritePage(SDDS_DATASET *SDDS_dataset);
#define SDDS_WriteTable(a) SDDS_WritePage(a)
epicsShareFuncSDDS extern int32_t SDDS_UpdatePage(SDDS_DATASET *SDDS_dataset, uint32_t mode);
#define FLUSH_TABLE 0x1UL
#define SDDS_UpdateTable(a) SDDS_UpdatePage(a, 0)
epicsShareFuncSDDS extern int32_t SDDS_SyncDataSet(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_SetColumn(SDDS_DATASET *SDDS_dataset, int32_t mode, void *data, int32_t rows, ...);
epicsShareFuncSDDS extern int32_t SDDS_SetColumnFromDoubles(SDDS_DATASET *SDDS_dataset, int32_t mode, double *data, int32_t rows, ...);
epicsShareFuncSDDS extern int32_t SDDS_SetColumnFromLongs(SDDS_DATASET *SDDS_dataset, int32_t mode, int32_t *data, int32_t rows, ...);
epicsShareFuncSDDS extern int32_t SDDS_SetParametersFromDoubles(SDDS_DATASET *SDDS_dataset, int32_t mode, ...);

#define SDDS_GET_BY_INDEX SDDS_BY_INDEX
#define SDDS_GET_BY_NAME SDDS_BY_NAME
epicsShareFuncSDDS extern int32_t SDDS_GetColumnInformation(SDDS_DATASET *SDDS_dataset, char *field_name, void *memory, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_GetParameterInformation(SDDS_DATASET *SDDS_dataset, char *field_name, void *memory, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_GetArrayInformation(SDDS_DATASET *SDDS_dataset, char *field_name, void *memory, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_GetAssociateInformation(SDDS_DATASET *SDDS_dataset, char *field_name, void *memory, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_ChangeColumnInformation(SDDS_DATASET *SDDS_dataset, char *field_name, void *memory, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_ChangeParameterInformation(SDDS_DATASET *SDDS_dataset, char *field_name, void *memory, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_ChangeArrayInformation(SDDS_DATASET *SDDS_dataset, char *field_name, void *memory, int32_t mode, ...);

epicsShareFuncSDDS extern void SDDS_SetReadRecoveryMode(int32_t mode);
epicsShareFuncSDDS extern int32_t SDDS_SetDefaultIOBufferSize(int32_t bufferSize);

/* prototypes for routines to read and use SDDS files  */
epicsShareFuncSDDS extern int32_t SDDS_InitializeInputFromSearchPath(SDDS_DATASET *SDDSin, char *file);
epicsShareFuncSDDS extern int32_t SDDS_InitializeInput(SDDS_DATASET *SDDS_dataset, char *filename);
epicsShareFuncSDDS extern int32_t SDDS_InitializeHeaderlessInput(SDDS_DATASET *SDDS_dataset, char *filename);
epicsShareFuncSDDS extern int32_t SDDS_GetRowLimit();
epicsShareFuncSDDS extern int32_t SDDS_SetRowLimit(int32_t limit);
epicsShareFuncSDDS extern int32_t SDDS_GotoPage(SDDS_DATASET *SDDS_dataset,int32_t page_number);
epicsShareFuncSDDS extern int32_t SDDS_ReadPage(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_ReadPageSparse(SDDS_DATASET *SDDS_dataset, uint32_t mode,
                                int32_t sparse_interval,
                                int32_t sparse_offset);
#define SDDS_ReadTable(a) SDDS_ReadPage(a)
epicsShareFuncSDDS extern int32_t SDDS_ReadAsciiPage(SDDS_DATASET *SDDS_dataset, int32_t sparse_interval,
                               int32_t sparse_offset);
epicsShareFuncSDDS extern int32_t SDDS_ReadRecoveryPossible(void);

epicsShareFuncSDDS extern int32_t SDDS_SetColumnFlags(SDDS_DATASET *SDDS_dataset, int32_t column_flag_value);
epicsShareFuncSDDS extern int32_t SDDS_SetRowFlags(SDDS_DATASET *SDDS_dataset, int32_t row_flag_value);
epicsShareFuncSDDS extern int32_t SDDS_GetRowFlag(SDDS_DATASET *SDDS_dataset, int32_t row);
epicsShareFuncSDDS extern int32_t SDDS_GetRowFlags(SDDS_DATASET *SDDS_dataset, int32_t *flag, int32_t rows);
epicsShareFuncSDDS extern int32_t SDDS_BufferedRead(void *target, size_t targetSize, FILE *fp, SDDS_FILEBUFFER *fBuffer);
#if defined(zLib)
epicsShareFuncSDDS extern int32_t SDDS_GZipBufferedRead(void *target, size_t targetSize, gzFile *gzfp, SDDS_FILEBUFFER *fBuffer);
#endif

/* logic flags for SDDS_AssertRowFlags and SDDS_AssertColumnFlags */
#define SDDS_FLAG_ARRAY  0x001UL
#define SDDS_INDEX_LIMITS 0x002UL
epicsShareFuncSDDS extern int32_t SDDS_AssertRowFlags(SDDS_DATASET *SDDS_dataset, uint32_t mode, ...);
/* modes for SDDS_SetColumnsOfInterest and SDDS_SetRowsOfInterest: */
#define SDDS_NAME_ARRAY 1
#define SDDS_NAMES_STRING 2
#define SDDS_NAME_STRINGS 3
#define SDDS_MATCH_STRING 4
#define SDDS_MATCH_EXCLUDE_STRING 5
#define SDDS_CI_NAME_ARRAY 6
#define SDDS_CI_NAMES_STRING 7
#define SDDS_CI_NAME_STRINGS 8
#define SDDS_CI_MATCH_STRING 9
#define SDDS_CI_MATCH_EXCLUDE_STRING 10

/* logic flags for SDDS_SetColumnsOfInterest, SDDS_SetRowsOfInterest, and SDDS_MatchRowsOfInterest: */
#define SDDS_AND               0x0001UL
#define SDDS_OR                0x0002UL
#define SDDS_NEGATE_MATCH      0x0004UL
#define SDDS_NEGATE_PREVIOUS   0x0008UL
#define SDDS_NEGATE_EXPRESSION 0x0010UL
#define SDDS_INDIRECT_MATCH    0x0020UL
#define SDDS_1_PREVIOUS        0x0040UL
#define SDDS_0_PREVIOUS        0x0080UL
/* used by MatchRowsOfInterest only at this point: */
#define SDDS_NOCASE_COMPARE    0x0100UL

epicsShareFuncSDDS extern int32_t SDDS_MatchColumns(SDDS_DATASET *SDDS_dataset, char ***match, int32_t matchMode, int32_t typeMode, ... );
epicsShareFuncSDDS extern int32_t SDDS_MatchParameters(SDDS_DATASET *SDDS_dataset, char ***match, int32_t matchMode, int32_t typeMode, ... );
epicsShareFuncSDDS extern int32_t SDDS_MatchArrays(SDDS_DATASET *SDDS_dataset, char ***match, int32_t matchMode, int32_t typeMode, ... );
epicsShareFuncSDDS extern int32_t SDDS_Logic(int32_t previous, int32_t match, uint32_t logic);
epicsShareFuncSDDS extern int32_t SDDS_SetColumnsOfInterest(SDDS_DATASET *SDDS_dataset, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_AssertColumnFlags(SDDS_DATASET *SDDS_dataset, uint32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_SetRowsOfInterest(SDDS_DATASET *SDDS_dataset, char *selection_column, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_MatchRowsOfInterest(SDDS_DATASET *SDDS_dataset, char *selection_column, char *label_to_match, int32_t logic);
epicsShareFuncSDDS extern int32_t SDDS_DeleteColumn(SDDS_DATASET *SDDS_dataset, char *column_name);
epicsShareFuncSDDS extern int32_t SDDS_DeleteParameter(SDDS_DATASET *SDDS_dataset, char *parameter_name);
epicsShareFuncSDDS extern int32_t SDDS_DeleteUnsetColumns(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_CountColumnsOfInterest(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_ColumnIsOfInterest(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern int32_t SDDS_ColumnCount(SDDS_DATASET *dataset);
epicsShareFuncSDDS extern int32_t SDDS_ParameterCount(SDDS_DATASET *dataset);
epicsShareFuncSDDS extern int32_t SDDS_ArrayCount(SDDS_DATASET *dataset);
epicsShareFuncSDDS extern int32_t SDDS_CountRowsOfInterest(SDDS_DATASET *SDDS_dataset);
#define SDDS_RowCount(SDDS_dataset) ((SDDS_dataset)->n_rows)
epicsShareFuncSDDS extern int32_t SDDS_DeleteUnsetRows(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_FilterRowsOfInterest(SDDS_DATASET *SDDS_dataset, char *filter_column, double lower, double upper, int32_t logic);
epicsShareFuncSDDS extern int32_t SDDS_ItemInsideWindow(void *data, int32_t index, int32_t type, double lower_limit, double upper_limit);
epicsShareFuncSDDS extern int32_t SDDS_FilterRowsByNumScan(SDDS_DATASET *SDDS_dataset, char *filter_column, uint32_t mode);
#define NUMSCANFILTER_INVERT 0x0001UL
#define NUMSCANFILTER_STRICT 0x0002UL

epicsShareFuncSDDS extern void *SDDS_GetColumn(SDDS_DATASET *SDDS_dataset, char *column_name);
epicsShareFuncSDDS extern void *SDDS_GetInternalColumn(SDDS_DATASET *SDDS_dataset, char *column_name);
epicsShareFuncSDDS extern double *SDDS_GetColumnInDoubles(SDDS_DATASET *SDDS_dataset, char *column_name);
epicsShareFuncSDDS extern int32_t *SDDS_GetColumnInLong(SDDS_DATASET *SDDS_dataset, char *column_name);
epicsShareFuncSDDS extern void *SDDS_GetNumericColumn(SDDS_DATASET *SDDS_dataset, char *column_name, int32_t desiredType);
epicsShareFuncSDDS extern void *SDDS_GetRow(SDDS_DATASET *SDDS_dataset, int32_t srow_index, void *memory);
epicsShareFuncSDDS extern void *SDDS_GetValue(SDDS_DATASET *SDDS_dataset, char *column_name, int32_t srow_index, void *memory);
epicsShareFuncSDDS extern void *SDDS_GetValueByIndex(SDDS_DATASET *SDDS_dataset, int32_t column_index, int32_t srow_index, void *memory);
epicsShareFuncSDDS extern void *SDDS_GetValueByAbsIndex(SDDS_DATASET *SDDS_dataset, int32_t column_index, int32_t srow_index, void *memory);
epicsShareFuncSDDS extern void *SDDS_GetParameter(SDDS_DATASET *SDDS_dataset, char *parameter_name, void *memory);
epicsShareFuncSDDS extern void *SDDS_GetParameterByIndex(SDDS_DATASET *SDDS_dataset, int32_t index, void *memory);
epicsShareFuncSDDS extern double *SDDS_GetParameterAsDouble(SDDS_DATASET *SDDS_dataset, char *parameter_name, double *data);
epicsShareFuncSDDS extern int32_t *SDDS_GetParameterAsLong(SDDS_DATASET *SDDS_dataset, char *parameter_name, int32_t *data);
epicsShareFuncSDDS extern char *SDDS_GetParameterAsString(SDDS_DATASET *SDDS_dataset, char *parameter_name, char **memory);
epicsShareFuncSDDS extern int32_t SDDS_GetParameters(SDDS_DATASET *SDDS_dataset, ...);
epicsShareFuncSDDS extern void *SDDS_GetFixedValueParameter(SDDS_DATASET *SDDS_dataset, char *parameter_name, void *memory);
epicsShareFuncSDDS extern int32_t SDDS_GetDescription(SDDS_DATASET *SDDS_dataset, char **text, char **contents);

epicsShareFuncSDDS extern void *SDDS_GetMatrixOfRows(SDDS_DATASET *SDDS_dataset, int32_t *n_rows);
epicsShareFuncSDDS extern void *SDDS_GetCastMatrixOfRows(SDDS_DATASET *SDDS_dataset, int32_t *n_rows, int32_t sddsType);
#define SDDS_ROW_MAJOR_DATA 1
#define SDDS_COLUMN_MAJOR_DATA 2
epicsShareFuncSDDS extern void *SDDS_GetMatrixFromColumn(SDDS_DATASET *SDDS_dataset, char *column_name, int32_t dimension1, int32_t dimension2, int32_t mode);
epicsShareFuncSDDS extern void *SDDS_GetDoubleMatrixFromColumn(SDDS_DATASET *SDDS_dataset, char *column_name, int32_t dimension1, int32_t dimension2, int32_t mode);

epicsShareFuncSDDS extern SDDS_ARRAY *SDDS_GetArray(SDDS_DATASET *SDDS_dataset, char *array_name, SDDS_ARRAY *memory);
#define SDDS_POINTER_ARRAY 1
#define SDDS_CONTIGUOUS_DATA 2
epicsShareFuncSDDS extern double *SDDS_GetArrayInDoubles(SDDS_DATASET *SDDS_dataset, char *array_name, int32_t *values);
epicsShareFuncSDDS extern int32_t SDDS_SetArrayVararg(SDDS_DATASET *SDDS_dataset, char *array_name, int32_t mode, void *data_pointer, ...);
epicsShareFuncSDDS extern int32_t SDDS_SetArray(SDDS_DATASET *SDDS_dataset, char *array_name, int32_t mode, void *data_pointer, int32_t *dimension);
epicsShareFuncSDDS extern int32_t SDDS_AppendToArrayVararg(SDDS_DATASET *SDDS_dataset, char *array_name, int32_t mode, void *data_pointer, int32_t elements, ...);


/* error-handling and utility routines: */
epicsShareFuncSDDS extern void *SDDS_Realloc(void *old_ptr, size_t new_size);
epicsShareFuncSDDS extern void *SDDS_Malloc(size_t size);
epicsShareFuncSDDS extern void SDDS_Free(void *mem);
epicsShareFuncSDDS extern void *SDDS_Calloc(size_t nelem, size_t elem_size);
epicsShareFuncSDDS extern int32_t SDDS_NumberOfErrors(void);
epicsShareFuncSDDS extern void SDDS_ClearErrors(void);
epicsShareFuncSDDS extern void SDDS_SetError(char *error_text);
epicsShareFuncSDDS extern void SDDS_Bomb(char *message);
epicsShareFuncSDDS extern void SDDS_Warning(char *message);
epicsShareFuncSDDS extern void SDDS_RegisterProgramName(const char *name);
#define SDDS_VERBOSE_PrintErrors 1
#define SDDS_EXIT_PrintErrors 2
epicsShareFuncSDDS extern void SDDS_PrintErrors(FILE *fp, int32_t mode);
#define SDDS_LAST_GetErrorMessages 0
#define SDDS_ALL_GetErrorMessages 1
epicsShareFuncSDDS extern char **SDDS_GetErrorMessages(int32_t *number, int32_t mode);

epicsShareFuncSDDS extern char **SDDS_GetColumnNames(SDDS_DATASET *SDDS_dataset, int32_t *number);
epicsShareFuncSDDS extern char **SDDS_GetParameterNames(SDDS_DATASET *SDDS_dataset, int32_t *number);
epicsShareFuncSDDS extern char **SDDS_GetAssociateNames(SDDS_DATASET *SDDS_dataset, int32_t *number);
epicsShareFuncSDDS extern char **SDDS_GetArrayNames(SDDS_DATASET *SDDS_dataset, int32_t *number);

epicsShareFuncSDDS extern COLUMN_DEFINITION *SDDS_GetColumnDefinition(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern COLUMN_DEFINITION *SDDS_CopyColumnDefinition(COLUMN_DEFINITION **target, COLUMN_DEFINITION *source);
epicsShareFuncSDDS extern int32_t SDDS_FreeColumnDefinition(COLUMN_DEFINITION *source);
epicsShareFuncSDDS extern int32_t SDDS_TransferColumnDefinition(SDDS_DATASET *target, SDDS_DATASET *source, char *name, char *newName);
epicsShareFuncSDDS extern int32_t SDDS_DefineColumnLikeParameter(SDDS_DATASET *target, SDDS_DATASET *source, char *name, char *newName);
epicsShareFuncSDDS extern int32_t SDDS_DefineColumnLikeArray(SDDS_DATASET *target, SDDS_DATASET *source, char *name, char *newName);
epicsShareFuncSDDS extern int32_t SDDS_TransferAllColumnDefinitions(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source,
                                              uint32_t mode);


epicsShareFuncSDDS extern PARAMETER_DEFINITION *SDDS_GetParameterDefinition(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern PARAMETER_DEFINITION *SDDS_CopyParameterDefinition(PARAMETER_DEFINITION **target, PARAMETER_DEFINITION *source);
epicsShareFuncSDDS extern int32_t SDDS_FreeParameterDefinition(PARAMETER_DEFINITION *source);
epicsShareFuncSDDS extern int32_t SDDS_TransferParameterDefinition(SDDS_DATASET *target, SDDS_DATASET *source, char *name, char *newName);
epicsShareFuncSDDS extern int32_t SDDS_DefineParameterLikeColumn(SDDS_DATASET *target, SDDS_DATASET *source, char *name, char *newName);
epicsShareFuncSDDS extern int32_t SDDS_DefineParameterLikeArray(SDDS_DATASET *target, SDDS_DATASET *source, char *name, char *newName);
#define SDDS_TRANSFER_KEEPOLD 0x01UL
#define SDDS_TRANSFER_OVERWRITE 0x02UL
epicsShareFuncSDDS extern int32_t SDDS_TransferAllParameterDefinitions(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source, 
                                                 uint32_t mode);

epicsShareFuncSDDS extern ARRAY_DEFINITION *SDDS_GetArrayDefinition(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern ARRAY_DEFINITION *SDDS_CopyArrayDefinition(ARRAY_DEFINITION **target, ARRAY_DEFINITION *source);
epicsShareFuncSDDS extern int32_t SDDS_FreeArrayDefinition(ARRAY_DEFINITION *source);
epicsShareFuncSDDS extern int32_t SDDS_TransferArrayDefinition(SDDS_DATASET *target, SDDS_DATASET *source, char *name, char *newName);
epicsShareFuncSDDS extern int32_t SDDS_TransferAllArrayDefinitions(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source,
                                             uint32_t mode);


epicsShareFuncSDDS extern ASSOCIATE_DEFINITION *SDDS_GetAssociateDefinition(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern ASSOCIATE_DEFINITION *SDDS_CopyAssociateDefinition(ASSOCIATE_DEFINITION **target, ASSOCIATE_DEFINITION *source);
epicsShareFuncSDDS extern int32_t SDDS_FreeAssociateDefinition(ASSOCIATE_DEFINITION *source);
epicsShareFuncSDDS extern int32_t SDDS_TransferAssociateDefinition(SDDS_DATASET *target, SDDS_DATASET *source, char *name, char *newName);

epicsShareFuncSDDS extern int32_t SDDS_GetColumnIndex(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern int32_t SDDS_GetParameterIndex(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern int32_t SDDS_GetArrayIndex(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern int32_t SDDS_GetAssociateIndex(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern int32_t SDDS_GetColumnType(SDDS_DATASET *SDDS_dataset, int32_t index);
epicsShareFuncSDDS extern int32_t SDDS_GetNamedColumnType(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern int32_t SDDS_GetParameterType(SDDS_DATASET *SDDS_dataset, int32_t index);
epicsShareFuncSDDS extern int32_t SDDS_GetNamedParameterType(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern int32_t SDDS_GetArrayType(SDDS_DATASET *SDDS_dataset, int32_t index);
epicsShareFuncSDDS extern int32_t SDDS_GetNamedArrayType(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern int32_t SDDS_GetTypeSize(int32_t type);
epicsShareFuncSDDS extern char *SDDS_GetTypeName(int32_t type);
epicsShareFuncSDDS extern int32_t SDDS_IdentifyType(char *typeName);

#define FIND_ANY_TYPE       0
#define FIND_SPECIFIED_TYPE 1
#define FIND_NUMERIC_TYPE   2
#define FIND_INTEGER_TYPE   3
#define FIND_FLOATING_TYPE  4

epicsShareFuncSDDS extern char *SDDS_FindColumn(SDDS_DATASET *SDDS_dataset,  int32_t mode, ...);
epicsShareFuncSDDS extern char *SDDS_FindParameter(SDDS_DATASET *SDDS_dataset,  int32_t mode, ...);
epicsShareFuncSDDS extern char *SDDS_FindArray(SDDS_DATASET *SDDS_dataset,  int32_t mode, ...);

epicsShareFuncSDDS extern int32_t SDDS_CheckColumn(SDDS_DATASET *SDDS_dataset, char *name, char *units, int32_t type, FILE *fp_message);
epicsShareFuncSDDS extern int32_t SDDS_CheckParameter(SDDS_DATASET *SDDS_dataset, char *name, char *units, int32_t type, FILE *fp_message);
epicsShareFuncSDDS extern int32_t SDDS_CheckArray(SDDS_DATASET *SDDS_dataset, char *name, char *units, int32_t type, FILE *fp_message);
epicsShareFuncSDDS extern int32_t SDDS_VerifyArrayExists(SDDS_DATASET *SDDS_dataset, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_VerifyColumnExists(SDDS_DATASET *SDDS_dataset, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_VerifyParameterExists(SDDS_DATASET *SDDS_dataset, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_PrintCheckText(FILE *fp, char *name, char *units, int32_t type, char *class_name, int32_t error_code);
#define SDDS_CHECK_OKAY 0
#define SDDS_CHECK_OK SDDS_CHECK_OKAY
#define SDDS_CHECK_NONEXISTENT 1
#define SDDS_CHECK_WRONGTYPE  2
#define SDDS_CHECK_WRONGUNITS  3

epicsShareFuncSDDS extern int32_t SDDS_IsActive(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_ForceInactive(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_LockFile(FILE *fp, char *filename, char *callerName);
epicsShareFuncSDDS extern int32_t SDDS_FileIsLocked(char *filename);
epicsShareFuncSDDS extern int32_t SDDS_BreakIntoLockedFile(char *filename);

epicsShareFuncSDDS extern int32_t SDDS_CopyString(char **target, char *source);
epicsShareFuncSDDS extern int32_t SDDS_CopyStringArray(char **target, char **source, int32_t n_strings);
epicsShareFuncSDDS extern int32_t SDDS_FreeStringArray(char **string, int32_t strings);
epicsShareFuncSDDS extern int32_t SDDS_VerifyPrintfFormat(char *format_string, int32_t type);
epicsShareFuncSDDS extern int32_t SDDS_HasWhitespace(char *string);
epicsShareFuncSDDS extern char *fgetsSkipComments(char *s, int32_t slen, FILE *fp, char skip_char);
epicsShareFuncSDDS extern char *fgetsSkipCommentsResize(char **s, int32_t *slen, FILE *fp, char skip_char);
#if defined(zLib)
epicsShareFuncSDDS extern char *fgetsGZipSkipComments(char *s, int32_t slen, gzFile *gzfp, char skip_char);
epicsShareFuncSDDS extern char *fgetsGZipSkipCommentsResize(char **s, int32_t *slen, gzFile *gzfp, char skip_char);
#endif
epicsShareFuncSDDS extern void SDDS_CutOutComments(char *s, char cc);
epicsShareFuncSDDS extern void SDDS_EscapeNewlines(char *s);
epicsShareFuncSDDS extern void SDDS_EscapeQuotes(char *s, char quote_char);
epicsShareFuncSDDS extern void SDDS_UnescapeQuotes(char *s, char quote_char);
epicsShareFuncSDDS extern int32_t SDDS_IsQuoted(char *string, char *position, char quotation_mark);
epicsShareFuncSDDS extern int32_t SDDS_GetToken(char *s, char *buffer, int32_t buflen);
epicsShareFuncSDDS extern int32_t SDDS_PadToLength(char *string, int32_t length);
epicsShareFuncSDDS extern void SDDS_EscapeCommentCharacters(char *string, char cc);
epicsShareFuncSDDS extern void SDDS_InterpretEscapes(char *s);

epicsShareFuncSDDS extern int32_t SDDS_ZeroMemory(void *mem, int32_t n_bytes);
epicsShareFuncSDDS extern int32_t SDDS_SetMemory(void *mem, int32_t n_elements, int32_t data_type, ...);
#define SDDS_PRINT_NOQUOTES 0x0001UL
epicsShareFuncSDDS extern int32_t SDDS_SprintTypedValue(void *data, int32_t index, int32_t type, char *format, char *buffer, uint32_t mode);
epicsShareFuncSDDS extern int32_t SDDS_PrintTypedValue(void *data, int32_t index, int32_t type, char *format, FILE *fp, uint32_t mode);
epicsShareFuncSDDS extern int32_t SDDS_WriteTypedValue(void *data, int32_t index, int32_t type, char *format, FILE *fp);
epicsShareFuncSDDS extern void *SDDS_CastValue(void *data, int32_t index, int32_t data_type, int32_t desired_type, void *memory);
epicsShareFuncSDDS extern void SDDS_RemovePadding(char *s);
epicsShareFuncSDDS extern int32_t SDDS_StringIsBlank(char *s);
epicsShareFuncSDDS extern void *SDDS_AllocateMatrix(int32_t size, int32_t dim1, int32_t dim2);
epicsShareFuncSDDS extern void SDDS_FreeMatrix(void **ptr, int32_t dim1);
epicsShareFuncSDDS extern void SDDS_FreeArray(SDDS_ARRAY *array);
epicsShareFuncSDDS extern void *SDDS_MakePointerArray(void *data, int32_t type, int32_t dimensions, int32_t *dimension);
epicsShareFuncSDDS extern int32_t SDDS_ApplyFactorToParameter(SDDS_DATASET *SDDS_dataset, char *name, double factor);
epicsShareFuncSDDS extern int32_t SDDS_ApplyFactorToColumn(SDDS_DATASET *SDDS_dataset, char *name, double factor);
epicsShareFuncSDDS extern int32_t SDDS_DeleteParameterFixedValues(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_SetDataMode(SDDS_DATASET *SDDS_dataset, int32_t newmode);
epicsShareFuncSDDS extern int32_t SDDS_CheckDataset(SDDS_DATASET *SDDS_dataset, const char *caller);
epicsShareFuncSDDS extern int32_t SDDS_CheckTabularData(SDDS_DATASET *SDDS_dataset, const char *caller);
epicsShareFuncSDDS extern int32_t SDDS_CheckDatasetStructureSize(int32_t size);
#define SDDS_CheckTableStructureSize(a) SDDS_CheckDatasetStructureSize(a)

#define TABULAR_DATA_CHECKS 0x0001UL
epicsShareFuncSDDS uint32_t SDDS_SetAutoCheckMode(uint32_t newMode);

epicsShareFuncSDDS extern int32_t SDDS_FlushBuffer(FILE *fp, SDDS_FILEBUFFER *fBuffer);
epicsShareFuncSDDS extern int32_t SDDS_BufferedWrite(void *target, size_t targetSize, FILE *fp, SDDS_FILEBUFFER *fBuffer);

epicsShareFuncSDDS extern int32_t SDDS_ScanData(char *string, int32_t type, int32_t field_length, void *data, int32_t index, int32_t is_parameter);

epicsShareFuncSDDS extern double SDDS_ConvertToDouble(int32_t type, void *data, int32_t index);
epicsShareFuncSDDS extern int32_t SDDS_ConvertToLong(int32_t type, void *data, int32_t index);

epicsShareFuncSDDS extern int32_t SDDS_WriteBinaryString(char *string, FILE *fp, SDDS_FILEBUFFER *fBuffer);
#if defined(zLib)
epicsShareFuncSDDS extern int32_t SDDS_GZipWriteBinaryString(char *string, gzFile *gzfp, SDDS_FILEBUFFER *fBuffer);
epicsShareFuncSDDS extern int32_t SDDS_GZipBufferedRead(void *target, size_t targetSize, gzFile *gzfp, SDDS_FILEBUFFER *fBuffer);
epicsShareFuncSDDS extern int32_t SDDS_GZipFlushBuffer(gzFile *gzfp, SDDS_FILEBUFFER *fBuffer);
epicsShareFuncSDDS extern int32_t SDDS_GZipBufferedWrite(void *target, size_t targetSize, gzFile *gzfp, SDDS_FILEBUFFER *fBuffer);
#endif


epicsShareFuncSDDS extern int32_t SDDS_CreateRpnMemory(char *name, short is_string);
epicsShareFuncSDDS extern int32_t SDDS_CreateRpnArray(char *name);

#if defined(RPN_SUPPORT)
epicsShareFuncSDDS extern int32_t SDDS_FilterRowsWithRpnTest(SDDS_DATASET *SDDS_dataset, char *rpn_test);
epicsShareFuncSDDS extern int32_t SDDS_StoreParametersInRpnMemories(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_StoreRowInRpnMemories(SDDS_DATASET *SDDS_dataset, int32_t row);
epicsShareFuncSDDS extern int32_t SDDS_StoreColumnsInRpnArrays(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_ComputeColumn(SDDS_DATASET *SDDS_dataset, int32_t column, char *equation);
epicsShareFuncSDDS extern int32_t SDDS_ComputeParameter(SDDS_DATASET *SDDS_dataset, int32_t column, char *equation);
#endif

#define SDDS_BIGENDIAN_SEEN      0x0001UL
#define SDDS_LITTLEENDIAN_SEEN   0x0002UL
#define SDDS_FIXED_ROWCOUNT_SEEN 0x0004UL
#define SDDS_BIGENDIAN         SDDS_BIGENDIAN_SEEN
#define SDDS_LITTLEENDIAN      SDDS_LITTLEENDIAN_SEEN
#define SDDS_FIXED_ROWCOUNT    SDDS_FIXED_ROWCOUNT_SEEN
epicsShareFuncSDDS extern int32_t SDDS_IsBigEndianMachine();
void SDDS_SwapShort(short *data);
epicsShareFuncSDDS extern void SDDS_SwapLong(int32_t *data);
void SDDS_SwapFloat(float *data);
void SDDS_SwapDouble(double *data);
epicsShareFuncSDDS extern int32_t SDDS_SwapEndsArrayData(SDDS_DATASET *SDDSin);
epicsShareFuncSDDS extern int32_t SDDS_SwapEndsParameterData(SDDS_DATASET *SDDSin) ;
epicsShareFuncSDDS extern int32_t SDDS_SwapEndsColumnData(SDDS_DATASET *SDDSin);






epicsShareFuncSDDS extern int32_t SDDS_ReadNonNativePage(SDDS_DATASET *SDDS_dataset);
int32_t SDDS_ReadNonNativePageSparse(SDDS_DATASET *SDDS_dataset, uint32_t mode,
                         int32_t sparse_interval,
                         int32_t sparse_offset);
int32_t SDDS_ReadNonNativeBinaryPage(SDDS_DATASET *SDDS_dataset, int32_t sparse_interval, int32_t sparse_offset);
int32_t SDDS_ReadNonNativeBinaryParameters(SDDS_DATASET *SDDS_dataset);
int32_t SDDS_ReadNonNativeBinaryArrays(SDDS_DATASET *SDDS_dataset);
int32_t SDDS_ReadNonNativeBinaryRow(SDDS_DATASET *SDDS_dataset, int32_t row, int32_t skip);
char *SDDS_ReadNonNativeBinaryString(FILE *fp, SDDS_FILEBUFFER *fBuffer, int32_t skip);
#if defined(zLib)
char *SDDS_ReadNonNativeGZipBinaryString(gzFile *gzfp, SDDS_FILEBUFFER *fBuffer, int32_t skip);
#endif



epicsShareFuncSDDS extern int32_t SDDS_WriteNonNativeBinaryPage(SDDS_DATASET *SDDS_dataset);
int32_t SDDS_WriteNonNativeBinaryParameters(SDDS_DATASET *SDDS_dataset);
int32_t SDDS_WriteNonNativeBinaryArrays(SDDS_DATASET *SDDS_dataset);
int32_t SDDS_WriteNonNativeBinaryRow(SDDS_DATASET *SDDS_dataset, int32_t row);

int32_t SDDS_WriteNonNativeBinaryString(char *string, FILE *fp, SDDS_FILEBUFFER *fBuffer);
#if defined(zLib)
int32_t SDDS_GZipWriteNonNativeBinaryString(char *string, gzFile *gzfp, SDDS_FILEBUFFER *fBuffer);
#endif


#ifdef __cplusplus
}
#endif

#endif
