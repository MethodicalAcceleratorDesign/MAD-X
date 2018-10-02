#include "madx.h"

// private globals

static int warn_numb  = 0;           /* Number of warnings */
static int warn_numbf = 0;           /* Number of warnings from fortran */
static int errorflag = 0;

// public interface

void
mad_err_getwarn(int* cwarn, int* fwarn)
{
  if (cwarn) *cwarn = warn_numb;
  if (fwarn) *fwarn = warn_numbf;
}

void
seterrorflagfort(int* errcode, const char* from, int *lf, const char* descr, int *ld)
{
/*Sets error flag - Used to comunicate occurance of an error.
  Mainly between c and fortran parts
  This one is called from Fortran
*/
  static const int maxlength = 400;
  char f[maxlength];
  char d[maxlength];
  int n = *lf , m = *ld;
  if (n >= maxlength ) n = maxlength - 1;
  if (m >= maxlength ) m = maxlength - 1;
  strncpy(f,from,n);
  strncpy(d,descr,m);
  f[n]=0;
  d[m]=0;
  seterrorflag(*errcode,f,d);
}

void
seterrorflag(int errcode, const char* from, const char* descr)
{
/*Sets error flag - Used to comunicate occurance of an error.
  Mainly between c and fortran parts*/
  errorflag = errcode;
  mad_error("seterrorflag","Errorcode: %d   Reported from %s:",errorflag,from);
  mad_error("seterrorflag","Description: %s",descr);
}

int
geterrorflag(void)
{
  /*returns errorflag status - used by fortran code to check if error occured*/
  return errorflag;
}

void
clearerrorflag(void)
{
  errorflag = 0;
}

#if 0 // not used...
char*
geterrrormessage()
{
  return errormessage;
}
#endif

void
warningnew(const char* t1, const char* fmt, ...)
{
/*prints warning on the standard error and accepts parameters printout with std C formatting*/
/*Piotr Skowronski CERN*/
  va_list         list;

  warn_numb++; /*I think that warnings should be counted even if the user does not want to see them*/
  fflush(0); /*flushes all the buffers -> so the warning appears in a correct place*/

  if (get_option("warn") == 0)
  {
    return;
  }

  va_start( list, fmt );

  fprintf(stdout,"++++++ warning: %s : ",t1); /*prints first part to the STDERR and +++....*/
  vfprintf(stdout, fmt, list); /*prints the second part and variables*/
  fprintf(stdout,"\n"); /*prints end of line*/
  fflush(stdout); /*flushes STDERR*/
  va_end(list);
}

void
mad_error(const char* t1, const char* fmt, ...)
{
/*prints warning on the standard error and accepts parameters printout with std C formatting*/
/*Piotr Skowronski CERN*/
  va_list         list;

  warn_numb++; /*I think that warnings should be counted even if the user does not want to see them*/
  fflush(0); /*flushes all the buffers -> so the warning appears in a correct place*/

  va_start( list, fmt );

  fprintf(stderr,"++++++ Error: %s : ",t1); /*prints first part to the STDERR and +++....*/
  vfprintf(stderr, fmt, list); /*prints the second part and variables*/
  fprintf(stderr,"\n"); /*prints end of line*/
  fflush(stderr); /*flushes STDERR*/
  va_end(list);
}

void
fatal_error(const char* t1, const char* t2)
  /*prints fatal error message, halts program */
{
  fprintf(stderr, "+=+=+= fatal: %s %s\n",t1,t2);
  if (get_option("no_fatal_stop ")==0) exit(1);
}

void
warning(const char* t1, const char* t2)
{
  if (get_option("warn")) {
    printf("++++++ warning: %s %s\n",t1,t2);
    warn_numb++;
  }
}

void
augmentfwarn(void)
{
/*increases counter of the fortran warnings*/
  warn_numbf++;
}

void
put_info(const char* t1, const char* t2)
{
  if (get_option("info") && get_option("warn"))
    printf("++++++ info: %s %s\n",t1,t2);
}

