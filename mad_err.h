#ifndef MAD_ERR_H
#define MAD_ERR_H

// interface

// used in C and Fortran
void  seterrorflagfort(int* errcode, const char* from, int *lf, const char* descr, int *ld);
void  seterrorflag(int errcode, const char* from, const char* descr);
int   geterrorflag(void);
void  clearerrorflag(void);
void  augmentfwarn(void);

void  put_info(char* t1, char* t2);
void  warning(char* t1, char* t2);
void  fatal_error(char* t1, char* t2);

// used in C only
void  warningnew(char* t1, char* fmt, ...);
void  error(char* t1, char* fmt, ...);

// used only in madx_finish
void  mad_err_getwarn(int* cwarn, int* fwarn);

#endif // MAD_ERR_H

