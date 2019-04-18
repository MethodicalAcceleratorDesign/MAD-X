#ifndef MAD_ERR_H
#define MAD_ERR_H

// interface

// used in C and Fortran
void  seterrorflagfort(int* errcode, const char* from, int *lf, const char* descr, int *ld);
void  seterrorflag    (int  errcode, const char* from,          const char* descr);
int   geterrorflag    (void);
void  clearerrorflag  (void);
void  augmentfwarn    (void);

void  put_info    (const char* t1, const char* t2);
void  warning     (const char* t1, const char* t2);
void  fatal_error (const char* t1, const char* t2);

// used in C only
void  warningnew  (const char* t1, const char* fmt, ...);
void  mad_error   (const char* t1, const char* fmt, ...);

// used only in madx_finish
void  mad_err_getwarn(int* cwarn, int* fwarn);

#endif // MAD_ERR_H
