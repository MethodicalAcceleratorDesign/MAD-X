#ifndef MAD_TRACK_H
#define MAD_TRACK_H

// types

struct in_cmd;
struct command;

// interface

void  pro_track(struct in_cmd*);
void  track_pteigen(double* eigen);
void  track_tables_dump(void);
void  track_tables_create(struct in_cmd*);
void  track_tables_delete(void);
void  track_start(struct command*);

/**
 * Used by copytrackstoarray
 * Used in madx_ptc_trackcavs.f90
 */
void  deletetrackstrarpositions(void);
int   getcurrentcmdname(char* string);
int   getnumberoftracks(void);
/**
 * Used in madx_ptc_trackcavs.f90
 */
int   gettrack(int* nt, double* x, double* px, double* y, double* py, double* t, double* pt);

int   next_start(double* x, double* px, double* y, double* py, double* t, double* deltae,
                 double* fx,double* phix, double* fy, double* phiy, double* ft,double* phit);

#endif // MAD_TRACK_H

