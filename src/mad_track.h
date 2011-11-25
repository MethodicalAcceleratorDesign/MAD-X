#ifndef MAD_TRACK_H
#define MAD_TRACK_H

// types

struct in_cmd;

// interface

void  pro_track(struct in_cmd* cmd);
void  track_observe(struct in_cmd* cmd);
void  track_run(struct in_cmd* cmd);
void  track_end(struct in_cmd* cmd);
void  track_pteigen(double* eigen);
void  track_track(struct in_cmd* cmd);
void  track_tables_dump(void);
void  track_tables_create(struct in_cmd* cmd);
void  track_start(struct command* comm);
void  track_ripple(struct in_cmd* cmd);

int next_start(double* x,double* px,double* y,double* py,double* t, double* deltae,double* fx,double* phix,double* fy,double* phiy, double* ft,double* phit);

const char* getcurrentelementname(void);
int   copytrackstoarray(void);
void  deletetrackstrarpositions(void);
int   getcurrentcmdname(char* string);
int   getnumberoftracks(void);

int   gettrack(int* nt, double* x,double* px,double* y,double* py,double* t,double* pt);

#endif // MAD_TRACK_H

