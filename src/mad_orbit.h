#ifndef MAD_ORBIT_H
#define MAD_ORBIT_H

// constants

enum Match_Mode{ kMatch_NoMatch = 0, kMatch_Std, kMatch_UseMacro, kMatch_PTCknobs };

// types

struct node;

struct val_mic {
  double before[2];
  double after[2];
};

struct id_mic {
  int   id_ttb;
  int   enable;
  struct val_mic val;
  struct node* p_node;
  struct id_mic *next;
  struct id_mic *previous;
};

struct id_mic2 {
  int   id_ttb[2];
  int   enable;
  struct val_mic val;
  struct node* p_node;
  struct node* p_node_s1;
  struct node* p_node_s2;
  struct id_mic2 *next;
  struct id_mic2 *previous;
};

struct orb_cor {
  double qx0;
  double qy0;
  double units;
  struct id_mic *cor_table;
  struct id_mic *mon_table;
};

struct orb_cor2 {
  double qx0;
  double qy0;
  double units;
  struct id_mic2 *cor_table;
  struct id_mic2 *mon_table;
};

// interface

void    pro_correct(struct in_cmd* cmd);
void    correct_correct(struct in_cmd* cmd);
void    correct_usemonitor(struct in_cmd* cmd);
void    correct_usekick(struct in_cmd* cmd);
void    correct_putorbit(struct in_cmd* cmd);
void    correct_getorbit(struct in_cmd* cmd); /* empty */
void    correct_option(struct in_cmd* cmd);
void    correct_readcorr(struct in_cmd* cmd);
void    correct_setcorr(struct in_cmd* cmd);

int     pro_correct_getactive(int ip, int *nm, int *nx, int *nc, double *corvec, double *monvec,char *conm);
void    pro_correct_write_results(double *monvec, double *resvec, double *corvec, int *nx, int *nc, int *nm, int imon, int icor, int ip);
void    pro_correct_fill_mon_table(int ip ,char *name, double old, double new_);
void    pro_correct_fill_corr_table(int ip ,char *name, double old, double new_);
void    pro_correct_make_mon_table(void);
void    pro_correct_make_corr_table(void);
double* pro_correct_response_line(int ip, int nc, int nm);
double* pro_correct_response_ring(int ip, int nc, int nm);
int     pro_correct_filter(int iplane, double sigcut);
void    pro_correct_write_cocu_table(void);
void    pro_correct_prtwiss(void);
int     pro_correct_getcorrs(struct in_cmd* cmd);
int     pro_correct_getorbit_ext(struct in_cmd* cmd);
int     pro_correct_getorbit(struct in_cmd* cmd);
int     pro_correct_gettables(int iplane, struct in_cmd* cmd);
int     pro_correct_getcommands(struct in_cmd* cmd);
void    pro_correct_option(struct in_cmd* cmd);

void    correct_correct1(struct in_cmd* cmd);

void    correct_correct2(struct in_cmd* cmd);
void    pro_correct2_fill_mon_table(int ip ,char *name, double old, double new_);
void    pro_correct2_fill_corr_table(int b, int ip ,char *name, double old, double new_);
void    pro_correct2_make_mon_table(void);
void    pro_correct2_make_corr_table(void);
double* pro_correct2_response_ring(int ip, int nc, int nm);
int     pro_correct2_getactive(int ip, int *nm, int *nx, int *nc, double *corvec, double *monvec,char *conm);
int     pro_correct2_getcorrs(struct in_cmd* cmd);
int     pro_correct2_getorbit(struct in_cmd* cmd);
int     pro_correct2_gettables(int iplane, struct in_cmd* cmd);
void    pro_correct2_write_results(double *monvec, double *resvec, double *corvec, int *nx, int *nc, int *nm, int imon, int icor, int ip);

void    fill_orbit_table(struct table* t_out, struct table* t_in);
void    store_orbit(struct command* comm, double* orbit);

// wrappers around orbf.f90 functions
int     c_micit(double *dmat,char *conm, double *monvec,double *corvec,double *resvec,int *nx,float rms,int imon,int icor,int niter);
void    c_haveit(double *dmat,double *monvec,double *corvec,double *resvec,int *nx,int imon,int icor);
int     c_svddec(double *dmat, int imon, int icor, int *sing, double *sngcut, double *sngval);
int     c_svdcorr(double *dmat, double *xin, double *cor, double *res, int *nx, int imon, int icor);

// from orbf.f90
void setupi_(int*, int*, int*, int*, int*, int*);
void primat_(int*, int*, int*);
void prdmat_(double*, int*, int*);

// for orbf.f90
uintptr_t locf_(char *iadr);
void      f_ctof(int *j, char *string, int *nel);

#endif // MAD_ORBIT_H

