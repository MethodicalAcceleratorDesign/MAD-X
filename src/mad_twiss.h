#ifndef MAD_TWISS_H
#define MAD_TWISS_H

// types

struct node;
struct table;
struct command;

// interface

void  pro_twiss(void);
void  store_beta0(struct in_cmd*);
void  store_savebeta(struct in_cmd*);
int   twiss_input(struct command*);

void  get_disp0(double* disp);
void  copy_twiss_data(double* twiss_data, int* offset, int* nval, int* interp_index);
void  complete_twiss_table(struct table*);
int   embedded_twiss(void);
void  print_eigenvectors_(double *eigenvectors);
#endif // MAD_TWISS_H

