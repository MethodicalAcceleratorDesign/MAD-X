#ifndef MAD_TWISS_H
#define MAD_TWISS_H

// types

struct node;
struct table;
struct command;

// interface

void  get_disp0(double* disp);
void  copy_twiss_data(double* twiss_data);
void  get_twiss_data(double* twiss_data);
void  complete_twiss_table(struct table* t);
void  exec_savebeta(void);
void  fill_beta0(struct command* beta0, struct node* node);
void  fill_twiss_header(struct table* t);
int   embedded_twiss(void);
void  pro_embedded_twiss(struct command* current_global_twiss);
void  pro_twiss(void);
void  set_twiss_deltas(struct command* comm);
void  store_beta0(struct in_cmd* cmd);
void  store_savebeta(struct in_cmd* cmd);
int   twiss_input(struct command* tw);

#endif // MAD_TWISS_H

