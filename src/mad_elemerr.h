#ifndef MAD_ELEMERR_H
#define MAD_ELEMERR_H

// types

struct in_cmd;

// interface

void  pro_error(struct in_cmd* cmd);
int   node_al_errors(double* errors);
int   node_fd_errors(double* errors);
int   node_rf_errors(double* errors, double *freq, double *harmon, double *lag );

#endif // MAD_ELEMERR_H

