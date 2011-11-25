#ifndef MAD_ELEMERR_H
#define MAD_ELEMERR_H

// types

struct in_cmd;

// interface

void  pro_error(struct in_cmd* cmd);
void  pro_error_make_efield_table(void);
void  error_eoption(struct in_cmd* cmd);
void  error_efield(struct in_cmd* cmd);
void  error_efcomp(struct in_cmd* cmd);
void  error_eprint(struct in_cmd* cmd);
void  error_ealign(struct in_cmd* cmd);
void  error_esave(struct in_cmd* cmd);
void  error_seterr(struct in_cmd* cmd);

int   node_al_errors(double* errors);
int   node_fd_errors(double* errors);

#endif // MAD_ELEMERR_H

