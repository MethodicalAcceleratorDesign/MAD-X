#ifndef MAD_MATCH_H
#define MAD_MATCH_H

// types

struct in_cmd;

// interface

void pro_match(struct in_cmd* cmd);

void match_action(struct in_cmd*);
void match_cell(struct in_cmd*);
void match_constraint(struct in_cmd*);
void match_couple(struct in_cmd*);
void match_end(struct in_cmd*);
void match_fix(struct in_cmd*);
void match_gweight(struct in_cmd*);
void match_global(struct in_cmd*);
int  match_input(struct command*);  /* OB 23.1.2002 */
void match_level(struct in_cmd*);
void match_match(struct in_cmd*);
void match_rmatrix(struct in_cmd*);
void match_tmatrix(struct in_cmd*);
void match_vary(struct in_cmd*);
void match_weight(struct in_cmd*);

// used by match.f90
void mtcond(int* print_flag, int* nf, double* fun_vec, int* stab_flag);

// used by matchjc.f90
int mtputconsname(char* noden, int* nodei , char* consn, int* consi);

#endif // MAD_MATCH_H


