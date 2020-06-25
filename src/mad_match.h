#ifndef MAD_MATCH_H
#define MAD_MATCH_H

// types

struct in_cmd;

// constants

enum Match_Mode{ kMatch_NoMatch = 0, kMatch_Std, kMatch_UseMacro, kMatch_PTCknobs };

// interface

void pro_match(struct in_cmd* cmd);

// used by match.f90
void mtcond(int* print_flag, int* nf, double* fun_vec, int* stab_flag);

// used by matchjc.f90
int mtputconsname(char* noden, int* nodei , char* consn, int* consi);

// used by match.f90, matchjc.f90, matchsa.f90
int next_vary(char* name, int* name_l, double* lower, double* upper, double* step, int* slope, double* opt);

#endif // MAD_MATCH_H


