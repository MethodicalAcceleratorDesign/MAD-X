#ifndef MAD_ORBIT_H
#define MAD_ORBIT_H

// types

struct in_cmd;
struct command;

// interface

void    pro_correct(struct in_cmd*);
void    store_orbit(struct command*, double* orbit);

// for orbf.f90
void    f_ctof(int *j, char *string, int *nel);

#endif // MAD_ORBIT_H

