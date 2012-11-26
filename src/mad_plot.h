#ifndef MAD_PLOT_H
#define MAD_PLOT_H

// types

struct in_cmd;

// interface

void    exec_plot(struct in_cmd* cmd);
double  plot_option(char* name);
void    get_title(char* tlt, int* l);
void    get_version(char* tlt, int* l);

#endif // MAD_PLOT_H

