#ifndef MAD_PLOT_H
#define MAD_PLOT_H

// types

struct in_cmd;

// interface

int     embedded_plot(void);
double  plot_option(char* name);
void    exec_plot(struct in_cmd* cmd);
void    gnuplot_append(char *gplfilename, char *psfilename);
void    get_version(char* tlt, int* l);
void    get_title(char* tlt, int* l);
int     interp_node(int *nint);
int     reset_interpolation(int *nint);

#endif // MAD_PLOT_H

