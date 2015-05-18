#ifndef MAD_APER_H
#define MAD_APER_H

// types

struct node;
struct in_cmd;

// interface

double  get_apertol(struct node* node, const char* par);    // used by mad_table.c
double  get_aperture(struct node* node, const char* par);   // used by mad_table.c
void    pro_aperture(struct in_cmd* cmd);             // used by mad_cmd.c

#endif // MAD_APER_H

