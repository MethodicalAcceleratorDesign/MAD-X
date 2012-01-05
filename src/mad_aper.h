#ifndef MAD_APER_H
#define MAD_APER_H

// types

struct node;
struct table;

struct aper_node            /* aperture limit node */
{
  char name[NAME_L];
  double n1;
  double s;
  char apertype[NAME_L];
  double aperture[4];
  double aper_tol[3];
};

struct aper_e_d             /* element displacement */
{
  char name[NAME_L];        /* element name */
  int curr;                 /* # of rows */
  double tab[E_D_MAX][3];   /* the table of read values */
};

// interface

double  get_apertol(struct node* node, char* par);    // used by mad_table.c
double  get_aperture(struct node* node, char* par);   // used by mad_table.c
void    pro_aperture(struct in_cmd* cmd);             // used by mad_cmd.c
struct aper_node* aperture(char *table, struct node* use_range[], struct table* tw_cp, int *tw_cnt);

#endif // MAD_APER_H

