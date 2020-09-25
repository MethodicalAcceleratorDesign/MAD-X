#ifndef MAD_ELEM_H
#define MAD_ELEM_H

// types
enum en_apertype{circle, ellipse, rectangle, lhcscreen, rectcircle, rectellipse, racetrack, octagon, custom, notdefined, custom_inter};
enum track_enums{non_existing, enum_other_bv, enum_lrad, enum_noise, enum_angle, enum_time_var};
struct node;
struct name_list;
struct command;
struct command_parameter;

struct element  /* each element is unique */
{
  char name[NAME_L];
  int def_type;                 /* 0 if defined separately,
                                   1 if inside sequence */
  int bv;                       /* bv: 0 false, 1 true (invert angle for
                                   sequence bv = -1) */
  double length;
  struct command* def;          /* pointer to defining command */
  struct element* parent;       /* pointer to parent of element */
                                /* *this for base_type elements (rbend etc.) */
  int stamp;
  struct element* base_type;    /* pointer to base_type of element */
                                /* *this for base_type elements (rbend etc.) */

  struct aperture* aper;
  double *tt_attrib;
  struct multipole* multip;
};

struct aperture
{
  enum en_apertype apertype;
  double *aper_offset;
  double *aperture;
  double *xlist;
  double *ylist;
  int length;
  int custom_inter;
};
struct multipole
{
  int nn;
  int ns;
  double *knl;
  double *ksl;
};

struct el_list /* contains list of element pointers sorted by name */
{
  int stamp;
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  struct name_list* list;       /* index list of names */
  struct element** elem;        /* element pointer list */
};

// interface

struct element* make_element(const char* name, const char* parent, struct command*, int flag);
struct element* clone_element(struct element*);
struct element* delete_element(struct element*);
void            update_element(struct element*, struct command* update);
void            update_element_children(struct element*, struct command* update);

void            dump_element(struct element*);
void            export_el_def(struct element*, char* string, int noexpr);
void            export_el_def_8(struct element*, char* string);

struct el_list* new_el_list(int length);
struct el_list* delete_el_list(struct el_list*);
struct element* find_element(const char* name, struct el_list*);
void            write_elems(struct el_list*, struct command_list*, FILE*, int noexpr);
void            write_elems_8(struct el_list*, struct command_list*, FILE*);

struct node*    new_elem_node(struct element*, int occ_cnt);
void            make_elem_node(struct element*, int occ_cnt);
char*           compound(char* e_name, int occ_cnt);

void    enter_element(struct in_cmd*);
void    element_name(char* name, int* l);
double  element_value(const struct node*, const char* par);
int     element_vector(const struct element*, const char* par, double* vector);

int     belongs_to_class(struct element*, const char*);
void    get_node_vector(const char* par, int* length, double* vector);
int     el_par_vector(int* total, double* vect);
double  el_par_value(const char* par, const struct element*);
double  el_par_value_recurse(const char* par, const struct element*);
void    fill_elem_var_list(struct element*, struct el_list*, struct var_list*);
void    add_to_el_list(struct element**, int inf, struct el_list*, int flag);
void    grow_el_list(struct el_list*);

void    set_aperture_element(struct element *el, struct command* def);
int     is_custom_set(void);
void    update_node_aperture(void);
void    check_for_update_in_seq(struct element* el, struct command* update, int nupdates);
// used by mad_mkthin.c
struct command_parameter* return_param(const char* par, const struct element*);
struct command_parameter* return_param_recurse(const char* par, const struct element*);

#endif // MAD_ELEM_H

