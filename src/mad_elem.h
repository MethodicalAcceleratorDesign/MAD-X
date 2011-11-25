#ifndef MAD_ELEM_H
#define MAD_ELEM_H

// types

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

struct element* make_element(char* name, char* parent, struct command* def, int flag);
void            make_elem_node(struct element* el, int occ_cnt);
struct element* clone_element(struct element* el);
struct element* new_element(char* name);
struct node*    new_elem_node(struct element* el, int occ_cnt);
struct el_list* new_el_list(int length);
struct element* delete_element(struct element* el);
struct el_list* delete_el_list(struct el_list* ell);
int     belongs_to_class(struct element* el, char* class);
char*   compound(char* e_name, int occ);
double  get_refpos(struct sequence* sequ);
double  element_value(struct node* node, char* par);
int     element_vector(struct element* el, char* par, double* vector);
void    get_node_vector(char* par, int* length, double* vector);
void    element_name(char* name, int* l);
int     el_par_vector(int* total, double* vect);
double  el_par_value(char* par, struct element* el);
double  el_par_value_recurse(char* par, struct element* elem);
struct command_parameter* return_param(char* par, struct element* elem);
struct command_parameter* return_param_recurse(char* par, struct element* elem);
void    enter_element(struct in_cmd* cmd);
void    enter_elm_reference(struct in_cmd* cmd, struct element* el, int flag);
void    fill_elem_var_list(struct element* el, struct el_list* ell, struct var_list* varl);
struct element* find_element(char* name, struct el_list* ell);
void    update_element(struct element* el, struct command* update);
void    add_to_el_list(struct element** el, int inf, struct el_list* ell, int flag);
void    grow_el_list(struct el_list* p);
void    dump_element(struct element* el);
void    dump_el_list(struct el_list* ell);
void    write_elems(struct el_list* ell, struct command_list* cl, FILE* file);
void    write_elems_8(struct el_list* ell, struct command_list* cl, FILE* file);
void    export_element(struct element* el, struct el_list* ell, FILE* file);
void    export_elem_8(struct element* el, struct el_list* ell, FILE* file);
void    export_el_def(struct element* el, char* string);
void    export_el_def_8(struct element* el, char* string);
void    export_el_par_8(struct command_parameter* par, char* string);
int     par_out_flag(char* base_name, char* par_name);

#endif // MAD_ELEM_H

