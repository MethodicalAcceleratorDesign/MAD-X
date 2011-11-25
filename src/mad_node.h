#ifndef MAD_NODE_H
#define MAD_NODE_H

// types

struct element;
struct command;
struct sequence;
struct expression;
struct constraint_list;
struct double_array;

struct node                /* the sequence is a linked list of nodes */
{
  char name[NAME_L];
  char* base_name;           /* basic type */
  struct node* previous;
  struct node* next;
  int share;               /* 0 normal, 1 if shared */
  int occ_cnt;             /* element occurrence count at node */
  int obs_point;           /* observation point number (tracking) */
  int sel_err;             /* error select flag */
  int sel_sector;          /* sectormap select flag */
  int con_cnt;             /* constraint counter */
  int enable;              /* flag for correctors and monitors: 0 off, 1 on */
  int moved;               /* temporary flag during sequence editing */
  int stamp;
  double position;         /* s position in sequence [m] */
  double at_value;
  double length;
  double dipole_bv;        /* +1 or -1 (if beam_bv AND element_bv) */
  double other_bv;         /* equal to beam_bv (+1 or -1) */
  double chkick;           /* calculated by orbit correction module */
  double cvkick;           /* calculated by orbit correction module */
  double match_data[74];   /* array for fast access to twiss data for match */
  struct expression* at_expr;
  char* from_name;
  struct element* p_elem;  /* pointer to element if any */
  struct sequence* p_sequ;  /* pointer to sequence if any */
  struct double_array* p_al_err; /* pointer to alignment error array */
  struct double_array* p_fd_err; /* pointer to field error array */
  struct command* savebeta; /* pointer to savebeta command if any */
  struct constraint_list* cl; /* pointer to constraint list during match */
  struct double_array* obs_orbit; /* for track observation point */
  struct double_array* orbit_ref; /* for threader orbit + cum. matrix */
};

struct node_list /* contains list of node pointers sorted by name */
{
  int stamp;
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  struct name_list* list;       /* index list of node (!) names */
                                /* node_name = el_name:occ_cnt */
  struct node** nodes;          /* node pointer list */
};

// interface

struct node*      new_node(char* name);
struct node_list* new_node_list(int length);
struct node*      clone_node(struct node* p, int flag);
struct node*      delete_node(struct node* p);
struct node*      delete_node_ring(struct node* start);
struct node_list* delete_node_list(struct node_list* l);
void    grow_node_list(struct node_list* p);
void    dump_node(struct node* node);
struct node* expand_node(struct node* node, struct sequence* top_sequ, struct sequence* sequ, double position);
double  node_value(char* par);
void    node_name(char* name, int* l);
double  get_node_pos(struct node* node, struct sequence* sequ); /*recursive */
double  hidden_node_pos(char* name, struct sequence* sequ);
void    link_in_front(struct node* new, struct node* el);
void    resequence_nodes(struct sequence* sequ);
void    store_node_value(char* par, double* value);
void    store_node_vector(char* par, int* length, double* vector);
void    add_to_node_list(struct node* node, int inf, struct node_list* nll);
int     count_nodes(struct sequence* sequ);
void    current_node_name(char* name, int* lg);
int     get_node_count(struct node* node);
int     advance_node(void);
double  line_nodes(struct char_p_array* flat);
void    node_string(char* key, char* string, int* l);
double  spec_node_value(char* par, int* number);
void    remove_from_node_list(struct node* node, struct node_list* nodes);
int     remove_one(struct node* node);
void    replace_one(struct node* node, struct element* el);
int     retreat_node(void);
void    set_node_bv(struct sequence* sequ);
void    set_new_position(struct sequence* sequ);

int type_ofCall advance_to_pos(char* table, int* t_pos);

#endif // MAD_NODE_H

