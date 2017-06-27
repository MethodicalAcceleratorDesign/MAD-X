#ifndef MAD_SEQ_H
#define MAD_SEQ_H

// types

struct el_list;
struct name_list;
struct vector_list;
struct constraint_list;
struct node;
struct node_list;
struct expression;
struct command;
struct table;
struct in_cmd;
struct element;

struct sequence
{
  /* original sequence */
  char name[NAME_L];
  char export_name[NAME_L];
  char* refpos;                 /* reference position for insertion */
  char* next_sequ;              /* name of following sequence (only survey) */
  int ref_flag;                 /* -1 for exit, 0 for centre, 1 for entry */
  int share;                    /* 0 normal, 1 if shared */
  int nested;                   /* 0 flat, 1 if nested */
  int con_cnt;                  /* constraint counter */
  int stamp;
  int line;                     /* set to 1 if origin is a line */
  int add_pass;                 /* number of additional passes */
  double length;                /* length as in declaration */
  struct expression* l_expr;    /* length expression as in declaration */
  struct node* start;           /* first node in sequence */
  struct node* end;             /* last node in sequence */
  struct node_list* nodes;      /* alphabetic list of nodes */
  struct el_list* cavities;     /* alphabetic list of cavities */
  struct el_list* crabcavities;     /* alphabetic list of crab cavities */
  struct command* beam;         /* pointer to beam attached */
  /* expanded sequence */
  int n_nodes;                  /* number of nodes when expanded */
  int start_node;               /* first node of current range in all_nodes */
  int pass_count;               /* number of executed passes */
  struct node* ex_start;        /* first node in expanded sequence */
  struct node* ex_end;          /* last node in expanded sequence */
  struct node* range_start;     /* first node of current range in sequence */
  struct node* range_end;       /* last node of current range in sequence */
  struct node** all_nodes;      /* sequential list of all nodes */
  struct node_list* ex_nodes;   /* alphabetic list of nodes (no drifts) */
  struct table* tw_table;       /* pointer to latest twiss table created */
  int           tw_valid;       /* true if current tw_table is valid */
  struct constraint_list* cl;   /* pointer to constraint list during match */
  struct vector_list* orbits;   /* pointer to list of stored orbits */
};

struct sequence_list /* contains list of sequence pointers sorted by name */
{
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  struct name_list* list;       /* index list of names */
  struct sequence** sequs;      /* sequence pointer list */
  int stamp;
};

// interface

struct node*     new_sequ_node(struct sequence*, int occ_cnt);
struct sequence* new_sequence(const char* name, int ref);
struct sequence* delete_sequence(struct sequence*);
struct sequence_list* new_sequence_list(int length);
struct sequence* find_sequence(const char* name, struct sequence_list*);

void    use_sequ(struct in_cmd*);
void    remove_from_sequ_list(struct sequence*, struct sequence_list*);
double  sequence_length(struct sequence*);
void    enter_sequence(struct in_cmd*);
int     aperture_count(struct sequence*);
void    enter_sequ_reference(struct in_cmd*, struct sequence*);
void    exec_dumpsequ(struct in_cmd*);
void    exec_save(struct in_cmd*);
void    exec_extract(struct in_cmd*);
void    expand_curr_sequ(int flag);
void    add_to_sequ_list(struct sequence*, struct sequence_list*);
void    reset_errors(struct sequence*);
void    reset_sector(struct sequence*, int val);
int     restart_sequ(void);
void    seq_edit_main(struct in_cmd*);
int     set_enable(const char* type, struct in_cmd*);
void    set_sequence(char* name);
int     set_cont_sequence(void);
int     sequ_check_valid_twiss(struct sequence*);

#endif // MAD_SEQ_H

