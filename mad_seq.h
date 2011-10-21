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

struct sequence
{
  /* original sequence */
  char name[NAME_L];
  char export_name[NAME_L];
  char* refpos;                 /* reference position for insertion */
  int ref_flag;                 /* -1 for exit, 0 for centre, 1 for entry */
  int share;                    /* 0 normal, 1 if shared */
  int nested;                   /* 0 flat, 1 if nested */
  int con_cnt;                  /* constraint counter */
  int stamp;
  int line;                     /* set to 1 if origin is a line */
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
  struct node* ex_start;        /* first node in expanded sequence */
  struct node* ex_end;          /* last node in expanded sequence */
  struct node* range_start;     /* first node of current range in sequence */
  struct node* range_end;       /* last node of current range in sequence */
  struct node** all_nodes;      /* sequential list of all nodes */
  struct node_list* ex_nodes;   /* alphabetic list of nodes (no drifts) */
  struct table* tw_table;       /* pointer to latest twiss table created */
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

struct node*          new_sequ_node(struct sequence* sequ, int occ_cnt);
struct sequence*      new_sequence(char* name, int ref);
struct sequence_list* new_sequence_list(int length);
struct sequence*      delete_sequence(struct sequence* sequ);
struct sequence_list* delete_sequence_list(struct sequence_list* sql);

struct node* install_one(struct element* el, char* from_name, double at_value, struct expression* at_expr, double position);
void    make_sequ_from_line(char* name);
void    insert_elem(struct sequence* sequ, struct node* node);
double  sequence_length(struct sequence* sequ);
void    remove_from_sequ_list(struct sequence* sequ, struct sequence_list* sql);
void    enter_sequence(struct in_cmd* cmd);
void    make_sequ_node(struct sequence* sequ, int occ_cnt);
void    dump_exp_sequ(struct sequence* sequ, int level);
void    dump_sequ(struct sequence* c_sequ, int level);
void    write_sequs(struct sequence_list* sql,struct command_list* cl, FILE* file);
void    make_occ_list(struct sequence* sequ);
int     type_ofCall aperture_count(struct sequence* sequ);
void    enter_sequ_reference(struct in_cmd* cmd, struct sequence* sequ);
void    all_node_pos(struct sequence* sequ);
void    exec_save(struct in_cmd* cmd);
void    exec_dumpsequ(struct in_cmd* cmd);
void    seq_flatten(struct sequence* sequ);
void    expand_sequence(struct sequence* sequ, int flag);
void    expand_curr_sequ(int flag);
struct sequence* extract_sequence(char* name, struct sequence* sequ, struct node* from, struct node* to, char* refpos);
void    fill_sequ_var_list(struct sequence_list* sql, struct el_list* ell, struct var_list* varl);
void    seq_edit_ex(struct sequence* seq);
void    seq_end_ex(void);
void    add_to_sequ_list(struct sequence* sequ, struct sequence_list* sql);
void    export_sequence(struct sequence* sequ, FILE* file);
void    export_sequ_8(struct sequence* sequ, struct command_list* cl, FILE* file);
void    reset_errors(struct sequence* sequ);
void    reset_sector(struct sequence* sequ, int val);
int     restart_sequ(void);
void    seq_cycle(struct in_cmd* cmd);
void    seq_install(struct in_cmd* cmd);
void    seq_end(struct in_cmd* cmd);
void    seq_edit(struct in_cmd* cmd);
void    seq_edit_main(struct in_cmd* cmd);
void    seq_move(struct in_cmd* cmd);
void    seq_reflect(struct in_cmd* cmd);
void    seq_remove(struct in_cmd* cmd);
void    seq_replace(struct in_cmd* cmd);
void    sequence_name(char* name, int* l);
int     set_enable(char* type, struct in_cmd* cmd);
void    use_sequ(struct in_cmd* cmd);

#endif // MAD_SEQ_H

