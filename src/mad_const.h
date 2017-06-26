#ifndef MAD_CONST_H
#define MAD_CONST_H

// types

struct node;
struct sequence;
struct expression;
struct command;

struct constraint /* contains one constraint */
{
  char name[NAME_L];
  int  type;                    /* 1 minimum */
                                /* 2 maximum */
                                /* 3 both 1 + 2 */
                                /* 4 value */
  int stamp;
  int n_pos;
  double value, c_min, c_max, weight, evaluated;
  struct expression *ex_value, *ex_c_min, *ex_c_max;
};

struct constraint_list /* contains list of constraints */
{
  int stamp;
  char name[NAME_L];
  int  max,                           /* max. pointer array size */
       curr;                          /* current occupation */
  struct constraint** constraints;    /* command pointer list */
};

// interface

struct constraint* clone_constraint(struct constraint*);
struct constraint* delete_constraint(struct constraint*);   // used by mad_match.c
void               dump_constraint(struct constraint*);     // used by mad_node.c

struct constraint_list* new_constraint_list(int length);
struct constraint_list* delete_constraint_list(struct constraint_list*);
//void                  dump_constraint_list(struct constraint_list*);

void fill_constraint_list(int type, struct command*, struct constraint_list*);
void update_sequ_constraints(struct sequence*, struct constraint_list*);
void update_node_constraints(struct node*, struct constraint_list*);

int  constraint_name(char* name, int* name_l, int* index);
int  next_constr_namepos(char* name);
int  next_constraint(char* name, int* name_l, int* type, double* value, double* c_min, double* c_max, double* weight, int* pos, double* evaluated, char* node_name, int* nn_len);
int  next_global(char* name, int* name_l, int* type, double* value, double* c_min, double* c_max, double* weight);

#endif // MAD_CONST_H

