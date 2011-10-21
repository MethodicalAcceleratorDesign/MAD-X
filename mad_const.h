#ifndef MAD_CONST_H
#define MAD_CONST_H

// types

struct node;
struct sequence;
struct expression;
struct command;
struct command_parameter;

struct constraint /* contains one constraint */
{
  char name[NAME_L];
  int  type;                    /* 1 minimum */
                                /* 2 maximum */
                                /* 3 both 1 + 2 */
                                /* 4 value */
  int stamp;
  double value, c_min, c_max, weight;
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

void add_to_constraint_list(struct constraint* cs, struct constraint_list* cl);
struct constraint* make_constraint(int type, struct command_parameter* par);
struct constraint* new_constraint(int type);
struct constraint_list* new_constraint_list(int length);
struct constraint* delete_constraint(struct constraint* cst);
struct constraint_list* delete_constraint_list(struct constraint_list* cl);
void  dump_constraint(struct constraint* c);
void  dump_constraint_list(struct constraint_list* cl);
void  grow_constraint_list(struct constraint_list* p);
void  fill_constraint_list(int type, struct command* cd, struct constraint_list* cl);
void  update_sequ_constraints(struct sequence* sequ, struct constraint_list* cl);
void  update_node_constraints(struct node* c_node, struct constraint_list* cl);

int   constraint_name(char* name, int* name_l, int* index);
int   next_constr_namepos(char* name);
int   next_constraint(char* name, int* name_l, int* type, double* value, double* c_min, double* c_max, double* weight);
int   next_global(char* name, int* name_l, int* type, double* value, double* c_min, double* c_max, double* weight);

#endif // MAD_CONST_H

