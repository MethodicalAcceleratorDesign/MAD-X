#ifndef MAD_VAR_H
#define MAD_VAR_H

// types

struct name_list;
struct command_list;
struct expression;
struct in_cmd;

struct constant
{
  char name[NAME_L];
  struct expression* exp;     /* pointer to defining expression (always) */
  double value;
  int stamp;
};

struct variable
{
  char name[NAME_L];
  int status;                 /* 0 value not evaluated, 1 evaluated */
  int type;                   /* 0 constant, 1 direct, 2 deferred, 3 string */
  int val_type;               /* 0 int 1 double (0..2) */
  char* string;               /* pointer to string if 3 */
  struct expression* expr;    /* pointer to defining expression (0..2) */
  double value;               /* (0..2) */
  int stamp;
};

struct var_list         /* contains list of variable pointers sorted by name */
{
  int stamp;
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  struct name_list* list;       /* index list of names */
  struct variable** vars;       /* variable pointer list */
};

// interface

struct variable* new_variable(char* name, double val, int val_type, int type, struct expression* expr, char* string);
struct var_list* new_var_list(int length);
struct var_list* clone_var_list(struct var_list* vl);
struct variable* delete_variable(struct variable* var);
struct var_list* delete_var_list(struct var_list* varl);
void    get_defined_constants(void);
void    grow_var_list(struct var_list* p);
void    dump_variable(struct variable* v);
void    write_vars(struct var_list* varl, struct command_list* cl, FILE* file);
void    write_vars_8(struct var_list* varl, struct command_list* cl, FILE* file);
char*   make_string_variable(char* string);
double  variable_value(struct variable* var);
double  get_variable(char* name);
char*   get_varstring(char* name);
void    enter_variable(struct in_cmd* cmd); /* stores variable contained in cmd */
struct variable* find_variable(char* name, struct var_list* varl);
void    set_sub_variable(char* comm, char* par, struct in_cmd* cmd);
void    add_to_var_list(struct variable* var, struct var_list* varl, int flag);
void    set_variable(char* name, double* value);
void    set_stringvar(char* name, char* string);
void    export_variable(struct variable* var, FILE* file);
void    export_var_8(struct variable* var, FILE* file);
int     predef_var(struct variable* var);
void    print_global(double delta);
// int     vary_name(char* name, int* name_l, int* index); // not used...

#endif // MAD_VAR_H

