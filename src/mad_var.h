#ifndef MAD_VAR_H
#define MAD_VAR_H

// types

struct in_cmd;
struct name_list;
struct command_list;
struct expression;

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

void             enter_variable(struct in_cmd*); /* stores variable contained in cmd */
struct variable* new_variable(const char* name, double val, int val_type, int type, struct expression*, char* string);
double           variable_value(struct variable*);

struct var_list* new_var_list(int length);
struct var_list* clone_var_list(struct var_list*);
struct var_list* delete_var_list(struct var_list*);
void             add_to_var_list(struct variable*, struct var_list*, int flag);
struct variable* find_variable(const char* name, struct var_list*);

char*   make_string_variable(char* string);
void    write_vars(struct var_list*, struct command_list*, FILE*, int noexpr);
void    write_vars_8(struct var_list*, struct command_list*, FILE*);
void    set_variable(const char* name, double* value);
void    set_stringvar(const char* name, char* string);
double  get_variable(const char* name);
char*   get_varstring(const char* name);
void    get_defined_constants(void);

#endif // MAD_VAR_H

