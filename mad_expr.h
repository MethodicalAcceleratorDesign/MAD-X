#ifndef MAD_EXPR_H
#define MAD_EXPR_H

// types

struct int_array;
struct double_array;
struct el_list;
struct var_list;
struct command_parameter;

struct expression
{
  char name[NAME_L];
  char* string;                 /* expression in string form */
  int status;                   /* status flag: 0 not evaluated
                                                1 evaluated */
  struct int_array* polish;     /* pointer to Polish notation, or NULL */
  double value;                 /* actual value */
  int stamp;
};

struct expr_list
{
  int stamp;
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  struct expression** list;     /* expression pointer list */
};

// interface

struct expression*  clone_expression(struct expression* p);
struct expression*  new_expression(char* in_string, struct int_array* polish);
struct expr_list*   new_expr_list(int length);
struct expr_list*   clone_expr_list(struct expr_list* p);
struct expression*  delete_expression(struct expression* expr);
struct expr_list*   delete_expr_list(struct expr_list* exprl);
void    grow_expr_list(struct expr_list* p);
void    dump_expression(struct expression* ex);
struct expression* make_expression(int n, char** toks);
double  expression_value(struct expression* expr, int flag); /* recursive */
void    fill_expr_list(char** toks, int s_start, int s_end, struct expr_list* p);
void    fill_expr_var_list(struct el_list* ell, struct expression* expr, struct var_list* varl);
double  double_from_expr(char** toks, int s_start, int s_end);
int     loc_expr(char** items, int nit, int start, int* end);
int     scan_expr(int c_item, char** item);
void    update_vector(struct expr_list* ell, struct double_array* da);
double  expr_combine(struct expression* exp1, double val1, char* oper, struct expression* exp2, double val2, struct expression** exp_comb);
double  combine_expr_expr(struct expression* exp1, char* oper, struct expression* exp2, struct expression** comb_exp);
double  combine_expr_val(struct expression* exp1, char* oper, double val2, struct expression** comb_exp);
double  combine_val_expr(double val1, char* oper, struct expression* exp2, struct expression** comb_exp);
struct expression*  compound_expr(struct expression* e1, double v1, char* oper, struct expression* e2, double v2);
struct expression*  scale_expr(struct expression* expr,double scale);
struct expression*  comb_param(struct command_parameter* param1, char* op, struct command_parameter* param2);

#endif // MAD_EXPR_H

