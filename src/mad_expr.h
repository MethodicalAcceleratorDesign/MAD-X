#ifndef MAD_EXPR_H
#define MAD_EXPR_H

// types

struct int_array;
struct double_array;
struct el_list;
struct var_list;

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

struct expression* make_expression(int n, char** toks);
struct expression* new_expression(const char* in_string, struct int_array*);
struct expression* clone_expression(struct expression*);
struct expression* delete_expression(struct expression*);
struct expression* scale_expr(struct expression*, double scale);
struct expression* compound_expr(struct expression*, double v1, const char* oper, struct expression*, double v2, int parentheses);
double             expr_combine(struct expression*, double v1, const char* oper, struct expression*, double v2, struct expression**);
double             expression_value(struct expression*, int flag); /* recursive */
void               dump_expression(struct expression*);

struct expr_list* new_expr_list(int length);
struct expr_list* clone_expr_list(struct expr_list*);
struct expr_list* delete_expr_list(struct expr_list*);
void              grow_expr_list(struct expr_list*);
void              fill_expr_list(char** toks, int s_start, int s_end, struct expr_list*);
void              fill_expr_var_list(struct el_list*, struct expression*, struct var_list*);
void              update_vector(struct expr_list*, struct double_array*);

double  double_from_expr(char** toks, int s_start, int s_end);
int     loc_expr(char** items, int nit, int start, int* end);
int     scan_expr(int c_item, char** item);

#endif // MAD_EXPR_H

