#ifndef MAD_EVAL_H
#define MAD_EVAL_H

// types

struct in_cmd;
struct int_array;

// interface

void    deco_init(void);
void    process(void);
void    pro_input(char* statement);
int     polish_expr(int c_item, char** item);   /* split input */
double  polish_value(struct int_array* deco, char* expr_string);
int     act_special(int type, char* statement);
void    print_value(struct in_cmd*);

#endif // MAD_EVAL_H

