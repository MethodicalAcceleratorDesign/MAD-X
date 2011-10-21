#ifndef MAD_RPN_H
#define MAD_RPN_H

void    deco_init(void);
int     polish_expr(int c_item, char** item);   /* split input */
double  polish_value(struct int_array* deco, char* expr_string);

#endif /// MAD_RPN_H

