#ifndef MAD_EVAL_H
#define MAD_EVAL_H

// types

struct name_list;

// interface

void    process(void);
void    pro_input(char* statement);
double  act_value(int pos, struct name_list* chunks);
int     act_special(int type, char* statement);

#endif // MAD_EVAL_H

