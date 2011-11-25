#ifndef MAD_LOGIC_H
#define MAD_LOGIC_H

// types

struct command;

// interface

int log_val(char* name, struct command* cmd);
int logic_expr(int nit, char* toks[]);
int simple_logic_expr(int nit, char* toks[]);

#endif // MAD_LOGIC_H

