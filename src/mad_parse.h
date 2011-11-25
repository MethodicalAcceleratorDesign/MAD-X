#ifndef MAD_PARSE_H
#define MAD_PARSE_H

// interface

int     v_length(char* form);
char*   v_format(char* string);
double  simple_double(char** toks, int start, int end);
int     get_val_num(char* in_string, int start, int end);
int     in_spec_list(char* string);
void    pre_split(char* inbuf, struct char_array* outbuf, int fill_flag);
char*   spec_join(char**, int);

// inline functions

#include <ctype.h>

static inline int
is_operand(char c) {
  return (isalnum(c) || c == '_' || c == '.');
}

static inline int
is_operator(char c) {
  return strchr("-+*/^", c) ? 1 : 0;
}

static inline int
is_expr_start(char c) {
  return strchr("-+(",c) || is_operand(c);
}

#endif // MAD_PARSE_H

