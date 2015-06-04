#ifndef MAD_PARSE_H
#define MAD_PARSE_H

// interface

char*   v_format(const char* string);
double  simple_double(char** toks, int start, int end);
int     in_spec_list(char* string);
void    pre_split(char* inbuf, struct char_array* outbuf, int fill_flag);
char*   spec_join(char**, int);

#endif // MAD_PARSE_H

