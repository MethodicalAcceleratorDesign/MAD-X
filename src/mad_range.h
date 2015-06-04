#ifndef MAD_RANGE_H
#define MAD_RANGE_H

void  remove_range(char* string, const char* s1, const char* s2);
void  remove_upto(char* string, const char* s1);

void  get_bracket_range(char* string, char lb, char rb, int* rs, int* re);
void  get_bracket_t_range(char* toks[], char lb, char rb, int start, int end, int* rs, int* re);

#endif // MAD_RANGE_H

