#ifndef MAD_OPTION_H
#define MAD_OPTION_H

int   get_option(char* str);
void  set_option(char* str, int* opt);
void  set_defaults(char* string); /* reset options, beam etc. to defaults */

#endif // MAD_OPTION_H
