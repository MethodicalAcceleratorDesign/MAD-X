#ifndef MAD_OPTION_H
#define MAD_OPTION_H

int   get_option   (const char* str);
void  set_option   (const char* str, int* opt);
void  set_defaults (const char* str); /* reset options, beam etc. to defaults */

#endif // MAD_OPTION_H
