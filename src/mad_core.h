#ifndef MAD_CORE_H
#define MAD_CORE_H

// constants

enum { CALL_LEVEL_ZERO }; // start call level of files in script

// interface

void madx_start(void);
void madx_input(int);
void madx_finish(void);

#endif // MAD_CORE_H

