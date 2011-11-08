#include <stdlib.h>
#include "mad_core.h"
#include "mad_wrap_f.h"

#define const // disable const for this module
#include "mad_main.h"
#undef  const

#if 0
// readonly global information about program's command line arguments
int     mad_argc;
char**  mad_argv;

// readonly global information about program's stack
void*   mad_stck_base;
s
int
main(int argc, char *argv[])
{
  mad_stck_base = &argc;
  
  mad_argc = argc;
  mad_argv = argv;
#endif

MAIN__()
{
  madx_start();
  madx_input(CALL_LEVEL_ZERO);
  madx_finish();
  
  return EXIT_SUCCESS;
}


