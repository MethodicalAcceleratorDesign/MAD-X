#include <stdlib.h>
#include "mad_core.h"
#include "mad_wrap_f.h"

#define const // disable const for this module
#include "mad_main.h"
#undef  const

// readonly global information about program's command line arguments
int     mad_argc;
char**  mad_argv;

// readonly global information about program's stack
void*   mad_stck_base;

#ifdef _GFORTRAN
void _gfortran_set_args (int argc, char *argv[]);
#endif

#ifdef _G95
void g95_runtime_start(int argc, char *argv[]);
void g95_runtime_stop (void);
#endif

#ifdef _LF95
// Lahey f95 specific (requires main to be MAIN__)
int
MAIN__()
{
  int a;
  mad_stck_base = &a;
  mad_argc = 0;
  mad_argv = 0;
#else
int
main(int argc, char *argv[])
{
  mad_stck_base = &argc;
  mad_argc = argc;
  mad_argv = argv;
#endif

#ifdef _GFORTRAN
// gfortran specific
  _gfortran_set_args(argc, argv);
#endif

#ifdef _G95
// g95 specific
  g95_runtime_start(argc, argv);
#endif

  madx_start();
  madx_input(CALL_LEVEL_ZERO);
  madx_finish();

#ifdef _G95
// g95 specific
  g95_runtime_stop();
#endif
  
  return EXIT_SUCCESS;
}


