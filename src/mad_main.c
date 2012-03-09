#include <stdlib.h>

#include "mad_extrn_f.h"
#include "mad_core.h"

#define const // disable const for this module
#include "mad_main.h"
#undef  const

// readonly global information about program's command line arguments
int     mad_argc;
char**  mad_argv;

// readonly global information about program's stack
void*   mad_stck_base;

#ifdef _GFORTRAN
void _gfortran_set_args    (int, char *[]);
void _gfortran_set_options (int, int   []);
#endif

#ifdef _NAGFOR
void f90_init   (int, char *[]);
void f90_finish (int);
#endif

#ifdef _G95
void g95_runtime_start (int, char *[]);
void g95_runtime_stop  (void);
#endif

#ifdef _LF95
// Lahey f95 specific (requires main to be MAIN__)
int MAIN__(void);
int MAIN__(void)
{
  int a = 0;
  mad_stck_base = &a;
  mad_argc = 0;
  mad_argv = 0;
#else

int
main(int argc, char *argv[])
{
  int a = 0;
  mad_stck_base = &a;
  mad_argc = argc;
  mad_argv = argv;
#endif // _LF95

#ifdef _GFORTRAN
  _gfortran_set_args(argc, argv);
  _gfortran_set_options(0, 0);
#endif

#ifdef _NAGFOR
  f90_init(argc, argv);
#endif

#ifdef _G95
  g95_runtime_start(argc, argv);
#endif

// madx main program
  madx_start();
  madx_input(CALL_LEVEL_ZERO);
  madx_finish();

#ifdef _NAGFOR
  f90_finish(EXIT_SUCCESS);
#endif

#ifdef _G95
  g95_runtime_stop();
#endif
  
  return EXIT_SUCCESS;
}


