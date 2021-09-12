#include "mad_extrn_f.h"
#include "mad_core.h"
#include "mad_err.h"

#define const // disable const for this module
#include "mad_main.h"
#undef  const

// readonly global information about program's command line arguments and stack base
#include <stddef.h>
int     mad_argc = 0;
char**  mad_argv = NULL;
void*   mad_stck = NULL;

#ifdef _GFORTRAN
#define _POSIX_C_SOURCE 200112L
#include <stdio.h>
#include <stdlib.h>

void _gfortran_set_args    (int, char *[]);
void _gfortran_set_options (int, int   []);
#else
#include <stdlib.h>
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
  mad_init(0, 0);

#else

int
main(int argc, char *argv[])
{
  mad_init(argc, argv);

#endif // _LF95

  mad_run ();
  mad_fini();

  return geterrorflag() ? EXIT_FAILURE : EXIT_SUCCESS;
}

void
mad_init(int argc, char *argv[])
{
  int a = 0;
  mad_stck = &a;
  mad_argc = argc;
  mad_argv = argv;

#ifdef _GFORTRAN
  _gfortran_set_args(argc, argv);
  _gfortran_set_options(0, 0);
#endif // _GFORTRAN

#ifdef _NAGFOR
  f90_init(argc, argv);
#endif

#ifdef _G95
  g95_runtime_start(argc, argv);
#endif

  madx_start();
}

void
mad_run(void)
{
  madx_input(CALL_LEVEL_ZERO);
}

void
mad_fini(void)
{
  madx_finish();

#ifdef _NAGFOR
  f90_finish(EXIT_SUCCESS);
#endif

#ifdef _G95
  g95_runtime_stop();
#endif
}

// Special environment setup for gfortran and I/O sync with C
// Check addr: nm madx64 | grep -E -w -e "_?init_env" -e _init | sort
// Check exec: export DYLD_PRINT_INITIALIZERS=1 ; ./madx64

#ifdef _GFORTRAN

#ifdef _DARWIN
// on Darwin, linking order ensures init_env to be called before _init
static void __attribute__((constructor))
#else
// on Linux|Win, linking order doesn't ensure any calling order, use priority
static void __attribute__((constructor(101)))
#endif
init_env (void)
{
  if (getenv("GFORTRAN_UNBUFFERED_PRECONNECTED") == 0) {
#ifndef _WIN32
    // known to leak memory (duplicate the string)!
    setenv("GFORTRAN_UNBUFFERED_PRECONNECTED", "y", 0);
#else
    // known to reference the string, use permanent buffer!
    int putenv(const char *string);
    putenv("GFORTRAN_UNBUFFERED_PRECONNECTED=y");
#endif
  }
}

#endif // _GFORTRAN
