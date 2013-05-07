#include <stdio.h>
#include "mad_extrn_f.h"
#include "mad_core.h"
#include "mad_err.h"
#include "mad_mem.h"

#define const // disable const for this module
#include "mad_main.h"
#undef  const

// readonly global information about program's command line arguments and stack base
int     mad_argc;
char**  mad_argv;
void*   mad_stck;

#ifdef _GFORTRAN
#define _POSIX_C_SOURCE 200112L
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

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

  // LD-2012: very ugly hack to make stdout unbuffered!!! any other idea?
  if (argc && getenv("GFORTRAN_UNBUFFERED_PRECONNECTED") == 0) {
#ifdef _WIN32
    int putenv(char *string);
    putenv("GFORTRAN_UNBUFFERED_PRECONNECTED=y");
#else
    setenv("GFORTRAN_UNBUFFERED_PRECONNECTED", "y", 0);
#endif
    execvp(argv[0], argv);
    // should never be reached...
    fprintf(stderr, "fatal error: unable to synchronize Fortran versus C I/O\n");
    exit(EXIT_FAILURE);
  }
#endif

#ifdef _NAGFOR
  f90_init(argc, argv);
#endif

#ifdef _G95
  g95_runtime_start(argc, argv);
#endif

#ifdef _USEGC
  GC_INIT();
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
