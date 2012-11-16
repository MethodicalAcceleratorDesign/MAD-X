#define _POSIX_C_SOURCE 200112L
#include "madx.h"

#ifndef _WIN32
#include <unistd.h>

int
intrac(void)
{
  /* returns non-zero if program is used interactively, else 0 */
  return isatty(fileno(in->input_files[0]));
}

#else // _WIN32
#include <io.h>

int
intrac(void)
{
  /* returns non-zero if program is used interactively, else 0 */
  return _isatty(_fileno(in->input_files[0]));
}

#endif

