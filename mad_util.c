#include "mad_extrn_f.h"

#ifndef _WIN32
#include <unistd.h>

int
intrac(void)
{
  /* returns non-zero inf program is used interactively, else 0 */
  return isatty(STDIN_FILENO);
}

#else // _WIN32

int
intrac(void)
{
  // always interactive on windows
  return 1;
}

#endif

