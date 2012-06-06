#include "mad_extrn_f.h"
#include "mad_util.h"

#ifndef _WIN32
#include <unistd.h>

int
intrac(void)
{
  /* returns non-zero inf program is used interactively, else 0 */
  return isatty(STDIN_FILENO);
}

#else // _WIN32
#include <io.h>

int
intrac(void)
{
  /* returns non-zero inf program is used interactively, else 0 */
  return _isatty(_fileno(stdin));
}

#endif

