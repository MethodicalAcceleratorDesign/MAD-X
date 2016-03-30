#define _POSIX_C_SOURCE 200112L
#include "madx.h"

//#ifndef _WIN32
#include <unistd.h>
#include <sys/stat.h>

int
intrac(void)
{
  /* returns non-zero if program is used interactively, else 0 */
  if (in->input_files[0] != stdin) return 0;

  struct stat stats;
  fstat(0, &stats);
  return (in->input_files[0] == stdin && S_ISFIFO(stats.st_mode)) || isatty(0);
}

//#else // _WIN32
//#include <io.h>
//
//int
//intrac(void)
//{
//  /* returns non-zero if program is used interactively, else 0 */
//  return _isatty(_fileno(in->input_files[0]));
//}
//
//#endif

