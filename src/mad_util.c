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

const char * 
strword(const char * str, const char * word)
{
  const char * p = str;
  for(;;) {
    p = strstr(p, word);
    if (p == NULL) break;

    if ((p==str) || !isalnum((unsigned char)p[-1])) {
      p += strlen(word);
      if (!isalnum((unsigned char)*p)) {
        break;  // found, quit
      }
    }
    // substring was found, but no word match, move by 1 char and retry
    p+=1;
  }
  return p;
}

