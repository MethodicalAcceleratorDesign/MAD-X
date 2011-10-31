#include "madx.h"

/*
void
cf77flush(void)
{
  fflush(stdout);
  fflush(stderr);
}
*/

void
write_nice(char* string, FILE* file)
{
  int n, pos, ssc;
  char *c = string;
  char k;
  supp_mul_char(' ', string);
  strcat(string, ";");
  n = strlen(string);
  while (n > LINE_FILL)
  {
    for (pos = LINE_FILL; pos > 10; pos--)
    {
      k = c[pos];
      if (strchr(" ,+-*/", k))  break;
    }
    c[pos] = '\0';
    fprintf(file, "%s\n", c);
    c[pos] = k;
    ssc = (uintptr_t) &c[pos] - (uintptr_t) c;
    n -= ssc;
    c = &c[pos];
  }
  fprintf(file, "%s\n", c);
}

void
write_nice_8(char* string, FILE* file)
{
  int n, pos, comma, ssc;
  char *c = string;
  char k;
  supp_mul_char(' ', string);
  strcat(string, ";");
  n = strlen(string);
  while (n > LINE_F_MAD8)
  {
    comma = 0;
    for (pos = LINE_F_MAD8; pos > 10; pos--)
    {
      k = c[pos];
      if (strchr(" ,+-*/", k))  break;
    }
    c[pos] = '\0';
    fprintf(file, "%s &\n", c);
    c[pos] = k;
    ssc = (uintptr_t) &c[pos] - (uintptr_t) c;
    n -= ssc;
    c = &c[pos];
  }
  fprintf(file, "%s\n", c);
}

