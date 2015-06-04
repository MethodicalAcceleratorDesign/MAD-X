#include "madx.h"

void
get_bracket_range(char* string, char lb, char rb, int* rs, int* re)
  /* find bracket range in string outside quotes (brackets are lb and rb) */
{
  int i, toggle = 0, level = 0, length = strlen(string);
  char quote = ' ';
  *rs = *re = -1;
  for (i = 0; i < length; i++)
  {
    if (toggle)
    {
      if (string[i] == quote)  toggle = 0;
    }
    else if(string[i] == '\'' || string[i] == '\"')
    {
      quote = string[i]; toggle = 1;
    }
    else if (string[i] == lb)
    {
      if (level++ == 0) *rs = i;
    }
    else if (string[i] == rb)
    {
      *re = i;
      if (--level == 0) return;
    }
  }
  *rs = -1;
}

void
get_bracket_t_range(char* toks[], char lb, char rb, int start, int end, int* rs, int* re)
  /* find bracket range in token list (brackets are lb and rb) */
{
  int i, level = 0;
  *rs = *re = start - 1;
  for (i = start; i <= end; i++)
  {
    if (*toks[i] == lb)
    {
      if (level++ == 0) *rs = i;
    }
    else if (*toks[i] == rb)
    {
      *re = i;
      if (--level == 0) return;
    }
  }
  *rs = start - 1;
}

void
remove_range(char* string, const char* s1, const char* s2)
  /* remove portion s1...s2 (included) in string */
{
  char* ps1 = strstr(string, s1);
  char* ps2 = strstr(string, s2);
  if (ps1 != NULL && ps2 != NULL && ps1 < ps2)
  {
    ps2++; ps2++;
    while (*ps2 != '\0')  *ps1++ = *ps2++;
    *ps1 = '\0';
  }
}

void
remove_upto(char* string, const char* s1)
  /* removes portion from start up to s2 (included) in string */
{
  char* ps1 = strstr(string, s1);
  if (ps1 != NULL)
  {
    ps1++; ps1++;
    while (*ps1 != '\0') *string++ = *ps1++;
    *string = '\0';
  }
}


