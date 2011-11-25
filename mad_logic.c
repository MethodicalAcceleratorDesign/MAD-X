#include "madx.h"

int
logic_expr(int nit, char* toks[])
{
  /* evaluates composite logical expression of type simple1 op simple2,
     where simple1 and simple2 are simple logical expressions, and
     op is either || or &&;
     returns -1 for invalid
     0 for false
     1 for true */
  int i = 0, k;
  char c = ' ';

  for (i = 0; i < nit; i++)
  {
    if (*toks[i] == '|' || *toks[i] == '&')
    {
      c = *toks[i]; break;
    }
  }
  if (i == nit) return simple_logic_expr(nit, toks);
  if (++i == nit || *toks[i] != c) return -1;
  k = simple_logic_expr(i-1, toks);
  if (c == '|')
  {
    if (k != 0) return 1;
    else if (++i == nit)  return -1;
    else return simple_logic_expr(nit-i, &toks[i]);
  }
  else
  {
    if (k == 0) return 0;
    else if (++i == nit)  return -1;
    else return simple_logic_expr(nit-i, &toks[i]);
  }
}

int
simple_logic_expr(int nit, char* toks[])
  /* evaluates a logical expression of type expr1 op expr2
     where op is one of ==, <>, <, >, <=, >=
     returns -1 if invalid
     0 if false
     1 if true */
{
  int i, t1, t2, l1_start, l1_end, l2_start, l2_end,
    logex = -1, brack = 0;
  double val1, val2;
  char c;

  for (i = 0; i < nit; i++)
  {
    c = *toks[i];
    switch (c)
    {
      case '(':
        brack++;
        break;
      case ')':
        brack--;
        break;
      case '<':
      case '>':
      case '=':
        if (brack == 0)  goto found;
    }
  }
  return -1;
  found:
  l1_start = 0;
  l2_start = (*toks[i+1] == '=' || *toks[i+1] == '>') ? i + 2 : i + 1;
  if ((t1 = loc_expr(toks, i, l1_start, &l1_end)) == 0) return -1;
  if (t1 == 2)
  {
    if (polish_expr(l1_end + 1 - l1_start, &toks[l1_start]) == 0)
      val1 = polish_value(deco, join(&toks[l1_start], l1_end + 1 - l1_start));
    else return -1;
  }
  else  val1 = simple_double(toks, l1_start, l1_end);
  if ((t2 = loc_expr(toks, nit, l2_start, &l2_end)) == 0) return -1;
  if (t2 == 2)
  {
    if (polish_expr(l2_end + 1 - l2_start, &toks[l2_start]) == 0)
      val2 = polish_value(deco, join(&toks[l2_start], l2_end + 1 - l2_start));
    else return -1;
  }
  else  val2 = simple_double(toks, l2_start, l2_end);
  if (c == '<')
  {
    if (*toks[i+1] == '=')     logex = val1 <= val2 ? 1 : 0;
    else if(*toks[i+1] == '>') logex = val1 != val2 ? 1 : 0;
    else                       logex = val1 <  val2 ? 1 : 0;
  }
  else if (c == '>')
  {
    if (*toks[i+1] == '=')  logex = val1 >= val2 ? 1 : 0;
    else                    logex = val1 >  val2 ? 1 : 0;
  }
  else if (c == '=' && *toks[i+1] == '=')  logex = val1 == val2 ? 1 : 0;
  return logex;
}

int
log_val(char* name, struct command* cmd)
  /* returns 0 = false, 1 = true for a logical command parameter */
{
  struct name_list* nl = cmd->par_names;
  struct command_parameter_list* pl = cmd->par;
  int pos = name_list_pos(name, nl);
  if (pos > -1 && nl->inform[pos]) /* "name" has beem read */
    return pl->parameters[pos]->double_value == zero ? 0 : 1;
  else return 0;
}


