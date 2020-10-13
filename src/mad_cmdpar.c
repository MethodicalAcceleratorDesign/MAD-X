#include "madx.h"

// public interface

struct command_parameter*
clone_command_parameter(const struct command_parameter* p)
{
  assert(p);
  struct command_parameter* clone = new_command_parameter(p->name, p->type);
  clone->call_def = p->call_def;
  switch (p->type)
  {
    case 4: // constraint
      clone->c_min = p->c_min;
      clone->c_max = p->c_max;
      clone->min_expr = clone_expression(p->min_expr);
      clone->max_expr = clone_expression(p->max_expr);
      /* FALLTHRU */

    case 0: // logical (not supported...)
    case 1: // integer (not supported...)
    case 2: // double
      clone->double_value = p->double_value;
      clone->expr = clone_expression(p->expr);
      break;
    case 3: // string
      clone->string = permbuff(p->string);
      clone->expr = NULL;
      break;
    case 11: // array of integers (not supported...)
    case 12: // array if doubles
      clone->double_array = clone_double_array(p->double_array);
      clone->expr_list = clone_expr_list(p->expr_list);
      break;
    case 13: // array of strings
      clone->m_string = clone_char_p_array(p->m_string);
  }
  return clone;
}

struct command_parameter*
new_command_parameter(const char* name, int type)
{
  const char *rout_name = "new_command_parameter";
  struct command_parameter* new = mycalloc(rout_name, 1, sizeof *new);
  strcpy(new->name, name); new->type = type;
  new->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", new->name);
  return new;
}

struct command_parameter_list*
new_command_parameter_list(int length)
{
  const char *rout_name = "new_command_parameter_list";
  struct command_parameter_list* il = mycalloc(rout_name, 1, sizeof *il);
  strcpy(il->name, "command_parameter_list");
  il->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", il->name);
  il->curr = 0;
  il->max = length;
  if (length > 0)
    il->parameters = mycalloc(rout_name, length, sizeof *il->parameters);
  return il;
}

struct command_parameter*
renew_command_parameter(struct command* cmd, const char* par)
/* Clone the given command parameter in the commands parameter list. This is
   useful for only replacing few writable parameters after clone_command_flat
   while sharing the other parameters. */
{
  int i = name_list_pos(par, cmd->par_names);
  if (i > -1) {
    return cmd->par->parameters[i] =
      clone_command_parameter(cmd->par->parameters[i]);
  }
  return NULL;
}

struct command_parameter*
delete_command_parameter(struct command_parameter* par)
{
  const char *rout_name = "delete_command_parameter";
  if (par == NULL) return NULL;
  if (stamp_flag && par->stamp != 123456)
    fprintf(stamp_file, "d_c_p double delete --> %s\n", par->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", par->name);
  if (par->expr != NULL)         delete_expression(par->expr);
  if (par->min_expr != NULL)     delete_expression(par->min_expr);
  if (par->max_expr != NULL)     delete_expression(par->max_expr);
  if (par->double_array != NULL) delete_double_array(par->double_array);
  if (par->expr_list != NULL)    delete_expr_list(par->expr_list);
  if (par->m_string != NULL)     delete_char_p_array(par->m_string, 0);
  myfree(rout_name, par);
  return NULL;
}

struct command_parameter_list*
delete_command_parameter_list(struct command_parameter_list* parl)
{
  const char *rout_name = "delete_command_parameter_list";
  if (parl == NULL) return NULL;
  if (stamp_flag && parl->stamp != 123456)
    fprintf(stamp_file, "d_c_p_l double delete --> %s\n", parl->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", parl->name);
  if (parl->parameters != NULL)
  {
    for (int i = 0; i < parl->curr; i++)
    {
      if (parl->parameters[i] != NULL)
      {
        parl->parameters[i] = delete_command_parameter(parl->parameters[i]);
      }
    }
    if (parl->parameters)  myfree(rout_name, parl->parameters);
  }
  myfree(rout_name, parl);
  return NULL;
}

void
dump_command_parameter(struct command_parameter* par)
{
  int i, k;
  char logic[2][8] = {"false", "true"};
  fprintf(prt_file, "parameter: %s   ", par->name);
  switch (par->type)
  {
    case 0:
      k = par->double_value;
      fprintf(prt_file, "logical: %s\n", logic[k]);
      break;
    case 1:
      if (par->expr != NULL)
      {
        dump_expression(par->expr);
        par->double_value = expression_value(par->expr, 2);
      }
      k = par->double_value;
      fprintf(prt_file, v_format("integer: %I\n"), k);
      break;
    case 2:
      if (par->expr != NULL)
      {
        dump_expression(par->expr);
        par->double_value = expression_value(par->expr, 2);
      }
      fprintf(prt_file, v_format("double value: %F\n"), par->double_value);
      break;
    case 11:
    case 12:
      if (par->double_array != NULL)
      {
        if (par->expr_list != NULL)
        {
          for (i = 0; i < par->double_array->curr; i++)
          {
            if (i < par->expr_list->curr && par->expr_list->list[i] != NULL)
              par->double_array->a[i]
                = expression_value(par->expr_list->list[i], 2);
          }
        }
        fprintf(prt_file, "double array: ");
        for (i = 0; i < par->double_array->curr; i++)
          fprintf(prt_file, v_format("%e "), par->double_array->a[i]);
        fprintf(prt_file, "\n");
      }
      break;
    case 3:
      fprintf(prt_file, "string: %s\n", par->string);
      break;
    case 13:
      dump_char_p_array(par->m_string);
  }
}

void
print_command_parameter(struct command_parameter* par)
{
  int i, k;
  char logic[2][8] = {"false", "true"};
  switch (par->type)
  {
    case 0:
      k = par->double_value;
      fprintf(prt_file, "%s = %s, ", par->name, logic[k]);
      break;
    case 1:
      k = par->double_value;
      fprintf(prt_file, v_format("%s = %I, "), par->name, k);
      break;
    case 2:
      fprintf(prt_file, v_format("%s = %F, "), par->name, par->double_value);
      break;
    case 11:
    case 12:
      if (par->double_array != NULL)
      {
        fprintf(prt_file, "double array: ");
        for (i = 0; i < par->double_array->curr; i++)
          fprintf(prt_file, v_format("%F, "), par->double_array->a[i]);
        fprintf(prt_file, "\n");
      }
      break;
    case 3:
      fprintf(prt_file, "%s = %s, ", par->name, par->string);
  }
}

void
grow_command_parameter_list(struct command_parameter_list* p)
{
  const char *rout_name = "grow_command_parameter_list";
  p->max *= 2;
  if (p->max == 0) p->max++;
  p->parameters = myrecalloc(rout_name, p->parameters, p->curr * sizeof *p->parameters, p->max * sizeof *p->parameters);
}

void
export_comm_par(struct command_parameter* par, char* string, int noexpr)
  /* exports a command parameter */
{
  int last;
  char num[2*NAME_L];
  strcat(string, ",");
  strcat(string, par->name);
  switch(par->type)
  {
    case 0:
      strcat(string, "=");
      if (par->double_value == zero) strcat(string, "false");
      else                           strcat(string, "true");
      break;
    case 1:
    case 2:
      if (par->expr != NULL && !noexpr) {
        strcat(string, ":=");
        strcat(string, par->expr->string);
      }
      else
      {
        if (par->expr != NULL)
          par->double_value = expression_value(par->expr, 2);

        // LD: use := for backward compatibility despite that it is a number...
        strcat(string, noexpr ? "=" : ":=");

        if (par->type == 1)
        {
          int k = par->double_value;
          sprintf(num, v_format("%I"), k);
        }
        else
          sprintf(num, v_format("%F"), par->double_value);
        strcat(string, supp_tb(num));
      }
      break;
    case 3:
      if (par->string != NULL)
      {
        strcat(string, "=");
        strcat(string, par->string);
      }
      break;
    case 11:
    case 12: {
      strcat(string, noexpr ? "=" : ":=");
      for (last = par->double_array->curr-1; last > 0; last--)
      {
        if (par->expr_list->list[last] != NULL)
        {
          if (zero_string(par->expr_list->list[last]->string) == 0) break;
        }
        else if (par->double_array->a[last] != zero) break;
      }

      strcat(string, "{");
      for (int i = 0; i <= last; i++) {

        if (i > 0) strcat(string, ",");

        if (par->expr_list->list[i] != NULL && noexpr)
          par->double_array->a[i] = expression_value(par->expr_list->list[i], 2);

        if (par->expr_list->list[i] != NULL && !noexpr) {
          strcat(string, par->expr_list->list[i]->string);
        }
        else {
          if (par->type == 11) {
            int k = par->double_array->a[i];
            sprintf(num, v_format("%I"), k);
          }
          else {
            sprintf(num, v_format("%F"), par->double_array->a[i]);
          }
          strcat(string, supp_tb(num));
        }
      }
      strcat(string, "}");
      break;
    }
    default: mad_error("export command param", "invalid type %d", par->type);
  }
}

int
command_par(const char* parameter, const struct command* cmd, struct command_parameter** cp)
  /* returns inform and the command parameter if found, else 0, NULL */
{
  if (cmd && cmd->par_names) {
    int i = name_list_pos(parameter, cmd->par_names);
    if (i > -1) {
      *cp = cmd->par->parameters[i];
      return cmd->par_names->inform[i];
    }
  }
  *cp = NULL;
  return 0;
}

struct expression*
command_par_expr(const char* parameter, struct command* cmd)
  /* returns a command parameter expression if found, else NULL */
{
  struct command_parameter* cp;
  command_par(parameter, cmd, &cp);
  return cp ? cp->expr : NULL;
}

double
command_par_special(const char* parameter, const struct element* el)
/* construct missing tilt from normal and skew  */
{
  double val = zero;

  if (strcmp(parameter, "tilt") == 0)
  {
    if ((val = command_par_value("tilt", el->def)) == zero)
    {
      val = zero;
    }
  }
  else val = command_par_value(parameter, el->def);
  return val;
}

char*
command_par_string(const char* parameter, const struct command* cmd)
  /* returns a command parameter string if found, else NULL */
{
  struct command_parameter* cp;
  command_par(parameter, cmd, &cp);
  return cp && cp->type == 3 ? cp->string : NULL;
}

char*
command_par_string_user(const char* parameter, const struct command* cmd)
  /* get command parameter string if explicitly set, call_def if passed
   * as bareword, else NULL. This will never return the DEFAULT value. */
{
  char* str;
  command_par_string_user2(parameter, cmd, &str);
  return str;
}

int
command_par_string_user2(const char* parameter, const struct command* cmd, char** str)
  /* get a command parameter string if explicitly set, call_def if passed
   * as bareword, else NULL. This will never return the DEFAULT value.
   * Returns `inform`, i.e. whether the parameter was specified explicitly.
   * This can be `1` even if `*str = NULL`, but `0` implies `*str=NULL` */
{
  *str = NULL;
  struct command_parameter* cp;
  int inf = command_par(parameter, cmd, &cp);
  if (inf && cp && cp->type == 3) {
    if (cp->string)        *str = cp->string;
    else if (cp->call_def) *str = cp->call_def->string;
  }
  return inf;
}

int
command_par_value_user2(const char* parameter, const struct command* cmd, double* val)
  /* returns a command parameter value val
     if found returns 1, else 0 */
{
  struct command_parameter* cp;
  int inf=command_par(parameter, cmd, &cp);
  if (inf && cp && cp->type < 3) {
    *val = cp->expr ? expression_value(cp->expr, 2) : cp->double_value;
  }
  else
   {
    *val = 0;
   }
  return inf;
}

int
command_par_string_or_calldef(const char* parameter, const struct command* cmd, char** str)
  /* returns command parameter string if explicitly set, otherwise call_def. */
  /* this function is only needed because mad_dict.c is not populated properly
   * and because the command decoding prefers
   *      `cp->string = DEFAULT`    over    `cp->string = NULL`
   * if the parameter is specified as bareword and DEFAULT is non-null
   * (which is IMO a bug).
   */
{
  *str = NULL;
  struct command_parameter* cp;
  int inf = command_par(parameter, cmd, &cp);
  if (cp && cp->type == 3) {
    if (inf && cp->string) *str = cp->string;
    else if (cp->call_def) *str = cp->call_def->string;
  }
  return inf;
}

void
add_cmd_parameter_clone(struct command* cmd, const struct command_parameter *param, const char* par_name, int inf)
/* hbu add an identical copy (clone) of param to cmd */
{
  if(param)
  {
    cmd->par->parameters[cmd->par->curr] = clone_command_parameter(param); /* set current to identical copy (clone) of param */
    add_to_name_list(par_name,inf,cmd->par_names);
    cmd->par->curr++;
  }
}

void
add_cmd_parameter_new(struct command* cmd, double par_value, const char* par_name, int inf)
/*hbu add a new param with one value to cmd */
{
  cmd->par->parameters[cmd->par->curr] = new_command_parameter(par_name, 2);
  cmd->par->parameters[cmd->par->curr]->double_value = par_value;
  add_to_name_list(par_name,inf,cmd->par_names);
  cmd->par->curr++;
}

double
command_par_value(const char* parameter, const struct command* cmd)
  /* returns a command parameter value if found, else zero */
{
  /*printf("command_par_value par = >>%s<<, c=%p\n",parameter,cmd);*/

  double value;
  command_par_value2(parameter, cmd, &value);
  return value;
}

int
command_par_value2(const char* parameter, const struct command* cmd, double* val)
  /* returns a command parameter value val
     if found returns 1, else 0 */
{
  struct command_parameter* cp;
  command_par(parameter, cmd, &cp);
  if (cp && cp->type < 3) {
    *val = cp->expr ? expression_value(cp->expr, 2) : cp->double_value;
    return 1;
  }
  *val = 0;
  return 0;
}


struct double_array*
command_par_array(const char* parameter, struct command* cmd)
  /* returns an updated command parameter array if found, else NULL */
{
  struct command_parameter* cp;
  command_par(parameter, cmd, &cp);
  if (cp && (cp->type == 11 || cp->type == 12)) {
    if (cp->expr_list)
      update_vector(cp->expr_list, cp->double_array);
    return cp->double_array;
  }
  return NULL;
}

int
command_par_vector(const char* parameter, struct command* cmd, double* vector)
  /* returns the length of, and an updated command parameter vector if not null
     if found, else 0 */

{
  struct double_array* array = command_par_array(parameter, cmd);
  if (!array) return 0;
  if (vector) copy_double(array->a, vector, array->curr);
  return array->curr;
}

void
set_command_par_value(const char* parameter, struct command* cmd, double val)
{
  struct command_parameter* cp;
  int i;
  if ((i = name_list_pos(parameter, cmd->par_names)) > -1)
  {
    cp = cmd->par->parameters[i];
    if (cp->type < 3)
    {
      cp->double_value = val;
      if (cp->expr != NULL) cp->expr = delete_expression(cp->expr);
      cmd->par_names->inform[i] = 1; /* mark as set */
    }
  }
}

struct command_parameter*
store_comm_par_def(char* toks[], int start, int end)
{
  struct command_parameter *pl[2];
  int i, j, jj, k, n, dummy, type, s_start, s_end, ss_end;
  char c = *toks[start];

  if (c == 'l')      type = 0;
  else if (c == 'i') type = 1;
  else if (c == 'r') type = 2;
  else if (c == 's') type = 3;
  else if (c == 'c') type = 4;
  else               return NULL; /* error */
  pl[0] = new_command_parameter(toks[start-3], type);
  pl[1] = pl[0]->call_def = new_command_parameter(toks[start-3], type);
  if (++start == end) return pl[0]; /* empty (r) etc. */
  start++; /* skip , */
  switch (type)
  {
    case 0:
      jj = 0;
      for (j = 0; j <= imin((end - start),2); j++)
      {
        if (strcmp(toks[start+j], "true") == 0
            )  pl[jj]->double_value = 1;
        jj++; j++; /* skip , */
      }
      break;
    case 1: /* int */
    case 2: /* double */
    case 4: /* constraint */
      for (i = 0; i < 2; i++)
      {
        if (start <= end)
        {
          get_bracket_t_range(toks, '{', '}', start, end, &s_start, &s_end);
          if (s_start < start) /* no curly bracket */
          {
            if (pl[0]->type > 10) /* no call defaults */
            {
              pl[1]->double_array = pl[0]->double_array;
              break;
            }
            if ((n = next_char(',', toks, start, end+1)) < 0) n = end+1;
            if ((dummy = loc_expr(toks, n, start, &k)) > 1)
            {
              if (polish_expr(k + 1 - start, &toks[start]) ==  0)
              {
                pl[i]->expr =
                  new_expression(join(&toks[start], k+1-start), deco);
                pl[i]->double_value = expression_value(pl[i]->expr, 2);
              }
            }
            else if (dummy > 0)
              pl[i]->double_value  = simple_double(toks, start, k);
            start = k + 2; /* skip , */
          }
          else
          {
            start = s_end + 1; s_start++; s_end--;
            pl[i]->double_array = new_double_array(s_end + 1 - s_start);
            pl[i]->expr_list = new_expr_list(s_end + 1 - s_start);
            fill_expr_list(toks, s_start, s_end, pl[i]->expr_list);
            update_vector(pl[i]->expr_list, pl[i]->double_array);
            pl[i]->type += 10;
          }
        }
      }
      break;
    case 3: /* string */
      if (start <= end)
      {
        get_bracket_t_range(toks, '{', '}', start, end, &s_start, &s_end);
        if (s_start < start) /* no curly bracket */
        {
          for (ss_end = start+1; ss_end < end; ss_end++)
          {
            if (*toks[ss_end] == ',') break;
          }
          if (strcmp(toks[start], none) != 0)
          {
            if (ss_end == start+1) pl[0]->string = permbuff(toks[start]);
            else pl[0]->string = permbuff(join(&toks[start], ss_end-start));
          }
          if (ss_end < end)
          {
            start = ss_end + 1; /* skip , */
            if (strcmp(toks[start], none) != 0)
              pl[1]->string = permbuff(toks[start]);
          }
        }
        else
        {
          start = s_end + 1; s_start++; s_end--;
          pl[0]->m_string = new_char_p_array(s_end + 1 - s_start);
          i = 0;
          for (j = 0; j < pl[0]->m_string->max; j++)
          {
            if (*toks[s_start+j] != ',' && strcmp(toks[s_start+j], none) != 0)
              pl[0]->m_string->p[i++] = permbuff(toks[s_start+j]);
          }
          pl[0]->m_string->curr = i;
          pl[0]->type += 10;
        }
      }
  }
  return pl[0];
}

void
store_comm_par_value(const char* parameter, double val, struct command* cmd)
{
  struct command_parameter* cp;
  command_par(parameter, cmd, &cp);
  if (cp) {
    cp->type = 2;
    if(cp->expr != NULL) cp->expr = delete_expression(cp->expr);
    cp->double_value = val;
  }
}

void
store_set(struct command* comm, int flag)
{
  char* p;
  char* name;
  current_format_f = clone_command(comm);
  struct command_parameter* cp;
  int i, lp, n = 0;
  if (command_par("format", comm, &cp) || !flag)
  {
    n++;
    for (i = 0; i < cp->m_string->curr; i++)
    {
      p = noquote(cp->m_string->p[i]);
      if (strchr(p, 's')) strcpy(string_format, p);
      else if (strpbrk(p, "id")) strcpy(int_format, p);
      else if (strpbrk(p, "feEgG")) strcpy(float_format, p);
      else if (strpbrk(p, "feEgGA")) strcpy(float_format, p);
    }
  }
  if (flag && command_par("sequence", comm, &cp))
  {
    n++;
    name = cp->string;
    if ((lp = name_list_pos(name, sequences->list)) > -1
        && sequences->sequs[lp]->ex_start != NULL)
      current_sequ = sequences->sequs[lp];
    else
    {
      warning("ignoring unknown or unused sequence:", name);
      return;
    }
  }
  if (n == 0)
  {
    warning("no parameter specified,", "ignored");
    return;
  }
}

void
fill_par_var_list(struct el_list* ell, struct command_parameter* par, struct var_list* varl)
  /* puts all variables an element parameter depends on, in a list */
{
  int i;
  switch (par->type)
  {
    case 1:
    case 2:
      if (par->expr != NULL) fill_expr_var_list(ell, par->expr, varl);
      break;
    case 11:
    case 12:
      for (i = 0; i < par->double_array->curr; i++)
        if (i < par->expr_list->curr && par->expr_list->list[i] != NULL)
          fill_expr_var_list(ell, par->expr_list->list[i], varl);
      break;
  }
}

int
par_present(const char* par, struct command* cmd)
  /* returns 1 if in cmd par is read, else returns 0 */
{
  struct command_parameter* cp;
  return command_par(par, cmd, &cp);
}

int
par_present_list(const char* par, struct command_list* c_list)
  /* returns 1 if in c_list par is read, else returns 0 */
{
  if (!c_list)
    return 0;
  for (int i = 0; i < c_list->curr; i++) {
    if (par_present(par, c_list->commands[i]))
      return 1;
  }
  return 0;
}

void
comm_para(const char* name, int* n_int, int* n_double, int* n_string,
  int* int_array, double* double_array, char* strings, int* string_lengths)
  /* returns the value for command parameter "name" being either
     one or several integers (including logicals),
     one or several doubles,
     one or several strings (packed in one, with length array)
     Input:
     name                  parameter name
     Output:
     n_int                 # integers
     n_double              # double
     n_string              # strings
     int_array             array for integers
     double_array          array for doubles
     strings               one string for all, packed
     string_lengths        length of each string in char

     ATTENTION: no check on sufficient array sizes
  */
{
  char buf[NAME_L];
  int i, l;
  struct command_parameter* cp;
  struct double_array* arr = NULL;
  *n_int = *n_double = *n_string = 0;

  mycpy(buf, name);
//  printf("comm_para: name(buf)='%s'\n", buf);

  if (this_cmd != NULL && this_cmd->clone != NULL) {
    command_par(buf, this_cmd->clone, &cp);
    if (cp) {

      switch (cp->type) {
      case 0:
	*n_int = 1;
	*int_array = cp->double_value;
	break;
      case 1:
	*n_int = 1;
	if (cp->expr == NULL) *int_array = cp->double_value;
	else *int_array = expression_value(cp->expr, 2);
	break;
      case 2:
	*n_double = 1;
	if (cp->expr == NULL) *double_array = cp->double_value;
	else *double_array = expression_value(cp->expr, 2);
	break;
      case 3:
	if (cp->string != NULL) {
	  *n_string = 1;
	  strfcpy(strings, cp->string, NAME_L);
          l = *string_lengths = imin(strlen(cp->string), NAME_L);
//          printf("comm_para3: strings='%s'[%d]\n", strings, l);
	}
	break;
      case 11:
      case 12:
	arr = cp->double_array;
	if (cp->expr_list != NULL) update_vector(cp->expr_list, arr);
	if (cp->type == 11) {
            for (i = 0; i < arr->curr; i++) int_array[i] = arr->a[i];
            *n_int = arr->curr;
	} else {
	  for (i = 0; i < arr->curr; i++) double_array[i] = arr->a[i];
	  *n_double = arr->curr;
	}
	break;
      case 13:
	for (i = 0; i < cp->m_string->curr; i++) {
	  strfcpy(strings, cp->m_string->p[i], NAME_L);
          l = string_lengths[i] = imin(strlen(cp->m_string->p[i]), NAME_L);
//          printf("comm_para13: strings='%s'[%d]\n", strings, l);
	  strings += l;
	}
	*n_string = cp->m_string->curr;
      }
    }
  }
}

void
store_comm_par_vector(const char* parameter, double* val, struct command* cmd)
{
  struct command_parameter* cp;
  command_par(parameter, cmd, &cp);
  if (cp && cp->double_array) {
      copy_double(val, cp->double_array->a, cp->double_array->curr);
      if (cp->expr_list != NULL)
        cp->expr_list = delete_expr_list(cp->expr_list);
  }
}


// LD: this function is horrible and probably buggy...
int
decode_par(struct in_cmd* cmd, int start, int number, int pos, int log)
{
  /* matches + stores in clone one keyword + following expression(s);
     returns the position of last token
     or -i where i is the start of an illegal item.
  */
  char** toks = cmd->tok_list->p;
  struct expression* expr = NULL;
  struct command_parameter* lp = cmd->cmd_def->par->parameters[pos];
  struct command_parameter* clp = cmd->clone->par->parameters[pos];
  int j, k, ks, i = start, e_type, ival, end, e_end, tot_end = 0, c_type = 0,
      val_type = 0, cnt = 0, con_flag = 0, t_num;
  double val = zero;

  if (lp->type < 10)
  {
    if (lp->type == 0)
    {
      start = i - log;
      if (log == 0)
      {
        if (i+2 < number && *toks[i+1] == '=')
        {
          if     (strncmp(toks[i+2], "true", 4) == 0)  ival = 1;
          else if(strncmp(toks[i+2], "false", 5) == 0) ival = 0;
          else return -i;
          end = i+2;
        }
        else
        {
          end = i;
          ival = lp->call_def->double_value;
        }
      }
      else
      {
        end = i;
        ival = 0;
      }
      tot_end = end;
      clp->double_value = ival;
    }
    else if (lp->type == 3)  /* string */
    {
      if (i+2 < number && *toks[i+1] == '=')
      {
        i += 2;
        for (j = i; j < number; j++)
          if (name_list_pos(alias(toks[j]), cmd->cmd_def->par_names) >= 0) break;
//        dirty quick fix for ticket #165
          if (*toks[j-1] == '-') j--;
          while (*toks[j-1] == ',') j--;
        tot_end = j - 1;
        clp->string = permbuff(noquote(join(&toks[i], j - i)));
      }
      else tot_end = i;
    }
    else if (lp->type < 3)  /* one int or double value */
    {
      if ((i+2 < number && *toks[i+1] == '=') ||
          (i+3 < number && *toks[i+1] == ':' && *toks[i+2] == '='))
      {
        val_type = (*toks[i+1] == ':' && *toks[i+2] == '=') ? 1 : 0;
        start = val_type + i + 2;
        for (t_num = start; t_num < number; t_num++) if(*toks[t_num] == ',')
          break;
        if ((e_type = loc_expr(toks, t_num, start, &end)) == 0) return -i;
        tot_end = end;
        if (e_type == 1) /* simple number */
        {
          val = simple_double(toks, start, end);
          clp->expr = NULL;
        }
        else /* expression */
        {
          if ((expr =
               make_expression(end + 1 - start, &toks[start])) == NULL)
            return -i;
          if (val_type) /* definition with ":=" */
          {
            val = zero;
            clp->expr = clone_expression(expr);
          }
          else
          {
            val = expression_value(expr, 2);
            clp->expr = NULL;
          }
          expr = delete_expression(expr);
        }
        clp->double_value = val;
      }
      else
      {
        clp->double_value = lp->call_def->double_value;
        clp->expr = NULL;
        tot_end = i;
      }
    }
    else if (lp->type == 4) { /* one constraint */
      if (i+1 < number && *toks[i+1] == ':')
        con_flag = 1, i++;
      else
        con_flag = 0; /* if != zero, := or :< or :> */

      if (i+2 < number) {
             if (*toks[i+1] == '=') c_type = 4;
        else if (*toks[i+1] == '>') c_type = 1;
        else if (*toks[i+1] == '<') c_type = 2;

        if (c_type) {
          start = i + 2;

          for (t_num = start; t_num < number; t_num++)
            if(*toks[t_num] == ',') break;

          if ((e_type = loc_expr(toks, t_num, start, &end)) == 0)
            return -i;

          tot_end = end;
          if (e_type == 1) { /* simple number */
            val = simple_double(toks, start, end);
            expr = NULL;
          }
          else { /* expression */
            if ((expr = make_expression(end + 1 - start, &toks[start])) == NULL)
              return -i;
            val = expression_value(expr, 2);
          }
        }
      }
      else {
        c_type = 4;
        val = zero;
        expr = NULL;
        tot_end = i;
      }
    }

         if (clp->c_type == 1 && c_type == 2) clp->c_type = 3; /* min */
    else if (clp->c_type == 2 && c_type == 1) clp->c_type = 3; /* max */
    else                                      clp->c_type = c_type;
    if (con_flag == 0) expr = NULL;

    switch(c_type)
    {
      case 1:
        clp->c_min = val;
        clp->min_expr = clone_expression(expr);
        break;
      case 2:
        clp->c_max = val;
        clp->max_expr = clone_expression(expr);
        break;
      case 4:
        clp->double_value = val;
        clp->expr = clone_expression(expr);
    }
  }
  else /* array */
  {
    if (lp->type == 13)  /* string array */
    {
      if (i+2 < number && *toks[i+1] == '=')
      {
        i += 2;
        get_bracket_t_range(toks, '{', '}', i, number-1, &start, &end);
        if (start >= i) /* {} found */
        {
          j = tot_end = end;
          start++;
        }
        else /* terminate with next keyword */
        {
          start = i;
          for (j = start; j < number; j++)
            if (name_list_pos(toks[j], cmd->cmd_def->par_names) > -1)
              break;
          tot_end = j - 1;
        }
        ks = start;
        for (k = start; k < j; k++)
        {
          if ((k+1 == j || *toks[k+1] == ','))
          {
            if (cnt == clp->m_string->max)
              grow_char_p_array(clp->m_string);
            clp->m_string->p[cnt++]
              = permbuff(noquote(join(&toks[ks], k+1-ks)));
            clp->m_string->curr = cnt;
            ks = k + 2;  k++;
          }
        }
      }
      else
      {
        clp->m_string->curr = 0;
        tot_end = i;
      }
    }
    else if ((i+2 < number && *toks[i+1] == '=') ||
             (i+3 < number && *toks[i+1] == ':' && *toks[i+2] == '='))
    {
      val_type = (*toks[i+1] == ':' && *toks[i+2] == '=') ? 1 : 0;
      i += val_type + 2;
      get_bracket_t_range(toks, '{', '}', i, number-1, &start, &end);
      if (start >= i) /* {} found */
      {
        j = tot_end = end;
        start++;
      }
      else /* terminate with next keyword */
      {
        start = i;
        for (j = start; j < number; j++)
          if (name_list_pos(toks[j], cmd->cmd_def->par_names) > -1)
            break;
        tot_end = j - 1;
      }
      if (clp->expr_list == NULL)
        clp->expr_list = new_expr_list(clp->double_array->max);
      while (start < j)
      {
        if ((end = next_char(',', toks, start, j)) < 0) end = j;
        e_type = loc_expr(toks, end, start, &e_end);
        if (e_end != end - 1) return -1;
        if (cnt == clp->double_array->max)
        {
          grow_double_array(clp->double_array);
          grow_expr_list(clp->expr_list);
        }
        if ((expr = clp->expr_list->list[cnt]) != NULL)
          delete_expression(expr);
        if ((expr =
             make_expression(end - start, &toks[start])) == NULL)
          return -i;
        if (val_type && e_type == 2)
        {
          clp->expr_list->list[cnt] = clone_expression(expr);
          val = zero;
        }
        else
        {
          clp->expr_list->list[cnt] = NULL;
          val = expression_value(expr, e_type);
        }
        expr = delete_expression(expr);
        clp->double_array->a[cnt++] = val;
        if (cnt > clp->double_array->curr)
          clp->double_array->curr = clp->expr_list->curr = cnt;
        start = end + 1;
      }
    }
    else return -i;
  }
  if (expr != NULL) delete_expression(expr);
  return tot_end;
}

const char*
alias(char* par_string) /* returns main parameter for alias */
{
  if (strcmp(par_string, "filename") == 0)  return file_string;
  else if(strcmp(par_string, "true_") == 0) return vrai;
  else if(strcmp(par_string, "false_") == 0) return faux;
  else return par_string;
}

int
log_val(const char* name, struct command* cmd)
  /* returns 0 = false, 1 = true for a logical command parameter */
{
  struct command_parameter* cp;
  return command_par(name, cmd, &cp) && cp->double_value != 0;
}

// public interface (used by Fortran)

double
get_value(const char* name, const char* par)
  /* this function is used by fortran to get the parameters values
     returns parameter value "par" for command or store "name" if present,
     else INVALID */
{
  int tmpi;
  double tmpd;

  mycpy(c_dum->c, name);
  mycpy(aux_buff->c, par);

  if (strcmp(c_dum->c, "beam") == 0)
    return command_par_value(aux_buff->c, current_beam);
  else if (strcmp(c_dum->c, "probe") == 0)
   {
    /*printf("Getting probe beam => %p\n", probe_beam);*/
    if (probe_beam == 0x0)
     {
       printf("\n\n get_value: PROBE IS NULL (name=%s, par=%s)!!!!!!\n\n\n",name,par);
       return 0.0;
     }
    return command_par_value(aux_buff->c, probe_beam);
   }
  else if (strcmp(c_dum->c, "survey") == 0)
  {
    return command_par_value(aux_buff->c, current_survey);
  }
  else if (strcmp(c_dum->c, "twiss") == 0)
  {
    return command_par_value(aux_buff->c, current_twiss);
  }
  else if (strcmp(c_dum->c, "sequence") == 0)
  {
    if (strcmp(aux_buff->c, "l") == 0) return sequence_length(current_sequ);
    else if (strcmp(aux_buff->c, "range_start") == 0)
      return (current_sequ->range_start->position
              - 0.5 * current_sequ->range_start->length);
    else if (strcmp(aux_buff->c, "add_pass") == 0)
      return ((double)current_sequ->add_pass);
    else return INVALID;
  }
  else
   {
   tmpi =  strcmp(c_dum->c, current_command->name);
   if (current_command != NULL &&  tmpi == 0)
     {
       tmpd = command_par_value(aux_buff->c, current_command);
       return tmpd;
     }
    else
     {
       return INVALID;
     }
   }
}

int
get_vector(const char* name, const char* par, double* vector)
  /* returns double "vector" for "par" of command or store "name";
     length is returned as function value (0 if not found) */
{
  mycpy(c_dum->c, name);
  mycpy(aux_buff->c, par);
  if (strcmp(c_dum->c, "threader") == 0)
    return command_par_vector(aux_buff->c, threader_par, vector);
  else return 0;
}

int
get_string(const char* name, const char* par, char* string)
  /* returns string for  value "par" of command or store "name" if present,
     length = string length, else length = 0 if not present */
{
  struct command* cmd;
  char* p;
  int length = 0;
  mycpy(c_dum->c, name);
  if (strcmp(c_dum->c, "beam") == 0)
  {
    mycpy(c_dum->c, par);
    if ((p = command_par_string(c_dum->c, current_beam)) != NULL)
    {
      strcpy(string, p); length = strlen(p);
    }
  }
  else if (strcmp(c_dum->c, "probe") == 0)
  {
    mycpy(c_dum->c, par);
    if ((p = command_par_string(c_dum->c, probe_beam)) != NULL)
    {
      strcpy(string, p); length = strlen(p);
    }
  }
  else if (strcmp(c_dum->c, "survey") == 0)
  {
    mycpy(c_dum->c, par);
      if ((p = command_par_string_user(c_dum->c, current_survey)) != NULL)
      {
        strcpy(string, p); length = strlen(p);
      }
  }
  /*else if (strcmp(c_dum->c, "ptc") == 0)
    {
    mycpy(c_dum->c, par);
    if (current_ptc != NULL) nl = current_ptc->par_names;
    if (nl != NULL )
    {
    if ((p = command_par_string(c_dum->c, current_ptc)) != NULL)
    {
    strcpy(string, p); length = strlen(p);
    }
    }
    }  */
  else if (strcmp(c_dum->c, "twiss") == 0)
  {
    mycpy(c_dum->c, par);
      if ((p = command_par_string_user(c_dum->c, current_twiss)) != NULL)
      {
        strcpy(string, p); length = strlen(p);
      }
  }
  else if (strcmp(c_dum->c, "sequence") == 0)
  {
    mycpy(c_dum->c, par);
    if (current_sequ != NULL && strcmp(c_dum->c, "name") == 0)
    {
      p = current_sequ->name;
      strcpy(string, p); length = strlen(p);
    }
  }
  else if (strcmp(c_dum->c, "element") == 0)
  {
    mycpy(c_dum->c, par);
    if (current_sequ != NULL && strcmp(c_dum->c, "name") == 0)
    {
      p = current_node->p_elem->name;
      strcpy(string, p); length = strlen(p);
    }
  }
  else
  {
/*     printf("<madxp.c: get_string>: Looking for command %s \n",c_dum->c);*/

    cmd = find_command(c_dum->c, stored_commands);
    if (cmd == NULL)
    {
      if (current_command != NULL)
      {
        if ( strcmp(c_dum->c, current_command->name) == 0)
        {
          cmd = current_command;
        }
      }
    }

    if (cmd != NULL)
    {
/*        printf("<madxp.c: get_string>: Found command %s \n",c_dum->c);*/
      mycpy(c_dum->c, par);
/*        printf("<madxp.c: get_string>: Looking for parameter %s \n",c_dum->c);*/
      if ((p = command_par_string(c_dum->c, cmd)) != NULL)
      {
        strcpy(string, p); length = strlen(p);
      }
      else
      {
        printf("<madxp.c: get_string>: Did not found parameter %s \n",c_dum->c);
      }
    }

    else
    {
      printf("<madxp.c: get_string>: Did not found command %s \n",c_dum->c);
    }
  }
  return length;
}

void
set_value(const char* name, const char* par, double* value)
  /* sets parameter value "par" for command or store "name" if present */
{
  mycpy(c_dum->c, name);
  mycpy(aux_buff->c, par);
  if (strcmp(c_dum->c, "beam") == 0)
    set_command_par_value(aux_buff->c, current_beam, *value);
  else if (strcmp(c_dum->c, "probe") == 0)
    set_command_par_value(aux_buff->c, probe_beam, *value);
  else if (strcmp(c_dum->c, "survey") == 0)
    set_command_par_value(aux_buff->c, current_survey, *value);
  else if (strcmp(c_dum->c, "twiss") == 0)
    set_command_par_value(aux_buff->c, current_twiss, *value);
  else if (current_command && strcmp(c_dum->c, current_command->name) == 0)
    set_command_par_value(aux_buff->c, current_command, *value);
}


