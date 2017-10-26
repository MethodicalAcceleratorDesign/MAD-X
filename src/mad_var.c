#include "madx.h"

static struct variable*
delete_variable(struct variable* var)
{
  const char *rout_name = "delete_variable";
  if (var == NULL)  return NULL;
  if (stamp_flag && var->stamp != 123456)
    fprintf(stamp_file, "d_v double delete --> %s\n", var->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", var->name);
  if (var->expr != NULL) delete_expression(var->expr);
  if (var->string != NULL) myfree(rout_name, var->string);
  myfree(rout_name, var);
  return NULL;
}

static void
grow_var_list(struct var_list* p)
{
  const char *rout_name = "grow_var_list";
  struct variable** v_loc = p->vars;
  int new = 2*p->max;

  p->max = new;
  p->vars = mycalloc(rout_name, new, sizeof *p->vars);
  for (int j = 0; j < p->curr; j++) p->vars[j] = v_loc[j];
  myfree(rout_name, v_loc);
}

static void
export_variable(struct variable* var, FILE* file, int noexpr)
  /* exports variable in mad-X format */
{
  if (var->status == 0) var->value = expression_value(var->expr, var->type);

  // LD: warning, the clear must be after the expression eval which uses c_dum->c
  //     ugly side effect...
  c_dum->c[0] = '\0';
  if (var->val_type == 0) strcat(c_dum->c, "int ");
  if (var->type == 0) strcat(c_dum->c, "const ");
  strcat(c_dum->c, var->name);

  if (var->type < 2 || noexpr) strcat(c_dum->c, " = ");
  else                         strcat(c_dum->c, " := ");

  if (var->expr != NULL && !noexpr) strcat(c_dum->c, var->expr->string);
  else if (var->val_type == 0)
  {
    int k = var->value;
    sprintf(c_join->c, "%d", k);
    strcat(c_dum->c, c_join->c);
  }
  else
  {
    sprintf(c_join->c, v_format("%F"), var->value);
    strcat(c_dum->c, supp_tb(c_join->c));
  }
  write_nice(c_dum->c, file);
}

static void
export_var_8(struct variable* var, FILE* file)
  /* exports variable in mad-8 format */
{
  int k;
  if (var->status == 0) var->value = expression_value(var->expr, var->type);

  // LD: warning, the clear must be after the expression eval which uses c_dum->c
  //     ugly side effect...
  c_dum->c[0] = '\0';
  if (var->type == 0)
  {
    strcat(c_dum->c, var->name);
    strcat(c_dum->c, ": constant = ");
  }
  else
  {
    strcat(c_dum->c, var->name);
    if (var->type < 2) strcat(c_dum->c, " = ");
    else               strcat(c_dum->c, " := ");
  }
  if (var->expr != NULL) strcat(c_dum->c, var->expr->string);
  else if (var->val_type == 0)
  {
    k = var->value; sprintf(c_join->c, v_format("%I"), k);
    strcat(c_dum->c, c_join->c);
  }
  else
  {
    sprintf(c_join->c, v_format("%F"), var->value);
    strcat(c_dum->c, supp_tb(c_join->c));
  }
  write_nice_8(c_dum->c, file);
}

static int
predef_var(struct variable* var)
  /* return 1 for predefined variable, else 0 */
{
  int pos = name_list_pos(var->name, variable_list->list);
  return (pos < start_var ? 1 : 0) ;
}

static void
set_sub_variable(char* comm, char* par, struct in_cmd* cmd)
{
  char* p;
  struct element* el;
  struct command *command, *keep_beam = current_beam;
  int end, start = cmd->decl_start, t_num, exp_type;
  double val = 0;

  for (t_num = start; t_num < cmd->tok_list->curr; t_num++)
    if (*(cmd->tok_list->p[t_num]) == ',') break;

  exp_type = loc_expr(cmd->tok_list->p, t_num, start, &end);

  if (exp_type == 1) /* literal constant */
    val = simple_double(cmd->tok_list->p, start, end);
  else if (polish_expr(end + 1 - start, &cmd->tok_list->p[start]) == 0)
    val = polish_value(deco, join(&cmd->tok_list->p[start], end + 1 - start));

  if (strncmp(comm, "beam", 4) == 0) {
    command = current_beam = find_command("default_beam", beam_list);
    if ((p = strchr(comm, '%')) != NULL) {
      if ((current_beam = find_command(++p, beam_list)) == NULL)
        current_beam = command;
    }
    set_command_par_value(par, current_beam, val);
  }
  else if ((el = find_element(comm, element_list)) != NULL)
    set_command_par_value(par, el->def, val);
  else if ((command = find_command(comm, stored_commands)) != NULL)
    set_command_par_value(par, command, val);
  else if ((command = find_command(comm, beta0_list)) != NULL)
    set_command_par_value(par, command, val);
  else if ((command = find_command(comm, defined_commands)) != NULL)
    set_command_par_value(par, command, val);
  current_beam = keep_beam;
}

// public interface

void
get_defined_constants(void)
{
  /* reads + stores the constants defined in madxdict.h */
  supp_char('\n', constant_def);
  pro_input(constant_def);
  start_var = variable_list->curr;
}

struct var_list*
delete_var_list(struct var_list* varl)
{
  const char *rout_name = "delete_var_list";
  if (varl == NULL) return NULL;
  if (stamp_flag && varl->stamp != 123456)
    fprintf(stamp_file, "d_v_l double delete --> %s\n", varl->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", varl->name);
  if (varl->list != NULL) delete_name_list(varl->list);
  if (varl->vars != NULL) myfree(rout_name, varl->vars);
  myfree(rout_name, varl);
  return NULL;
}

struct variable*
find_variable(const char* name, struct var_list* varl)
{
  int pos;
  if ((pos = name_list_pos(name, varl->list)) < 0)
    return NULL;
  return varl->vars[pos];
}

double
variable_value(struct variable* var)
{
  int k;
  double val = zero;
  if (var->type < 2 && var->status > 0) val = var->value;
  else if(var->expr == NULL) val = var->value;
  else
  {
    var->value = val = expression_value(var->expr, var->type);
    var->status = 1;
  }
  if (var->val_type == 0) /* int */
  {
    k = val; val = k;
  }
  return val;
}

struct var_list*
clone_var_list(struct var_list* vl)
{
  int i, l = vl->curr > 0 ? vl->curr : 1;
  struct var_list* clone;
  clone = new_var_list(l);
  strcpy(clone->name, vl->name);
  clone->list = clone_name_list(vl->list);
  for (i = 0; i < vl->curr; i++) clone->vars[i] = vl->vars[i];
  clone->curr = vl->curr;
  return clone;
}

struct variable*
new_variable(const char* name, double val, int val_type, int type, struct expression* expr, char* string)
{
  const char *rout_name = "new_variable";
  struct variable* var = mycalloc(rout_name, 1, sizeof *var);
  strcpy(var->name, name);
  var->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", var->name);
  var->value = val;
  var->type = type;
  var->val_type = val_type;
  if ((var->expr = expr) == NULL) var->status = 1;
  if (string != NULL) var->string = tmpbuff(string);
  return var;
}

struct var_list*
new_var_list(int length)
{
  const char *rout_name = "new_var_list";
  struct var_list* var = mycalloc(rout_name, 1, sizeof *var);
  strcpy(var->name, "var_list");
  var->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", var->name);
  var->list = new_name_list(var->name, length);
  var->vars = mycalloc(rout_name, length, sizeof *var->vars);
  var->max = length;
  return var;
}

char*
get_varstring(const char* name)
{
  struct variable* var;
  char *ret; // *p, not used
  ret = NULL;
  mycpy(c_dum->c, name);
  if (strstr(c_dum->c, "->") == NULL) /* variable */ // (p = not used
    if ((var = find_variable(c_dum->c, variable_list)) != NULL)
      ret = var->string;

  return ret;
}

void
write_vars(struct var_list* varl, struct command_list* cl, FILE* file, int noexpr)
{
  for (int i = 0; i < varl->curr; i++) {
    if (predef_var(varl->vars[i]) == 0 && pass_select_list_str(varl->vars[i]->name, cl))
      export_variable(varl->vars[i], file, noexpr);
  }
}

void
write_vars_8(struct var_list* varl, struct command_list* cl, FILE* file)
{
  int i;
  for (i = 0; i < varl->curr; i++)
  {
    if (predef_var(varl->vars[i]) == 0
        && pass_select_list_str(varl->vars[i]->name, cl))
      export_var_8(varl->vars[i], file);
  }
}

char*
make_string_variable(char* string)
  /* creates + stores a variable containing a character string */
{
  char* name = get_new_name();
  struct variable* var = new_variable(name, zero, 3, 0, NULL, string);
  add_to_var_list(var, variable_list, 0);
  return var->name;
}

void
enter_variable(struct in_cmd* cmd) /* stores variable contained in cmd */
{
  struct variable* var;
  struct expression* expr = NULL;
  int k, end, type = 0, name_pos = 0, start = cmd->decl_start, val_type = 0;
  double val = 0;
  char comm[NAME_L];
  char par[NAME_L];
  char* name;
  char *p, *n, *q = comm;
  int exp_type = loc_expr(cmd->tok_list->p, cmd->tok_list->curr,
                          start, &end);
  switch (cmd->sub_type)
  {
    case 2:
    case 3:
      val_type = 0;
      type = 0;
      name_pos = 2;
      break;
    case 4:
    case 5:
      val_type = 1;
      type = 0;
      name_pos = 2;
      break;
    case 6:
      val_type = 0;
      type = 1;
      name_pos = 1;
      break;
    case 7:
      val_type = 0;
      type = 2;
      name_pos = 1;
      break;
    case 8:
      val_type = 1;
      type = 1;
      name_pos = 1;
      break;
    case 9:
      val_type = 1;
      type = 2;
      name_pos = 1;
      break;
    case 10:
    case 11:
      val_type = 1;
      type = 0;
      name_pos = 1;
      break;
    case 12:
      val_type = 1;
      type = 1;
      name_pos = 0;
      break;
    case 13:
      val_type = 1;
      type = 2;
      name_pos = 0;
      break;
    default:
      fatal_error("illegal command sub_type in:",
                  join(cmd->tok_list->p, cmd->tok_list->curr));
  }
  n = name = permbuff(cmd->tok_list->p[name_pos]);
  if (exp_type == 0)
  {
    warning("illegal expression set to 0 in:",
            join_b(cmd->tok_list->p, cmd->tok_list->curr));
  }
  else
  {
    if ((p = strstr(name,"->")) != NULL) /* element or command parameter */
    {
      while (n < p)  *(q++) = *(n++);
      *q = '\0';
      q = par; n++; n++;
      while (*n != '\0')  *(q++) = *(n++);
      *q = '\0';
      set_sub_variable(comm, par, cmd);
    }
    else if ((var = find_variable(name, variable_list)) != 0
             && var->type == 0)
    {
      warning("ignored: attempt to redefine constant:", var->name);
    }
    else if (exp_type == 1) /* literal constant */
    {
      val = simple_double(cmd->tok_list->p, start, end);
      var = new_variable(name, val, val_type, type, NULL, NULL);
      add_to_var_list(var, variable_list, 1);
    }
    else
    {
      if (polish_expr(end + 1 - start, &cmd->tok_list->p[start]) == 0)
      {
        if (type == 2) /* deferred: expression kept */
        {
          expr = new_expression(join(&cmd->tok_list->p[start], end + 1 - start), deco);
          val = 0; // LD 2012.10.16: drop warning due to expression_value(expr, type);
        }
        else
        {
          expr = NULL;
          val = polish_value(deco, join(&cmd->tok_list->p[start], end + 1 - start));
        }
        if (val_type == 0)
        {
          if (fabs(val) < 2.e9)  k = val;
          else                   k = 0;
          val = k;
        }
        var = new_variable(name, val, val_type, type, expr, NULL);
        add_to_var_list(var, variable_list, 1);
      }
      else
      {
        warning("illegal expression set to 0 in:",
                join_b(cmd->tok_list->p, cmd->tok_list->curr));
      }
    }
  }
}

void
add_to_var_list( /* adds variable to alphabetic variable list */
  struct variable* var, struct var_list* varl, int flag)
  /* flag = 0: undefined reference in expression, 1: definition
     2: separate list, do not drop variable */
{
  int pos; // , j; not used

  if ((pos = name_list_pos(var->name, varl->list)) > -1)
  {
    if (flag == 1)
    {
      if (varl->list->inform[pos] == 1)
      {
        put_info(var->name, "redefined");
/*
  printf("Old Value:\n");
  export_variable(varl->vars[pos], stdout);
  printf("New Value:\n");
  export_variable(var, stdout);
  printf("File %s line %d\n",filenames[in->curr], currentline[in->curr] );
*/
      }
      else varl->list->inform[pos] = flag;
    }
    if (flag < 2) delete_variable(varl->vars[pos]);
    varl->vars[pos] = var;
  }
  else
  {
    if (varl->curr == varl->max) grow_var_list(varl);
    add_to_name_list(permbuff(var->name), flag, varl->list); // j = not used
    varl->vars[varl->curr++] = var;
  }
}

void
set_stringvar(const char* name, char* string)
{
  /* sets variable name->string to string */
//  char* p;
  struct variable* var;
  mycpy(c_dum->c, name);
  if (strstr(c_dum->c, "->") == NULL) /* variable */ // (p = not used
  {
    if ((var = find_variable(c_dum->c, variable_list)) != NULL)
    {
      if (var->type == 3) var->string = string;
    }
    else
    {
      var = new_variable(c_dum->c, zero, 0, 3, NULL, string);
      add_to_var_list(var, variable_list, 1);
    }
  }
}

// public interface (used by Fortran)

double
get_variable(const char* name)
{
  char comm[NAME_L];
  char par[NAME_L];
  double val = zero;
  struct variable* var;
  struct element* el;
  struct command* cmd;
  char *p, *n = c_dum->c, *q = comm;
  mycpy(c_dum->c, name);
  if ((p = strstr(c_dum->c, "->")) == NULL) /* variable */
  {
    if ((var = find_variable(c_dum->c, variable_list)) != NULL)
      val = variable_value(var);
  }
  else /* element or command parameter */
  {
    while (n < p)  *(q++) = *(n++);
    *q = '\0';
    q = par; n++; n++;
    while (*n != '\0')  *(q++) = *(n++);
    *q = '\0';
    if (((el = find_element(comm, element_list)) && (cmd = el->def))
            || (cmd = find_command(comm, stored_commands))
            || (cmd = find_command(comm, beta0_list))
            || (cmd = find_command(comm, defined_commands))) {
      val = command_par_value(par, cmd);
    }
  }
  return val;
}

void
set_variable(const char* name, double* value)
{
  /* sets variable name to value */
  char comm[NAME_L];
  char par[NAME_L];
  struct variable* var;
  double val = *value;
  struct element* el;
  struct command* cmd;
  char *p, *n = c_dum->c, *q = comm;
  mycpy(c_dum->c, name);
  if ((p = strstr(c_dum->c, "->")) == NULL) /* variable */
  {
    if ((var = find_variable(c_dum->c, variable_list)) != NULL)
    {
      if (var->type == 0)
        warning("ignored: attempt to redefine constant:", var->name);
      else if (var->type < 3)
      {
        var->value = val;
        var->type = 1;
        if (var->expr != NULL)  var->expr = delete_expression(var->expr);
      }
    }
    else
    {
      var = new_variable(c_dum->c, val, 1, 1, NULL, NULL);
      add_to_var_list(var, variable_list, 1);
    }
  }
  else /* element or command parameter */
  {
    while (n < p)  *(q++) = *(n++);
    *q = '\0';
    q = par; n++; n++;
    while (*n != '\0')  *(q++) = *(n++);
    *q = '\0';
    if (((el = find_element(comm, element_list)) && (cmd = el->def))
            || (cmd = find_command(comm, stored_commands))
            || (cmd = find_command(comm, beta0_list))
            || (cmd = find_command(comm, defined_commands))) {
      set_command_par_value(par, cmd, val);
    }
  }
}

