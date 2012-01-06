#include "madx.h"

static void
exec_delete_sequ(char* name)
{
  struct sequence* keep = current_sequ;
  int spos;
  if ((spos = name_list_pos(name, sequences->list)) >= 0)
  {
    current_sequ = sequences->sequs[spos];
    if (current_sequ->ex_start != NULL) /* delete expanded */
    {
      current_sequ->ex_nodes = delete_node_list(current_sequ->ex_nodes);
      current_sequ->ex_start = delete_node_ring(current_sequ->ex_start);
      current_sequ->orbits = delete_vector_list(current_sequ->orbits);
    }
    sequences->sequs[spos] = delete_sequence(current_sequ);
    remove_from_sequ_list(current_sequ, sequences);
    current_sequ = keep;  
  }
  else warning("sequence to be deleted does not exist:", name);
}

static void
exec_delete_table(char* name)
{
  struct table_list* tl;
  int j, k, pos;
  for (j = 0; j < all_table_lists->curr; j++)
  {
    tl = all_table_lists->table_lists[j];
    if ((pos = name_list_pos(name, tl->names)) >= 0)
    {
      tl->tables[pos] = delete_table(tl->tables[pos]);
      k = remove_from_name_list(name, tl->names);
      tl->tables[k] = tl->tables[--tl->curr];
      return;
    }
  }
}

// public interface

void
exec_option(void)
{
  if (get_option("reset")) set_defaults("option");
  if (get_option("tell")) print_command(options);

}

void
exec_help(struct in_cmd* cmd)
  /* prints list of commands */
{
  char** toks = cmd->tok_list->p;
  int i, k = 0, pos, n = cmd->tok_list->curr;
  if (n == 1)
  {
    while (special_comm_cnt[k] > 0) k++;
    puts("special commands - no further help:");
    puts(" ");
    for (i = 0; i < k-1; i++)
    {
      if (strchr(special_comm_desc[i], '(') != NULL)
        fprintf(prt_file, "%s<condition>){<statements(s)>}\n",
                &special_comm_desc[i][0]);
      else if (strchr(special_comm_desc[i], '{') != NULL)
        fprintf(prt_file, "%s<statements(s)>}\n",
                &special_comm_desc[i][0]);
      else fprintf(prt_file, "%s{<statements(s)>}\n",
                   &special_comm_desc[i][0]);
    }
    fprintf(prt_file, "<name>:line(...);\n");
    puts(" ");
    puts("normal commands or predefined particles:");
    dump_name_list(defined_commands->list);
  }
  else
  {
    for (i = 1; i < n; i++)
    {
      if ((pos = name_list_pos(toks[i], defined_commands->list)) > -1)
        dump_command(defined_commands->commands[pos]);
      else puts("no help for this command - try help; (no arguments)");
    }
  }
}

void
exec_assign(struct in_cmd* cmd)
  /* executes output unit assignment */
{
  char* p;
  char tmp[FNAME_L];
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int pos = name_list_pos("echo", nl);
  if (prt_file != stdout)  fclose(prt_file);
  if (nl->inform[pos])
  {
    p = pl->parameters[pos]->string; strcpy(tmp, p);
    if (strcmp(stolower(tmp), "terminal") == 0)  prt_file = stdout;
    else
    {
      if (assign_start == 0)
      {
        assign_start = 1; prt_file = fopen(p, "w");
      }
      else prt_file = fopen(p, "a");
    }
  }
  else prt_file = stdout;
}

void
exec_call(struct in_cmd* cmd)
  /* handles calling external files */
{
  struct command_parameter_list* pl = cmd->clone->par;
  struct name_list* nl = cmd->clone->par_names;
  int pos = name_list_pos("file", nl);
  int top = in->curr;
  if (nl->inform[pos])
  {
    if (down_unit(pl->parameters[pos]->string)) madx_input(top);
  }
  else warning("call without filename:", "ignored");
}

void
exec_cmd_delete(struct in_cmd* cmd)
/* handles all delete request through "delete" command */
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int pos;
  char* name;
  pos = name_list_pos("sequence", nl);
  if (nl->inform[pos])
  {
    name = pl->parameters[pos]->string;
    exec_delete_sequ(name);
  }
  pos = name_list_pos("table", nl);
  if (nl->inform[pos])
  {
    name = pl->parameters[pos]->string;
    exec_delete_table(name);
  }
}

void
exec_show(struct in_cmd* cmd)
  /* executes "show" command */
{
  struct element* el;
  struct variable* var;
  char** toks = cmd->tok_list->p;
  int i, pos, n = cmd->tok_list->curr;
  for (i = 1; i < n; i++)
  {
    if (strcmp(toks[i],","))
    {
      if (strncmp(toks[i], "beam", 4) == 0) show_beam(toks[i]);
      else if ((pos = name_list_pos(toks[i], defined_commands->list)) > -1)
      {
        if (strcmp(toks[i], "option") == 0) dump_command(options);
        else if (strcmp(toks[i], "eoption") == 0 && current_eopt != NULL)
          dump_command(current_eopt);
        else dump_command(defined_commands->commands[pos]);
      }
      else if ((pos = name_list_pos(toks[i], beta0_list->list)) > -1)
        dump_command(beta0_list->commands[pos]);
      else if ((el = find_element(toks[i], element_list)) != NULL)
        dump_element(el);
      else if ((var = find_variable(toks[i], variable_list)))
      {
        if (var->expr)  fprintf(prt_file, "%s := %s ;\n", toks[i], var->expr->string);
        else fprintf(prt_file, v_format("%s = %F ;\n"), toks[i], var->value);
      }
      else fprintf(prt_file, "%s not found\n;", toks[i]);
    }
  }
}

void
exec_create_table(struct in_cmd* cmd)
  /* makes a user defined table */
{
  char rout_name[] = "exec_create_table";
  struct table* t;
  int* t_types;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  struct char_p_array* m;
  char** t_c;
  int j, pos = name_list_pos("table", nl);
  char* name = NULL;
  int  ncols = 0;  /*number of columns*/

  if (nl->inform[pos] == 0)
  {
    warning("no table name:", "ignored");
    return;
  }
  if ((name = pl->parameters[pos]->string) == NULL)
  {
    warning("no table name: ", "ignored");
    return;
  }
  if ((pos = name_list_pos(name, table_register->names)) > -1)
  {
    warning("table already exists: ", "ignored");
    return;
  }

  pos = name_list_pos("column", nl);
  if (nl->inform[pos] == 0)
  {
    warning("table without columns: ", "ignored");
    return;
  }
  m = pl->parameters[pos]->m_string;
  ncols = m->curr;
  /* now make table */
  t_types = mymalloc(rout_name, ncols*sizeof(int));
  t_c = mymalloc(rout_name, (ncols+1)*sizeof(char*));

  for (j = 0; j < m->curr; j++)
  {
    if (*m->p[j] == '_')
    {
      t_types[j] = 3; /* type string */
      t_c[j] = permbuff(&m->p[j][1]);
    }
    else
    {
      t_types[j] = 2; /* type double */
      t_c[j] = permbuff(m->p[j]);
    }
  }
  t_c[ncols] = blank;
  t = make_table(name, "user", t_c, t_types, USER_TABLE_LENGTH);
  t->org_cols = 0;  /* all entries are "added" */
  add_to_table_list(t, table_register);
  myfree(rout_name, t_c); myfree(rout_name, t_types);
  t->dynamic = 1;
}

void
exec_store_coguess(struct in_cmd* cmd)
  /* stores the initial orbit guess of the user */
{
  struct name_list* nl = cmd->clone->par_names;
  int pos = name_list_pos("tolerance", nl);
  double tol;
  if (nl->inform[pos])
  {
    tol = command_par_value("tolerance", cmd->clone);
    set_variable("twiss_tol", &tol);
  }
  store_orbit(cmd->clone, guess_orbit);
  guess_flag = 1;
}

void
exec_dump(struct in_cmd* cmd)
  /* write a table out */
{
  struct table* t;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int pos = name_list_pos("table", nl);
  char* name = NULL;
  char *f, filename[FNAME_L];
  if (nl->inform[pos] == 0)
  {
    warning("dump without table name:", "ignored");
    return;
  }
  if ((name = pl->parameters[pos]->string) == NULL)
  {
    warning("dump without table name:", "ignored");
    return;
  }
  pos = name_list_pos("file", nl);
  if (nl->inform[pos] == 0) strcpy(filename, "terminal");
  else if ((f = pl->parameters[pos]->string) == NULL
           || *f == '\0') strcpy(filename, name);
  else strcpy(filename,f);
  if ((pos = name_list_pos(name, table_register->names)) > -1)
  {
    t = table_register->tables[pos];
    out_table(name, t, filename);
  }
  else
  {
    warning("table name not found:", "ignored");
  }
}

void
exec_fill_table(struct in_cmd* cmd)
  /* adds variables to a table */
{
  struct table* t;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int pos = name_list_pos("table", nl);
  char* name = NULL;
  int row,curr;
  if (nl->inform[pos] == 0)
  {
    warning("no table name:", "ignored");
    return;
  }
  if ((name = pl->parameters[pos]->string) == NULL)
  {
    warning("no table name: ", "ignored");
    return;
  }
  pos=name_list_pos("row", nl);
  row=(int) pl->parameters[pos]->double_value;
  if ((pos = name_list_pos(name, table_register->names)) > -1)
  {
    t = table_register->tables[pos];
    if (row<0) {
      add_vars_to_table(t);
      if (++t->curr == t->max) grow_table(t);
    } else {
      row--;
      curr=t->curr;
      if (row < t->curr) {
        t->curr=row;}
      else {
        t->curr--;
      }
      add_vars_to_table(t);
      t->curr=curr;
    }
  }
  else warning("table not found: ", "ignored");
  return;
}

void
exec_setvars_table(struct in_cmd* cmd)
  /* set variables from a table */
{
  struct table* t;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int pos = name_list_pos("table", nl);
  char* name = NULL;
  int row,curr;
  if (nl->inform[pos] == 0)
  {
    warning("no table name:", "ignored");
    return;
  }
  if ((name = pl->parameters[pos]->string) == NULL)
  {
    warning("no table name: ", "ignored");
    return;
  }
  current_node = NULL; /* to distinguish from other table fills */
  pos=name_list_pos("row", nl);
  row=(int) pl->parameters[pos]->double_value;
  if ((pos = name_list_pos(name, table_register->names)) > -1)
  {
    t = table_register->tables[pos];
    row--;
    curr=t->curr;
    if ((row < t->curr) && (row >-1)) {
      t->curr=row;}
    else {
      t->curr--;
    }
    set_vars_from_table(t);
    t->curr=curr;
  }
  else warning("table not found: ", "ignored");
  return;
}

void
exec_print(struct in_cmd* cmd)
  /* prints text from "print" command to current output unit */
{
  struct command_parameter_list* pl = cmd->clone->par;
  struct name_list* nl = cmd->clone->par_names;
  int pos = name_list_pos("text", nl);
  if (nl->inform[pos]) fprintf(prt_file,"%s\n", pl->parameters[pos]->string);
}

