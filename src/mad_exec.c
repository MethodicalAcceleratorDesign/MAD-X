#include "madx.h"

static void
exec_delete_sequ(char* name)
{
  struct sequence* keep = current_sequ;
  int spos;
  if ((spos = name_list_pos(name, sequences->list)) >= 0) {
    current_sequ = sequences->sequs[spos];
    if (current_sequ->ex_start != NULL) { /* delete expanded */
      current_sequ->ex_nodes = delete_node_list(current_sequ->ex_nodes);
      current_sequ->ex_start = delete_node_ring(current_sequ->ex_start);
      current_sequ->orbits = delete_vector_list(current_sequ->orbits);
    }
    sequences->sequs[spos] = delete_sequence(current_sequ);
    remove_from_sequ_list(current_sequ, sequences);
    current_sequ = keep;
  }
  else warning("sequence to be deleted does not exist:", name);
  return;
}

void
exec_delete_table(const char* name)
{
  struct table_list* tl;
  int j, k, pos;
  for (j = 0; j < all_table_lists->curr; j++) {
    tl = all_table_lists->table_lists[j];
    if ((pos = name_list_pos(name, tl->names)) >= 0) {
      tl->tables[pos] = delete_table(tl->tables[pos]);
      k = remove_from_name_list(name, tl->names);
      tl->tables[k] = tl->tables[--tl->curr];
      return;
    }
  }
}

// public interface

void
exec_option(struct in_cmd* cmd)
{
  if (get_option("reset")) set_defaults("option");
  if (get_option("tell"))  print_command(options);

  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int pos = name_list_pos("rand", nl);
  if (nl->inform[pos]) {
    const char *kind = pl->parameters[pos]->string;
    pos = name_list_pos("randid", nl);
    int rng_id = pl->parameters[pos]->double_value;
    setrand(kind, rng_id);
  }
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
  int cut = name_list_pos("truncate", nl);

  if (prt_file != stdout)  fclose(prt_file);
  if (nl->inform[pos]) {
    p = pl->parameters[pos]->string; strcpy(tmp, p);
    if (strcmp(stolower(tmp), "terminal") == 0)
      prt_file = stdout;
    else {
      p = str2path(p);
      if (assign_start == 0) {
        assign_start = 1;
        prt_file = fopen(p, "w");
      }
      else if (!nl->inform[cut] || !pl->parameters[cut]->double_value)
        prt_file = fopen(p, "a");
      else
        prt_file = fopen(p, "w");

      if (!prt_file) {
        warning("unable to open assigned file: ", p);
        prt_file = stdout;
      }
    }
  }
  else prt_file = stdout;
}

void
exec_removefile(struct in_cmd* cmd)
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int pos = name_list_pos("file", nl);

  if (nl->inform[pos]) {
    char *src = str2path(pl->parameters[pos]->string);
    if (remove(src))
      warning("unable to remove file: ", pl->parameters[pos]->string);
  }
}

void
exec_renamefile(struct in_cmd* cmd)
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int pos = name_list_pos("file", nl);
  int new = name_list_pos("to", nl);

  if (nl->inform[pos] && nl->inform[new]) {
    char *src = str2path(pl->parameters[pos]->string);
    char *dst = str2path(pl->parameters[new]->string);
    if (rename(src,dst)) warning("unable to rename file: ", src);
  }
}

void
exec_copyfile(struct in_cmd* cmd)
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int pos = name_list_pos("file", nl);
  int new = name_list_pos("to", nl);
  int flg = name_list_pos("append", nl);

  if (nl->inform[pos] && nl->inform[new]) {
    char *src_s = str2path(pl->parameters[pos]->string);
    char *dst_s = str2path(pl->parameters[new]->string);

    FILE *src = fopen(src_s, "r");
    if (!src) {
      warning("unable to open in read mode file: ", src_s);
      return;
    }

    const char *mode = "w";
    if (nl->inform[flg] && pl->parameters[flg]->double_value)
      mode = "a";

    FILE *dst = fopen(dst_s, mode);
    if (!dst) {
      warning("unable to open in write mode file: ", dst_s);
      fclose(src);
      return;
    }

    int c;
    while ((c = fgetc(src)) != EOF) fputc(c, dst);

    if (!feof(src))
      warning("unable to copy entirely file: ", src_s);

    fclose(src);
    fclose(dst);
  }
}

void
exec_call(struct in_cmd* cmd)
  /* handles calling external files */
{
  struct command_parameter_list* pl = cmd->clone->par;
  struct name_list* nl = cmd->clone->par_names;
  int pos = name_list_pos("file", nl);
  int top = in->curr;

  if (nl->inform[pos]) {
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
  if (nl->inform[pos]) {
    name = pl->parameters[pos]->string;
    exec_delete_sequ(name);
  }

  pos = name_list_pos("table", nl);
  if (nl->inform[pos]) {
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

  for (i = 1; i < n; i++) {
    if (strcmp(toks[i],",")) {
      if (strncmp(toks[i], "beam", 4) == 0) show_beam(toks[i]);
      else if ((pos = name_list_pos(toks[i], defined_commands->list)) > -1) {
        if (strcmp(toks[i], "option") == 0) dump_command(options);
        else if (strcmp(toks[i], "eoption") == 0 && current_eopt != NULL)
          dump_command(current_eopt);
        else dump_command(defined_commands->commands[pos]);
      }
      else if ((pos = name_list_pos(toks[i], beta0_list->list)) > -1)
        dump_command(beta0_list->commands[pos]);
      else if ((el = find_element(toks[i], element_list)) != NULL)
        dump_element(el);
      else if ((var = find_variable(toks[i], variable_list))) {
        if (var->expr)
          fprintf(prt_file, v_format("%S := %S ;\n"), toks[i], var->expr->string);
        else
          fprintf(prt_file, v_format("%S  = %F ;\n"), toks[i], var->value);
      }
      else fprintf(prt_file, "%s not found\n", toks[i]);
    }
  }
  return;
}

void
exec_create_table(struct in_cmd* cmd)
  /* makes a user defined table */
{
  const char *rout_name = "exec_create_table";
  struct table* t;
  int* t_types;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  struct char_p_array* m;
  const char**t_c;
  int j, pos = name_list_pos("table", nl);
  char* name = NULL;
  int  ncols = 0;  /*number of columns*/

  if (nl->inform[pos] == 0) {
    warning("no table name:", "ignored");
    return;
  }

  if ((name = pl->parameters[pos]->string) == NULL) {
    warning("no table name: ", "ignored");
    return;
  }

  if ((pos = name_list_pos(name, table_register->names)) > -1) {
    warning("table already exists: ", "ignored");
    return;
  }

  pos = name_list_pos("column", nl);
  if (nl->inform[pos] == 0) {
    warning("table without columns: ", "ignored");
    return;
  }

  m = pl->parameters[pos]->m_string;
  ncols = m->curr;

  /* now make table */
  t_types = mymalloc_atomic(rout_name, ncols * sizeof *t_types);
  t_c = mymalloc(rout_name, (ncols+1) * sizeof *t_c);

  for (j = 0; j < m->curr; j++) {
    if (*m->p[j] == '_') {
      t_types[j] = 3; /* type string */
      t_c[j] = permbuff(&m->p[j][1]);
    }
    else {
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

  return;
}

void
exec_store_coguess(struct in_cmd* cmd)
  /* stores the initial orbit guess of the user */
{
  struct name_list* nl = cmd->clone->par_names;
  double tol, toldefault=1.e-6;

  int pos = name_list_pos("tolerance", nl);
  if (nl->inform[pos])  {
    tol = command_par_value("tolerance", cmd->clone);
    set_variable("twiss_tol", &tol);
  }
  store_orbit(cmd->clone, guess_orbit);
  guess_flag = 1;

  /* 2014-May-30  13:55:50  ghislain: clear option added to cancel coguess */
  if (log_val("clear", cmd->clone)) {
    set_variable("twiss_tol",&toldefault);
    zero_double(guess_orbit, 6);
    guess_flag=0;
  }

  return;
}

void
exec_dump(struct in_cmd* cmd)
  /* write a table out */
{
  struct command_parameter_list* pl = cmd->clone->par;
  struct name_list* nl = cmd->clone->par_names;
  char *f, filename[FNAME_L], *name = NULL;
  int pos;

  // get "table" command parameter
  if ((pos = name_list_pos("table", nl) < 0) || nl->inform[pos] == 0 ||
      (name = pl->parameters[pos]->string) == NULL) {
    warning("dump without table name:", "ignored");
    return;
  }

  // get "file" command parameter
  ;
  if ((pos = name_list_pos("file", nl)) < 0 || nl->inform[pos] == 0)
    strcpy(filename, "terminal"); // write to console
  else if ((f = pl->parameters[pos]->string) == NULL || *f == '\0')
    strcpy(filename, name); // write to file with same name as table
  else
    strcpy(filename,f);

  // get table from registered tables
  if ((pos = name_list_pos(name, table_register->names)) < 0) {
    warning("table not found:", "ignored");
    return;
  }

  struct table* t = table_register->tables[pos];
  out_table(name, t, filename);

  return;
}

void
exec_shrink_table(struct in_cmd* cmd)
  /* removes rows from a table */
{
  struct table* t;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int pos = name_list_pos("table", nl);
  char* name = NULL;
  int row;

  if (nl->inform[pos] == 0) {
    warning("no table name:", "ignored");
    return;
  }

  if ((name = pl->parameters[pos]->string) == NULL) {
    warning("no table name: ", "ignored");
    return;
  }

  if ((pos = name_list_pos(name, table_register->names)) < 0) {
    warning("table name not found:", "ignored");
    return;
  }

  t = table_register->tables[pos];
  pos = name_list_pos("row", nl);
  row = pos >= 0 ? pl->parameters[pos]->double_value : t->curr - 1;

  if (row < 0) row = t->curr + row;
  if (row < 0 || row > t->curr) {
    warning("row index out of bounds:", " ignored");
    return;
  }
  t->curr = row;
}

void
exec_fill_table(struct in_cmd* cmd)
  /* adds variables to a table */
{
  struct command_parameter_list* pl = cmd->clone->par;
  struct name_list* nl = cmd->clone->par_names;
  char* name = NULL;
  int pos, row;
  double scale;

  if ((pos = name_list_pos("table", nl)) < 0 || nl->inform[pos] == 0 ||
      (name = pl->parameters[pos]->string) == NULL) {
    warning("no table name:", "ignored");
    return;
  }

  if ((pos = name_list_pos(name, table_register->names)) < 0) {
    warning("table not found:", "ignored");
    return;
  }

  struct table* t = table_register->tables[pos];

  if ((pos = name_list_pos("row", nl)) < 0)
    row = t->curr+1;
  else {
    row = pl->parameters[pos]->double_value;
    if (row < 1) row = t->curr+row+1;   // reflect negative index
    if (row < 1 || row > t->curr+1) {   // bounds check
      warning("row index out of bounds:", " ignored");
      return;
    }
  }
  pos =  name_list_pos("scale", nl);
  scale =  pl->parameters[pos]->double_value;

  int cols = t->org_cols, curr = t->curr;
  t->org_cols = 0;    t->curr = row - 1;
  add_vars_to_table(t,scale);
  t->org_cols = cols; t->curr = curr;

//   printf("table fill: %s [%d/%d]\n", name, t->curr, t->max);

  if (row == t->curr+1) // enlarge if needed
    if (++t->curr == t->max) grow_table(t);
}

void
exec_fill_knob_table(struct in_cmd* cmd)
  /* add knob to variables with weights from a table*/
{
  struct command_parameter_list* pl = cmd->clone->par;
  struct name_list* nl = cmd->clone->par_names;
  const char* name = NULL, *knob = NULL;
  int pos, row;
  double scale;

  if ((pos = name_list_pos("table", nl)) < 0 || nl->inform[pos] == 0 ||
      (name = pl->parameters[pos]->string) == NULL) {
    warning("no table name:", "ignored");
    return;
  }

  if ((pos = name_list_pos(name, table_register->names)) < 0) {
    warning("table not found:", "ignored");
    return;
  }

  struct table* t = table_register->tables[pos];

  if ((pos = name_list_pos("row", nl)) < 0)
    row = t->curr+1;
  else {
    row = pl->parameters[pos]->double_value;
    if (row < 1) row = t->curr+row+1;   // reflect negative index
    if (row < 1 || row > t->curr+1) {   // bounds check
      warning("row index out of bounds:", " ignored");
      return;
    }
  }

  pos = name_list_pos("knob", nl);
  knob = pos >= 0 ? pl->parameters[pos]->string : NULL;

  if (knob == NULL) {
    warning("invalid knob, not found:", " ignored");
    return;
  }

  pos = name_list_pos("scale", nl);
  scale =  pl->parameters[pos]->double_value;

  double varvalue[t->num_cols];
  for (int i = 0; i < t->num_cols; i++) {
    if (t->columns->inform[i] < 3) {
      varvalue[i] = get_variable(t->columns->names[i]);
      //printf("--%s--%g--\n",t->columns->names[i],varvalue[i]);
    }
  }
  //printf("--%s--%g--\n",knob,knobvalue);

  double knobvalue = get_variable(knob)+1;
  set_variable(knob, &knobvalue);

  int cols = t->org_cols, curr = t->curr;
  t->org_cols = 0;    t->curr = row - 1;
  for (int i = 0; i < t->num_cols; i++) {
    if (t->columns->inform[i] < 3) {
      t->d_cols[i][row-1] = scale *
                   (get_variable(t->columns->names[i]) - varvalue[i]);
      //printf("%d\n",row-1);
    }
  }
  t->org_cols = cols; t->curr = curr;
  knobvalue -= 1;
  set_variable(knob, &knobvalue);

  if (row == t->curr+1) // enlarge if needed
    if (++t->curr == t->max) grow_table(t);
}

void
exec_setvars_table(struct in_cmd* cmd)
  /* set variables from a table */
{
  struct command_parameter_list* pl = cmd->clone->par;
  struct name_list* nl = cmd->clone->par_names;
  const char* name = NULL;
  int pos, row;

  if ((pos = name_list_pos("table", nl)) < 0 || nl->inform[pos] == 0 ||
      (name = pl->parameters[pos]->string) == NULL) {
    warning("no table name:", "ignored");
    return;
  }

  if ((pos = name_list_pos(name, table_register->names)) < 0) {
    warning("table not found:", "ignored");
    return;
  }

  struct table* t = table_register->tables[pos];

  if ((pos = name_list_pos("row", nl)) < 0)
    row = t->curr;
  else {
    row = (int) pl->parameters[pos]->double_value;
    if (row < 1) row = t->curr+row+1;  // reflect negative index
    if (row < 1 || row > t->curr) {    // bounds check
      warning("row index out of bounds:", " ignored");
      return;
    }
  }

  current_node = NULL; /* to distinguish from other table fills, remanent! */
  int curr = t->curr;
  t->curr = row - 1;
  set_vars_from_table(t);
  t->curr = curr;
}


void
exec_setvars_lin_table(struct in_cmd* cmd)
  /* set variables from a table by linear interpolation between values in two rows */
{
  struct command_parameter_list* pl = cmd->clone->par;
  struct name_list* nl = cmd->clone->par_names;
  const char* name = NULL, *param = NULL;
  char expr[10*NAME_L];
  int pos, row1, row2;

  if ((pos = name_list_pos("table", nl)) < 0 || nl->inform[pos] == 0 ||
      (name = pl->parameters[pos]->string) == NULL) {
    warning("no table name:", "ignored");
    return;
  }

  if ((pos = name_list_pos(name, table_register->names)) < 0) {
    warning("table not found:", "ignored");
    return;
  }

  struct table* t = table_register->tables[pos];

  pos  = name_list_pos("row1", nl);
  row1 = pos >= 0 ? (int) pl->parameters[pos]->double_value : t->curr;
  pos  = name_list_pos("row2", nl);
  row2 = pos >= 0 ? (int) pl->parameters[pos]->double_value : t->curr;
  pos  = name_list_pos("param", nl);
  param = pos >= 0 ? pl->parameters[pos]->string : "interp";

  if (abs(row1) > t->curr || row1 == 0){
    warning("row1 index out of bounds:", " ignored");
    return;
  }
  if (abs(row2) > t->curr || row2 == 0){
    warning("row2 index out of bounds:", " ignored");
    return;
  }

  /* negative row numbers are counting backwards from last row */
  /* transform into positive values */
  if (row1 < 0) row1 = t->curr+row1+1;
  if (row2 < 0) row2 = t->curr+row2+1;

  current_node = NULL; /* to distinguish from other table fills, remanent! */

  for (int i = 0; i < t->num_cols; i++) {
    if (t->columns->inform[i] < 3) {
      const char *colname = t->columns->names[i];
      double val1 = t->d_cols[i][row1-1];
      double val2 = t->d_cols[i][row2-1];
      // 2014-Aug-18  17:15:08  ghislain:
      // value := val1*param + val2*(1-param) ;
      // sprintf(expr,"%s:=%10.16g*(%s)%+10.16g*(1-(%s));", colname,val1,param,val2,param);
      // is counterintuitve for interpolation between val1 and val2 and should instead be
      // value := val1 + param*(val2-val1) = val1*(1-param) + val2*param;
      sprintf(expr, "%s:=%10.16g*(1-(%s))%+10.16g*(%s);", colname, val1, param, val2, param);
      pro_input(expr);
    } else if (t->columns->inform[i] == 3) {
      set_stringvar(t->columns->names[i],t->s_cols[i][row1-1]) ;
    }
  }
}

void
exec_setvars_knob_table(struct in_cmd* cmd)
  /* add knob to variables with weights from a table*/
{
  struct command_parameter_list* pl = cmd->clone->par;
  struct name_list* nl = cmd->clone->par_names;
  const char* name = NULL, *knob = NULL;
  char expr[10000];
  char subexpr[100];
  int pos, row;
  int noappend;

  if ((pos = name_list_pos("table", nl)) < 0 || nl->inform[pos] == 0 ||
      (name = pl->parameters[pos]->string) == NULL) {
    warning("no table name:", "ignored");
    return;
  }

  if ((pos = name_list_pos(name, table_register->names)) < 0) {
    warning("table not found:", "ignored");
    return;
  }

  struct table* t = table_register->tables[pos];

  pos  = name_list_pos("row", nl);
  row  = pos >= 0 ? (int) pl->parameters[pos]->double_value : t->curr;
  pos  = name_list_pos("knob", nl);
  knob = pos >= 0 ? pl->parameters[pos]->string : NULL;
  pos  = name_list_pos("noappend",nl);
  noappend= pl->parameters[pos]->double_value;

  if (abs(row) > t->curr || row == 0) {
    warning("row index out of bounds:", " ignored");
    return;
  }

  if (knob == NULL) {
    warning("invalid knob, not found:", " ignored");
    return;
  }

  /* negative row numbers are counting backwards from last row */
  /* transform into positive values */
  if (row < 0) row = t->curr + 1 + row;

  current_node = NULL; /* to distinguish from other table fills, remament! */

  struct variable* var = NULL;
  for (int i = 0; i < t->num_cols; i++) {
    if (t->columns->inform[i] < 3) {
      const char *colname = t->columns->names[i];
      double val = t->d_cols[i][row-1];
      sprintf(subexpr,"%+24.16g*%s", val, knob);
      if ((noappend==0) && (var = find_variable(colname, variable_list))) {
        if (var->expr)
          sprintf(expr, "%s := %s %s;", colname, var->expr->string, subexpr);
        else
          if (var->value==0){
            sprintf(expr, "%s := %s;", colname, subexpr);
          } else {
            sprintf(expr, "%s := %+24.16g %s;", colname, var->value, subexpr);
          };
      } else
          sprintf(expr, "%s := %s;", colname, subexpr);
      pro_input(expr);
    }
  }
}

void
exec_setvars_const_table(struct in_cmd* cmd)
  /* add knob to variables with weights from a table*/
{
  struct command_parameter_list* pl = cmd->clone->par;
  struct name_list* nl = cmd->clone->par_names;
  const char* name = NULL;
  int pos;

  if ((pos = name_list_pos("table", nl)) < 0 || nl->inform[pos] == 0 ||
      (name = pl->parameters[pos]->string) == NULL) {
    warning("no table name:", "ignored");
    return;
  }

  if ((pos = name_list_pos(name, table_register->names)) < 0) {
    warning("table not found:", "ignored");
    return;
  }

  struct table* t = table_register->tables[pos];

  pos = name_list_pos("const", nl);
  double constant = pl->parameters[pos]->double_value;

  current_node = NULL; /* to distinguish from other table fills, remament! */

  for (int i = 0; i < t->num_cols; i++) {
    if (t->columns->inform[i] < 3) {
      const char *colname = t->columns->names[i];
      set_variable(colname, &constant);
    }
  }
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

void // this function extend print_value from mad_eval.c
exec_printf(struct in_cmd* cmd)
{
  struct command_parameter_list* pl = cmd->clone->par;
  struct name_list* nl = cmd->clone->par_names;

  // retrieve output format from text=""
  int txt_pos = name_list_pos("text", nl);
  if (!nl->inform[txt_pos]) { warning("missing text=:",""); return; }
  char *txt_str = v_format(pl->parameters[txt_pos]->string);

  // check for value=...
  int val_pos = name_list_pos("value", nl);
  if (!nl->inform[val_pos]) { warning("missing value=:",""); return; }

  // retrieve vector of values from value=...
  int val_n = command_par_vector("value", cmd->clone, NULL);
  if (val_n < 100) val_n = 100;
  double val[val_n];
  command_par_vector("value", cmd->clone, val);

  // enough to print a full twiss row, anyway C limits is 127, and var_form is long enough
  fprintf(prt_file, txt_str,
    val[ 0], val[ 1], val[ 2], val[ 3], val[ 4], val[ 5], val[ 6], val[ 7], val[ 8], val[ 9],
    val[10], val[11], val[12], val[13], val[14], val[15], val[16], val[17], val[18], val[19],
    val[20], val[21], val[22], val[23], val[24], val[25], val[26], val[27], val[28], val[29],
    val[30], val[31], val[32], val[33], val[34], val[35], val[36], val[37], val[38], val[39],
    val[40], val[41], val[42], val[43], val[44], val[45], val[46], val[47], val[48], val[49],
    val[50], val[51], val[52], val[53], val[54], val[55], val[56], val[57], val[58], val[59],
    val[60], val[61], val[62], val[63], val[64], val[65], val[66], val[67], val[68], val[69],
    val[70], val[71], val[72], val[73], val[74], val[75], val[76], val[77], val[78], val[79],
    val[80], val[81], val[82], val[83], val[84], val[85], val[86], val[87], val[88], val[89],
    val[90], val[91], val[92], val[93], val[94], val[95], val[96], val[97], val[98], val[99]);
  fprintf(prt_file, "\n");
}

