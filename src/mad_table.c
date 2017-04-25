#include "madx.h"

#ifndef _WIN32
#include <sys/utsname.h> // for uname
#endif

// private functions

static char*
get_table_string(char* left, char* right)
{
/* for command tabstring(table,column,row) where table = table name, */
/* column = name of a column containing strings, row = integer row number */
/* starting at 0, returns the string found in that column and row, else NULL */
  int col, ntok, pos, row;
  char** toks;
  struct table* table;
  *right = '\0';
  strcpy(c_dum->c, ++left);
  supp_char(',', c_dum->c);
  mysplit(c_dum->c, tmp_p_array);
  toks = tmp_p_array->p; ntok = tmp_p_array->curr;
  if (ntok == 3 && (pos = name_list_pos(toks[0], table_register->names)) > -1)
  {
    table = table_register->tables[pos];
    if ((col = name_list_pos(toks[1], table->columns)) > -1)
    {
      row = atoi(toks[2]);
      if(row > 0 && row <= table->curr && table->s_cols[col])
        return table->s_cols[col][row-1];
    }
  }
  return NULL;
}

static char*
get_table_index(char* left, char* right)
{
// for command tabindex(table,column,row,name) or
//             tabindex(table,column,name) where
// table = table name,
// column = name of a column containing strings,
// row is the starting row to search for (default 1),
// name = start of the name to search.
  int ntok, pos;
  char** toks;

  *right = '\0';
  strcpy(c_dum->c, ++left);
  supp_char(',', c_dum->c);
  mysplit(c_dum->c, tmp_p_array);
  toks = tmp_p_array->p; ntok = tmp_p_array->curr;
  if ((ntok == 3 || ntok == 4) &&
      (pos = name_list_pos(toks[0], table_register->names)) > -1)
  {
    struct table* table = table_register->tables[pos];
    char *name = toks[3 - (ntok == 3)]; // name is in position 2 or 3
    int len = strlen(name);
    int col = name_list_pos(toks[1], table->columns);
    int row = ntok == 4 ? atoi(toks[2]) : 1; // row may or may not be present
    if (col > -1 && row > 0 && table->s_cols[col]) {
      for (; row <= table->curr; row++) {
        if (mystrnicmp(table->s_cols[col][row-1], name, len) == 0) {
          sprintf(c_dum->c, "%d", row);
          return permbuff(c_dum->c);
        }
      }
    }
  }
  return NULL;
}

static int
tab_name_code(const char* name, const char* t_name)
  /* returns 1 if name corresponds to t_name, else 0 */
{
  char tmp[2*NAME_L];
  char *p, *n = one_string;
  strcpy(tmp, name);
  if ((p = strstr(tmp, "->")) != NULL)
  {
    *p = '\0'; p = strstr(name, "->"); p++; n = ++p;
  }
  if (strchr(t_name, ':'))
  {
    strcat(tmp, ":"); strcat(tmp, n);
  }
  return (strcmp(tmp, t_name) == 0 ? 1 : 0);
}

static int
table_row(struct table* table, const char* name)
{
  int i, j, ret = -1;
  for (i = 0; i < table->num_cols; i++)
  {
    if(table->columns->inform[i] == 3) {
      if (debuglevel > 2)
        printf("table_row: Column %d named <<%s>> is of strings. We use it to find the name.\n",
               i,table->columns->names[i]);
      break;
    }
  }

  if (i < table->num_cols) {
    for (j = 0; j < table->curr; j++)
    {
      if (debuglevel > 2) printf("table_row: Comparing <<%s>> <<%s>>\n",name, table->s_cols[i][j]);
      if (tab_name_code(name, table->s_cols[i][j])) break;
    }
    if (j < table->curr) ret = j;
  }
  else
  {
    if (debuglevel > 1) printf("Can not find a column to search for row containing %s\n",name);
  }
/*  if(ret==-1) fatal_error("Name of row not found", name);*/
  if(ret==-1) warning("table_row: Name of row not found:",name);
  return ret;
}

static void
add_table_vars(struct name_list* cols, struct command_list* select)
  /* 1: adds user selected variables to table - always type 2 = double
     2: adds aperture variables apertype (string) + aper_1, aper_2 etc. */
{
  int i, j, k, n, pos;
  char* var_name;
  char tmp[12];
  struct name_list* nl;
  struct command_parameter_list* pl;
  for (i = 0; i < select->curr; i++)
  {
    nl = select->commands[i]->par_names;
    pl = select->commands[i]->par;
    pos = name_list_pos("column", nl);
    if (nl->inform[pos])
    {
      for (j = 0; j < pl->parameters[pos]->m_string->curr; j++)
      {
        var_name = pl->parameters[pos]->m_string->p[j];
        if (strcmp(var_name, "apertype") == 0)
        {
          if ((n = aperture_count(current_sequ)) > 0)
          {
            add_to_name_list(permbuff("apertype"), 3, cols);
            for (k = 0; k < n; k++)
            {
              sprintf(tmp, "aper_%d", k+1);
              add_to_name_list(permbuff(tmp), 2, cols);
            }
          }
        }
        else if (name_list_pos(var_name, cols) < 0) /* not yet in list */
          add_to_name_list(permbuff(var_name), 2, cols);
      }
    }
  }
}

static void
grow_table_list(struct table_list* tl)
{
  const char *rout_name = "grow_table_list";
  struct table** t_loc = tl->tables;
  int new = 2*tl->max;

  grow_name_list(tl->names);
  tl->max = new;
  tl->tables = mycalloc(rout_name, new, sizeof *tl->tables);
  for (int j = 0; j < tl->curr; j++) tl->tables[j] = t_loc[j];
  myfree(rout_name, t_loc);
}

static void
grow_table_list_list(struct table_list_list* tll)
{
  const char *rout_name = "grow_table_list_list";
  struct table_list** t_loc = tll->table_lists;
  int new = 2*tll->max;

  tll->max = new;
  tll->table_lists = mycalloc(rout_name, new, sizeof *tll->table_lists);
  for (int j = 0; j < tll->curr; j++) tll->table_lists[j] = t_loc[j];
  myfree(rout_name, t_loc);
}

static void
add_to_table_list_list(struct table_list* table_list, struct table_list_list* tll)
  /* adds a table_list to a list of table_lists */
{
  int j;
  for (j = 0; j < tll->curr; j++)
    if (tll->table_lists[j] == table_list) return;
  if (tll->curr == tll->max) grow_table_list_list(tll);
  tll->table_lists[tll->curr++] = table_list;
}

static void
write_table(struct table* t, const char* filename)
  /* writes rows with columns listed in row and col */
{
  assert(t && filename);
  char l_name[NAME_L];
  struct int_array* col = t->col_out;
  struct int_array* row = t->row_out;
  int i, j, k, tmp, n;
  time_t now;
  struct tm* tm;
#if 0
  char sys_name[200];
#ifndef _WIN32
  struct utsname u;
  i = uname(&u); /* get system name */
  strcpy(sys_name, u.sysname);
#else // _WIN32
  strcpy(sys_name, "Win32");
#endif
#endif

  time(&now);    /* get system time */
  tm = localtime(&now); /* split system time */
  if (strcmp(filename, "terminal") == 0) out_file = stdout;
  else if ((out_file = fopen(filename, "w")) == NULL)
  {
    warning("cannot open output file:", filename); return;
  }
  if (t != NULL)
  {
    strcpy(l_name, t->name);
    n = strlen(t->name);
    fprintf(out_file,
            "@ NAME             %%%02ds \"%s\"\n", n,
            stoupper(l_name));

    strcpy(l_name, t->type);
    n = strlen(t->type);
    fprintf(out_file,
            "@ TYPE             %%%02ds \"%s\"\n", n,
            stoupper(l_name));

    if (t->header != NULL)
    {
      for (j = 0; j < t->header->curr; j++)
        fprintf(out_file, "%s\n", t->header->p[j]);
    }
    if (title != NULL)
    {
      n = strlen(title);
      fprintf(out_file,
              "@ TITLE            %%%02ds \"%s\"\n", n, title);
    }

    n = strlen(version_name)+strlen(version_ostype)+strlen(version_arch)+2;
    fprintf(out_file,
            "@ ORIGIN           %%%02ds \"%s %s %s\"\n",
            n, version_name, version_ostype, version_arch);

    fprintf(out_file,
            "@ DATE             %%08s \"%02d/%02d/%02d\"\n",
            tm->tm_mday, tm->tm_mon+1, tm->tm_year%100);

    fprintf(out_file,
            "@ TIME             %%08s \"%02d.%02d.%02d\"\n",
            tm->tm_hour, tm->tm_min, tm->tm_sec);
    fprintf(out_file, "* ");

    for (i = 0; i < col->curr; i++)
    {
      strcpy(l_name, t->columns->names[col->i[i]]);
//      printf("name: col %d->%d:%s\n",i, col->i[i], l_name);

      if (t->columns->inform[col->i[i]] == 1)
        fprintf(out_file, v_format("%NIs "), stoupper(l_name));
      else if (t->columns->inform[col->i[i]] == 2)
        fprintf(out_file, v_format("%NFs "), stoupper(l_name));
      else if (t->columns->inform[col->i[i]] == 3)
        fprintf(out_file, v_format("%S "), stoupper(l_name));
    }
    fprintf(out_file, "\n");

    fprintf(out_file, "$ ");
    for (i = 0; i < col->curr; i++)
    {
//      printf("type: col %d->%d:%d\n",i, col->i[i], t->columns->inform[col->i[i]]);

      if (t->columns->inform[col->i[i]] == 1)
        fprintf(out_file, v_format("%NIs "),"%d");
      else if (t->columns->inform[col->i[i]] == 2)
        fprintf(out_file, v_format("%NFs "),"%le");
      else if (t->columns->inform[col->i[i]] == 3)
        fprintf(out_file, v_format("%S "),"%s");
    }
    fprintf(out_file, "\n");

    for (j = 0; j < row->curr; j++)
    {
      if (row->i[j])
      {
        if (t->l_head[j] != NULL)
        {
          for (k = 0; k < t->l_head[j]->curr; k++)
            fprintf(out_file, "%s\n", t->l_head[j]->p[k]);
        }
        for (i = 0; i < col->curr; i++)
        {
//          printf("row %d col %d datatype %d \n",j,i, t->columns->inform[col->i[i]] );
          if (t->columns->inform[col->i[i]] == 1) {
            tmp = t->d_cols[col->i[i]][j];
            fprintf(out_file, v_format(" %I"), tmp);
          }
          else if (t->columns->inform[col->i[i]] == 2) {
            fprintf(out_file, v_format(" %F"), t->d_cols[col->i[i]][j]);
//          printf("%s[%2d,%2d]=%+8.5f    ",t->name,col->i[i],j,t->d_cols[col->i[i]][j]);
          }
          else if (t->columns->inform[col->i[i]] == 3) {
            char *pc = c_dum->c;
            pc[0]='"', k=1;
            if (t->s_cols[col->i[i]][j] != NULL) {
//              printf("%s[%2d,%2d]=%s    ",t->name,col->i[i],j, t->s_cols[col->i[i]][j]);
              strcpy(&c_dum->c[1], t->s_cols[col->i[i]][j]);
              stoupper(c_dum->c);
              pc = strip(c_dum->c); /* remove :<occ_count> */
              k = strlen(pc);
            }
            pc[k++] = '"'; pc[k] = '\0';
            fprintf(out_file, v_format(" %S "), pc);
          }
        }
        fprintf(out_file, "\n");
      }
    }
    if (strcmp(filename, "terminal") != 0) fclose(out_file);
  }
}

#if 0 // not used...
static int
result_from_normal(char* name_var, int* order, double* val)
  /* returns value of table normal_results corresponding to the given variable name
     and to the given orders
     function value return:
     0  OK
     -1 table  does not exist
     -2 column does not exist
     -3 row    does not exist
  */
{
  int row,k,found,pos;
  char string[AUX_LG],n_var[AUX_LG];
  double d_val=zero;
  struct table* t;

  pos = name_list_pos("normal_results", table_register->names);
  t = table_register->tables[pos];

  *val = zero;
  found = 0;
  mycpy(n_var, name_var);
  for (row = 1; row <= t->curr; row++)
  {
    k = string_from_table_row("normal_results","name", &row, string);
    if (k != 0) return k;
    if (strcmp(string,n_var) == 0)
    {
      found = 1;
      k = double_from_table_row("normal_results","order1", &row, &d_val);
      if ((int)d_val != order[0]) found = 0;
      k = double_from_table_row("normal_results","order2", &row, &d_val);
      if ((int)d_val != order[1]) found = 0;
      k = double_from_table_row("normal_results","order3", &row, &d_val);
      if ((int)d_val != order[2]) found = 0;
      k = double_from_table_row("normal_results","order4", &row, &d_val);
      if ((int)d_val != order[3]) found = 0;
    }
    if (found == 1) break;
  }
  if (found == 1)
    k = double_from_table_row("normal_results","value", &row, &d_val);
  *val = d_val;
  return 0;
}
#endif

#if 0 // not used...
static struct table*
read_his_table(struct in_cmd* cmd)
  /* reads and stores TFS table */
{
  struct table* t = NULL;
  struct char_p_array* tcpa = NULL;
  struct name_list* tnl = NULL;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int pos = name_list_pos("file", nl);
  int i, k, error = 0;
  char *cc, *filename, *type = NULL, *tmp, *name;

  if(nl->inform[pos] && (filename = pl->parameters[pos]->string) != NULL)
  {
    if ((tab_file = fopen(filename, "r")) == NULL)
    {
      warning("cannot open file:", filename); return NULL;
    }
  }
  else
  {
    warning("no filename,","ignored"); return NULL;
  }
  while (fgets(aux_buff->c, aux_buff->max, tab_file))
  {
    cc = strtok(aux_buff->c, " \"\n");
    if (*cc == '@')
    {
      if ((tmp = strtok(NULL, " \"\n")) != NULL
          && strcmp(tmp, "TYPE") == 0)
      {
        if ((name = strtok(NULL, " \"\n")) != NULL) /* skip format */
        {
          if ((name = strtok(NULL, " \"\n")) != NULL)
            type = permbuff(stolower(name));
        }
      }
    }
    else if (*cc == '*' && tnl == NULL)
    {
      tnl = new_name_list("table_names", 20);
      while ((tmp = strtok(NULL, " \"\n")) != NULL)
        add_to_name_list(permbuff(stolower(tmp)), 0, tnl);
    }
    else if (*cc == '$' && tcpa == NULL)
    {
      if (tnl == NULL)
      {
        warning("formats before names","skipped"); return NULL;
      }
      tcpa = new_char_p_array(20);
      while ((tmp = strtok(NULL, " \"\n")) != NULL)
      {
        if (tcpa->curr == tcpa->max) grow_char_p_array(tcpa);
        if (strcmp(tmp, "%s") == 0)       tnl->inform[tcpa->curr] = 3;
        else if (strcmp(tmp, "%hd") == 0) tnl->inform[tcpa->curr] = 1;
        else                              tnl->inform[tcpa->curr] = 2;
        tcpa->p[tcpa->curr++] = permbuff(tmp);
      }
    }
    else
    {
      if(t == NULL)
      {
        if (type == NULL)
        {
          warning("TFS table without type,","skipped"); error = 1;
        }
        else if (tcpa == NULL)
        {
          warning("TFS table without formats,","skipped"); error = 1;
        }
        else if (tnl == NULL)
        {
          warning("TFS table without column names,","skipped"); error = 1;
        }
        else if (tnl->curr == 0)
        {
          warning("TFS table: empty column name list,","skipped");
          error = 1;
        }
        else if (tnl->curr != tcpa->curr)
        {
          warning("TFS table: number of names and formats differ,",
                  "skipped");
          error = 1;
        }
        if (error)
        {
          delete_name_list(tnl); return NULL;
        }
        t = new_table(type, "input", 500, tnl);
      }
      for (i = 0; i < tnl->curr; i++)
      {
        if (t->curr == t->max) grow_table(t);
        tmp = tcpa->p[i];
        if (strcmp(tmp,"%s") == 0)
          t->s_cols[i][t->curr] = tmpbuff(stolower(cc));
        else if (strcmp(tmp,"%d") == 0 || strcmp(tmp,"%hd") == 0)
         {
          sscanf(cc, tmp, &k); t->d_cols[i][t->curr] = k;
         }
        else sscanf(cc, tmp, &t->d_cols[i][t->curr]);
        if (i+1 < tnl->curr)
        {
          if ((cc =strtok(NULL, " \"\n")) == NULL)
          {
            warning("read_his_table: incomplete table line starting with:", aux_buff->c);
            return NULL;
          }
        }
      }
      t->curr++;
    }
  }
  fclose(tab_file);
  t->origin = 1;
  add_to_table_list(t, table_register);
  return NULL;
}
#endif

static void
set_selected_rows(struct table* t, struct command_list* select, struct command_list* deselect)
{
  int i, j;

  if (!current_sequ) {
    warning("No current selection available, skipping select", t->name);
    return;
  }

  c_range_start = get_node_count(current_sequ->range_start);
  c_range_end = get_node_count(current_sequ->range_end);

  get_select_t_ranges(select, deselect, t);
  if (select != NULL) {
    for (j = 0; j < t->curr; j++)  t->row_out->i[j] = 0;
    for (i = 0; i < select->curr; i++) {
      for (j = s_range->i[i]; j <= e_range->i[i]; j++) {
        if (t->row_out->i[j] == 0) {
          if (!t->s_cols[0])
            warning("Invalid column type (string expected)", t->name);
          else
            t->row_out->i[j] = pass_select(t->s_cols[0][j], select->commands[i]);
        }
      }
    }
  }
  if (deselect != NULL) {
    for (i = 0; i < deselect->curr; i++) {
      for (j = sd_range->i[i]; j <= ed_range->i[i]; j++) {
        if (t->row_out->i[j] == 1) {
          if (!t->s_cols[0])
            warning("Invalid column type (string expected)", t->name);
          else
            t->row_out->i[j] = 1 - pass_select(t->s_cols[0][j], deselect->commands[i]);
        }
      }
    }
  }
}

// public interface

struct table*
new_table(const char* name, const char* type, int rows, struct name_list* cols)
{
  const char *rout_name = "new_table";
  int i, n = cols->curr;
  struct table* t = mycalloc(rout_name, 1, sizeof *t);

  strcpy(t->name, name);
  strcpy(t->type, type);
  t->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", "table");
  t->columns = cols;
  t->num_cols = t->org_cols = n;
  t->s_cols = mycalloc(rout_name, n, sizeof *t->s_cols);
  t->d_cols = mycalloc(rout_name, n, sizeof *t->d_cols);
  t->max = ++rows; /* +1 because of separate augment_count */
  for (i = 0; i < n; i++) {
    if (cols->inform[i] < 3)
      t->d_cols[i] = mycalloc_atomic(rout_name, rows, sizeof *t->d_cols[0]);
    else if (cols->inform[i] == 3)
      t->s_cols[i] = mycalloc(rout_name, rows, sizeof *t->s_cols[0]);
  }
  t->row_out = new_int_array(rows);
  t->col_out = new_int_array(n);
  t->node_nm = new_char_p_array(rows);
  t->p_nodes = mycalloc(rout_name, rows, sizeof *t->p_nodes);
  t->l_head  = mycalloc(rout_name, rows, sizeof *t->l_head);
  return t;
}

struct table_list*
new_table_list(int size)
{
  const char *rout_name = "new_table_list";
  struct table_list* tl = mycalloc(rout_name, 1, sizeof *tl);
  strcpy(tl->name, "table_list");
  tl->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", tl->name);
  tl->max = size;
  tl->curr = 0;
  tl->names = new_name_list(tl->name, size);
  tl->tables = mycalloc(rout_name, size, sizeof *tl->tables);
  add_to_table_list_list(tl, all_table_lists);
  return tl;
}

struct table_list_list*
new_table_list_list(int size)
{
  const char *rout_name = "new_table_list_list";
  struct table_list_list* tll = mycalloc(rout_name, 1, sizeof *tll);
  strcpy(tll->name, "table_list_list");
  tll->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", tll->name);
  tll->max = size;
  tll->curr = 0;
  tll->table_lists = mycalloc(rout_name, size, sizeof *tll->table_lists);
  return tll;
}

void
check_table(char* string)
  /* replaces argument of "table" if any by a string variable */
{
  char *pa, *pb, *pt, *pl, *pr, *sv;
  pa = string;
  while ((pb = strstr(pa, "table")) != NULL)
  {
    if (is_token(pb, string, 5))
    {
      if (quote_level(pa, pb) == 0)
      {
        mystrcpy(c_join, pa);                             // global var
        pt = strstr(c_join->c, "table");
        if ((pl = strchr(pt, '(')) == NULL) return;
        if ((pr = strchr(pl, ')')) == NULL) return;
        *pl = '\0';
        *pr = '\0';
        sv = make_string_variable(++pl);
        *pa ='\0';
        strcat(string, c_join->c);
        strcat(string, " ( ");
        strcat(string, sv);
        strcat(string, " ) ");
        strcat(string, ++pr);
      }
    }
    pa = ++pb;
  }
}

void
check_tabindex(char* string)
  /* replaces argument of "tabindex(tab_name, col_name, row_start, string)" if any by the row index
     return 0 if not found
   */
{
  char *pa, *pb, *pt, *pl, *pr, *sv;
  pa = string;
  while ((pb = strstr(pa, "tabindex")) != NULL)
  {
    if (is_token(pb, string, strlen("tabindex")))
    {
      if (quote_level(pa, pb) == 0)
      {
        mystrcpy(c_join, pa);                             // global var
        pt = strstr(c_join->c, "tabindex");
        if ((pl = strchr(pt, '(')) == NULL) return;
        if ((pr = strchr(pl, ')')) == NULL) return;
        if ((sv = get_table_index(pl,pr)) == NULL) sv = permbuff("0");
        *pt = '\0';
        *pa ='\0';
        strcat(string, c_join->c);
        strcat(string, sv);
        strcat(string, ++pr);
      }
    }
    pa = ++pb;
  }
}

void
check_tabstring(char* string)
  /* replaces tabstring(tab_name, col_name, row_number) by the
     string found in that column/row of table tab_name, or by "_void_"
     if not found */
{
  char *pa, *pb, *pt, *pl, *pr, *sv;
  pa = string;
  while ((pb = strstr(pa, "tabstring")) != NULL)
  {
    if (is_token(pb, string, 9))
    {
      if (quote_level(pa, pb) == 0)
      {
        mystrcpy(c_join, pa);
        pt = strstr(c_join->c, "tabstring");
        if ((pl = strchr(pt, '(')) == NULL) return;
        if ((pr = strchr(pl, ')')) == NULL) return;
        if ((sv = get_table_string(pl,pr)) == NULL) sv = permbuff("_void_");
        *pt = '\0';
        *pa ='\0';
        strcat(string, c_join->c);
        strcat(string, sv);
        strcat(string, ++pr);
      }
    }
    pa = ++pb;
  }
}

double
table_value(void)
{
  double val = zero;
  int ntok, pos, col, row;
  char** toks;
  struct table* table;
  char temp[NAME_L];

  if (current_variable != NULL && current_variable->string != NULL) {
    strcpy(c_dum->c, current_variable->string);
    supp_char(',', c_dum->c);
    mysplit(c_dum->c, tmp_p_array);
    toks = tmp_p_array->p; ntok = tmp_p_array->curr;
    /* toks[0] -> table name
       toks[1] -> row name
       toks[2] -> col name
    */
    if (ntok > 1)
     {
      if ((pos = name_list_pos(toks[0], table_register->names)) > -1)
       {
        table = table_register->tables[pos];
        if ((col = name_list_pos(toks[ntok-1], table->columns)) > -1)
         {
          if (ntok > 2)
           { /* find row - else current (dynamic), or 0 */
            /* start mod - HG 26.3.2011 */
            if (ntok > 5)
             { /* check for [ count ] and convert to ->count */
               if (*toks[2] == '[' && *toks[4] == ']')
                {
                   strcat(toks[1], "->");
                   strcat(toks[1], toks[3]);
                 }
             }
            /* end mod - HG 26.3.2011 */
            row = table_row(table, toks[1]);
           }
          else
           {
             if (table->dynamic)
              {
                row = table->curr;
              }
             else
              {
                row = 0;
              }
           }

          if (row >= 0) /*in case table_row(table, toks[1]); returns  -1 == row not found*/
           {
             if (col >= table->num_cols)
              {
                printf("trying to get column %d out of range %d\n",col,table->num_cols);
                if (get_option("no_fatal_stop ")==0) exit(1);
                return val;
              }

             if (row >= table->max)
              {
                printf("trying to get row %d of range %d\n",row,table->max);
                if (get_option("no_fatal_stop ")==0) exit(1);
                return val;
              }

             /*printf("val %f         col %d of %d >>>>  row %d of %d ",val, col, table->num_cols, row, table->max);  */
             val = table->d_cols[col][row];
           }

        }/*column found*/
        else if ((ntok == 3) && ((col = name_list_pos(toks[1], table->columns)) > -1))
         {
          row = atoi(toks[2])-1;
          if(row < table->curr)
           {
            val = table->d_cols[col][row];
           }
        }
        else if(ntok == 2)
         {
          strncpy(temp, toks[1], NAME_L);
          if (mystricmp(temp, "tablelength") == 0)
           {
             val = table->curr;
           }
         }
      } /*pos > -1, table name found in the list*/
    } /*ntok > 0*/
  }/*current variable*/

  return val;
}

struct column_info
table_get_column(char* table_name, char* column_name)
{
  struct column_info info={NULL,0,'V',0};
  int pos, col; // not used , i;
  struct table* table;
  if ((pos = name_list_pos(table_name, table_register->names)) > -1) {
    table = table_register->tables[pos];
    if ((col = name_list_pos(column_name, table->columns)) > -1) {
      //printf("col: n %d type %d\n",col,table->columns->inform[col]);
      info.length = table->curr;
      if (table->columns->inform[col]==1) {
        info.data=table->d_cols[col];
        info.datatype='i';
        info.datasize=sizeof(double);
      } else if (table->columns->inform[col]==2) {
        info.data=table->d_cols[col];
        info.datatype='d';
        info.datasize=sizeof(double);
      } else if (table->columns->inform[col]==3) {
        info.data=table->s_cols[col];
        info.datasize=NAME_L;
        info.datatype='S';
      };
    }
  }
  return info;
}

struct char_p_array *
table_get_header(char* table_name)
{
  int pos;
  if ((pos = name_list_pos(table_name, table_register->names)) > -1)
    return table_register->tables[pos]->header;
  // table was not found, we return 0 pointer..
  return NULL;
}

void
augment_count(const char* table) /* increase table occ. by 1, fill missing */
{
  int pos;
  struct table* t;
  mycpy(c_dum->c, table);
  if ((pos = name_list_pos(c_dum->c, table_register->names)) > -1)
    t = table_register->tables[pos];
  else {
    warning("Can not find table",table);
    return;
  }

  if (strcmp(t->type, "twiss") == 0) complete_twiss_table(t);

  if (t->num_cols > t->org_cols)  add_vars_to_table(t,1);

  if (t->p_nodes != NULL) t->p_nodes[t->curr] = current_node;

  if (t->node_nm != NULL)
  {
    t->node_nm->p[t->curr] = current_node->name;
    t->node_nm->curr = t->curr;
  }
  if (++t->curr == t->max) grow_table(t);
}

void
augmentcountonly(const char* table) /* increase table occ. by 1 */
{
  int pos;
  struct table* t;
  mycpy(c_dum->c, table);
  if ((pos = name_list_pos(c_dum->c, table_register->names)) > -1)
    t = table_register->tables[pos];
  else
  {
    warning("Can not find table",table);
    return;
  }

  if (t->num_cols > t->org_cols)  add_vars_to_table(t,1);

  if (++t->curr == t->max) grow_table(t);
}

void
add_to_table_list(struct table* t, struct table_list* tl)
  /* adds table t to table list tl */
{
  int pos; //, j; not used
  if ((pos = name_list_pos(t->name, tl->names)) < 0)
  {
    if (tl->curr == tl->max) grow_table_list(tl);
    add_to_name_list(tmpbuff(t->name), 0, tl->names); // j = not used
    tl->tables[tl->curr++] = t;
  }
  else
  {
    tl->tables[pos] = delete_table(tl->tables[pos]);
    tl->tables[pos] = t;
  }
}

struct table *
detach_table_from_table_list(const char *name, struct table_list* tl)
  /* detaches table named 'name' from table list 'tl' */
  /* returns NULL if the table isn't found in table list 'tl' */
{
  struct table *retval = NULL;
  int pos;
  if ((pos = name_list_pos(name, tl->names)) > -1)
  {
    int k = remove_from_name_list(tl->tables[pos]->name,
                                  tl->names);

    retval = tl->tables[pos];

    tl->tables[k] = tl->tables[--tl->curr];
  }
  return retval;
}

int
remove_table_from_table_list(const char *name, struct table_list* tl)
  /* removes table named 'name' from table list 'tl' */
  /* returns -1 if the table isn't found in table list 'tl' */
{
  struct table *t = NULL;
  if ((t = detach_table_from_table_list(name, tl))) {
    delete_table(t);
  }
  return t ? 0 : -1;
}

void
add_vars_to_table(struct table* t, double scale)
  /* fills user-defined variables into current table_row) */
{

  for (int i = t->org_cols; i < t->num_cols; i++)
  {
    if (t->columns->inform[i] < 3)
    {
      if (strstr(t->columns->names[i], "aper_"))
        t->d_cols[i][t->curr] = get_aperture(current_node, t->columns->names[i]);
      else if (strstr(t->columns->names[i], "aptol_"))
        t->d_cols[i][t->curr] = get_apertol(current_node, t->columns->names[i]);
      else {
        t->d_cols[i][t->curr] = get_variable(t->columns->names[i])*scale;
      }
    }
    else if (current_node)
     {
       char* p = command_par_string(t->columns->names[i], current_node->p_elem->def) ;
       t->s_cols[i][t->curr] = tmpbuff(p ? p : "none");
     }
    else
     {
       t->s_cols[i][t->curr] = tmpbuff(get_varstring(t->columns->names[i]));
     }
  }
}

void
set_vars_from_table(struct table* t)
  /* set variables from current table_row) */
{
  int i;

  for (i = 0; i < t->num_cols; i++)
  {
    if (t->columns->inform[i] ==2)
    {
      set_variable(t->columns->names[i],&t->d_cols[i][t->curr]) ;
    }
    else if (t->columns->inform[i] ==3)
    {
      set_stringvar(t->columns->names[i],t->s_cols[i][t->curr]) ;
    }
  }
}

struct table*
delete_table(struct table* t)
{
  const char *rout_name = "delete_table";
  int i, j;
  if (t == NULL) return NULL;
  if (stamp_flag && t->stamp != 123456)
    fprintf(stamp_file, "d_t double delete --> %s\n", t->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", "table");
  if (t->header != NULL) t->header = delete_char_p_array(t->header, 1);
  if (t->col_out != NULL) t->col_out = delete_int_array(t->col_out);
  if (t->row_out != NULL) t->row_out = delete_int_array(t->row_out);
  if (t->node_nm != NULL) t->node_nm = delete_char_p_array(t->node_nm, 0);
  for (i = 0; i < t->curr; i++) {
    if (t->l_head[i] != NULL)
      t->l_head[i] = delete_char_p_array(t->l_head[i], 1);
  }
  if (t->l_head)  myfree(rout_name, t->l_head);
  if (t->p_nodes) myfree(rout_name, t->p_nodes);

  if (t->d_cols) {
    for (i = 0; i < t->num_cols; i++)
      if (t->columns->inform[i] < 3 && t->d_cols[i])
        myfree(rout_name, t->d_cols[i]);
    myfree(rout_name, t->d_cols);
  }

  if (t->s_cols) {
    for (i = 0; i < t->num_cols; i++)
     {
       if (t->columns->inform[i] == 3 && t->s_cols[i])
       {
         for (j = 0; j < t->curr; j++)
           {
              /*printf("%d %d %s %#x\n",i,j,t->s_cols[i][j], t->s_cols[i][j]); */
              myfree(rout_name, t->s_cols[i][j]);
           }
         /*printf("Deleting column  %d %#x\n",i,t->s_cols[i]);*/

         myfree(rout_name, t->s_cols[i]);
       }
    }
    /*printf("\n\n%d %s %#x\n\n\n",i,t->s_cols[i], t->s_cols[i]);   */
    myfree(rout_name, t->s_cols);
  }
  t->columns = delete_name_list(t->columns);
  myfree(rout_name, t);
  return NULL;
}

void
double_table(char* table)
{
  int pos;
  struct table* t;

  mycpy(c_dum->c, table);
  if ((pos = name_list_pos(c_dum->c, table_register->names)) > -1)
    t = table_register->tables[pos];
  else return;
  grow_table(t);
}

void
grow_table(struct table* t) /* doubles number of rows */
{
  const char *rout_name = "grow_table";
  int i, j, new = 2*t->max;
  char** s_loc;
  struct char_p_array* t_loc = t->node_nm;
  double* d_loc;
  struct node** p_loc = t->p_nodes;
  struct char_p_array** pa_loc = t->l_head;

  t->max = new;
  t->p_nodes = mycalloc(rout_name, new, sizeof *t->p_nodes);
  t->l_head  = mycalloc(rout_name, new, sizeof *t->l_head);
  t->node_nm = new_char_p_array(new);

  for (i = 0; i < t->curr; i++) {
    t->node_nm->p[i] = t_loc->p[i];
    t->p_nodes[i] = p_loc[i];
    t->l_head[i] = pa_loc[i];
  }
  delete_char_p_array(t_loc, 0);
  myfree(rout_name, pa_loc);
  t->node_nm->curr = t->curr; myfree(rout_name, p_loc);
  for (j = 0; j < t->num_cols; j++) {
    if ((s_loc = t->s_cols[j]) != NULL) {
      t->s_cols[j] = mycalloc(rout_name, new, sizeof *t->s_cols[0]);
      for (i = 0; i < t->curr; i++) t->s_cols[j][i] = s_loc[i];
      myfree(rout_name, s_loc);
    }
  }
  for (j = 0; j < t->num_cols; j++) {
    if ((d_loc = t->d_cols[j]) != NULL) {
      t->d_cols[j] = mycalloc_atomic(rout_name, new, sizeof *t->d_cols[0]);
      for (i = 0; i < t->curr; i++) t->d_cols[j][i] = d_loc[i];
      myfree(rout_name, d_loc);
    }
  }
}

void
print_table(struct table* t)
{
  int i, j, k, l, n, tmp, wpl = 4;
  if (t != NULL)
  {
    fprintf(prt_file, "\n");
    fprintf(prt_file, "++++++ table: %s\n", t->name);
    l = (t->num_cols-1) / wpl + 1;
    for (k = 0; k < l; k++)
    {
      n = wpl*(k+1) > t->num_cols ? t->num_cols : wpl*(k+1);
      fprintf(prt_file, "\n");
      for (i = wpl*k; i < n; i++)
      {
        if (t->columns->inform[i] == 1)
          fprintf(prt_file, v_format("%NIs "), t->columns->names[i]);
        else if (t->columns->inform[i] == 2)
          fprintf(prt_file, v_format("%NFs "), t->columns->names[i]);
        else if (t->columns->inform[i] == 3)
          fprintf(prt_file, v_format("%S "), t->columns->names[i]);
      }
      fprintf(prt_file, "\n");
      for (j = 0; j < t->curr; j++)
      {
        for (i = wpl*k; i < n; i++)
        {
          if (t->columns->inform[i] == 1)
          {
            tmp = t->d_cols[i][j];
            fprintf(prt_file, v_format("%I "), tmp);
          }
          else if (t->columns->inform[i] == 2)
            fprintf(prt_file, v_format("%F "), t->d_cols[i][j]);
          else if (t->columns->inform[i] == 3)
            fprintf(prt_file, v_format("%S "), t->s_cols[i][j]);
        }
        fprintf(prt_file, "\n");
      }
    }
  }
}

void
make_map_table(int* map_table_max_rows)
{
  assert(map_table_max_rows);
  assert(table_register->names);

  remove_table_from_table_list("map_table", table_register);

  /* initialise table */
  map_table = make_table("map_table", "map_tab", map_tab_cols,
                         map_tab_types, *map_table_max_rows);

  assert(map_table);
  assert(table_register);

  add_to_table_list(map_table, table_register);
  map_table->dynamic = 1;
  reset_count("map_table");

  /*
  printf("Creating map table  Done \n");
  pos = name_list_pos("map_table", table_register->names);
  printf("Checking position of the table: pos %d \n", pos);
  */

}

struct table*
make_table(const char* name, const char* type, const char* const* table_cols, const int* table_types, int rows)
{
  struct table* t;
  struct name_list *cols;
  struct command_list* scl;
  int i, n = 0;
  while (*table_cols[n] != ' ')
  {
/*     printf("make table %s col %d %s\n",name, n, table_cols[n]);*/
    n++;
  }
  cols = new_name_list("columns", n);
  for (i = 0; i < n; i++)
    add_to_name_list(table_cols[i], table_types[i], cols);
  if ((scl = find_command_list(name, table_select)) != NULL && scl->curr > 0)
    add_table_vars(cols, scl);
  t = new_table(name, type, rows, cols);
  t->org_cols = n;
  return t;
}

void
reset_count(const char* table) /* resets table counter to zero */
{
  int pos;
  struct table* t;
  mycpy(c_dum->c, table);
  if ((pos = name_list_pos(c_dum->c, table_register->names)) > -1)
    t = table_register->tables[pos];
  else return;
  t->curr = 0;
}

void
sector_out(char* sector_table_name, double* pos, double* kick, double* rmatrix, double* tmatrix)
{
  int i;
  int j;
  int k;
  int index;

  /* the name is not \0 finished and displays ugly in C */
  /* but well, this is how table names are handled... */

  char * elementName = current_node->p_elem->name;

  string_to_table_curr( sector_table_name, "name", elementName );
  double_to_table_curr( sector_table_name, "pos", pos );

  /* 6 kicks */
  for (i=0; i<6; i++){
    char kickStr[2+1];
    sprintf(kickStr,"k%i",i+1);
    kickStr[2]='\0';
    index = i;
    double_to_table_curr( sector_table_name, kickStr,&kick[index]);
  }
  /* 36 R-matrix terms */
  for (j=0; j<6; j++){
    for (i=0; i<6; i++){
      char rStr[3+1];
      sprintf(rStr,"r%i%i",i+1,j+1);
      rStr[3]='\0';
      index = i+j*6;
      double_to_table_curr( sector_table_name, rStr,&rmatrix[index]);
    }
  }
  /* 216 T-matrix terms */
  for (k=0; k<6; k++){
    for (j=0; j<6; j++){
      for (i=0;i<6; i++){
	char tStr[4+1];
	sprintf(tStr,"t%i%i%i",i+1,j+1,k+1);
	tStr[4]='\0';
	index = i+j*6+k*36;
	double_to_table_curr( sector_table_name, tStr,&tmatrix[index]);
      }
    }
  }

  augment_count( sector_table_name ); /* move to next record */
}

void
out_table(const char* tname, struct table* t, const char* filename)
  /* output of a table */
{
  struct command_list* scl = find_command_list(tname, table_select);
  struct command_list* dscl = find_command_list(tname, table_deselect);

  while (t->num_cols > t->col_out->max)
    grow_int_array(t->col_out);

  while (t->curr > t->row_out->max)
    grow_int_array(t->row_out);

  t->row_out->curr = t->curr;
  if (par_present("full", NULL, scl))
    put_info("obsolete option 'full'"," ignored on 'select'");

  for (int j = 0; j < t->curr    ; j++) t->row_out->i[j] = 1;
  for (int j = 0; j < t->num_cols; j++) t->col_out->i[j] = j;

  t->col_out->curr = t->num_cols;
  if ((scl != NULL && scl->curr > 0) || (dscl != NULL && dscl->curr > 0))
  {
    set_selected_columns(t, scl);
    set_selected_rows(t, scl, dscl);
  }
  write_table(t, filename);
}

struct table*
read_table(struct in_cmd* cmd)
  /* reads and stores TFS table */
{
  struct table* t = NULL;
  struct char_p_array* tcpa = NULL;
  struct name_list* tnl = NULL;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int pos = name_list_pos("file", nl);
  short sk;
  int i, k, error = 0;
  char *cc, *filename, *type = NULL, *tmp, *name;

  char* namtab;

  if ((namtab = command_par_string("table",cmd->clone)) != NULL) {
    printf("Want to make named table: %s\n",namtab);
  }
  else
  {
    if (get_option("debug")) {
      printf("No table name requested\n");
      printf("Use default name (i.e. name from file) \n");
    }
    namtab = NULL;
  }

  if(nl->inform[pos] && (filename = pl->parameters[pos]->string) != NULL)
  {
    if ((tab_file = fopen(filename, "r")) == NULL)
    {
      fatal_error("cannot open file:", filename); return NULL;
    }
  }
  else
  {
    warning("no filename,","ignored"); return NULL;
  }
  while (fgets(aux_buff->c, aux_buff->max, tab_file))
  {
    supp_char('\r', aux_buff->c);
    cc = strtok(aux_buff->c, " \"\n");
    if (*cc == '@')
    {
      if ((tmp = strtok(NULL, " \"\n")) != NULL
          && strcmp(tmp, "TYPE") == 0)
      {
        if ((name = strtok(NULL, " \"\n")) != NULL) /* skip format */
        {
          if ((name = strtok(NULL, " \"\n")) != NULL)
            type = permbuff(stolower(name));
        }
      }
      else if (strcmp(tmp, "NAME") == 0)
      {
        if ((name = strtok(NULL, " \"\n")) != NULL) /* skip format */
        {
          if ((name = strtok(NULL, " \"\n")) != NULL)
            namtab = permbuff(stolower(name));
        }
      }
    }
    else if (*cc == '*' && tnl == NULL)
    {
      tnl = new_name_list("table_names", 20);
      while ((tmp = strtok(NULL, " \"\n")) != NULL)
        add_to_name_list(permbuff(stolower(tmp)), 0, tnl);
    }
    else if (*cc == '$' && tcpa == NULL)
    {
      if (tnl == NULL)
      {
        warning("formats before names","skipped"); return NULL;
      }
      tcpa = new_char_p_array(20);
      while ((tmp = strtok(NULL, " \"\n")) != NULL)
      {
        if (tcpa->curr == tcpa->max) grow_char_p_array(tcpa);
        if (strcmp(tmp, "%s") == 0)       tnl->inform[tcpa->curr] = 3;
        else if (strcmp(tmp, "%d") == 0)  tnl->inform[tcpa->curr] = 1;
        else if (strcmp(tmp, "%hd") == 0) tnl->inform[tcpa->curr] = 1;
        else                              tnl->inform[tcpa->curr] = 2;
        tcpa->p[tcpa->curr++] = permbuff(tmp);
      }
    }
    else
    {
      if(t == NULL)
      {
        if (type == NULL)
        {
          warning("TFS table without type,","skipped"); error = 1;
        }
        else if (tcpa == NULL)
        {
          warning("TFS table without formats,","skipped"); error = 1;
        }
        else if (tnl == NULL)
        {
          warning("TFS table without column names,","skipped"); error = 1;
        }
        else if (tnl->curr == 0)
        {
          warning("TFS table: empty column name list,","skipped");
          error = 1;
        }
        else if (tnl->curr != tcpa->curr)
        {
          warning("TFS table: number of names and formats differ,",
                  "skipped");
          error = 1;
        }
        if (error)
        {
          delete_name_list(tnl); return NULL;
        }
        if(namtab != NULL) {
          t = new_table(namtab, type,    500, tnl);
        }
        else
        {
          t = new_table(type, type,    500, tnl);
        }
      }
      for (i = 0; i < tnl->curr; i++)
      {
        if (t->curr == t->max) grow_table(t);
        tmp = tcpa->p[i];
        if (strcmp(tmp,"%s") == 0) t->s_cols[i][t->curr] = stolower(tmpbuff(cc));
        else if (strcmp(tmp,"%d") == 0)
        {
          sscanf(cc, tmp, &k); t->d_cols[i][t->curr] = k;
        }
        else if (strcmp(tmp,"%hd") == 0)
        {
          sscanf(cc, tmp, &sk); t->d_cols[i][t->curr] = sk;
        }
        else sscanf(cc, tmp, &t->d_cols[i][t->curr]);
        if (i+1 < tnl->curr)
        {
          if ((cc =strtok(NULL, " \"\n")) == NULL)
          {
            warning("read_table: incomplete table line starting with:", aux_buff->c);
            return NULL;
          }
        }
      }
      t->curr++;
    }
  }
  fclose(tab_file);
  if ((tab_file = fopen(filename, "r")) == NULL)
  {
    warning("cannot open file:", filename); return NULL;
  }
/* read & store table header */
  t->header = new_char_p_array(50);
  while (fgets(aux_buff->c, aux_buff->max, tab_file))
  {
    supp_char('\r', aux_buff->c);
    if ((*aux_buff->c != ' ') &&
        ((*aux_buff->c == '@') || (*aux_buff->c == '*')))
    {
      if (t->header->curr == t->header->max) grow_char_p_array(t->header);
      int len = strlen(aux_buff->c)+1;
      t->header->p[t->header->curr] = mymalloc_atomic("read_table", len * sizeof *t->header->p[0]);
      strcpy(t->header->p[t->header->curr], aux_buff->c);
      t->header->curr++;
    }
  }
  fclose(tab_file);
  t->origin = 1;
  add_to_table_list(t, table_register);
  return NULL;
}

int
get_table_range(char* range, struct table* table, int* rows)
  /* returns start and end row (rows[0] and rows[1])
     of a range in a table; 0 if not found, 1 (1 row) or 2 (> 1) */
{
  char tmp[2*NAME_L], buf[5*NAME_L], *c[2], *p;
  int i, n;

  if (table == NULL) return 0;

  n = 1;
  rows[1] = rows[0] = 0;

  stolower(strcpy(buf, range));

  c[0] = buf;
  if ((p = strchr(buf, '/')) != NULL) {
    *p = 0; // cut the string in two pieces (i.e. the ranges)
    c[1] = p+1;
    n = 2;
  }

//  fprintf(stderr, "get_table_range: range='%s', c[0]='%s', c[1]='%s'\n", range, c[0], c[1]);

  for (i = 0; i < n; i++) {
    if (*c[i] == '#') {
      if (strncmp(c[i], "#s", 2) == 0) rows[i] = 0;
      else if (strncmp(c[i], "#e", 2) == 0) rows[i] = table->curr - 1;
      else {
        warning("illegal table range ignored:", range);
        return 0;
      }
    }
    else {
      strcpy(tmp, c[i]);
      if (square_to_colon(tmp) == 0) {
        warning("illegal table range ignored:", range);
        return 0;
      }
      if ((rows[i] = char_p_pos(tmp, table->node_nm)) < 0) {
        warning("illegal table range ignored:", range);
        return 0;
      }
    }
  }

  if (n == 1) rows[1] = rows[0];

  return n;
}

void
table_range(char* table, char* range, int* rows)
  /* returns first and last row numbers (start=1) in rows
     or 0 if table or range invalid */
{
  int pos;
  struct table* t;
  char buf[5*NAME_L];

  rows[0] = rows[1] = 0;
  stolower(mycpy(buf, table));
  if ((pos = name_list_pos(buf, table_register->names)) > -1) {
    t = table_register->tables[pos];
    mycpy(buf, range);
    get_table_range(buf, t, rows);
    rows[0]++, rows[1]++;
  } else {
    warning("invalid table name, range ignored (invalid results may occur!) for table", buf);
  }

//  fprintf(stderr, "table_range: row[0]='%d', row[1]='%d', table->curr=%d\n", rows[0], rows[1], t->curr);
}

struct table*
read_my_table(struct in_cmd* cmd)
     /* reads and stores TFS table */
{
  struct table* t = NULL;
  struct char_p_array* tcpa = NULL;
  struct name_list* tnl = NULL;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int pos = name_list_pos("file", nl);
  int i, k, error = 0;
  short  sk;
  char *cc, *filename, *type = NULL, *tmp, *name;
  double tmpd;
  char* namtab;

  if ((namtab = command_par_string("table",cmd->clone)) != NULL) {
       printf("Want to make named table: %s\n",namtab);
  } else {
       if (get_option("debug")) {
         printf("No table name requested\n");
         printf("Use default name (i.e. name from file) \n");
       }
       namtab = NULL;
  }

  if(nl->inform[pos] && (filename = pl->parameters[pos]->string) != NULL)
    {
     if ((tab_file = fopen(filename, "r")) == NULL)
       {
         fatal_error("cannot open file:", filename); return NULL; /* frs: to avoid unwanted results */
       }
    }
  else
    {
     warning("no filename,","ignored"); return NULL;
    }
  while (fgets(aux_buff->c, aux_buff->max, tab_file))
    {
     cc = strtok(aux_buff->c, " \"\n");
     if (*cc == '@')
       {
       if ((tmp = strtok(NULL, " \"\n")) != NULL
              && strcmp(tmp, "TYPE") == 0)
        {
         if ((name = strtok(NULL, " \"\n")) != NULL) /* skip format */
           {
            if ((name = strtok(NULL, " \"\n")) != NULL)
                  type = permbuff(stolower(name));
           }
        }
       }
     else if (*cc == '*' && tnl == NULL)
       {
      tnl = new_name_list("table_names", 20);
        while ((tmp = strtok(NULL, " \"\n")) != NULL)
            add_to_name_list(permbuff(stolower(tmp)), 0, tnl);
       }
     else if (*cc == '$' && tcpa == NULL)
       {
      if (tnl == NULL)
        {
         warning("formats before names","skipped"); return NULL;
        }
      tcpa = new_char_p_array(20);
        while ((tmp = strtok(NULL, " \"\n")) != NULL)
        {
         if (tcpa->curr == tcpa->max) grow_char_p_array(tcpa);
           if (strcmp(tmp, "%s") == 0)       tnl->inform[tcpa->curr] = 3;
           else if (strcmp(tmp, "%hd") == 0) tnl->inform[tcpa->curr] = 1;
           else if (strcmp(tmp, "%d") == 0)  tnl->inform[tcpa->curr] = 1;
           else                              tnl->inform[tcpa->curr] = 2;
           tcpa->p[tcpa->curr++] = permbuff(tmp);
        }
       }
     else
       {
        if(t == NULL)
          {
         if (type == NULL)
           {
            warning("TFS table without type,","skipped"); error = 1;
           }
         else if (tcpa == NULL)
           {
            warning("TFS table without formats,","skipped"); error = 1;
           }
         else if (tnl == NULL)
           {
            warning("TFS table without column names,","skipped"); error = 1;
           }
         else if (tnl->curr == 0)
           {
            warning("TFS table: empty column name list,","skipped");
              error = 1;
           }
         else if (tnl->curr != tcpa->curr)
           {
            warning("TFS table: number of names and formats differ,",
                       "skipped");
              error = 1;
           }
           if (error)
           {
            delete_name_list(tnl); return NULL;
           }
           if(namtab != NULL) {
             t = new_table(namtab, type,    500, tnl);
           } else {
             t = new_table(type, type,    500, tnl);
           }
        }
      for (i = 0; i < tnl->curr; i++)
        {
         if (t->curr == t->max) grow_table(t);
         tmp = tcpa->p[i];

        /*
         printf("read_my_table %d <%s> <%s> \n", i, tcpa->p[i], tnl->names[i] );
        */

         if (strcmp(tmp,"%s") == 0)
            {
             /* printf("reading format %s \n", tmp); */
             t->s_cols[i][t->curr] = stolower(tmpbuff(cc));
             /* printf("read %d %d = %s \n", i,t->curr, t->s_cols[i][t->curr]); */

              /*printf("read_my_table coln [%d %d]=%s\n",i,t->curr,t->s_cols[i][t->curr]);*/
            }
           else if (strcmp(tmp,"%d") == 0 )
           {
             /* printf("reading format %s \n", tmp); */
             sscanf(cc, tmp, &k);
             /* printf("read %d %d = %d \n", i,t->curr, k); */
             t->d_cols[i][t->curr] = k;
           }
           else if (strcmp(tmp,"%hd") == 0 )
           {
              /* printf("reading format %s \n", tmp); */
              sscanf(cc, tmp, &sk);
              /* printf("read %d %d = %d \n", i,t->curr, sk); */
              t->d_cols[i][t->curr] = sk;
           }
           else
           {
              /* printf("reading format %s \n", tmp); */
              sscanf(cc, tmp, &tmpd);
              /* printf("read %d %d = %f \n", i,t->curr, tmpd); */
              t->d_cols[i][t->curr] = tmpd;
           }


           if (i+1 < tnl->curr)
           {
             if ((cc =strtok(NULL, " \"\n")) == NULL)
              {
               warning("read_my_table: incomplete table line starting with:", aux_buff->c);
                 return NULL;
              }
           }
        }
        t->curr++;
       }
    }
  fclose(tab_file);
  t->origin = 1;
  add_to_table_list(t, table_register);
  return NULL;
}

void
set_selected_columns(struct table* t, struct command_list* select)
{
  int i, j, pos, k, n = 0;
  char* p;
  struct name_list* nl;
  struct command_parameter_list* pl;
  if (select && par_present("column", NULL, select))
  {
    for (j = 0; j < t->num_cols; j++)  /* deselect all columns */
      t->col_out->i[j] = 0;
    t->col_out->curr = 0;
    for (i = 0; i < select->curr; i++)
    {
      nl = select->commands[i]->par_names;
      pl = select->commands[i]->par;
      pos = name_list_pos("column", nl);
      if (nl->inform[pos])
      {
        for (j = 0; j < pl->parameters[pos]->m_string->curr; j++)
        {
          if (strcmp(pl->parameters[pos]->m_string->p[j], "re") == 0)
          {
            for (k = 0; k < t->num_cols; k++)
            {
              if (strncmp("re", t->columns->names[k], 2) == 0)
              {
                if (k <  t->num_cols
                    && int_in_array(k, n, t->col_out->i) == 0)
                  t->col_out->i[n++] = k;
              }
            }
          }
          else if (strcmp(pl->parameters[pos]->m_string->p[j], "eign") == 0)
          {
            for (k = 0; k < t->num_cols; k++)
            {
              if (strncmp("eign", t->columns->names[k], 2) == 0)
              {
                if (k <  t->num_cols
                    && int_in_array(k, n, t->col_out->i) == 0)
                  t->col_out->i[n++] = k;
              }
            }
          }
          else if (strcmp(pl->parameters[pos]->m_string->p[j],
                          "apertype") == 0)
          {
            for (k = 0; k < t->num_cols; k++)
            {
              if (strncmp("aper", t->columns->names[k], 4) == 0)
              {
                if (k <  t->num_cols
                    && int_in_array(k, n, t->col_out->i) == 0)
                  t->col_out->i[n++] = k;
              }
            }
          }
          else
          {
            p = pl->parameters[pos]->m_string->p[j];
            if ((k = name_list_pos(p, t->columns)) > -1)
            {
              if (k <  t->num_cols
                  && int_in_array(k, n, t->col_out->i) == 0)
                t->col_out->i[n++] = k;
            }
          }
        }
      }
    }
    t->col_out->curr = n;
  }
}

/*
  Grouping accessors
*/

#if 0
int
str_from_table(char* table, char* name, int* row, char* val)
     /* WH 22.06.2004, corrected from: char_from_table */
     /* returns val at position row in column with name "name".
        function value return:
        0  OK
        -1 table  does not exist
        -2 column does not exist
        -3 row    does not exist
     */
{
  int pos;
  struct table* t;

  strcpy(val,"No-Name");
  mycpy(c_dum->c, table);
  if ((pos = name_list_pos(c_dum->c, table_register->names)) > -1)
    t = table_register->tables[pos];
  else return -1;
  mycpy(c_dum->c, name);
  if ((pos = name_list_pos(c_dum->c, t->columns)) < 0) return -2;
  if (*row > t->curr)  return -3;
   strncpy(val,t->s_cols[pos][*row-1],NAME_L);
  // completely buggy and slow
  // while (strlen(val)<NAME_L) val[strlen(val)]=' ';
  val[NAME_L-1] = '\0';
  return 0;
}

int
str_from_tablet(struct table *t, char* name, int* row, char* val)
     /* WH 22.06.2004, corrected from: char_from_table */
     /* returns val at position row in column with name "name".
        function value return:
        0  OK
        -1 table  does not exist
        -2 column does not exist
        -3 row    does not exist
     */
{
  int pos;

  strcpy(val,"No-Name");
  mycpy(c_dum->c, name);
  if ((pos = name_list_pos(c_dum->c, t->columns)) < 0) return -2;
  if (*row > t->curr)  return -3;
   strncpy(val,t->s_cols[pos][*row-1],NAME_L);
  // competely buggy, useless and slow
  // while (strlen(val)<NAME_L) val[strlen(val)]=' ';
  val[NAME_L-1] = '\0';
  return 0;
}

// dangerous function that uses table->node_nm sometimes corrupted or unmaintained
int
nodename_from_table_row(const char* table, const int* row, char* string)
  /* returns NODE NAME at position row (name is DISCARDED).
     function value return:
     0  OK
     -1 table  does not exist
     -2 column does not exist // not used...
     -3 row    does not exist
  */
{
  char buf[NAME_L];
  struct table* tbl;
  int pos;

  *string = '\0';

  mycpy(buf, table);
  if ((pos = name_list_pos(buf, table_register->names)) < 0 ||
     !(tbl = table_register->tables[pos])) {
    warning("nodename_from_table_row: name of table not found:" , buf);
    return -1;
  }
  if (*row < 1 || *row > tbl->curr) {
    warning("nodename_from_table_row: row out of range", "");
    return -3;
  }

  strcpy(string, tbl->node_nm->p[*row-1]);
  return 0;
}

#endif

int
table_length(const char* table)
  /* returns no. of rows in table */
{
  char tbl_s[NAME_L];
  struct table *tbl;
  int pos;

  mycpy(tbl_s, table);
  if ((pos = name_list_pos(tbl_s, table_register->names)) < 0 ||
     !(tbl = table_register->tables[pos])) {
    warning("table_length: table not found:", tbl_s);
    return 0;
  }
  return tbl->curr;
}

int
table_exists(const char* table)
  /* returns no. of rows in table */
{
  char tbl_s[NAME_L];
  int pos;

  mycpy(tbl_s, table);
  if ((pos = name_list_pos(tbl_s, table_register->names)) < 0 ||
     !table_register->tables[pos])
    return 0;

  return 1;
}

int
table_column_exists(const char* table, const char *name)
{
  char tbl_s[NAME_L], col_s[NAME_L];
  struct table *tbl;
  int pos;

  mycpy(tbl_s, table);
  if ((pos = name_list_pos(tbl_s, table_register->names)) < 0 ||
     !(tbl = table_register->tables[pos]))
    return 0;

  mycpy(col_s, name);
  if (name_list_pos(col_s, tbl->columns) < 0)
    return 0;

  return 1;
}

int
table_cell_exists(const char* table, const char* name, const int* row)
{
  char tbl_s[NAME_L], col_s[NAME_L];
  struct table *tbl;
  int pos;

  mycpy(tbl_s, table);
  if ((pos = name_list_pos(tbl_s, table_register->names)) < 0 ||
     !(tbl = table_register->tables[pos]))
    return 0;

  mycpy(col_s, name);
  if (name_list_pos(col_s, tbl->columns) < 0)
    return 0;

  return *row >= 1 && *row <= tbl->curr;
}

int
table_header_exists(const char* table, const char *name)
{
  char tbl_s[NAME_L], hdr_s[NAME_L], buf[256];
  struct table *tbl;
  int pos, hdr;
  char *p;

  mycpy(tbl_s, table);
  if ((pos = name_list_pos(tbl_s, table_register->names)) < 0 ||
     !(tbl = table_register->tables[pos]))
    return 0;

  mycpy(hdr_s, name);
  if (tbl->header) {
    for (hdr = 0; hdr < tbl->header->curr; hdr++) {
      strcpy(buf, &tbl->header->p[hdr][1]);
      if ((p=strtok(buf, " \"\n")) && mystricmp(p, hdr_s) == 0)
        return 1;
    }
  }
  return 0;
}

int
double_from_table_header(const char* table, const char* name, double* val)
  /* returns val from table header at position "name", if present.
     function value return:
     0  OK
     -1 table does not exist
     -2 header or parameter does not exist
     -3 parameter value does not exist or is not a number
  */
{
  char tbl_s[NAME_L], hdr_s[NAME_L], buf[256];
  struct table *tbl;
  int pos, hdr;
  char *p;

  *val = 0.0;

  mycpy(tbl_s, table);
  if ((pos = name_list_pos(tbl_s, table_register->names)) < 0 ||
     !(tbl = table_register->tables[pos])) {
    warning("double_from_table_header: table not found:", tbl_s);
    return -1;
  }

  mycpy(hdr_s, name);
  if (tbl->header) {
    for (hdr = 0; hdr < tbl->header->curr; hdr++) {
      strcpy(buf, &tbl->header->p[hdr][1]);
      if ((p=strtok(buf, " \"\n")) && mystricmp(p, hdr_s) == 0) {
        if (strstr(strtok(NULL, " \"\n"), "%le") == NULL) {
          warning("double_from_table_header: parameter without value in table header:", (sprintf(buf,"%s->%s",tbl_s,hdr_s),buf));
          return -3;
        }
        if (sscanf(strtok(NULL, " \"\n"), "%le", val) != 1) {
          warning("double_from_table_header: invalid parameter value in table header:", (sprintf(buf,"%s->%s",tbl_s,hdr_s),buf));
          return -3;
        }
        return 0;
      }
    }
    warning("double_from_table_header: parameter not found in table header:", (sprintf(buf,"%s->%s",tbl_s,hdr_s),buf));
    return -2;
  } else {
    warning("double_from_table_header: table has no header:", tbl_s);
    return -2;
  }
}

int
double_from_table_row(const char* table, const char* name, const int* row, double* val)
  /* returns val at position row in column with name "name".
     function value return:
     0  OK
     -1 table  does not exist
     -2 column does not exist
     -3 row    does not exist
  */
{
  char tbl_s[NAME_L], col_s[NAME_L], buf[5*NAME_L];
  struct table *tbl;
  int pos, col;

  *val = 0.0;

  mycpy(tbl_s, table);
  if ((pos = name_list_pos(tbl_s, table_register->names)) < 0 ||
     !(tbl = table_register->tables[pos])) {
    warning("double_from_table_row: table not found:", tbl_s);
    return -1;
  }

  mycpy(col_s, name);
  if ((col = name_list_pos(col_s, tbl->columns)) < 0) {
    warning("double_from_table_row: column not found:", (sprintf(buf,"%s->%s",tbl_s,col_s),buf));
    return -2;
  }
  if (tbl->columns->inform[col] >= 3) {
    warning("double_from_table_row: invalid column type:", (sprintf(buf,"%s->%s",tbl_s,col_s),buf));
    return -2;
  }
  if (*row < 1 || *row > tbl->curr) {
    warning("double_from_table_row: row out of range:", (sprintf(buf,"%s->%s[1>=%d<=%d]",tbl_s,col_s,*row,tbl->curr),buf));
    return -3;
  }

  *val = tbl->d_cols[col][*row-1];
  return 0;
}

int
string_from_table_row(const char* table, const char* name, const int* row, char* string)
  /* returns string at position row in column with name "name".
     assumes string to be long enough...
     function value return:
     0  OK
     -1 table  does not exist
     -2 column does not exist
     -3 row    does not exist
  */
{
  char tbl_s[NAME_L], col_s[NAME_L], buf[5*NAME_L];
  struct table* tbl;
  int pos, col;

  *string = '\0';

  mycpy(tbl_s, table);
  if ((pos = name_list_pos(tbl_s, table_register->names)) < 0 ||
     !(tbl = table_register->tables[pos])) {
    warning("string_from_table_row: table not found:", tbl_s);
    return -1;
  }
  mycpy(col_s, name);
  if ((col = name_list_pos(col_s, tbl->columns)) < 0) {
    warning("string_from_table_row: column not found:", (sprintf(buf,"%s->%s",tbl_s,col_s),buf));
    return -2;
  }
  if (tbl->columns->inform[col] != 3) {
    warning("string_from_table_row: invalid column type:", (sprintf(buf,"%s->%s",tbl_s,col_s),buf));
    return -2;
  }
  if (*row < 1 || *row > tbl->curr) {
    warning("string_from_table_row: row out of range:", (sprintf(buf,"%s->%s[1>=%d<=%d]",tbl_s,col_s,*row,tbl->curr),buf));
    return -3;
  }

  strcpy(string, tbl->s_cols[col][*row-1]);
  return 0;
}

int
double_to_table_row(const char* table, const char* name, const int* row, const double* val)
  /* puts val at row position in column with name "name".
     0  OK
     -1 table  does not exist
     -2 column does not exist
     -3 row    does not exist
  */
{
  char tbl_s[NAME_L], col_s[NAME_L], buf[5*NAME_L];
  struct table* tbl;
  int pos, col;

  mycpy(tbl_s, table);
  if ((pos = name_list_pos(tbl_s, table_register->names)) < 0 ||
     !(tbl = table_register->tables[pos])) {
    warning("double_to_table_row: table not found:", tbl_s);
    return -1;
  }
  mycpy(col_s, name);
  if ((col = name_list_pos(col_s, tbl->columns)) < 0) {
    warning("double_to_table_row: column not found:", (sprintf(buf,"%s->%s",tbl_s,col_s),buf));
    return -2;
  }
  if (tbl->columns->inform[col] >= 3) {
    warning("double_to_table_row: invalid column type:", (sprintf(buf,"%s->%s",tbl_s,col_s),buf));
    return -2;
  }
  if (*row < 1 || *row > tbl->curr) {
    warning("double_to_table_row: row out of range:", (sprintf(buf,"%s->%s[1>=%d<=%d]",tbl_s,col_s,*row,tbl->curr),buf));
    return -3;
  }

  tbl->d_cols[col][*row-1] = *val;
  return 0;
}

int
string_to_table_row(const char* table, const char* name, const int *row, const char* string)
  /* puts string at row position in column with name "name".
     0  OK
     -1 table  does not exist
     -2 column does not exist
     -3 row    does not exist
  */
{
  char tbl_s[NAME_L], col_s[NAME_L], buf[5*NAME_L];
  struct table* tbl;
  int pos, col;

  mycpy(tbl_s, table);
  if ((pos = name_list_pos(tbl_s, table_register->names)) < 0 ||
     !(tbl = table_register->tables[pos])) {
    warning("string_to_table_row: table not found:", tbl_s);
    return -1;
  }
  mycpy(col_s, name);
  if ((col = name_list_pos(col_s, tbl->columns)) < 0) {
    warning("string_to_table_row: column not found:", (sprintf(buf,"%s->%s",tbl_s,col_s),buf));
    return -2;
  }
  if (tbl->columns->inform[col] != 3) {
    warning("string_to_table_row: invalid column type:", (sprintf(buf,"%s->%s",tbl_s,col_s),buf));
    return -2;
  }
  if (*row < 1 || *row > tbl->curr) {
    warning("string_to_table_row: row out of range:", (sprintf(buf,"%s->%s[1>=%d<=%d]",tbl_s,col_s,*row,tbl->curr),buf));
    return -3;
  }

  if (tbl->s_cols[col][*row-1])
    myfree("string_to_table_row", tbl->s_cols[col][*row-1]);

  mycpy(buf, string);
  if (strcmp(buf, "name") == 0)
    tbl->s_cols[col][*row-1] = tmpbuff(current_node->name);
  else if (strcmp(buf, "base_name") == 0)
    tbl->s_cols[col][*row-1] = tmpbuff(current_node->base_name);
  else
    tbl->s_cols[col][*row-1] = tmpbuff(buf);
  return 0;
}

int
double_to_table_curr(const char* table, const char* name, const double* val)
  /* puts val at current position in column with name "name".
     The table count is increased separately with "augment_count"
     0  OK
     -1 table  does not exist
     -2 column does not exist
     -3 row    does not exist (need expansion)
  */
{
  char tbl_s[NAME_L], col_s[NAME_L], buf[5*NAME_L];
  struct table* tbl;
  int pos, col;

  mycpy(tbl_s, table);
  if ((pos = name_list_pos(tbl_s, table_register->names)) < 0 ||
     !(tbl = table_register->tables[pos])) {
    warning("double_to_table_curr: table not found:", tbl_s);
    return -1;
  }
  mycpy(col_s, name);
  if ((col = name_list_pos(col_s, tbl->columns)) < 0) {
    warning("double_to_table_curr: column not found:", (sprintf(buf,"%s->%s",tbl_s,col_s),buf));
    return -2;
  }
  if (tbl->columns->inform[col] >= 3) {
    warning("double_to_table_curr: invalid column type:", (sprintf(buf,"%s->%s",tbl_s,col_s),buf));
    return -2;
  }
  if (tbl->curr >= tbl->max) {
    warning("double_to_table_curr: row out of range (need expansion):", (sprintf(buf,"%s->%s[%d<%d]",tbl_s,col_s,tbl->curr,tbl->max),buf));
    return -3;
  }

  tbl->d_cols[col][tbl->curr] = *val;
  return 0;
}

int
vector_to_table_curr(const char* table, const char* name, const double* vals, const int* nval)
  /* puts nval values of array vals at the current line into columns starting with column name
     The table count is increased separately with "augment_count"
     0  OK
     -1 table  does not exist
     -2 column does not exist
     -3 row    does not exist (need expansion)
  */
{
  assert(table);
  assert(name);
  assert(vals);
  assert(nval);

  char tbl_s[NAME_L], col_s[NAME_L], buf[5*NAME_L];
  struct table* tbl;
  int pos, col, last, j;

  mycpy(tbl_s, table);
  if ((pos = name_list_pos(tbl_s, table_register->names)) < 0 ||
     !(tbl = table_register->tables[pos])) {
    warning("vector_to_table_curr: table not found:", tbl_s);
    return -1;
  }


  mycpy(col_s, name);
  if ((col = name_list_pos(col_s, tbl->columns)) < 0) {
    warning("vector_to_table_curr: column not found: ", (sprintf(buf,"%s->%s",tbl_s,col_s),buf));
    return -2;
  }
  if (tbl->curr >= tbl->max) {
    warning("vector_to_table_curr: row out of range (need expansion):",
      (sprintf(buf,"%s->%s[%d<%d]",tbl_s,col_s,tbl->curr,tbl->max),buf));
    return -3;
  }

  if (col + *nval > tbl->num_cols) {
    warning("vector_to_table_curr: too many values provided - vector truncated:",
      (sprintf(buf,"%s->%s[%d<=%d]",tbl_s,col_s,col + *nval,tbl->num_cols),buf));
    last = tbl->num_cols;
  } else
    last = col + *nval;

  for (j = col; j < last; j++) {
    if (tbl->columns->inform[j] >= 3)
      warning("vector_to_table_curr: invalid column type - value skipped:",
        (sprintf(buf,"%s->%s",tbl_s,tbl->columns->names[j]),buf));
    else
      tbl->d_cols[j][tbl->curr] = vals[j-col];
  }

  return 0;
}

int
string_to_table_curr(const char* table, const char* name, const char* string)
  /* puts string at current position in column with name "name".
     The table count is increased separately with "augment_count"
     0  OK
     -1 table  does not exist
     -2 column does not exist
     -3 row    does not exist (need expansion)
  */
{
  char tbl_s[NAME_L], col_s[NAME_L], buf[5*NAME_L];
  struct table* tbl;
  int pos, col;

  mycpy(tbl_s, table);
  if ((pos = name_list_pos(tbl_s, table_register->names)) < 0 ||
     !(tbl = table_register->tables[pos])) {
    warning("string_to_table_curr: table not found:", tbl_s);
    return -1;
  }


  mycpy(col_s, name);
  if ((col = name_list_pos(col_s, tbl->columns)) < 0) {
    warning("string_to_table_curr: column not found:", (sprintf(buf,"%s->%s",tbl_s,col_s),buf));
    return -2;
  }
  if (tbl->columns->inform[col] != 3) {
    warning("string_to_table_curr: invalid column type:", (sprintf(buf,"%s->%s",tbl_s,col_s),buf));
    return -2;
  }
  if (tbl->curr >= tbl->max) {
    warning("string_to_table_curr: row out of range (need expansion):", (sprintf(buf,"%s->%s[%d<%d]",tbl_s,col_s,tbl->curr,tbl->max),buf));
    return -3;
  }

  if (tbl->s_cols[col][tbl->curr])
    myfree("string_to_table_curr", tbl->s_cols[col][tbl->curr]);

  mycpy(buf, string);
  if (strcmp(buf, "name") == 0)
    tbl->s_cols[col][tbl->curr] = tmpbuff(current_node->name);
  else if (strcmp(buf, "base_name") == 0)
    tbl->s_cols[col][tbl->curr] = tmpbuff(current_node->base_name);
  else
    tbl->s_cols[col][tbl->curr] = tmpbuff(buf);
  return 0;
}

int
comment_to_table_curr(const char* table, const char* comment, const int* length)
  /* Saves the comment string at the current line.
     This comment is then printed in front of this line.
     Several calls to the same current line are possible.
     0  OK
     -1 table does not exist
  */
{
  char tbl_s[NAME_L];
  struct table* tbl;
  int pos;

  mycpy(tbl_s, table);
  if ((pos = name_list_pos(tbl_s, table_register->names)) < 0 ||
     !(tbl = table_register->tables[pos])) {
    warning("comment_to_table_curr: table not found:" , tbl_s);
    return -1;
  }

  strncpy(c_dum->c, comment, *length); c_dum->c[*length] = '\0';
  if (tbl->l_head[tbl->curr] == NULL)
    tbl->l_head[tbl->curr] = new_char_p_array(2);
  else if (tbl->l_head[tbl->curr]->curr == tbl->l_head[tbl->curr]->max)
    grow_char_p_array(tbl->l_head[tbl->curr]);
  tbl->l_head[tbl->curr]->p[tbl->l_head[tbl->curr]->curr++] = tmpbuff(c_dum->c);
  return 0;
}

void
rename_table(struct table *tbl, const char *name )
{
  strncpy(tbl->name, name, NAME_L-1);
  tbl->name[NAME_L-1] = '\0';
}

void
table_add_header(struct table* t, const char* format, ...)
{
    va_list args;
    va_start(args, format);
    vsprintf(c_dum->c, v_format(format), args);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);
    va_end(args);
}

#if 0 // not used...
/*
  LD: 2012.11.29
  - These "slow" accessors were added for stable access of "muting" tables
    that is tables where columns and rows may change during data processing
  - As for other access-by-name, it looks for the first column of strings
  - The row name used for searching is the mangled name as stored in the
    table, that is with the trailing :# (count number), if any.
*/
static int
get_table_row(const struct table* tbl, const char* name)
{
  int col, row = tbl->curr;

  for (col = 0; col < tbl->num_cols; col++)
    if(tbl->columns->inform[col] == 3) break;

  if (col < tbl->num_cols)
    for (row = 0; row < tbl->curr; row++)
      if (!strcmp(name, tbl->s_cols[col][row])) break;

  return row == tbl->curr ? -1 : row;
}

double
get_table_value(const char* tbl_s, const char *row_s, const char *col_s)
{
  int pos, row, col;

  if ((pos = name_list_pos(tbl_s, table_register->names)) > -1) {
    const struct table *tbl = table_register->tables[pos];
    if ((col = name_list_pos(col_s, tbl->columns)) > -1) {
      if ((row = get_table_row(tbl, row_s)) > -1)
        return tbl->d_cols[col][row];

      else warning("get_table_value: name of row not found:"   , row_s);
    } else warning("get_table_value: name of column not found:", col_s);
  }   else warning("get_table_value: name of table not found:" , tbl_s);

  return 0;
}

void
set_table_value(const char* tbl_s, const char *row_s, const char *col_s, double *val)
{
  int pos, row, col;

  if ((pos = name_list_pos(tbl_s, table_register->names)) > -1) {
    const struct table *tbl = table_register->tables[pos];
    if ((col = name_list_pos(col_s, tbl->columns)) > -1) {
      if ((row = get_table_row(tbl, row_s)) > -1)
        tbl->d_cols[col][row] = *val;

      else warning("get_table_value: name of row not found:"   , row_s);
    } else warning("get_table_value: name of column not found:", col_s);
  }   else warning("get_table_value: name of table not found:" , tbl_s);
}

#endif
