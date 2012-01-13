#include "madx.h"

#ifndef _WIN32
#include <sys/utsname.h> // for uname
#endif

// private functions

#if 0 // not used...
static int
table_org(char* table)
  /* returns origin: 0  this job, 1 read or unknown */
{
  int pos;
  int org = 1;
  mycpy(c_dum->c, table);
  if ((pos = name_list_pos(c_dum->c, table_register->names)) > -1)
    org = table_register->tables[pos]->origin;
  return org;
}
#endif

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

static int
tab_name_code(char* name, char* t_name)
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
table_row(struct table* table, char* name)
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
  char rout_name[] = "grow_table_list";
  struct table** t_loc = tl->tables;
  int j, new = 2*tl->max;

  grow_name_list(tl->names);
  tl->max = new;
  tl->tables = mycalloc(rout_name,new, sizeof(struct table*));
  for (j = 0; j < tl->curr; j++) tl->tables[j] = t_loc[j];
  myfree(rout_name, t_loc);
}

static void
grow_table_list_list(struct table_list_list* tll)
{
  char rout_name[] = "grow_table_list_list";
  struct table_list** t_loc = tll->table_lists;
  int j, new = 2*tll->max;

  tll->max = new;
  tll->table_lists = mycalloc(rout_name,new, sizeof(struct table_list*));
  for (j = 0; j < tll->curr; j++) tll->table_lists[j] = t_loc[j];
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
write_table(struct table* t, char* filename)
  /* writes rows with columns listed in row and col */
{
  char l_name[NAME_L];
  char sys_name[200], t_pc[2*NAME_L];
  char* pc = t_pc;
  struct int_array* col = t->col_out;
  struct int_array* row = t->row_out;
  int i, j, k, tmp, n;
  time_t now;
  struct tm* tm;
#ifndef _WIN32
  struct utsname u;
  i = uname(&u); /* get system name */
  strcpy(sys_name, u.sysname);
#else // _WIN32
  strcpy(sys_name, "Win32");
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

    n = strlen(version_name)+strlen(sys_name)+1;
    fprintf(out_file,
            "@ ORIGIN           %%%02ds \"%s %s\"\n",
            n, version_name, sys_name);

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
/*          printf("row %d col %d datatype %d \n",j,i, t->columns->inform[col->i[i]] );*/
          if (t->columns->inform[col->i[i]] == 1)
          {
            tmp = t->d_cols[col->i[i]][j];
            fprintf(out_file, v_format(" %I"), tmp);
          }
          else
            if (t->columns->inform[col->i[i]] == 2)
            {
              fprintf(out_file, v_format(" %F"), t->d_cols[col->i[i]][j]);
              /*printf("%s[%2d,%2d]=%+8.5f    ",t->name,col->i[i],j,t->d_cols[col->i[i]][j]);*/
            }
            else if (t->columns->inform[col->i[i]] == 3)
            {
              pc[0] = c_dum->c[0] = '\"';
              if (t->s_cols[col->i[i]][j] != NULL)
              {
                strcpy(&c_dum->c[1], t->s_cols[col->i[i]][j]);
                stoupper(c_dum->c);
                pc = strip(c_dum->c); /* remove :<occ_count> */
                k = strlen(pc);
              }
              else k = 1;
              pc[k++] = '\"'; pc[k] = '\0';
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
static void
string_to_table_row(char* table, char* name, int* row, char* string)
  /* puts string at current position in column with name "name".
     The table count is increased separately with "augment_count" */
{
  int pos;
  struct table* t;

  mycpy(c_dum->c, table);
  if ((pos = name_list_pos(c_dum->c, table_register->names)) > -1)
    t = table_register->tables[pos];
  else return;
  mycpy(c_dum->c, name);
  if ((pos = name_list_pos(c_dum->c, t->columns)) >= 0
      && t->columns->inform[pos] == 3)
  {
    mycpy(c_dum->c, string);
    t->s_cols[pos][*row-1] = tmpbuff(c_dum->c);
  }
}
#endif

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
    k = string_from_table("normal_results","name", &row, string);
    if (k != 0) return k;
    if (strcmp(string,n_var) == 0)
    {
      found = 1;
      k = double_from_table("normal_results","order1", &row, &d_val);
      if ((int)d_val != order[0]) found = 0;
      k = double_from_table("normal_results","order2", &row, &d_val);
      if ((int)d_val != order[1]) found = 0;
      k = double_from_table("normal_results","order3", &row, &d_val);
      if ((int)d_val != order[2]) found = 0;
      k = double_from_table("normal_results","order4", &row, &d_val);
      if ((int)d_val != order[3]) found = 0;
    }
    if (found == 1) break;
  }
  if (found == 1)
    k = double_from_table("normal_results","value", &row, &d_val);
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
            warning("incomplete table line starting with:", aux_buff->c);
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
  c_range_start = get_node_count(current_sequ->range_start);
  c_range_end = get_node_count(current_sequ->range_end);
  get_select_t_ranges(select, deselect, t);
  if (select != 0)
  {
    for (j = 0; j < t->curr; j++)  t->row_out->i[j] = 0;
    for (i = 0; i < select->curr; i++)
    {
      for (j = s_range->i[i]; j <= e_range->i[i]; j++)
      {
        if (t->row_out->i[j] == 0) t->row_out->i[j]
                                     = pass_select(t->s_cols[0][j], select->commands[i]);
      }
    }
  }
  if (deselect != NULL)
  {
    for (i = 0; i < deselect->curr; i++)
    {
      for (j = sd_range->i[i]; j <= ed_range->i[i]; j++)
      {
        if (t->row_out->i[j] == 1) t->row_out->i[j]
                                     = 1 - pass_select(t->s_cols[0][j], deselect->commands[i]);
      }
    }
  }
}

// public interface

struct table*
new_table(char* name, char* type, int rows, struct name_list* cols)
{
  char rout_name[] = "new_table";
  int i, n = cols->curr;
  struct table* t
    = (struct table*) mycalloc(rout_name,1, sizeof(struct table));

  strcpy(t->name, name);
  strcpy(t->type, type);
  t->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", "table");
  t->columns = cols;
  t->num_cols = t->org_cols = n;
  t->s_cols = (char***) mycalloc(rout_name,n, sizeof(char**));
  t->d_cols = (double**) mycalloc(rout_name,n, sizeof(double*));
  t->max = ++rows; /* +1 because of separate augment_count */
  for (i = 0; i < n; i++)
  {
    if (cols->inform[i] < 3)
      t->d_cols[i] = (double*) mycalloc(rout_name,rows, sizeof(double));
    else if (cols->inform[i] == 3)
      t->s_cols[i] = (char**) mycalloc(rout_name,rows, sizeof(char*));
  }
  t->row_out = new_int_array(rows);
  t->col_out = new_int_array(n);
  t->node_nm = new_char_p_array(rows);
  t->p_nodes = (struct node**) mycalloc(rout_name,rows, sizeof(struct nodes*));
  t->l_head
    = (struct char_p_array**)
    mycalloc(rout_name,rows, sizeof(struct char_p_array*));
  return t;
}

struct table_list*
new_table_list(int size)
{
  char rout_name[] = "new_table_list";
  struct table_list* tl
    = (struct table_list*) mycalloc(rout_name,1, sizeof(struct table_list));
  strcpy(tl->name, "table_list");
  tl->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", tl->name);
  tl->max = size;
  tl->curr = 0;
  tl->names = new_name_list(tl->name, size);
  tl->tables
    = (struct table**) mycalloc(rout_name,size, sizeof(struct table*));
  add_to_table_list_list(tl, all_table_lists);
  return tl;
}

struct table_list_list*
new_table_list_list(int size)
{
  char rout_name[] = "new_table_list_list";
  struct table_list_list* tll
    = (struct table_list_list*) 
    mycalloc(rout_name,1, sizeof(struct table_list_list));
  strcpy(tll->name, "table_list_list");
  tll->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", tll->name);
  tll->max = size;
  tll->curr = 0;
  tll->table_lists
    = (struct table_list**) 
    mycalloc(rout_name,size, sizeof(struct table_list*));
  return tll;
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
    if (ntok > 1) {
      if ((pos = name_list_pos(toks[0], table_register->names)) > -1) {
        table = table_register->tables[pos];
        if ((col = name_list_pos(toks[ntok-1], table->columns)) > -1) {
          if (ntok > 2) { /* find row - else current (dynamic), or 0 */
            /* start mod - HG 26.3.2011 */
            if (ntok > 5) { /* check for [ count ] and convert to ->count */
              if (*toks[2] == '[' && *toks[4] == ']') {
                strcat(toks[1], "->");
                strcat(toks[1], toks[3]);
              }
	          }
	          /* end mod - HG 26.3.2011 */
            row = table_row(table, toks[1]);
          }
          else if (table->dynamic) row = table->curr;
          else row = 0;
          val = table->d_cols[col][row];
        }
        else if ((ntok == 3) && ((col = name_list_pos(toks[1], table->columns)) > -1)) {
          row = atoi(toks[2])-1;
          if(row < table->curr) val = table->d_cols[col][row];
        }
        else if(ntok == 2) {
          strncpy(temp, toks[1], NAME_L);
          if (strcmp(stolower(temp), "tablelength") == 0) val = table->curr;
	      }
      }
    }
  }
  return val;
}

void
augment_count(char* table) /* increase table occ. by 1, fill missing */
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

  if (t->num_cols > t->org_cols)  add_vars_to_table(t);

  if (t->p_nodes != NULL) t->p_nodes[t->curr] = current_node;

  if (t->node_nm != NULL)
  {
    t->node_nm->p[t->curr] = current_node->name;
    t->node_nm->curr = t->curr;
  }
  if (++t->curr == t->max) grow_table(t);
}

void
augmentcountonly(char* table) /* increase table occ. by 1 */
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

  if (t->num_cols > t->org_cols)  add_vars_to_table(t);

  if (++t->curr == t->max) grow_table(t);
}

int
char_from_table(char* table, char* name, int* row, char* val)
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
  strncpy(val,t->node_nm->p[*row-1],NAME_L);
  while (strlen(val)<=NAME_L) val[strlen(val)]=' ';
  return 0;
}

void
comment_to_table(char* table, char* comment, int* length)
  /* Saves the comment string at the current line.
     This comment is then printed in front of this line.
     Several calls to the same current line are possible. */
{
  int pos;
  struct table* t;
  mycpy(c_dum->c, table);
  if ((pos = name_list_pos(c_dum->c, table_register->names)) > -1)
    t = table_register->tables[pos];
  else return;
  strncpy(c_dum->c, comment, *length); c_dum->c[*length] = '\0';
  if (t->l_head[t->curr] == NULL)
    t->l_head[t->curr] = new_char_p_array(2);
  else if (t->l_head[t->curr]->curr == t->l_head[t->curr]->max)
    grow_char_p_array(t->l_head[t->curr]);
  t->l_head[t->curr]->p[t->l_head[t->curr]->curr++] = tmpbuff(c_dum->c);
}

void
add_to_table_list(struct table* t, struct table_list* tl)
  /* adds table t to table list tl */
{
  int pos, j;
  if ((pos = name_list_pos(t->name, tl->names)) < 0)
  {
    if (tl->curr == tl->max) grow_table_list(tl);
    j = add_to_name_list(tmpbuff(t->name), 0, tl->names);
    tl->tables[tl->curr++] = t;
  }
  else
  {
    tl->tables[pos] = delete_table(tl->tables[pos]);
    tl->tables[pos] = t;
  }
}

void
add_vars_to_table(struct table* t)
  /* fills user-defined variables into current table_row) */
{
  int i;
  char* p;

  for (i = t->org_cols; i < t->num_cols; i++)
  {
    if (t->columns->inform[i] < 3)
    {
      if (strstr(t->columns->names[i], "aper_"))
        t->d_cols[i][t->curr]
          = get_aperture(current_node, t->columns->names[i]);
      else if (strstr(t->columns->names[i], "aptol_"))
        t->d_cols[i][t->curr]
          = get_apertol(current_node, t->columns->names[i]);
      else t->d_cols[i][t->curr] = get_variable(t->columns->names[i]);
    }
    else if (current_node)
    {
      if ((p = command_par_string(t->columns->names[i],
                                  current_node->p_elem->def)) == NULL)
        t->s_cols[i][t->curr] = tmpbuff("none");
      else t->s_cols[i][t->curr] = tmpbuff(p);
    }
    else t->s_cols[i][t->curr] = get_varstring(t->columns->names[i]);
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
  char rout_name[] = "delete_table";
  int i, j;
  if (t == NULL) return NULL;
  if (stamp_flag && t->stamp != 123456)
    fprintf(stamp_file, "d_t double delete --> %s\n", t->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", "table");
  if (t->header != NULL) t->header = delete_char_p_array(t->header, 1);
  if (t->col_out != NULL) t->col_out = delete_int_array(t->col_out);
  if (t->row_out != NULL) t->row_out = delete_int_array(t->row_out);
  if (t->node_nm != NULL) t->node_nm = delete_char_p_array(t->node_nm, 0);
  for (i = 0; i < t->curr; i++)
  {
    if (t->l_head[i] != NULL)
      t->l_head[i] = delete_char_p_array(t->l_head[i], 1);
  }
  if (t->l_head)  myfree(rout_name, t->l_head);
  if (t->p_nodes) myfree(rout_name, t->p_nodes);
  if (t->d_cols)
  {
    for (i = 0; i < t->num_cols; i++)
      if (t->columns->inform[i] < 3 && t->d_cols[i]) myfree(rout_name, t->d_cols[i]);
    myfree(rout_name, t->d_cols);
  }
  if (t->s_cols)
  {
    for (i = 0; i < t->num_cols; i++)
    {
      if (t->columns->inform[i] == 3 && t->s_cols[i])
      {
        for (j = 0; j < t->curr; j++)
          if (t->s_cols[i][j]) myfree(rout_name, t->s_cols[i][j]);
        myfree(rout_name, t->s_cols[i]);
      }
    }
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
  char rout_name[] = "grow_table";
  int i, j, new = 2*t->max;
  char** s_loc;
  struct char_p_array* t_loc = t->node_nm;
  double* d_loc;
  struct node** p_loc = t->p_nodes;
  struct char_p_array** pa_loc = t->l_head;

  t->max = new;
  t->p_nodes = (struct node**) mycalloc(rout_name,new, sizeof(struct node*));
  t->l_head
    = (struct char_p_array**)
    mycalloc(rout_name,new, sizeof(struct char_p_array*));
  t->node_nm = new_char_p_array(new);

  for (i = 0; i < t->curr; i++)
  {
    t->node_nm->p[i] = t_loc->p[i];
    t->p_nodes[i] = p_loc[i];
    t->l_head[i] = pa_loc[i];
  }
  delete_char_p_array(t_loc, 0);
  myfree(rout_name, pa_loc);
  t->node_nm->curr = t->curr; myfree(rout_name, p_loc);
  for (j = 0; j < t->num_cols; j++)
  {
    if ((s_loc = t->s_cols[j]) != NULL)
    {
      t->s_cols[j] = (char**) mycalloc(rout_name,new, sizeof(char*));
      for (i = 0; i < t->curr; i++) t->s_cols[j][i] = s_loc[i];
      myfree(rout_name, s_loc);
    }
  }
  for (j = 0; j < t->num_cols; j++)
  {
    if ((d_loc = t->d_cols[j]) != NULL)
    {
      t->d_cols[j] = (double*) mycalloc(rout_name,new, sizeof(double));
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
double_to_table(char* table, char* name, double* val)
  /* puts val at current position in column with name "name".
     The table count is increased separately with "augment_count" */
{
  int pos;
  struct table* t;

  /*  printf("double_to_table <%s> <%s> <%f>\n",table,name,*val);*/

  mycpy(c_dum->c, table);
  if ((pos = name_list_pos(c_dum->c, table_register->names)) > -1)
  {
    t = table_register->tables[pos];
  }
  else
  {
    printf("Can not find table %s\n",table);
    return;
  }
  mycpy(c_dum->c, name);
  if ((pos = name_list_pos(c_dum->c, t->columns)) >= 0 && t->columns->inform[pos] < 3)
  {
    t->d_cols[pos][t->curr] = *val;
  }
  else
  {
    printf("Position of column %s is %d\n",name,pos);
  }
}

void
double_to_table_row(char* table, char* name, int* row, double* val)
  /* puts val at row position in column with name "name".
     The table count is increased separately with "augment_count" */
{
  int pos;
  struct table* t;

  mycpy(c_dum->c, table);
  if ((pos = name_list_pos(c_dum->c, table_register->names)) > -1)
    t = table_register->tables[pos];
  else return;
  mycpy(c_dum->c, name);
  if ((pos = name_list_pos(c_dum->c, t->columns)) >= 0
      && t->columns->inform[pos] < 3) t->d_cols[pos][*row-1] = *val;
}

int
double_from_table(char* table, char* name, int* row, double* val)
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

  *val = zero;
  mycpy(c_dum->c, table);
  if ((pos = name_list_pos(c_dum->c, table_register->names)) > -1)
    t = table_register->tables[pos];
  else return -1;
  mycpy(c_dum->c, name);
  if ((pos = name_list_pos(c_dum->c, t->columns)) < 0) return -2;
  if (*row > t->curr)  return -3;
  *val = t->d_cols[pos][*row-1];
  return 0;
}

int
string_from_table(char* table, char* name, int* row, char* string)
  /* returns val at position row in column with name "name".
     function value return:
     0  OK
     -1 table  does not exist
     -2 column does not exist
     -3 row    does not exist
     struct command_parameter* cp;
     struct double_array* arr = NULL;
  */
{
  int pos,l;
  struct table* t;

  mycpy(c_dum->c, table);
  if ((pos = name_list_pos(c_dum->c, table_register->names)) > -1)
    t = table_register->tables[pos];
  else return -1;
  mycpy(c_dum->c, name);
  if ((pos = name_list_pos(c_dum->c, t->columns)) < 0) return -2;
  if (*row > t->curr)  return -3;
  l = strlen(t->s_cols[pos][*row-1]);
  mycpy(string, t->s_cols[pos][*row-1]);
  return 0;
}

void
make_map_table(int* map_table_max_rows)
{
  int k, pos;
  if ((pos = name_list_pos("map_table", table_register->names)) > -1)
  {
    delete_table(table_register->tables[pos]);
    k = remove_from_name_list(table_register->tables[pos]->name,
                              table_register->names);
    table_register->tables[k] = table_register->tables[--table_register->curr];
  }
  /* initialise table */
  map_table = make_table("map_table", "map_tab", map_tab_cols,
                         map_tab_types, *map_table_max_rows);
  add_to_table_list(map_table, table_register);
  map_table->dynamic = 1;
  reset_count("map_table");
}

struct table*
make_table(char* name, char* type, char** table_cols, int* table_types, int rows)
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

int
get_table_range(char* range, struct table* table, int* rows)
  /* returns start and end row (rows[0] and rows[1])
     of a range in a table; 0 if not found, 1 (1 row) or 2 ( > 1) */
{
  int i, n;
  char* c[2];
  char tmp[NAME_L], dumtex[3*NAME_L];;
  rows[0] = rows[1] = 0;
  mycpy(c_dum->c, range); stolower(c_dum->c); strcpy(dumtex, c_dum->c);
  c[0] = strtok(c_dum->c, "/");
  if ((c[1] = strtok(NULL,"/")) == NULL) /* only one element given */
    n = 1;
  else n = 2;
  for (i = 0; i < n; i++)
  {
    if (*c[i] == '#')
    {
      if (strncmp(c[i], "#s", 2) == 0) rows[i] = 0;
      else if (strncmp(c[i], "#e", 2) == 0) rows[i] = table->curr - 1;
      else
      {
        warning("illegal table range ignored:", dumtex);
        return 0;
      }
    }
    else
    {
      strcpy(tmp, c[i]);
      if (square_to_colon(tmp) == 0)
      {
        warning("illegal table range ignored:", dumtex);
        return 0;
      }
      if ((rows[i] = char_p_pos(tmp, table->node_nm)) < 0)
      {
        warning("illegal table range ignored:", dumtex);
        return 0;
      }
    }
  }
  if (n == 1) rows[1] = rows[0];
  return n;
}

void
reset_count(char* table) /* resets table counter to zero */
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

  string_to_table( sector_table_name, "name", elementName );
  double_to_table( sector_table_name, "pos", pos );

  /* 6 kicks */
  for (i=0; i<6; i++){
    char kickStr[2+1];
    sprintf(kickStr,"k%i",i+1);
    kickStr[2]='\0';
    index = i;
    double_to_table( sector_table_name, kickStr,&kick[index]);
  }
  /* 36 R-matrix terms */
  for (j=0; j<6; j++){
    for (i=0; i<6; i++){
      char rStr[3+1];
      sprintf(rStr,"r%i%i",i+1,j+1);
      rStr[3]='\0';
      index = i+j*6;
      double_to_table( sector_table_name, rStr,&rmatrix[index]);
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
	double_to_table( sector_table_name, tStr,&tmatrix[index]);
      }
    }
  }

  augment_count( sector_table_name ); /* move to next record */
}

void
headvalue(char* table_name, char* par, double* value)
/* returns the value of header parameter par from table table_name if present,
   else 10^12 */
{
  int i, pos;
  char lpar[NAME_L], ltab[NAME_L];
  char* tp;
  struct table* tab;
  *value = ten_p_12;
  mycpy(ltab, table_name);
  stolower(ltab);
  if ((pos = name_list_pos(ltab, table_register->names)) > -1)
  {
    tab = table_register->tables[pos];
    mycpy(lpar, par);
    if (tab->header)
    {
      for (i = 0; i < tab->header->curr; i++)
      {
        strcpy(aux_buff->c, &tab->header->p[i][1]);
        if ((tp =strtok(aux_buff->c, " \"\n")) &&
            compare_no_case(tp, lpar) == 0)
        {
          if (strstr(strtok(NULL, " \"\n"), "%le") != NULL)
          {
            sscanf(strtok(NULL, " \"\n"), "%le", value);
            break;
          }
        }
      }
    }
  }
  return;
}

void
out_table(char* tname, struct table* t, char* filename)
  /* output of a table */
{
  int j;

  struct command_list* scl = find_command_list(tname, table_select);
  struct command_list* dscl = find_command_list(tname, table_deselect);
  while (t->num_cols > t->col_out->max)
    grow_int_array(t->col_out);
  while (t->curr > t->row_out->max)
    grow_int_array(t->row_out);
  t->row_out->curr = t->curr;
  if (par_present("full", NULL, scl))
    put_info("obsolete option 'full'"," ignored on 'select'");
  for (j = 0; j < t->curr; j++) t->row_out->i[j] = 1;
  for (j = 0; j < t->num_cols; j++) t->col_out->i[j] = j;
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
            warning("incomplete table line starting with:", aux_buff->c);
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
      t->header->p[t->header->curr]
        = (char*) mymalloc("read_table", strlen(aux_buff->c)+1);
      strcpy(t->header->p[t->header->curr], aux_buff->c);
      t->header->curr++;
    }
  }
  fclose(tab_file);
  t->origin = 1;
  add_to_table_list(t, table_register);
  return NULL;
}

void
string_to_table(char* table, char* name, char* string)
  /* buffers + puts "string"
     at current position in column with name "name".
     The table count is increased separately with "augment_count" */
{
  int pos;
  struct table* t;

  mycpy(c_dum->c, table);
  if ((pos = name_list_pos(c_dum->c, table_register->names)) > -1)
    t = table_register->tables[pos];
  else return;
  mycpy(c_dum->c, name);
  if ((pos = name_list_pos(c_dum->c, t->columns)) >= 0
      && t->columns->inform[pos] == 3)
  {
    mycpy(c_dum->c, string);
    if (strcmp(c_dum->c, "name") == 0)
      t->s_cols[pos][t->curr] = tmpbuff(current_node->name);
    else t->s_cols[pos][t->curr] = tmpbuff(c_dum->c);
  }
}

int
table_length(char* table)
  /* returns no. of rows in table */
{
  int pos;
  int length = 0;
  mycpy(c_dum->c, table);
  if ((pos = name_list_pos(c_dum->c, table_register->names)) > -1)
    length = table_register->tables[pos]->curr;
  return length;
}

void
table_range(char* table, char* range, int* rows)
  /* returns first and last row numbers (start=1) in rows
     or 0 if table or range invalid */
{
  int pos;
  struct table* t;

  rows[0] = rows[1] = 0;
  mycpy(c_dum->c, table);
  if ((pos = name_list_pos(c_dum->c, table_register->names)) > -1)
  {
    t = table_register->tables[pos];
    get_table_range(range, t, rows);
    rows[0]++; rows[1]++;
  }
}

void
vector_to_table(char* table, char* col, int* nval, double* vals)
  /* puts nval values of array vals at the current line into columns
     starting with column whose name is in "col";
     The table count is increased separately with "augment_count" */
{
  int j, pos, c_pos, last = 0;
  struct table* t;

  mycpy(c_dum->c, table);
  if ((pos = name_list_pos(c_dum->c, table_register->names)) > -1)
    t = table_register->tables[pos];
  else return;
  mycpy(c_dum->c, col);
  if ((c_pos = name_list_pos(c_dum->c, t->columns)) > -1)
    last = mymin(c_pos + *nval, t->num_cols);
  for (j = c_pos; j < last; j++)
    if (t->columns->inform[j] < 3) t->d_cols[j][t->curr] = vals[j-c_pos];
}

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
  while (strlen(val)<NAME_L) val[strlen(val)]=' ';
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
  while (strlen(val)<NAME_L) val[strlen(val)]=' ';
  val[NAME_L-1] = '\0';
  return 0;
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
           if (strcmp(tmp,"%s") == 0) t->s_cols[i][t->curr] = stolower(tmpbuff(cc));
           else if (strcmp(tmp,"%d") == 0 )
           {
            sscanf(cc, tmp, &k); t->d_cols[i][t->curr] = k;
           }
           else if (strcmp(tmp,"%hd") == 0 )
           {
            sscanf(cc, tmp, &sk); t->d_cols[i][t->curr] = sk;
           }
           else sscanf(cc, tmp, &t->d_cols[i][t->curr]);
           if (i+1 < tnl->curr)
           {
              if ((cc =strtok(NULL, " \"\n")) == NULL)
              {
               warning("incomplete table line starting with:", aux_buff->c);
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

