/* FS & TdA 15.03.2004: fix missing myfree of p_loc in routine grow_table */
/* TdA 15.03.2004: 2 new routines for plotting "interp_node" and "reset_interpolation" */
#ifdef _WRAP_FORTRAN_CALLS
#include "fortran_wrappers.h"
#endif

#ifdef _WIN32
#ifndef _UINTPTR_T_
#define _UINTPTR_T_
#define uintptr_t unsigned int	/* 32 bytes-long (should be 64 on WIN64) */
#endif
#else
#include <stdint.h>		/* uintptr_t, to fit pointers into integers of correct size */
#endif

int add_drifts(struct node* c_node, struct node* end)
{
  struct node *d1;
  struct element *drift;
  double pos, dl, el2;
  int cnt = 0;
  pos = c_node->position
    - c_node->length / two;
  while (c_node != NULL)
  {
    cnt++;
    el2 = c_node->length / two;
    dl = c_node->position - el2 - pos;
    if (dl + ten_m_6 < zero)
    {
      sprintf(c_dum->c, "%s, length %e", c_node->name, dl);
      fatal_error("negative drift in front of", c_dum->c);
    }
    else if (dl > ten_m_6)
    {
      cnt++;
      drift = get_drift(dl);
      d1 = new_elem_node(drift, 0);
      link_in_front(d1, c_node);
      d1->position = pos + dl / two;
    }
    pos = c_node->position + el2;
    if (c_node == end) break;
    c_node = c_node->next;
  }
  return cnt;
}

void add_table_vars(struct name_list* cols, struct command_list* select)
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

void add_to_command_list(char* label, struct command* comm,
                         struct command_list* cl, int flag)
  /* adds command comm to the command list cl */
  /* flag for printing a warning */
{
  int pos, j;
  if ((pos = name_list_pos(label, cl->list)) > -1)
  {
    if (flag) put_info(label, "redefined");
    if (cl != defined_commands && cl != stored_commands)
      delete_command(cl->commands[pos]);
    cl->commands[pos] = comm;
  }
  else
  {
    if (cl->curr == cl->max) grow_command_list(cl);
    j = add_to_name_list(permbuff(label), 0, cl->list);
    cl->commands[cl->curr++] = comm;
  }
}

void add_to_command_list_list(char* label, struct command_list* cl,
                              struct command_list_list* sl)
  /* adds command list cl to command-list list sl */
{
  int pos, j;
  if ((pos = name_list_pos(label, sl->list)) > -1)
  {
    delete_command_list(sl->command_lists[pos]);
    sl->command_lists[pos] = cl;
  }
  else
  {
    if (sl->curr == sl->max) grow_command_list_list(sl);
    j = add_to_name_list(permbuff(label), 0, sl->list);
    sl->command_lists[sl->curr++] = cl;
  }
}

void add_to_constraint_list(struct constraint* cs, struct constraint_list* cl)
  /* add constraint cs to the constraint list cl */
{
  if (cl->curr == cl->max) grow_constraint_list(cl);
  cl->constraints[cl->curr++] = cs;
}

void add_to_el_list( /* adds element to alphabetic element list */
  struct element** el, int inf, struct el_list* ell, int flag)
  /* inf is entered in the namelist */
  /*  flag < 0: do not delete if already present, do not warn */
  /*       = 0: delete, but do not warn */
  /*       = 1: delete & warn */
  /*       = 2: warn and ignore if already present - resets *el to old */
{
  int pos, j;
  struct node* p_node;
  if ((pos = name_list_pos((*el)->name, ell->list)) > -1)
  {
    if (flag > 1)
    {
      warning("implicit element re-definition ignored:", (*el)->name);
      *el = ell->elem[pos];
    }
    else
    {
      if (flag > 0)
      {
        put_info("element redefined:", (*el)->name);
/*
 *         printf("File %s line %d\n",filenames[in->curr], currentline[in->curr] );
 *         printf("Old Definition:\n");
 *         dump_element(ell->elem[pos]);
 *         printf("New Definition:\n");
 *         dump_element(*el);
 */

      }
      if (flag >= 0 && ell == element_list)
      {
        for (j = 0; j < ell->curr; j++) /* set existing pointers to new */
        {
          if (ell->elem[j] != ell->elem[pos]
              && ell->elem[j]->parent == ell->elem[pos])
            ell->elem[j]->parent = *el;
        }
        for (j = 0; j < sequences->curr; j++)
        {
          p_node = sequences->sequs[j]->start;
          while (p_node && p_node != sequences->sequs[j]->end)
          {
            if (p_node->p_elem == ell->elem[pos]) p_node->p_elem = *el;
            p_node = p_node->next;
          }
          if (strcmp((*el)->base_type->name, "rfcavity") == 0 &&
	      find_element((*el)->name, sequences->sequs[j]->cavities) != NULL)
	    sequences->sequs[j]->cavities->elem[name_list_pos((*el)->name,
	    sequences->sequs[j]->cavities->list)] = *el;
        }
        delete_element(ell->elem[pos]);
      }
      ell->elem[pos] = *el;
    }
  }
  else
  {
    if (ell->curr == ell->max) grow_el_list(ell);
    j = add_to_name_list(permbuff((*el)->name), inf, ell->list);
    ell->elem[ell->curr++] = *el;
  }
}

void add_to_macro_list( /* adds macro to alphabetic macro list */
  struct macro* macro, struct macro_list* nll)
{
  int pos, j;
  if ((pos = name_list_pos(macro->name, nll->list)) > -1)
  {
    warning("macro redefined:", macro->name);
    delete_macro(nll->macros[pos]);
    nll->macros[pos] = macro;
  }
  else
  {
    if (nll->curr == nll->max) grow_macro_list(nll);
    j = add_to_name_list(permbuff(macro->name), 0, nll->list);
    nll->macros[nll->curr++] = macro;
  }
  /* RDM new matching*/
  if (match_is_on==2)
  {

    for(j=0; j < MAX_MATCH_MACRO;j++)
    {
      if (match2_macro_name[j]==NULL)
      {
        break;
      }
    }

    if (j  >= MAX_MATCH_MACRO  )
    {
      printf("Max number of match macros reached. Augmenting.\n");
      match2_augmentnmacros();
      j = MAX_MATCH_MACRO - 1;
    }

    match2_macro_name[j]=macro->name;

  }
  /* RDM new matching*/
}

int add_to_name_list(char* name, int inf, struct name_list* vlist)
  /* adds name to alphabetic name list vlist */
  /* inf is an integer kept with name */
{
  int j, num, low = 0, mid, high = vlist->curr - 1, pos = 0, ret;

  if (name == NULL) return -1;

  ret = name_list_pos(name, vlist);
  if ( ret < 0)
  {
    while (low <= high)
    {
      mid = (low + high) / 2;
      if ((num = strcmp(name, vlist->names[vlist->index[mid]])) < 0)
      {
        high = mid - 1; pos = mid;
      }
      else if (num > 0) {
        low  = mid + 1; pos = low;
      }
    }
    ret = vlist->curr;
    if (vlist->curr == vlist->max) grow_name_list(vlist);
    for (j = vlist->curr; j > pos; j--) vlist->index[j] = vlist->index[j-1];
    vlist->index[pos] = vlist->curr;
    vlist->inform[vlist->curr] = inf;
    vlist->names[vlist->curr++] = name;
  }
  else  vlist->inform[ret] = inf;
  return ret;
}

void add_to_node_list( /* adds node to alphabetic node list */
  struct node* node, int inf, struct node_list* nll)
{
  int pos, j;
  if ((pos = name_list_pos(node->name, nll->list)) < 0)
  {
    if (nll->curr == nll->max) grow_node_list(nll);
    j = add_to_name_list(node->name, inf, nll->list);
    nll->nodes[nll->curr++] = node;
  }
}

void add_to_sequ_list(struct sequence* sequ, struct sequence_list* sql)
  /* adds sequence sequ to sequence list sql */
{
  int i;
  int firstfreeslot = -1;
  for (i = 0; i < sql->curr; i++) if (sql->sequs[i] == sequ)  return;

  /*printf("add_to_sequ_list %s \n",sequ->name);*/

  for (i = 0; i < sql->curr; i++)
  {
    if (sql->sequs[i] == 0x0)
     {
       /*printf("add_to_sequ_list: pos %d is NULL\n",i);*/
       firstfreeslot = i;
       continue;
     }
   
    if (strcmp(sql->sequs[i]->name, sequ->name) == 0)
    {
      /*printf("add_to_sequ_list sequence with this name is already in: %s \n",sequ->name);*/
      sql->sequs[i] = sequ;
      sql->list->names[i] = sequ->name;
      return;
    }
  }
  
  if (firstfreeslot >= 0)
   {/*This protects agains problems sequence redefinition*/
     /*printf("add_to_sequ_list: adding at found free slot\n");*/
     sql->sequs[firstfreeslot] = sequ;
   }
  else
   { 
     /*printf("add_to_sequ_list: adding at new slot\n");*/
     if (sql->curr == sql->max) grow_sequence_list(sql);
     sql->sequs[sql->curr++] = sequ;
   }  

  add_to_name_list(sequ->name, 0, sql->list);
  
}

void add_to_table_list(struct table* t, struct table_list* tl)
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

void add_to_table_list_list(struct table_list* table_list,
                            struct table_list_list* tll)
  /* adds a table_list to a list of table_lists */
{
  int j;
  for (j = 0; j < tll->curr; j++) 
     if (tll->table_lists[j] == table_list) return;
  if (tll->curr == tll->max) grow_table_list_list(tll);
  tll->table_lists[tll->curr++] = table_list;
}

void add_to_var_list( /* adds variable to alphabetic variable list */
  struct variable* var, struct var_list* varl, int flag)
  /* flag = 0: undefined reference in expression, 1: definition
     2: separate list, do not drop variable */
{
  int pos, j;

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
    j = add_to_name_list(permbuff(var->name), flag, varl->list);
    varl->vars[varl->curr++] = var;
  }
}

void add_vars_to_table(struct table* t)
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

void set_vars_from_table(struct table* t)
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

void set_variable(char* name, double* value)
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
    if ((el = find_element(comm, element_list)) != NULL)
      set_command_par_value(par, el->def, val);
    else if ((cmd = find_command(comm, stored_commands)) != NULL)
      set_command_par_value(par, cmd, val);
    else if ((cmd = find_command(comm, beta0_list)) != NULL)
      set_command_par_value(par, cmd, val);
    else if ((cmd = find_command(comm, defined_commands)) != NULL)
      set_command_par_value(par, cmd, val);
  }
}

void set_stringvar(char* name, char* string)
{
  /* sets variable name->string to string */
  char* p;
  struct variable* var;
  mycpy(c_dum->c, name);
  if ((p = strstr(c_dum->c, "->")) == NULL) /* variable */
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

int char_cnt(char c, char* string)
  /* returns number of occurrences of character c in string */
{
  int i, k = 0;
  for (i = 0; i < strlen(string); i++) if(string[i] == c) k++;
  return k;
}

int char_p_pos(char* name, struct char_p_array* p)
  /* returns the position of name in character pointer array p,
     or -1 if not found */
{
  int i;
  for (i = 0; i < p->curr; i++) if (strcmp(name, p->p[i]) == 0) return i;
  return -1;
}

struct char_p_array* clone_char_p_array(struct char_p_array* p)
{
  int i;
  struct char_p_array* clone = new_char_p_array(p->max);
  for (i = 0; i < p->curr; i++) clone->p[i] = permbuff(p->p[i]);
  clone->curr = p->curr;
  return clone;
}

struct command* clone_command(struct command* p)
{
  int i;
  struct command* clone = new_command(p->name, 0, p->par->curr,
                                      p->module, p->group, p->link_type,
                                      p->mad8_type);
  copy_name_list(clone->par_names, p->par_names);
  clone->par->curr = p->par->curr;
  for (i = 0; i < p->par->curr; i++)
    clone->par->parameters[i] =
      clone_command_parameter(p->par->parameters[i]);
  return clone;
}

struct command_parameter* clone_command_parameter(struct command_parameter* p)
{
  struct command_parameter* clone = new_command_parameter(p->name, p->type);
  clone->call_def = p->call_def;
  switch (p->type)
  {
    case 4:
      clone->c_min = p->c_min;
      clone->c_max = p->c_max;
      clone->min_expr = clone_expression(p->min_expr);
      clone->max_expr = clone_expression(p->max_expr);
    case 0:
    case 1:
    case 2:
      clone->double_value = p->double_value;
      clone->expr = clone_expression(p->expr);
      break;
    case 3:
      clone->string = p->string;
      clone->expr = NULL;
      break;
    case 11:
    case 12:
      clone->double_array = clone_double_array(p->double_array);
      clone->expr_list = clone_expr_list(p->expr_list);
      break;
    case 13:
      clone->m_string = clone_char_p_array(p->m_string);
  }
  return clone;
}

struct double_array* clone_double_array(struct double_array* p)
{
  int i;
  struct double_array* clone = new_double_array(p->curr);
  clone->curr = p->curr;
  for (i = 0; i < p->curr; i++) clone->a[i] = p->a[i];
  return clone;
}

struct element* clone_element(struct element* el)
{
  struct element* clone = new_element(el->name);
  clone->length = el->length;
  clone->bv = el->bv;
  clone->def = el->def;
  clone->parent = el;
  clone->base_type = el->base_type;
  return clone;
}

struct expression* clone_expression(struct expression* p)
{
  struct expression* clone;
  if (p == NULL) return NULL;
  clone = new_expression(p->string, p->polish);
  clone->status = p->status;
  clone->value = p->value;
  return clone;
}

struct expr_list* clone_expr_list(struct expr_list* p)
{
  int i;
  struct expr_list* clone;
  if (p == NULL)  return NULL;
  clone = new_expr_list(p->curr);
  for (i = 0; i < p->curr; i++) clone->list[i] = clone_expression(p->list[i]);
  clone->curr = p->curr;
  return clone;
}

struct int_array* clone_int_array(struct int_array* p)
{
  int i;
  struct int_array* clone = new_int_array(p->curr);
  clone->curr = p->curr;
  for (i = 0; i < p->curr; i++) clone->i[i] = p->i[i];
  return clone;
}

struct macro* clone_macro(struct macro* org)
{
  int i;
  struct macro* clone
    = new_macro(org->n_formal, org->body->curr, org->tokens->curr);
  if (org->body->curr > 0) strcpy(clone->body->c, org->body->c);
  clone->body->curr = org->body->curr;
  for (i = 0; i < org->tokens->curr; i++)
    clone->tokens->p[i] = org->tokens->p[i];
  clone->tokens->curr = org->tokens->curr;
  for (i = 0; i < org->n_formal; i++)
    clone->formal->p[i] = org->formal->p[i];
  clone->n_formal = org->n_formal;
  return clone;
}

struct name_list* clone_name_list(struct name_list* p)
{
  int i, l = p->curr > 0 ? p->curr : 1;
  char name[2*NAME_L];
  struct name_list* clone;
  strcpy(name, p->name); strcat(name, "_clone");
  clone = new_name_list(name, l);
  for (i = 0; i < p->curr; i++) clone->index[i] = p->index[i];
  for (i = 0; i < p->curr; i++) clone->inform[i] = p->inform[i];
  for (i = 0; i < p->curr; i++) clone->names[i] = p->names[i];
  clone->curr = p->curr;
  return clone;
}

struct var_list* clone_var_list(struct var_list* vl)
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

struct node* clone_node(struct node* p, int flag)
{
  /* Transfers errors from original nodes if flag != 0;
     this is needed for SXF input  */
  struct node* clone = new_node(p->name);
  strcpy(clone->name,p->name);
  clone->base_name = p->base_name;
  clone->occ_cnt = p->occ_cnt;
  clone->sel_err = p->sel_err;
  clone->position = p->position;
  clone->at_value = p->at_value;
  clone->length = p->length;
  clone->at_expr = p->at_expr;
  clone->from_name = p->from_name;
  clone->p_elem = p->p_elem;
  clone->p_sequ = p->p_sequ;
  clone->savebeta = p->savebeta;
  if (flag)
  {
    clone->p_al_err = p->p_al_err;
    clone->p_fd_err = p->p_fd_err;
  }
  return clone;
}
double combine_expr_expr(struct expression* exp1, char* oper, 
                         struct expression* exp2, struct expression** comb_exp)
{
  strcpy(c_dum->c, exp1->string);
  strcat(c_dum->c, oper);
  strcat(c_dum->c, exp2->string);
  mysplit(c_dum->c, tmp_p_array);
  *comb_exp = make_expression(tmp_p_array->curr, tmp_p_array->p);
  return expression_value(*comb_exp, 2);
}

double combine_expr_val(struct expression* exp1, char* oper, 
                        double val2, struct expression** comb_exp)
{
  strcpy(c_dum->c, exp1->string);
  sprintf(aux_buff->c, "%.12g", val2);
  strcat(c_dum->c, oper);
  strcat(c_dum->c, aux_buff->c);
  mysplit(c_dum->c, tmp_p_array);
  *comb_exp = make_expression(tmp_p_array->curr, tmp_p_array->p);
  return expression_value(*comb_exp, 2);
}

double combine_val_expr(double val1, char* oper, 
                        struct expression* exp2, struct expression** comb_exp)

{
  sprintf(c_dum->c, "%.12g", val1);
  strcat(c_dum->c, oper);
  strcat(c_dum->c, exp2->string);
  mysplit(c_dum->c, tmp_p_array);
  *comb_exp = make_expression(tmp_p_array->curr, tmp_p_array->p);
  return expression_value(*comb_exp, 2);
}

void conv_char(char* string, struct int_array* tint)
  /*converts character string to integer array, using ascii code */
{
  int i, l = strlen(string),
    n = (l < tint->max-1) ? l : tint->max-1;
  tint->i[0] = n;
  for (i = 0; i < n; i++)  tint->i[i+1] = (int) string[i];
}

void copy_double(double* source, double* target, int n)
  /* copies n double precision values from source to target */
{
  int j;
  for (j = 0; j < n; j++)  target[j] = source[j];
}

void copy_name_list(struct name_list* out, struct name_list* in)
  /* copies namelist in to namelist out */
{
  int i, l = in->curr > 0 ? in->curr : 1;
  while (out->max < l) grow_name_list(out);
  for (i = 0; i < in->curr; i++) out->index[i]  = in->index[i];
  for (i = 0; i < in->curr; i++) out->inform[i] = in->inform[i];
  for (i = 0; i < in->curr; i++) out->names[i]  = in->names[i];
  out->curr = in->curr;
}

struct char_array* delete_char_array(struct char_array* pa)
{
  char rout_name[] = "delete_char_array";
  if (pa == NULL)  return NULL;
  if (pa->c != NULL)  myfree(rout_name, pa->c);
  myfree(rout_name, pa);
  return NULL;
}

struct char_p_array* delete_char_p_array(struct char_p_array* pa, int flag)
  /* flag = 0: delete only pointer array, = 1: delete char. buffers, too */
{
  char rout_name[] = "delete_char_p_array";
  int i;
  if (pa == NULL)  return NULL;
  if (stamp_flag && pa->stamp != 123456)
    fprintf(stamp_file, "d_c_p_a double delete --> %s\n", pa->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", pa->name);
  if (flag)
  {
    for (i = 0; i < pa->curr; i++)  myfree(rout_name, pa->p[i]);
  }
  if (pa->p != NULL)  myfree(rout_name, pa->p);
  myfree(rout_name, pa);
  return NULL;
}

struct command* delete_command(struct command* cmd)
{
  char rout_name[] = "delete_command";
  if (cmd == NULL) return NULL;
  if (stamp_flag && cmd->stamp != 123456)
    fprintf(stamp_file, "d_c double delete --> %s\n", cmd->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", cmd->name);
  if (cmd->par != NULL)  delete_command_parameter_list(cmd->par);
  if (cmd->par_names != NULL) delete_name_list(cmd->par_names);
  myfree(rout_name, cmd);
  return NULL;
}

struct command_list* delete_command_list(struct command_list* cl)
{
  char rout_name[] = "delete_command_list";
  int i;
  if (cl == NULL) return NULL;
  if (stamp_flag && cl->stamp != 123456)
    fprintf(stamp_file, "d_c_l double delete --> %s\n", cl->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", cl->name);
  if (cl->list != NULL) delete_name_list(cl->list);
  for (i = 0; i < cl->curr; i++) delete_command(cl->commands[i]);
  if (cl->commands) myfree(rout_name, cl->commands);
  myfree(rout_name, cl);
  return NULL;
}

struct command_parameter*
delete_command_parameter(struct command_parameter* par)
{
  char rout_name[] = "delete_command_parameter";
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
  char rout_name[] = "delete_command_parameter_list";
  int i;
  if (parl == NULL) return NULL;
  if (stamp_flag && parl->stamp != 123456)
    fprintf(stamp_file, "d_c_p_l double delete --> %s\n", parl->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", parl->name);
  if (parl->parameters != NULL)
  {
    for (i = 0; i < parl->curr; i++)
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

int compare_no_case(char* string_1, char* string_2)
/* like strcmp, but ignoring case */
{
  int ret;
  char rout_name[] = "compare_no_case";
  char* s1 = mymalloc(rout_name, strlen(string_1)+1);
  char* s2 = mymalloc(rout_name, strlen(string_2)+1);
  strcpy(s1, string_1); stolower(s1);
  strcpy(s2, string_2); stolower(s2);
  ret = strcmp(s1, s2);
  myfree(rout_name, s1); myfree(rout_name, s2);
  return ret;
}

struct constraint* delete_constraint(struct constraint* cst)
{
  char rout_name[] = "delete_constraint";
  if (cst == NULL)  return NULL;
  if (stamp_flag && cst->stamp != 123456)
    fprintf(stamp_file, "d_c double delete --> %s\n", cst->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", "constraint");
  myfree(rout_name, cst);
  return NULL;
}

struct constraint_list* delete_constraint_list(struct constraint_list* cl)
{
  char rout_name[] = "delete_constraint_list";
  if (cl == NULL)  return NULL;
  if (stamp_flag && cl->stamp != 123456)
    fprintf(stamp_file, "d_c_l double delete --> %s\n", cl->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", "constraint_list");
  myfree(rout_name, cl);
  return NULL;
}

struct double_array* delete_double_array(struct double_array* a)
{
  char rout_name[] = "delete_double_array";
  if (a != NULL)
  {
    if (a->a != NULL) myfree(rout_name, a->a);
    myfree(rout_name, a);
  }
  return NULL;
}

struct element* delete_element(struct element* el)
{
  char rout_name[] = "delete_element";
  if (el == NULL)  return NULL;
  if (stamp_flag && el->stamp != 123456)
    fprintf(stamp_file, "d_e double delete --> %s\n", el->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", el->name);
  myfree(rout_name, el);
  return NULL;
}

struct el_list* delete_el_list(struct el_list* ell)
{
  char rout_name[] = "delete_el_list";
  if (ell->list == NULL) return NULL;
  if (stamp_flag && ell->stamp != 123456)
    fprintf(stamp_file, "d_e_l double delete --> %s\n", ell->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", ell->name);
  delete_name_list(ell->list);
  if (ell->elem != NULL) myfree(rout_name, ell->elem);
  myfree(rout_name, ell);
  return NULL;
}

struct expression* delete_expression(struct expression* expr)
{
  char rout_name[] = "delete_expression";
  if (expr == NULL) return NULL;
  if (stamp_flag && expr->stamp != 123456)
    fprintf(stamp_file, "d_ex double delete --> %s\n", expr->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", expr->name);
  if (expr->polish != NULL) expr->polish = delete_int_array(expr->polish);
  if (expr->string != NULL) myfree(rout_name, expr->string);
  myfree(rout_name, expr);
  return NULL;
}

struct expr_list* delete_expr_list(struct expr_list* exprl)
{
  char rout_name[] = "delete_expr_list";
  int i;
  if (exprl == NULL) return NULL;
  if (stamp_flag && exprl->stamp != 123456)
    fprintf(stamp_file, "d_ex_l double delete --> %s\n", exprl->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", exprl->name);
  if (exprl->list != NULL)
  {
    for (i = 0; i < exprl->curr; i++)
      if (exprl->list[i] != NULL)  delete_expression(exprl->list[i]);
    myfree(rout_name, exprl->list);
  }
  myfree(rout_name, exprl);
  return NULL;
}

struct in_cmd* delete_in_cmd(struct in_cmd* cmd)
{
  char rout_name[] = "delete_in_cmd";
  if (cmd == NULL) return NULL;
  if (stamp_flag && cmd->stamp != 123456)
    fprintf(stamp_file, "d_i_c double delete --> %s\n", cmd->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", cmd->name);
  if (cmd->tok_list != NULL)
    cmd->tok_list = delete_char_p_array(cmd->tok_list, 0);
  myfree(rout_name, cmd);
  return NULL;
}

struct int_array* delete_int_array(struct int_array* i)
{
  char rout_name[] = "delete_int_array";
  if (i == NULL)  return NULL;
  if (stamp_flag && i->stamp != 123456)
    fprintf(stamp_file, "d_i_a double delete --> %s\n", i->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", i->name);
  if (i->i != NULL) myfree(rout_name, i->i);
  myfree(rout_name, i);
  return NULL;
}

struct macro* delete_macro(struct macro* macro)
{
  char rout_name[] = "delete_macro";
  if (macro == NULL)  return NULL;
  if (stamp_flag && macro->stamp != 123456)
    fprintf(stamp_file, "d_m double delete --> %s\n", macro->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", macro->name);
  if (macro->formal != NULL) delete_char_p_array(macro->formal, 0);
  if (macro->tokens != NULL) delete_char_p_array(macro->tokens, 0);
  if (macro->body != NULL) delete_char_array(macro->body);
  myfree(rout_name, macro);
  return NULL;
}

struct name_list* delete_name_list(struct name_list* l)
{
  char rout_name[] = "delete_name_list";
  if (l == NULL) return NULL;
  if (stamp_flag && l->stamp != 123456)
    fprintf(stamp_file, "d_n_l double delete --> %s\n", l->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", l->name);
  if (l->index != NULL)  myfree(rout_name, l->index);
  if (l->inform != NULL)  myfree(rout_name, l->inform);
  if (l->names != NULL)   myfree(rout_name, l->names);
  myfree(rout_name, l);
  return NULL;
}

struct node* delete_node(struct node* p)
{
  char rout_name[] = "delete_node";
  if (p == NULL) return NULL;
  if (stamp_flag && p->stamp != 123456)
    fprintf(stamp_file, "d_n double delete --> %s\n", p->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", p->name);
  if (p->p_al_err) p->p_al_err = delete_double_array(p->p_al_err);
  if (p->p_fd_err) p->p_fd_err = delete_double_array(p->p_fd_err);
  myfree(rout_name, p);
  return NULL;
}

struct node* delete_node_ring(struct node* start)
{
  struct node *p, *q;
  if (start == NULL) return NULL;
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", "node_ring");
  q = start->next;
  while (q != NULL && q != start)
  {
    p = q; q = q->next; delete_node(p);
  }
  delete_node(start);
  return NULL;
}

struct node_list* delete_node_list(struct node_list* l)
{
  char rout_name[] = "delete_node_list";
  if (l == NULL)  return NULL;
  if (stamp_flag && l->stamp != 123456)
    fprintf(stamp_file, "d_no_l double delete --> %s\n", l->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", l->name);
  if (l->nodes != NULL)  myfree(rout_name, l->nodes);
  if (l->list != NULL)  delete_name_list(l->list);
  myfree(rout_name, l);
  return NULL;
}

struct sequence* delete_sequence(struct sequence* sequ)
{
  char rout_name[] = "delete_sequence";
  if (sequ->ex_start != NULL)
  {
    sequ->ex_nodes = delete_node_list(sequ->ex_nodes);
    sequ->ex_start = delete_node_ring(sequ->ex_start);
    sequ->orbits = delete_vector_list(sequ->orbits);
    myfree(rout_name, sequ->all_nodes);
  }
  if (sequ->l_expr) sequ->l_expr = delete_expression(sequ->l_expr);
  sequ->nodes = delete_node_list(sequ->nodes);
  sequ->start = delete_node_ring(sequ->start);
  if (sequ->cavities) sequ->cavities = delete_el_list(sequ->cavities);
  myfree(rout_name, sequ);
  return NULL;
}

struct sequence_list* delete_sequence_list(struct sequence_list* sql)
{
  char rout_name[] = "delete_sequence_list";
  if (sql == NULL) return NULL;
  if (stamp_flag && sql->stamp != 123456)
    fprintf(stamp_file, "d_s_l double delete --> %s\n", sql->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", sql->name);
  if (sql->list != NULL) delete_name_list(sql->list);
  if (sql->sequs != NULL) myfree(rout_name, sql->sequs);
  myfree(rout_name, sql);
  return NULL;
}

struct table* delete_table(struct table* t)
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

struct variable* delete_variable(struct variable* var)
{
  char rout_name[] = "delete_variable";
  if (var == NULL)  return NULL;
  if (stamp_flag && var->stamp != 123456)
    fprintf(stamp_file, "d_v double delete --> %s\n", var->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", var->name);
  if (var->expr != NULL) delete_expression(var->expr);
  if (var->string != NULL) myfree(rout_name, var->string);
  myfree(rout_name, var);
  return NULL;
}

struct var_list* delete_var_list(struct var_list* varl)
{
  char rout_name[] = "delete_var_list";
  if (varl == NULL) return NULL;
  if (stamp_flag && varl->stamp != 123456)
    fprintf(stamp_file, "d_v_l double delete --> %s\n", varl->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", varl->name);
  if (varl->list != NULL) delete_name_list(varl->list);
  if (varl->vars != NULL) myfree(rout_name, varl->vars);
  myfree(rout_name, varl);
  return NULL;
}

struct vector_list* delete_vector_list(struct vector_list* vector)
{
  char rout_name[] = "delete_vector_list";
  int j;
  if (vector == NULL) return NULL;
  if (vector->names != NULL)
  {
    for (j = 0; j < vector->names->curr; j++)
      if (vector->vectors[j]) delete_double_array(vector->vectors[j]);
    delete_name_list(vector->names);
  }
  if (vector->vectors != NULL) myfree(rout_name, vector->vectors);
  myfree(rout_name, vector);
  return NULL;
}

void disable_line( /* prevents line from further expansion by "use" */
  char* name, struct macro_list* nll)
{
  int pos;
  if ((pos = name_list_pos(name, nll->list)) > -1) nll->macros[pos]->dead = 1;
}

void dump_char_array(struct char_array* a)
{
  char* c = a->c;
  int n = 0, l_cnt = 60, k;
  while (n < a->curr)
  {
    k = a->curr - n; if (k > l_cnt) k = l_cnt;
    strncpy(c_dum->c, c, k);
    c += k; n += k;
    c_dum->c[k] = '\0';
    fprintf(prt_file, "%s\n", c_dum->c);
  }
}

void dump_char_p_array(struct char_p_array* p)
{
  int i;
  for (i = 0; i < p->curr; i++) fprintf(prt_file, "%s\n", p->p[i]);
}

void dump_command(struct command* cmd)
{
  int i;
  fprintf(prt_file, "command: %s  module: %s\n",
          cmd->name, cmd->module);
  for (i = 0; i < cmd->par->curr; i++)
    dump_command_parameter(cmd->par->parameters[i]);
}

void dump_command_parameter(struct command_parameter* par)
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

void dump_constraint(struct constraint* c)
{
  fprintf(prt_file,
          v_format("name: %s type: %I value: %F min: %F max: %F weight: %F\n"),
          c->name, c->type, c->value, c->c_min, c->c_max, c->weight);
}

void dump_constraint_list(struct constraint_list* cl)
{
  int i;
  for (i = 0; i < cl->curr; i++)
  {
    if (cl->constraints[i]) dump_constraint(cl->constraints[i]);
  }
}

void dump_element(struct element* el)
{
  fprintf(prt_file, v_format("+++ dumping element %S  parent %S\n"),
          el->name, el->parent->name);
  dump_command(el->def);
}

void dump_el_list(struct el_list* ell)
{
  int i;
  for (i = 0; i < ell->curr; i++) dump_element(ell->elem[i]);
}

void dump_expression(struct expression* ex)
{
  ex->value = expression_value(ex, 2);
  fprintf(prt_file, v_format("expression: %s :: value: %F\n"),
          ex->string, ex->value);
}

void dump_exp_sequ(struct sequence* sequ, int level)
  /* executes the command dumpsequ, ... */
{
  struct node* c_node;
  int j;
  double suml = zero;
  puts("+++++++++ dump expanded sequence +++++++++");
  c_node = sequ->ex_start;
  while(c_node != NULL)
  {
    suml += c_node->length;
    if (level > 2)
    {
      dump_node(c_node);
      if (c_node->p_al_err != NULL)
      {
        puts("alignment errors:");
        for (j = 0; j < c_node->p_al_err->curr; j++)
          printf(v_format("%F "), c_node->p_al_err->a[j]);
        printf("\n");
      }
      if (c_node->p_fd_err != NULL)
      {
        puts("field errors:");
        for (j = 0; j < c_node->p_fd_err->curr; j++)
          printf(v_format("%e "), c_node->p_fd_err->a[j]);
        printf("\n");
      }
      if (level > 3 && c_node->p_elem != NULL)  dump_element(c_node->p_elem);
    }
    else if (level > 0 && strcmp(c_node->base_name, "drift") != 0)
      fprintf(prt_file, v_format("%S: at = %F  flag = %I\n"), c_node->name,
              c_node->position, c_node->enable);
    if (c_node == sequ->ex_end)  break;
    c_node = c_node->next;
  }
  fprintf(prt_file, v_format("=== sum of node length: %F\n"), suml);
}

void dump_in_cmd(struct in_cmd* p_inp)
{
  fprintf(prt_file, "%s: type =%d, sub_type = %d, decl_start = %d\n",
          p_inp->label, p_inp->type, p_inp->sub_type, p_inp->decl_start);
  if (p_inp->cmd_def != NULL)
  {
    fprintf(prt_file, "defining command: %s\n", p_inp->cmd_def->name);
    /* dump_command(p_inp->cmd_def); */
  }
}

void dump_int_array(struct int_array* ia)
{
  int i;
  fprintf(prt_file, "dump integer array, length: %d\n", ia->curr);
  for (i = 0; i < ia->curr; i++)
  {
    fprintf(prt_file, v_format("%d "), ia->i[i]);
    if ((i+1)%10 == 0) fprintf(prt_file, "\n");
  }
  if (ia->curr%10 != 0) fprintf(prt_file, "\n");
}

void dump_macro(struct macro* m)
{
  fprintf(prt_file, "name: %s\n", m->name);
  if (m->formal != NULL) dump_char_p_array(m->formal);
  dump_char_array(m->body);
  if (m->tokens != NULL) dump_char_p_array(m->tokens);
}

void dump_macro_list(struct macro_list* ml)
{
  int i;
  puts("++++++ dump of macro list");
  for (i = 0; i < ml->curr; i++) dump_macro(ml->macros[i]);
}

void dump_name_list(struct name_list* nl)
{
  int i;
  puts(" ");
  for (i = 0; i < nl->curr; i++)
  {
    fprintf(prt_file, v_format("%S %I\n"),
            nl->names[nl->index[i]], nl->inform[nl->index[i]]);
  }
}

void dump_node(struct node* node)
{
  int i;
  char pname[NAME_L] = "NULL", nname[NAME_L] = "NULL", 
       from_name[NAME_L] = "NULL";
  if (node->previous != NULL) strcpy(pname, node->previous->name);
  if (node->next != NULL) strcpy(nname, node->next->name);
  if (node->from_name != NULL) strcpy(from_name, node->from_name);
  fprintf(prt_file, 
  v_format("name: %S  occ: %I base: %S  from_name: %S at_value: %F  position: %F\n"),
          node->name, node->occ_cnt, node->base_name, from_name, 
          node->at_value, node->position);
  fprintf(prt_file, v_format("  names of - previous: %S  next: %S\n"),
          pname, nname);
  if (node->cl != NULL)  for (i = 0; i < node->cl->curr; i++)
    dump_constraint(node->cl->constraints[i]);
}

void dump_sequ(struct sequence* c_sequ, int level)
{
  struct node* c_node;
  double suml = zero;
  fprintf(prt_file, v_format("+++ dump sequence: %S\n"), c_sequ->name);
  c_node = c_sequ->start;
  while(c_node != NULL)
  {
    suml += c_node->length;
    if (level > 2)
    {
      dump_node(c_node);
      if (level > 3 && c_node->p_elem != NULL)  dump_element(c_node->p_elem);
    }
    else if (level > 0 && strcmp(c_node->base_name, "drift") != 0)
      fprintf(prt_file, v_format("%S: at = %F\n"),
              c_node->name, c_node->position);
    if (c_node == c_sequ->end)  break;
    c_node = c_node->next;
  }
  fprintf(prt_file, v_format("=== sum of node length: %F\n"), suml);
}

void dump_variable(struct variable* v)
{
  fprintf(prt_file, "=== dumping variable %s\n", v->name);
}

void exec_delete_sequ(char* name)
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

void exec_delete_table(char* name)
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

void export_element(struct element* el, struct el_list* ell, FILE* file)
  /* recursive to have parents always in front for MAD-8 */
{
  int pos = name_list_pos(el->name, ell->list);
  char out[AUX_LG];
  if (pos >= 0)
  {
    if (ell->list->inform[pos] == 0)  /* not yet written */
    {
      export_element(el->parent, ell, file);
      strcpy(out, el->name);
      strcat(out, ": ");
      strcat(out, el->parent->name);
      export_el_def(el, out);
      write_nice(out, file);
      ell->list->inform[pos] = 1;
    }
  }
}

void export_comm_par(struct command_parameter* par, char* string)
  /* exports a command parameter */
{
  int i, k, last;
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
      strcat(string, ":=");
      if (par->expr != NULL) strcat(string, par->expr->string);
      else
      {
        if (par->type == 1)
        {
          k = par->double_value; sprintf(num, v_format("%I"), k);
        }
        else sprintf(num, v_format("%F"), par->double_value);
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
    case 12:
      strcat(string, ":=");
      for (last = par->double_array->curr-1; last > 0; last--)
      {
        if (par->expr_list->list[last] != NULL)
        {
          if (zero_string(par->expr_list->list[last]->string) == 0) break;
        }
        else if (par->double_array->a[last] != zero) break;
      }
      strcat(string, "{");
      for (i = 0; i <= last; i++)
      {
        if (i > 0) strcat(string, ",");
        if (par->expr_list->list[i] != NULL)
          strcat(string, par->expr_list->list[i]->string);
        else
        {
          if (par->type == 11)
          {
            k = par->double_array->a[i]; sprintf(num, v_format("%I"), k);
          }
          else sprintf(num, v_format("%F"), par->double_array->a[i]);
          strcat(string, supp_tb(num));
        }
      }
      strcat(string, "}");
  }
}

void export_elem_8(struct element* el, struct el_list* ell, FILE* file)
  /* exports an element in mad-8 format */
  /* recursive to have parents always in front for MAD-8 */
{
  int pos = name_list_pos(el->name, ell->list);
  char out[AUX_LG];
  if (pos >= 0)
  {
    if (ell->list->inform[pos] == 0)  /* not yet written */
    {
      export_elem_8(el->parent, ell, file);
      strcpy(out, el->name);
      strcat(out, ": ");
      strcat(out, el->parent->name);
      export_el_def_8(el, out);
      write_nice_8(out, file);
      ell->list->inform[pos] = 1;
    }
  }
}

void export_el_def(struct element* el, char* string)
  /* exports an element definition in mad-X format */
{
  int i;
  struct command* def = el->def;
  struct command_parameter* par;
  for (i = 0; i < def->par->curr; i++)
  {
    par = def->par->parameters[i];
    if (def->par_names->inform[i]
        && par_out_flag(el->base_type->name, par->name))
      export_comm_par(par, string);
  }
}

void export_el_def_8(struct element* el, char* string)
  /* exports an element definition in mad-8 format */
{
  struct command* def = el->def;
  struct command_parameter* par;
  int i, div = 1;
  double val[] = {0, 0, 0};
  char* base_name = el->base_type->name;
  char num[2*NAME_L];
  /* special treatment for tilt */
  for (i = 0; i < def->par->curr; i++)
  {
    par = def->par->parameters[i];
    if (def->par_names->inform[i])
    {
      if (strcmp(base_name, "quadrupole") == 0)
      {
        div = 2;
        if (strcmp(par->name,"k1") == 0) val[0] = command_par_special("k1", el);
        else if (strcmp(par->name,"k1s") == 0)
          val[1] = command_par_special("k1s", el);
        else if (strcmp(par->name,"tilt") == 0)
          val[2] = command_par_special("tilt", el);
        else if(par_out_flag(el->base_type->name, par->name))
          export_el_par_8(par, string);
      }
      else if (strcmp(base_name, "sextupole") == 0)
      {
        div = 3;
        if (strcmp(par->name,"k2") == 0) val[0] = command_par_special("k2", el);
        else if (strcmp(par->name,"k2s") == 0)
          val[1] = command_par_special("k2s", el);
        else if (strcmp(par->name,"tilt") == 0)
          val[2] = command_par_special("tilt", el);
        else if(par_out_flag(el->base_type->name, par->name))
          export_el_par_8(par, string);
      }
      else if (strcmp(base_name, "octupole") == 0)
      {
        div = 4;
        if (strcmp(par->name,"k3") == 0) val[0] = command_par_special("k3", el);
        else if (strcmp(par->name,"k3s") == 0)
          val[1] = command_par_special("k3s", el);
        else if (strcmp(par->name,"tilt") == 0)
          val[2] = command_par_special("tilt", el);
        else if(par_out_flag(el->base_type->name, par->name))
          export_el_par_8(par, string);
      }
      else if (strcmp(base_name, "elseparator") == 0)
      {
        if (strcmp(par->name,"ex") == 0) val[0] = command_par_special("ex", el);
        else if (strcmp(par->name,"ey") == 0)
          val[1] = command_par_special("ey", el);
        else if(par_out_flag(el->base_type->name, par->name))
          export_el_par_8(par, string);
      }
      else if(par_out_flag(el->base_type->name, par->name))
        export_el_par_8(par, string);
    }
  }
  if (val[1] != zero)
    val[2] = atan2(val[1], val[0]) / div;
  if (val[2] != zero) val[0] = sqrt(val[0]*val[0]+val[1]*val[1]);
  if (val[0] != zero)
  {
    strcat(string, ",");
    if (strcmp(base_name, "quadrupole") == 0)
    {
      strcat(string, "k1 =");
    }
    else if (strcmp(base_name, "sextupole") == 0)
    {
      strcat(string, "k2 =");
    }
    else if (strcmp(base_name, "octupole") == 0)
    {
      strcat(string, "k3 =");
    }
    else if (strcmp(base_name, "elseparator") == 0)
    {
      strcat(string, "e =");
    }
    sprintf(num, v_format("%F"), val[0]);
    strcat(string, supp_tb(num));
    if (val[2] != zero)
    {
      strcat(string, ",tilt =");
      sprintf(num, v_format("%F"), val[2]);
      strcat(string, supp_tb(num));
    }
  }
}

void export_el_par_8(struct command_parameter* par, char* string)
  /* exports an element parameter in mad-8 format */
{
  int i, k, last, vtilt = 0;
  char num[2*NAME_L], tmp[8], tmpt[8];
  switch(par->type)
  {
    case 0:
      strcat(string, ",");
      strcat(string, par->name);
      strcat(string, " =");
      if (par->double_value == zero) strcat(string, "false");
      else                           strcat(string, "true");
      break;
    case 1:
    case 2:
      strcat(string, ",");
      strcat(string, par->name);
      strcat(string, "=");
      if (par->expr != NULL && strcmp(par->name, "harmon") != 0)
        strcat(string, par->expr->string);
      else
      {
        if (par->type == 1)
        {
          k = par->double_value; sprintf(num, v_format("%I"), k);
        }
        else sprintf(num, v_format("%F"), par->double_value);
        strcat(string, supp_tb(num));
      }
      break;
    case 3:
      if (par->string)
      {
        strcat(string, ",");
        strcat(string, par->name);
        strcat(string, "=");
        strcat(string, par->string);
      }
      break;
    case 11:
    case 12:
      vtilt = strcmp(par->name, "ks") == 0 ? 1 : 0;
      for (last = par->double_array->curr-1; last > 0; last--)
      {
        if (par->expr_list->list[last] != NULL)
        {
          if (zero_string(par->expr_list->list[last]->string) == 0) break;
        }
        else if (par->double_array->a[last] != zero) break;
      }
      for (i = 0; i <= last; i++)
      {
        if (par->expr_list->list[i] != NULL
            && !zero_string(par->expr_list->list[i]->string))
        {
          strcat(string, ",");
          sprintf(tmp, " k%dl =", i);
          sprintf(tmpt, ", t%d", i);
          strcat(string, tmp);
          strcat(string, par->expr_list->list[i]->string);
          if (vtilt) strcat(string, tmpt);
        }
        else if (par->double_array->a[i] != zero)
        {
          strcat(string, ",");
          sprintf(tmp, " k%dl =", i);
          sprintf(tmpt, ", t%d", i);
          if (par->type == 11)
          {
            k = par->double_array->a[i]; sprintf(num, "%d", k);
          }
          else sprintf(num, v_format("%F"), par->double_array->a[i]);
          strcat(string, tmp);
          strcat(string, supp_tb(num));
          if (vtilt) strcat(string, tmpt);
        }
      }
  }
}

void export_sequence(struct sequence* sequ, FILE* file)
  /* exports sequence in mad-X format */
{
  char num[2*NAME_L];
  struct element* el;
  struct sequence* sq;
  struct node* c_node = sequ->start;
  int exp_par_flag;
  char rpos[3][6] = {"exit", "centre", "entry"};
  *c_dum->c = '\0';
  if (sequ->share) strcat(c_dum->c, "shared ");
  strcat(c_dum->c, sequ->export_name);
  strcat(c_dum->c, ": sequence");
  if (sequ->ref_flag)
  {
    strcat(c_dum->c, ", refer = ");
    strcat(c_dum->c, rpos[sequ->ref_flag+1]);
  }
  if (sequ->refpos != NULL)
  {
    strcat(c_dum->c, ", refpos = ");
    strcat(c_dum->c, sequ->refpos);
  }
  strcat(c_dum->c, ", l = ");
  if (sequ->l_expr != NULL) strcat(c_dum->c, sequ->l_expr->string);
  else
  {
    sprintf(num, v_format("%F"), sequence_length(sequ));
    strcat(c_dum->c, supp_tb(num));
  }
  write_nice(c_dum->c, file);
  while(c_node != NULL)
  {
    exp_par_flag = 0;
    *c_dum->c = '\0';
    if (strchr(c_node->name, '$') == NULL
        && strstr(c_node->name, "_p_") == NULL
        && strcmp(c_node->base_name, "drift") != 0)
    {
      if ((el = c_node->p_elem) != NULL)
      {
        if (c_node->p_elem->def_type)
        {
          strcat(c_dum->c, el->name);
          strcat(c_dum->c, ": ");
          strcat(c_dum->c, el->parent->name);
          exp_par_flag = 1;
        }
        else strcat(c_dum->c, el->name);
      }
      else if ((sq = c_node->p_sequ) != NULL) strcat(c_dum->c, sq->name);
      else fatal_error("save error: node without link:", c_node->name);
      strcat(c_dum->c, ", at = ");
      if (c_node->at_expr != NULL) strcat(c_dum->c, c_node->at_expr->string);
      else
      {
        sprintf(num, v_format("%F"), c_node->at_value);
        strcat(c_dum->c, supp_tb(num));
      }
      if (c_node->from_name != NULL)
      {
        strcat(c_dum->c, ", from = ");
        strcat(c_dum->c, c_node->from_name);
      }
      if (exp_par_flag) export_el_def(c_node->p_elem, c_dum->c);
      write_nice(c_dum->c, file);
    }
    if (c_node == sequ->end)  break;
    c_node = c_node->next;
  }
  strcpy(c_dum->c, "endsequence");
  write_nice(c_dum->c, file);
}

void export_sequ_8(struct sequence* sequ, struct command_list* cl, FILE* file)
  /* exports sequence in mad-8 format */
{
  char num[2*NAME_L];
  int exp_par_flag;
  struct element* el;
  struct sequence* sq;
  struct node* c_node = sequ->start;
  if (pass_select_list(sequ->name, cl) == 0)  return;
  *c_dum->c = '\0';
  strcat(c_dum->c, sequ->export_name);
  strcat(c_dum->c, ": sequence");
  if (sequ->ref_flag ==1)  strcat(c_dum->c, ", refer=entry");
  write_nice_8(c_dum->c, file);
  while(c_node != NULL)
  {
    exp_par_flag = 0;
    *c_dum->c = '\0';
    if (strchr(c_node->name, '$') == NULL
        && strcmp(c_node->base_name, "drift") != 0)
    {
      if ((el = c_node->p_elem) != NULL)
      {
        if (c_node->p_elem->def_type)
        {
          strcat(c_dum->c, el->name);
          strcat(c_dum->c, ": ");
          strcat(c_dum->c, el->parent->name);
          exp_par_flag = 1;
        }
        else strcat(c_dum->c, el->name);
      }
      else if ((sq = c_node->p_sequ) != NULL) strcat(c_dum->c, sq->name);
      else fatal_error("save error: node without link:", c_node->name);
      strcat(c_dum->c, ", at = ");
      if (c_node->at_expr != NULL) strcat(c_dum->c, c_node->at_expr->string);
      else
      {
        sprintf(num, v_format("%F"), c_node->at_value);
        strcat(c_dum->c, supp_tb(num));
      }
      if (c_node->from_name != NULL)
      {
        strcat(c_dum->c, ", from = ");
        strcat(c_dum->c, c_node->from_name);
      }
      if (exp_par_flag)  export_el_def_8(c_node->p_elem, c_dum->c);
      write_nice_8(c_dum->c, file);
    }
    if (c_node == sequ->end)  break;
    c_node = c_node->next;
  }
  strcpy(c_dum->c, sequ->name);
  strcat(c_dum->c, "_end: marker, at = ");
  sprintf(num, v_format("%F"), sequence_length(sequ));
  strcat(c_dum->c,num);
  write_nice_8(c_dum->c, file);
  strcpy(c_dum->c, "endsequence");
  write_nice_8(c_dum->c, file);
}

void export_variable(struct variable* var, FILE* file)
  /* exports variable in mad-X format */
{
  int k;
  *c_dum->c = '\0';
  if (var->status == 0) var->value = expression_value(var->expr, var->type);
  if (var->val_type == 0) strcat(c_dum->c, "int ");
  if (var->type == 0) strcat(c_dum->c, "const ");
  strcat(c_dum->c, var->name);
  if (var->type < 2) strcat(c_dum->c, " = ");
  else               strcat(c_dum->c, " := ");
  if (var->expr != NULL) strcat(c_dum->c, var->expr->string);
  else if (var->val_type == 0)
  {
    k = var->value; sprintf(c_join->c, "%d", k); strcat(c_dum->c, c_join->c);
  }
  else
  {
    sprintf(c_join->c, v_format("%F"), var->value);
    strcat(c_dum->c, supp_tb(c_join->c));
  }
  write_nice(c_dum->c, file);
}

void export_var_8(struct variable* var, FILE* file)
  /* exports variable in mad-8 format */
{
  int k;
  *c_dum->c = '\0';
  if (var->status == 0) var->value = expression_value(var->expr, var->type);
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

double find_value(char* name, int ntok, char** toks)
  /* returns value found in construct "name = value", or INVALID */
{
  double val = INVALID;
  int j;
  for (j = 0; j < ntok; j++)
  {
    if (strcmp(toks[j], name) == 0)
    {
      if (j+2 < ntok && *toks[j+1] == '=')
      {
        sscanf(toks[j+2], "%lf", &val);
        break;
      }
    }
  }
  return val;
}

double frndm()
  /* returns random number r with 0 <= r < 1 from flat distribution */
{
  const double one = 1;
  double scale = one / MAX_RAND;
  if (next_rand == NR_RAND)  irngen();
  return scale*irn_rand[next_rand++];
}

void ftoi_array(struct double_array* da, struct int_array* ia)
  /* converts and copies double array into integer array */
{
  int i, l = da->curr;
  while (l >= ia->max)  grow_int_array(ia);
  for (i = 0; i < l; i++) ia->i[i] = da->a[i];
  ia->curr = l;
}

void grow_char_array( /* doubles array size */
  struct char_array* p)
{
  char rout_name[] = "grow_char_array";
  char* p_loc = p->c;
  int j, new = 2*p->max;

  p->max = new;
  p->c = (char*) mymalloc(rout_name, new);
  for (j = 0; j < p->curr; j++) p->c[j] = p_loc[j];
  myfree(rout_name, p_loc);
}

void grow_char_array_list( /* doubles list size */
  struct char_array_list* p)
{
  char rout_name[] = "grow_char_array_list";
  struct char_array** c_loc = p->ca;
  int j, new = 2*p->max;

  p->max = new;
  p->ca
    = (struct char_array**) mycalloc(rout_name,new,
                                     sizeof(struct char_array*));
  for (j = 0; j < p->curr; j++) p->ca[j] = c_loc[j];
  myfree(rout_name, c_loc);
}

void grow_char_p_array( /* doubles array size */
  struct char_p_array* p)
{
  char rout_name[] = "grow_char_p_array";
  char** p_loc = p->p;
  int j, new = 2*p->max;

  p->max = new;
  p->p = (char**) mycalloc(rout_name,new, sizeof(char*));
  for (j = 0; j < p->curr; j++) p->p[j] = p_loc[j];
  myfree(rout_name, p_loc);
}

void grow_command_list( /* doubles list size */
  struct command_list* p)
{
  char rout_name[] = "grow_command_list";
  struct command** c_loc = p->commands;
  int j, new = 2*p->max;

  p->max = new;
  p->commands
    = (struct command**) mycalloc(rout_name,new, sizeof(struct command*));
  for (j = 0; j < p->curr; j++) p->commands[j] = c_loc[j];
  myfree(rout_name, c_loc);
}

void grow_command_list_list( /* doubles list size */
  struct command_list_list* p)
{
  char rout_name[] = "grow_command_list_list";
  struct command_list** c_loc = p->command_lists;
  int j, new = 2*p->max;

  p->max = new;
  p->command_lists = (struct command_list**)
    mycalloc(rout_name,new, sizeof(struct command_list*));
  for (j = 0; j < p->curr; j++) p->command_lists[j] = c_loc[j];
  myfree(rout_name, c_loc);
}

void grow_command_parameter_list( /* doubles list size */
  struct command_parameter_list* p)
{
  char rout_name[] = "grow_command_parameter_list";
  struct command_parameter** c_loc = p->parameters;
  int j, new = 2*p->max;

  p->max = new;
  p->parameters = (struct command_parameter**)
    mycalloc(rout_name,new, sizeof(struct command_parameter*));
  for (j = 0; j < p->curr; j++) p->parameters[j] = c_loc[j];
  myfree(rout_name, c_loc);
}

void grow_constraint_list( /* doubles list size */
  struct constraint_list* p)
{
  char rout_name[] = "grow_constraint_list";
  struct constraint** c_loc = p->constraints;
  int j, new = 2*p->max;

  p->max = new;
  p->constraints = (struct constraint**)
    mycalloc(rout_name,new, sizeof(struct constraint*));
  for (j = 0; j < p->curr; j++) p->constraints[j] = c_loc[j];
  myfree(rout_name, c_loc);
}

void grow_double_array( /* doubles array size */
  struct double_array* p)
{
  char rout_name[] = "grow_double_array";
  double* a_loc = p->a;
  int j, new = 2*p->max;

  p->max = new;
  p->a = (double*) mymalloc(rout_name,new * sizeof(double));
  for (j = 0; j < p->curr; j++) p->a[j] = a_loc[j];
  myfree(rout_name, a_loc);
}

void grow_el_list( /* doubles list size */
  struct el_list* p)
{
  char rout_name[] = "grow_el_list";
  struct element** e_loc = p->elem;
  int j, new = 2*p->max;
  p->max = new;
  p->elem
    = (struct element**) mycalloc(rout_name,new, sizeof(struct element*));
  for (j = 0; j < p->curr; j++) p->elem[j] = e_loc[j];
  myfree(rout_name, e_loc);
}

void grow_expr_list( /* doubles list size */
  struct expr_list* p)
{
  char rout_name[] = "grow_expr_list";
  struct expression** e_loc = p->list;
  int j, new = 2*p->max;
  p->max = new;
  p->list
    = (struct expression**) mycalloc(rout_name,new, sizeof(struct expression*));
  for (j = 0; j < p->curr; j++) p->list[j] = e_loc[j];
  myfree(rout_name, e_loc);
}

void grow_in_buff_list( /* doubles list size */
  struct in_buff_list* p)
{
  char rout_name[] = "grow_in_buff_list";
  struct in_buffer** e_loc = p->buffers;
  FILE** f_loc = p->input_files;
  int j, new = 2*p->max;
  p->max = new;
  p->buffers
    = (struct in_buffer**) mycalloc(rout_name,new, sizeof(struct in_buffer*));
  for (j = 0; j < p->curr; j++) p->buffers[j] = e_loc[j];
  myfree(rout_name, e_loc);
  p->input_files = mycalloc(rout_name, new, sizeof(FILE*));
  for (j = 0; j < p->curr; j++) p->input_files[j] = f_loc[j];
  myfree(rout_name, f_loc);
}

void grow_in_cmd_list( /* doubles list size */
  struct in_cmd_list* p)
{
  char rout_name[] = "grow_in_cmd_list";
  struct in_cmd** c_loc = p->in_cmds;
  int j, new = 2*p->max;

  p->max = new;
  p->in_cmds
    = (struct in_cmd**) mycalloc(rout_name,new, sizeof(struct in_cmd*));
  for (j = 0; j < p->curr; j++) p->in_cmds[j] = c_loc[j];
  myfree(rout_name, c_loc);
}

void grow_int_array( /* doubles array size */
  struct int_array* p)
{
  char rout_name[] = "grow_int_array";
  int* i_loc = p->i;
  int j, new = 2*p->max;

  p->max = new;
  p->i = (int*) mymalloc(rout_name,new * sizeof(int));
  for (j = 0; j < p->curr; j++) p->i[j] = i_loc[j];
  myfree(rout_name, i_loc);
}

void grow_macro_list( /* doubles list size */
  struct macro_list* p)
{
  char rout_name[] = "grow_macro_list";
  struct macro** n_loc = p->macros;
  int j, new = 2*p->max;
  p->max = new;
  p->macros = (struct macro**) mycalloc(rout_name,new, sizeof(struct macro*));
  for (j = 0; j < p->curr; j++) p->macros[j] = n_loc[j];
  myfree(rout_name, n_loc);
}

void grow_name_list( /* doubles list size */
  struct name_list* p)
{
  char rout_name[] = "grow_name_list";
  char** n_loc = p->names;
  int* l_ind = p->index;
  int* l_inf = p->inform;
  int j, new = 2*p->max;

  p->max = new;
  p->names = (char**) mycalloc(rout_name,new, sizeof(char*));
  p->index = (int*) mycalloc(rout_name,new, sizeof(int));
  p->inform = (int*) mycalloc(rout_name,new, sizeof(int));
  for (j = 0; j < p->curr; j++)
  {
    p->names[j] = n_loc[j];
    p->index[j] = l_ind[j];
    p->inform[j] = l_inf[j];
  }
  myfree(rout_name, n_loc);
  myfree(rout_name, l_ind);
  myfree(rout_name, l_inf);
}

void grow_node_list( /* doubles list size */
  struct node_list* p)
{
  char rout_name[] = "grow_node_list";
  struct node** n_loc = p->nodes;
  int j, new = 2*p->max;
  p->max = new;
  p->nodes = (struct node**) mycalloc(rout_name,new, sizeof(struct node*));
  for (j = 0; j < p->curr; j++) p->nodes[j] = n_loc[j];
  myfree(rout_name, n_loc);
}

void grow_sequence_list(struct sequence_list* l)
{
  char rout_name[] = "grow_sequence_list";
  struct sequence** sloc = l->sequs;
  int j, new = 2*l->max;
  l->max = new;
  l->sequs
    = (struct sequence**) mycalloc(rout_name,new, sizeof(struct sequence*));
  for (j = 0; j < l->curr; j++) l->sequs[j] = sloc[j];
  myfree(rout_name, sloc);
}

void double_table(char* table)
{
  int pos;
  struct table* t;

  mycpy(c_dum->c, table);
  if ((pos = name_list_pos(c_dum->c, table_register->names)) > -1)
    t = table_register->tables[pos];
  else return;
  grow_table(t);
}
void grow_table(struct table* t) /* doubles number of rows */
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

void grow_table_list(struct table_list* tl)
{
  char rout_name[] = "grow_table_list";
  struct table** t_loc = tl->tables;
  int j, new = 2*tl->max;

  grow_name_list(tl->names);
  tl->max = new;
  tl->tables = (struct table**) mycalloc(rout_name,new, sizeof(struct table*));
  for (j = 0; j < tl->curr; j++) tl->tables[j] = t_loc[j];
  myfree(rout_name, t_loc);
}

void grow_table_list_list(struct table_list_list* tll)
{
  char rout_name[] = "grow_table_list_list";
  struct table_list** t_loc = tll->table_lists;
  int j, new = 2*tll->max;

  tll->max = new;
  tll->table_lists = (struct table_list**) 
      mycalloc(rout_name,new, sizeof(struct table_list*));
  for (j = 0; j < tll->curr; j++) tll->table_lists[j] = t_loc[j];
  myfree(rout_name, t_loc);
}

void grow_var_list( /* doubles list size */
  struct var_list* p)
{
  char rout_name[] = "grow_var_list";
  struct variable** v_loc = p->vars;
  int j, new = 2*p->max;

  p->max = new;
  p->vars
    = (struct variable**) mycalloc(rout_name,new, sizeof(struct variable*));
  for (j = 0; j < p->curr; j++) p->vars[j] = v_loc[j];
  myfree(rout_name, v_loc);
}

void grow_vector_list( /* doubles list size */
  struct vector_list* p)
{
  char rout_name[] = "grow_vector_list";
  struct double_array** v_loc = p->vectors;
  int j, new = 2*p->max;

  p->max = new;
  p->vectors
    = (struct double_array**) mycalloc(rout_name,new,
                                       sizeof(struct double_array*));
  for (j = 0; j < p->curr; j++) p->vectors[j] = v_loc[j];
  myfree(rout_name, v_loc);
}

double grndm()
  /* returns random number x from normal distribution */
{
  double xi1 = 2*frndm()-one, xi2=2*frndm()-one, zzr;
  while ((zzr = xi1*xi1+xi2*xi2) > one)
  {
    xi1 = 2*frndm()-one; xi2=2*frndm()-one;
  }
  zzr = sqrt(-2*log(zzr)/zzr);
  return xi1*zzr;
}

int inbounds(char* p, int n, char**p_list)
  /* checks whether pointer p is inside one of the ranges in p_list,
     returns 1 for yes, 0 for no */
{
  int j;
  for (j = 0; j < n; j++)
    if (p > p_list[j] && p < p_list[++j]) return 1;
  return 0;
}

void init55(int seed)
  /* initializes random number algorithm */
{
  int i, ii, k = 1, j = abs(seed)%MAX_RAND;
  irn_rand[NR_RAND-1] = j;
  for (i = 0; i < NR_RAND-1; i++)
  {
    ii = (ND_RAND*(i+1))%NR_RAND;
    irn_rand[ii-1] = k;
    if ((k = j - k) < 0) k += MAX_RAND;
    j = irn_rand[ii-1];
  }
  /* warm up */
  for (i = 0; i < 3; i++) irngen();
}

int intrac()
  /* returns non-zero inf program is used interactively, else 0 */
{
  return ((int) isatty(0));
}

void irngen()
  /* creates random number for frndm() */
{
  int i, j;
  for (i = 0; i < NJ_RAND; i++)
  {
    if ((j = irn_rand[i] - irn_rand[i+NR_RAND-NJ_RAND]) < 0) j += MAX_RAND;
    irn_rand[i] = j;
  }
  for (i = NJ_RAND; i < NR_RAND; i++)
  {
    if ((j = irn_rand[i] - irn_rand[i-NJ_RAND]) < 0) j += MAX_RAND;
    irn_rand[i] = j;
  }
  next_rand = 0;
}

int is_token(char* pb, char* string, int slen)
{
  char* pbl = pb;
  if ((pbl == string || *(--pbl) == ' ')
      && (*(pb+slen) == '\0' || *(pb+slen) == ' '))  return 1;
  else return 0;
}

char* join(char** it_list, int n)
  /* joins n character strings into one */
{
  int j;
  *c_join->c = '\0';
  for (j = 0; j < n; j++) strcat(c_join->c, it_list[j]);
  return c_join->c;
}

char* join_b(char** it_list, int n)
  /* joins n character strings into one, blank separated */
{
  char* target;
  int j, k = 0;
  target = c_join->c;
  for (j = 0; j < n; j++)
  {
    strcpy(&target[k], it_list[j]);
    k += strlen(it_list[j]);
    target[k++] = ' ';
  }
  target[k] = '\0';
  return target;
}

void* mycalloc(char* caller, size_t nelem, size_t size)
#ifdef _MEM_LEAKS
{
  /* calls calloc, checks for memory granted */
  void* p;
  int* i_p;
  size_t l_size = nelem*size + sizeof(double);
  if ((p = calloc(1, l_size)) == NULL)
   {
#ifdef _DONOTCATCHOVERFLOW
      warning("mycalloc: memory overflow, called from routine:", caller);
      warning("mycalloc: Program will crash now","");
#else
      fatal_error("memory overflow, called from routine:", caller);
#endif      
   } 
  mtable[-item_no]=p;
  i_p = (int*) p;
  *i_p++ = item_no;
  *i_p = l_size;
  fprintf(stderr,"ALLOCATE called by %s \n",caller);
  fprintf(stderr,"[Allocated item %i (size %i)] \n",item_no,l_size);
  if (item_no == -MTABLE_SIZE)
  {
    fatal_error("Too many allocs!!!", "MTABLE_SIZE");
  }
  item_no = item_no - 1;
  return (void *)((char*)p+sizeof(double));
}
#endif
#ifndef _MEM_LEAKS
{
  /* calls calloc, checks for memory granted */
  void* p;
  int* i_p;
  size_t l_size = nelem*size + sizeof(double);
  if ((p = calloc(1, l_size)) == NULL)
   {
#ifdef _DONOTCATCHOVERFLOW
      warning("mycalloc: memory overflow, called from routine:", caller);
      warning("mycalloc: Program will crash now","");
#else
      fatal_error("memory overflow, called from routine:", caller);
#endif      
   } 
  i_p = (int*) p; *i_p = FREECODE;
  return ((char*)p+sizeof(double));
}
#endif
void mycpy(char* sout, char* sin)
  /* copies string, ends at any non-ascii character including 0 */
{
  char *p, *q;
  int l = 1;

  p = sin;  q = sout;
  while (*p > ' ' && *p <= '~' && l < 2*NAME_L)
  {
    *q++ = *p++;  l++;
  }
  *q = '\0';
}

void myfree(char* rout_name, void* p)
#ifdef _MEM_LEAKS
{
  int my_size,my_item_no,myend,old_item_no;
  char* l_p = (char*)p - sizeof(double);
  int* i_p = (int*) l_p;
  myfree_caller = rout_name;
/* Look for the integer address (backwards) in mtable */
  myend=-item_no-1;
  old_item_no=0;
  while (myend > 0)
  {
    if (mtable[myend] == i_p)
    {
      old_item_no=-myend;
      mtable[myend]=NULL;
      break;
    }
    else
    {
      myend=myend-1;
    }
  }
  if ( old_item_no == 0)
  {
    fatal_error("Free memory error!!!, called from routine:", myfree_caller);
  }
  my_item_no = *i_p++;
  my_size = *i_p;
  if (my_item_no != old_item_no)
  {
    fatal_error("Memory item number discrepancy!!!:", myfree_caller);
  }
  fprintf(stderr,"DEALLOCATE called by %s \n",myfree_caller);
  fprintf(stderr,"[Deallocated item %i (size %i)] \n",my_item_no,my_size);
  free(l_p);
  myfree_caller = none;
}
#endif
#ifndef _MEM_LEAKS
{
  char* l_p = (char*)p - sizeof(double);
  int* i_p = (int*) l_p;
  myfree_caller = rout_name;
  if (*i_p == FREECODE)
  {
    *i_p = 0; free(l_p);
  }
  myfree_caller = none;
}
#endif
void* mymalloc(char* caller, size_t size)
#ifdef _MEM_LEAKS
{
  /* calls malloc, checks for memory granted */
  void* p;
  int* i_p;
  size_t l_size = size + sizeof(double);
  if ((p = malloc(l_size)) == NULL)
   {
#ifdef _DONOTCATCHOVERFLOW
      warning("mymalloc: memory overflow, called from routine:", caller);
      warning("mymalloc:","Program will crash now");
#else
      fatal_error("memory overflow, called from routine:", caller);
#endif      
   } 
  i_p = (int*) p;
  mtable[-item_no]=i_p;
  *i_p++ = item_no;
  *i_p = l_size;
  fprintf(stderr,"ALLOCATE called by %s \n",caller);
  fprintf(stderr,"[Allocated item %i (size %i)] \n",item_no,l_size);
  if (item_no == -MTABLE_SIZE)
  {
    fatal_error("Too many allocs!!!", "MTABLE_SIZE");
  }
  item_no = item_no - 1;
  return (void *)((char*)p+sizeof(double));
}
#endif
#ifndef _MEM_LEAKS
{
  /* calls malloc, checks for memory granted */
  void* p;
  int* i_p;
  size_t l_size = size + sizeof(double)+2;
/*  printf("xxxx %d xxxx\n",l_size);*/
  if ((p = malloc(l_size)) == NULL)
   {
#ifdef _DONOTCATCHOVERFLOW
      warning("mymalloc: memory overflow, called from routine:", caller);
      warning("mymalloc:","Program will crash now");
#else
      fatal_error("memory overflow, called from routine:", caller);
#endif      
   } 
  i_p = (int*) p; *i_p = FREECODE;
  return (void *)((char*)p+sizeof(double));
}
#endif
char* mystrchr(char* string, char c)
  /* returns strchr for character c, but only outside strings included
     in single or double quotes */
{
  char quote = ' '; /* only for the compiler */
  int toggle = 0;
  while (*string != '\0')
  {
    if (toggle)
    {
      if (*string == quote) toggle = 0;
    }
    else if(*string == '\'' || *string == '\"')
    {
      quote = *string; toggle = 1;
    }
    else if (*string == c) return string;
    string++;
  }
  return NULL;
}

void mystrcpy(struct char_array* target, char* source)
{
  /* string copy to char_array with size adjustment */
  while (strlen(source) > target->max) grow_char_array(target);
  strcpy(target->c, source);
}

char* mystrstr(char* string, char* s)
  /* returns strstr for s, but only outside strings included
     in single or double quotes */
{
  char quote = ' '; /* only for the compiler */
  int toggle = 0, n = strlen(s);
  if (n == 0)  return NULL;
  while (*string != '\0')
  {
    if (toggle)
    {
      if (*string == quote) toggle = 0;
    }
    else if(*string == '\'' || *string == '\"')
    {
      quote = *string; toggle = 1;
    }
    else if (strncmp(string, s, n) == 0) return string;
    string++;
  }
  return NULL;
}

void myrepl(char* in, char* out, char* string_in, char* string_out)
  /* replaces all occurrences of "in" in string_in by "out"
     in output string string_out */
{
  int n, add, l_in = strlen(in), l_out = strlen(out);
  char* cp;
  char tmp[8];
  while ((cp = strstr(string_in, in)) != NULL)
  {
    while (string_in != cp) *string_out++ = *string_in++;
    string_in += l_in;
    if (*out == '$')
    {
      n = get_variable(&out[1]);
      sprintf(tmp,"%d", n); add = strlen(tmp);
      strncpy(string_out, tmp, add);
      string_out += add;
    }
    else
    {
      strncpy(string_out, out, l_out);
      string_out += l_out;
    }
  }
  strcpy(string_out, string_in);
}

int name_list_pos(char* p, struct name_list* vlist)
{
  int num, mid, low = 0, high = vlist->curr - 1;
  while (low <= high)
  {
    mid = (low + high) / 2;
    if ((num=strcmp(p, vlist->names[vlist->index[mid]])) < 0)  high = mid - 1;
    else if ( num > 0) low  = mid + 1;
    else               return vlist->index[mid];
  }
  return -1;
}

struct char_array* new_char_array(int length)
{
  char rout_name[] = "new_char_array";
  struct char_array* il =
    (struct char_array*) mycalloc(rout_name,1, sizeof(struct char_array));
  il->stamp = 123456;
  il->curr = 0;
  il->max = length;
  il->c = (char*) mymalloc(rout_name,length);
  return il;
}

struct char_array_list* new_char_array_list(int size)
{
  char rout_name[] = "new_char_array_list";
  struct char_array_list* tl = (struct char_array_list*)
    mycalloc(rout_name,1, sizeof(struct char_array_list));
  strcpy(tl->name, "char_array_list");
  tl->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", tl->name);
  tl->ca
    = (struct char_array**) mycalloc(rout_name,size, sizeof(struct char_array*));
  tl->max = size;
  return tl;
}

struct char_p_array* new_char_p_array(int length)
{
  char rout_name[] = "new_char_p_array";
  struct char_p_array* il
    = (struct char_p_array*) mycalloc(rout_name,1, sizeof(struct char_p_array));
  strcpy(il->name, "char_p_array");
  il->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", il->name);
  il->curr = 0;
  il->max = length;
  il->p = (char**) mycalloc(rout_name,length, sizeof(char*));
  memset(il->p,0,length*sizeof(char*));
  return il;
}

struct command* new_command(char* name, int nl_length, int pl_length,
                            char* module, char* group, int link, int mad_8)
{
  char rout_name[] = "new_command";
  char loc_name[2*NAME_L];
  struct command* new
    = (struct command*) mycalloc(rout_name,1, sizeof(struct command));
  strcpy(loc_name, name); strcat(loc_name, "_param");
  new->stamp = 123456;
  strcpy(new->name, name);
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", loc_name);
  strcpy(new->module, module);
  strcpy(new->group, group);
  new->link_type = link;
  new->mad8_type = mad_8;
  if (nl_length == 0) nl_length = 1;
  new->par_names = new_name_list(loc_name, nl_length);
  new->par = new_command_parameter_list(pl_length);
  return new;
}

struct command_list* new_command_list(char* l_name, int length)
{
  char rout_name[] = "new_command_list";
  struct command_list* il =
    (struct command_list*) mycalloc(rout_name,1, sizeof(struct command_list));
  strcpy(il->name, l_name);
  il->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", il->name);
  il->curr = 0;
  il->max = length;
  il->list = new_name_list(il->name, length);
  il->commands
    = (struct command**) mycalloc(rout_name,length, sizeof(struct command*));
  return il;
}

struct command_list_list* new_command_list_list(int length)
{
  char rout_name[] = "new_command_list_list";
  struct command_list_list* il =
    (struct command_list_list*)
    mycalloc(rout_name,1, sizeof(struct command_list_list));
  strcpy(il->name, "command_list_list");
  il->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", il->name);
  il->curr = 0;
  il->max = length;
  il->list = new_name_list(il->name, length);
  il->command_lists
    = (struct command_list**)
    mycalloc(rout_name,length, sizeof(struct command_list*));
  return il;
}


struct command_list_list* delete_command_list_list( struct command_list_list* ll)
{
  char rout_name[] = "delete_command_list_list";
  
  int i;
  
  if (ll == 0x0) return 0x0;
  
  for (i = 0; i < ll->curr; i++)
   {
     delete_command_list(ll->command_lists[i]);
   }

  myfree(rout_name, ll->command_lists );
  
  delete_name_list(ll->list);

  myfree(rout_name, ll );

  return 0x0;
}


struct command_parameter* new_command_parameter(char* name, int type)
{
  char rout_name[] = "new_command_parameter";
  struct command_parameter* new
    = (struct command_parameter*)
    mycalloc(rout_name,1, sizeof(struct command_parameter));
  strcpy(new->name, name); new->type = type;
  new->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", new->name);
  return new;
}

struct command_parameter_list* new_command_parameter_list(int length)
{
  int i;
  char rout_name[] = "new_command_parameter_list";
  struct command_parameter_list* il =
    (struct command_parameter_list*)
    mycalloc(rout_name,1, sizeof(struct command_parameter_list));
  strcpy(il->name, "command_parameter_list");
  il->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", il->name);
  il->curr = 0;
  il->max = length;
  if (length > 0)
  {
    il->parameters = (struct command_parameter**)
      mycalloc(rout_name,length, sizeof(struct command_parameter*));
    for (i = 0; i < length; i++) il->parameters[i] = NULL;
  }
  return il;
}

struct constraint* new_constraint(int type)
{
  char rout_name[] = "new_constraint";
  struct constraint* new = (struct constraint*)
    mycalloc(rout_name,1, sizeof(struct constraint));
  strcpy(new->name, "constraint");
  new->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", new->name);
  new->type = type;
  return new;
}

struct constraint_list* new_constraint_list(int length)
{
  char rout_name[] = "new_constraint_list";
  struct constraint_list* il
    = (struct constraint_list*)
    mycalloc(rout_name,1, sizeof(struct constraint_list));
  strcpy(il->name, "constraint_list");
  il->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", il->name);
  il->curr = 0;
  il->max = length;
  il->constraints = (struct constraint**)
    mycalloc(rout_name,length, sizeof(struct constraint*));
  return il;
}

struct double_array* new_double_array(int length)
{
  char rout_name[] = "new_double_array";
  struct double_array* il
    = (struct double_array*)
    mycalloc(rout_name,1, sizeof(struct double_array));
  il->stamp = 123456;
  il->curr = 0;
  il->max = length;
  il->a = (double*) mycalloc(rout_name,length, sizeof(double));
  return il;
}

struct element* new_element(char* name)
{
  char rout_name[] = "new_element";
  struct element* el
    = (struct element*) mycalloc(rout_name,1, sizeof(struct element));
  strcpy(el->name, name);
  el->stamp = 123456;
  el->def = 0x0;
  el->parent = 0x0;
  el->base_type = 0x0;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", el->name);
  return el;
}

struct el_list* new_el_list(int length)
{
  char rout_name[] = "new_el_list";
  struct el_list* ell
    = (struct el_list*) mycalloc(rout_name,1, sizeof(struct el_list));
  strcpy(ell->name, "el_list");
  ell->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", ell->name);
  ell->list = new_name_list(ell->name, length);
  ell->elem
    = (struct element**) mycalloc(rout_name,length, sizeof(struct element*));
  ell->max = length;
  return ell;
}

double expr_combine(struct expression* exp1, double val1, char* oper, 
                    struct expression*  exp2, double val2, 
                    struct expression** exp_comb)
{
  double val = 0;

  if (exp1 == NULL && exp2 == NULL)
  {
   *exp_comb = NULL;
   switch(oper[1])
   {
    case '+':
	val = val1 + val2;
        break;
    case '-':
        val = val1 - val2;
   }
  }
  else if(exp1 == NULL) val = combine_val_expr(val1, oper, exp2, exp_comb);
  else if(exp2 == NULL) val = combine_expr_val(exp1, oper, val2, exp_comb);
  else                  val = combine_expr_expr(exp1, oper, exp2, exp_comb);
  return val;
}

struct node* new_elem_node(struct element* el, int occ_cnt)
{
  struct node* p;
  p = new_node(compound(el->name, occ_cnt));
  p->p_elem = el;
  p->length = el->length;
  p->base_name = el->base_type->name;
  p->occ_cnt = occ_cnt;
  return p;
}

struct expression* new_expression(char* in_string, struct int_array* polish)
{
  char rout_name[] = "new_expression";
  int j;
  struct expression* ex =
    (struct expression*) mycalloc(rout_name,1,sizeof(struct expression));
  strcpy(ex->name, "expression");
  ex->stamp = 123456;
  ex->string = (char*) mymalloc(rout_name,strlen(in_string)+1);
  strcpy(ex->string, in_string);
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", ex->name);
  if (polish != NULL)
  {
    ex->polish = new_int_array(polish->curr);
    ex->polish->curr = polish->curr;
    for (j = 0; j < polish->curr; j++) ex->polish->i[j] = polish->i[j];
  }
  return ex;
}

struct expr_list* new_expr_list(int length)
{
  char rout_name[] = "new_expr_list";
  struct expr_list* ell =
    (struct expr_list*) mycalloc(rout_name,1, sizeof(struct expr_list));
  strcpy(ell->name, "expr_list");
  ell->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", ell->name);
  ell->list  =
    (struct expression**) mycalloc(rout_name,length, sizeof(struct expression*));
  ell->max = length;
  return ell;
}

struct in_buffer* new_in_buffer(int length)
{
  char rout_name[] = "new_in_buffer";
  struct in_buffer* new =
    (struct in_buffer*) mycalloc(rout_name,1, sizeof(struct in_buffer));
  strcpy(new->name, "in_buffer");
  new->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", new->name);
  new->c_a = new_char_array(length);
  new->flag = -1;
  return new;
}

struct in_buff_list* new_in_buff_list(int length)
{
  char rout_name[] = "new_inbuf_list";
  struct in_buff_list* bll =
    (struct in_buff_list*) mycalloc(rout_name,1, sizeof(struct in_buff_list));
  strcpy(bll->name, "in_buff_list");
  bll->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", bll->name);
  bll->buffers =
    (struct in_buffer**) mycalloc(rout_name,length, sizeof(struct in_buffer*));
  bll->input_files =
    (FILE**) mycalloc(rout_name,length, sizeof(FILE*));
  bll->max = length;
  return bll;
}

struct in_cmd* new_in_cmd(int length)
{
  char rout_name[] = "new_in_cmd";
  struct in_cmd* new
    = (struct in_cmd*) mycalloc(rout_name,1, sizeof(struct in_cmd));
  strcpy(new->name, "in_cmd");
  new->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", new->name);
  new->tok_list = new_char_p_array(length);
  return new;
}

struct in_cmd_list* new_in_cmd_list(int length)
{
  char rout_name[] = "new_in_cmd_list";
  struct in_cmd_list* il =
    (struct in_cmd_list*) mycalloc(rout_name,1, sizeof(struct in_cmd_list));
  strcpy(il->name, "in_cmd_list");
  il->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", il->name);
  il->curr = 0;
  il->max = length;
  il->labels = new_name_list(il->name, length);
  il->in_cmds
    = (struct in_cmd**) mycalloc(rout_name,length, sizeof(struct in_cmd*));
  return il;
}

struct int_array* new_int_array(int length)
{
  char rout_name[] = "new_int_array";
  struct int_array* il =
    (struct int_array*) mycalloc(rout_name,1, sizeof(struct int_array));
  strcpy(il->name, "int_array");
  il->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", il->name);
  il->curr = 0;
  il->max = length;
  il->i = (int*) mymalloc(rout_name,length * sizeof(int));
  return il;
}

struct macro* new_macro(int n_formal, int length, int p_length)
{
  char rout_name[] = "new_macro";
  struct macro* m
    = (struct macro*) mycalloc(rout_name,1, sizeof(struct macro));
  strcpy(m->name, "macro");
  m->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", m->name);
  if ((m->n_formal  = n_formal) > 0) m->formal = new_char_p_array(n_formal);
  if (p_length > 0) m->tokens = new_char_p_array(p_length);
  ++length;
  m->body = new_char_array(length + 1);
  return m;
}

struct macro_list* new_macro_list(int length)
{
  char rout_name[] = "new_macro_list";
  struct macro_list* nll =
    (struct macro_list*) mycalloc(rout_name,1, sizeof(struct macro_list));
  strcpy(nll->name, "macro_list");
  nll->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", nll->name);
  nll->list = new_name_list(nll->name, length);
  nll->macros
    = (struct macro**) mycalloc(rout_name,length, sizeof(struct macro*));
  nll->max = length;
  return nll;
}

struct name_list* new_name_list(char* list_name, int length)
{
  char rout_name[] = "new_name_list";
  struct name_list* il =
    (struct name_list*) mycalloc(rout_name,1, sizeof(struct name_list));
  strcpy(il->name, list_name);
  il->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", il->name);
  il->names = (char**) mycalloc(rout_name,length, sizeof(char*));
  il->index = (int*) mycalloc(rout_name,length, sizeof(int));
  il->inform = (int*) mycalloc(rout_name,length, sizeof(int));
  il->max = length;
  return il;
}

struct node* new_node(char* name)
{
  char rout_name[] = "new_node";
  struct node* p = (struct node*) mycalloc(rout_name,1, sizeof(struct node));
  strcpy(p->name, name);
  p->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", p->name);
  return p;
}

struct node_list* new_node_list(int length)
{
  char rout_name[] = "new_node_list";
  struct node_list* nll =
    (struct node_list*) mycalloc(rout_name,1, sizeof(struct node_list));
  strcpy(nll->name, "node_list");
  nll->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", nll->name);
  nll->list = new_name_list(nll->name, length);
  nll->nodes
    = (struct node**) mycalloc(rout_name,length, sizeof(struct node*));
  nll->max = length;
  return nll;
}

struct sequence* new_sequence(char* name, int ref)
{
  char rout_name[] = "new_sequence";
  struct sequence* s = mycalloc(rout_name,1, sizeof(struct sequence));
  strcpy(s->name, name);
  s->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", s->name);
  s->ref_flag = ref;
  s->nodes = new_node_list(10000);
  return s;
}

struct sequence_list* new_sequence_list(int length)
{
  char rout_name[] = "new_sequence_list";
  struct sequence_list* s
    = mycalloc(rout_name,length, sizeof(struct sequence_list));
  strcpy(s->name, "sequence_list");
  s->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", s->name);
  s->max = length;
  s->list = new_name_list(s->name, length);
  s->sequs
    = (struct sequence**) mycalloc(rout_name,length, sizeof(struct sequence*));
  return s;
}

struct node* new_sequ_node(struct sequence* sequ, int occ_cnt)
{
  struct node* p;
  p = new_node(compound(sequ->name, occ_cnt));
  p->p_sequ = sequ;
  p->length = sequence_length(sequ);
  p->base_name = permbuff("sequence");
  return p;
}

struct table* new_table(char* name, char* type, int rows,
                        struct name_list* cols)
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

struct table_list* new_table_list(int size)
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

struct table_list_list* new_table_list_list(int size)
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

struct variable* new_variable(char* name, double val, int val_type,
                              int type, struct expression* expr, char* string)
{
  char rout_name[] = "new_variable";
  struct variable* var =
    (struct variable*) mycalloc(rout_name,1, sizeof(struct variable));
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

struct var_list* new_var_list(int length)
{
  char rout_name[] = "new_var_list";
  struct var_list* var
    = (struct var_list*) mycalloc(rout_name,1, sizeof(struct var_list));
  strcpy(var->name, "var_list");
  var->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", var->name);
  var->list = new_name_list(var->name, length);
  var->vars
    = (struct variable**) mycalloc(rout_name,length, sizeof(struct variable*));
  var->max = length;
  return var;
}

struct vector_list* new_vector_list(int length)
  /* creates a name list and pointer list
     for double arrays with initial length "length".
  */
{
  char rout_name[] = "new_vector_list";
  struct vector_list* vector
    = (struct vector_list*) mycalloc(rout_name,1, sizeof(struct vector_list));
  vector->max = length;
  vector->names = new_name_list("vector_list", length);
  vector->vectors
    = (struct double_array**) mycalloc(rout_name, length,
                                       sizeof(struct double_array*));
  return vector;
}

char next_non_blank(char* string)
  /* returns next non-blank in string outside quotes, else blank */
{
  int i, toggle = 0, l = strlen(string);
  char quote = ' ';
  for (i = 0; i < l; i++)
  {
    if (toggle)
    {
      if (string[i] == quote)  toggle = 0;
    }
    else if (string[i] == '\'' || string[i] == '\"')
    {
      quote = string[i]; toggle = 1;
    }
    else if (string[i] != ' ')  return string[i];
  }
  return ' ';
}

int next_non_blank_pos(char* string)
  /* returns position of next non-blank in string outside quotes, else -1 */
{
  int i, toggle = 0, l = strlen(string);
  char quote = ' ';
  for (i = 0; i < l; i++)
  {
    if (toggle)
    {
      if (string[i] == quote)  toggle = 0;
    }
    else if (string[i] == '\'' || string[i] == '\"')
    {
      quote = string[i]; toggle = 1;
    }
    else if (string[i] != ' ')  return i;
  }
  return -1;
}

char* noquote(char* string)
{
  char* c = string;
  char* d = c;
  char k;
  if (string != NULL)
  {
    k = *c;
    if (k == '\"' || k == '\'')
    {
      d++;
      while (*d != k) *c++ = *d++;
      *c = '\0';
    }
  }
  return string;
}

int par_out_flag(char* base_name, char* par_name)
{
  /* marks the element parameters that are to be written on "save" */
  if (strcmp(par_name,"at") == 0 || strcmp(par_name,"from") == 0) return 0;
  if (strcmp(base_name, "multipole") == 0
      && strcmp(par_name,"l") == 0) return 0;
  if (strcmp(base_name, "rcollimator") == 0
      && strcmp(par_name,"lrad") == 0) return 0;
  if (strcmp(base_name, "ecollimator") == 0
      && strcmp(par_name,"lrad") == 0) return 0;
  return 1;
}

int predef_var(struct variable* var)
  /* return 1 for predefined variable, else 0 */
{
  int pos = name_list_pos(var->name, variable_list->list);
  return (pos < start_var ? 1 : 0) ;
}

void print_command(struct command* cmd)
{
  int i;
  fprintf(prt_file, "command: %s\n", cmd->name);
  for (i = 0; i < cmd->par->curr; i++)
  {
    print_command_parameter(cmd->par->parameters[i]);
    if ((i+1)%3 == 0) fprintf(prt_file, "\n");
  }
  if (i%3 != 0) fprintf(prt_file, "\n");
}

void print_command_parameter(struct command_parameter* par)
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

void print_global(double delta)
{
  char tmp[NAME_L], trad[4];
  double alfa = get_value("probe", "alfa");
  double freq0 = get_value("probe", "freq0");
  double gamma = get_value("probe", "gamma");
  double beta = get_value("probe", "beta");
  double circ = get_value("probe", "circ");
  double bcurrent = get_value("probe", "bcurrent");
  double npart = get_value("probe", "npart");
  double energy = get_value("probe", "energy");
  int kbunch = get_value("probe", "kbunch");
  int rad = get_value("probe", "radiate");
  double gamtr = zero, t0 = zero, eta;
  get_string("probe", "particle", tmp);
  if (rad) strcpy(trad, "T");
  else     strcpy(trad, "F");
  if (alfa > zero) gamtr = sqrt(one / alfa);
  else if (alfa < zero) gamtr = sqrt(-one / alfa);
  if (freq0 > zero) t0 = one / freq0;
  eta = alfa - one / (gamma*gamma);
  puts(" ");
  printf(" Global parameters for %ss, radiate = %s:\n\n",
         tmp, trad);
  printf(v_format(" C         %F m          f0        %F MHz\n"),circ, freq0);
  printf(v_format(" T0        %F musecs     alfa      %F \n"), t0, alfa);
  printf(v_format(" eta       %F            gamma(tr) %F \n"), eta, gamtr);
  printf(v_format(" Bcurrent  %F A/bunch    Kbunch    %I \n"),
         bcurrent, kbunch);
  printf(v_format(" Npart     %F /bunch     Energy    %F GeV \n"),
         npart,energy);
  printf(v_format(" gamma     %F            beta      %F\n"),
         gamma, beta);
}

void print_rfc()
  /* prints the rf cavities present */
{
  double freq0, harmon, freq;
  int i, n = current_sequ->cavities->curr;
  struct element* el;
  if (n == 0)  return;
  freq0 = command_par_value("freq0", probe_beam);
  printf("\n RF system: \n");
  printf(v_format(" %S %NFs %NFs %NFs %NFs %NFs\n"),
         "Cavity","length[m]","voltage[MV]","lag","freq[MHz]","harmon");
  for (i = 0; i < n; i++)
  {
    el = current_sequ->cavities->elem[i];
    if ((harmon = el_par_value("harmon", el)) > zero)
    {
      freq = freq0 * harmon;
      printf(v_format(" %S %F %F %F %F %F\n"),
             el->name, el->length, el_par_value("volt", el),
             el_par_value("lag", el), freq, harmon);
    }
  }
}

void print_table(struct table* t)
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

void print_value(struct in_cmd* cmd)
{
  char** toks = &cmd->tok_list->p[cmd->decl_start];
  int n = cmd->tok_list->curr - cmd->decl_start;
  int j, s_start = 0, end, type, nitem;
  while (s_start < n)
  {
    for (j = s_start; j < n; j++) if (*toks[j] == ',') break;
    if ((type = loc_expr(toks, j, s_start, &end)) > 0)
    {
      nitem = end + 1 - s_start;
      if (polish_expr(nitem, &toks[s_start]) == 0)
        fprintf(prt_file, v_format("%s = %F ;\n"),
                spec_join(&toks[s_start], nitem), 
                polish_value(deco, join(&toks[s_start], nitem)));
      else
      {
        warning("invalid expression:", spec_join(&toks[s_start], nitem));
        return;
      }
      s_start = end+1;
      if (s_start < n-1 && *toks[s_start] == ',') s_start++;
    }
    else
    {
      warning("invalid expression:", spec_join(&toks[s_start], n));
      return;
    }
  }
}

int quote_level(char* string, char* send)
{
/* returns the level count of quotation marks " and ' inside string between */
/* start of string and send */
  int level = 0;
  char* p;
  char c = ' ';
  for (p = string; p < send; p++)
  {
    if (level == 0)
    {
      if (*p == '\"' || *p == '\'')
      {
        c = *p; level++;
      }
    }
    else if(*p == c) level--;
  }
  return level;
}

int remove_colon(char** toks, int number, int start)
  /* removes colon behind declarative part for MAD-8 compatibility */
{
  int i, k = start;
  for (i = start; i < number; i++)  if (*toks[i] != ':') toks[k++] = toks[i];
  return k;
}

void remove_from_macro_list( /* removes macro from alphabetic macro list */
  struct macro* macro, struct macro_list* nll)
{
  int pos;
  if ((pos = name_list_pos(macro->name, nll->list)) > -1)
  {
    remove_from_name_list(macro->name, nll->list);
    delete_macro(nll->macros[pos]);
    nll->macros[pos] = nll->macros[nll->curr--];
  }
}

int remove_from_name_list(char* name, struct name_list* nl)
{
  int j, i, k = -1;
  for (i = 0; i < nl->curr; i++)
    if (strcmp(nl->names[nl->index[i]], name) == 0) break;
  
  
  if (i < nl->curr)
  {
   
    k = nl->index[i];
    for (j = i+1; j < nl->curr; j++) nl->index[j-1] = nl->index[j];
    
    for (j = 0; j < nl->curr-1; j++)
      if(nl->index[j] == nl->curr-1) break;
      
    nl->index[j] = k;
    nl->inform[k] = nl->inform[nl->curr-1];
    nl->names[k] = nl->names[--nl->curr];
  }
  
  return k;
}

void remove_from_sequ_list(struct sequence* sequ, struct sequence_list* sql)
  /* removes sequence sequ from sequence list sql */
{
  int i, pos;
  for (i = 0; i < sql->curr; i++) if (sql->sequs[i] == sequ)  break;
  pos = name_list_pos(sequ->name, sql->list);
  remove_from_name_list(sequ->name, sql->list);
  for (i = pos+1; i < sql->curr; i++) sql->sequs[i-1] = sql->sequs[i];
  sql->curr--;  
  return;
}

void replace(char* buf, char in, char out)
  /* replaces character in by character out in string buf */
{
  int j, l = strlen(buf);
  for (j = 0; j < l; j++)  if (buf[j] == in)  buf[j] = out;
}

double sequence_length(struct sequence* sequ)
{
  double val = 0;
  if (sequ)
  {
   if (sequ->l_expr)  
   val = sequ->length = expression_value(sequ->l_expr,2);
   else val = sequ->length;
  }
  return val;
}

int square_to_colon(char* string)
  /* sets occurrence count behind colon, possibly replacing [] */
{
  char* t;
  int k = strlen(string);
  if ((t = strchr(string, '[')) == NULL)
  {
    string[k++] = ':'; string[k++] = '1'; string[k] = '\0';
  }
  else
  {
    *t = ':';
    if ((t = strchr(string, ']')) == NULL)  return 0;
    else *t = '\0';
  }
  return strlen(string);
}

char* stolower(char* s)  /* converts string to lower in place */
{
  char *c = s;
  int j;
  for (j = 0; j < strlen(s); j++)
  {
    *c = (char) tolower((int) *c); c++;
  }
  return s;
}

void stolower_nq(char* s)
  /* converts string to lower in place outside quotes */
{
  char *c = s;
  int j, toggle = 0;
  char quote = ' '; /* just to suit the compiler */
  for (j = 0; j < strlen(s); j++)
  {
    if (toggle)
    {
      if (*c == quote) toggle = 0;
    }
    else if (*c == '\"' || *c == '\'')
    {
      toggle = 1; quote = *c;
    }
    else 
     {
       *c = (char) tolower((int) *c);
     }  
    c++;
  }
}

char* stoupper(char* s)  /* converts string to upper in place */
{
  char *c = s;
  int j;
  for (j = 0; j < strlen(s); j++)
  {
    *c = (char) toupper((int) *c); c++;
  }
  return s;
}

int string_cnt(char c, int n, char* toks[])
  /* returns number of strings in toks starting with character c */
{
  int i, k = 0;
  for (i = 0; i < n; i++) if(*toks[i] == c) k++;
  return k;
}

char* strip(char* name)
  /* strip ':' and following off */
{
  char* p;
  strcpy(tmp_key, name);
  if ((p = strchr(tmp_key, ':')) != NULL) *p = '\0';
  return tmp_key;
}

void supp_char(char c, char* string)
  /* suppresses character c in string */
{
  char* cp = string;
  while (*string != '\0')
  {
    if (*string != c)  *cp++ = *string;
    string++;
  }
  *cp = '\0';
}

int supp_lt(char* inbuf, int flag)
  /* suppress leading, trailing blanks and replace some special char.s*/
{
  int l = strlen(inbuf), i, j;
  replace(inbuf, '\x9', ' '); /* tab */
  replace(inbuf, '\xd', ' '); /* Windows e-o-l */
  if (flag == 0)  replace(inbuf, '\n', ' '); /* e-o-l */
  supp_tb(inbuf); /* suppress trailing blanks */
  if ((l = strlen(inbuf)) > 0)
  {
    for (j = 0; j < l; j++) if (inbuf[j] != ' ') break; /* leading blanks */
    if (j > 0)
    {
      for (i = 0; i < l - j; i++) inbuf[i] = inbuf[i+j];
      inbuf[i] = '\0';
    }
  }
  return strlen(inbuf);
}

void supp_mul_char(char c, char* string)
  /* reduces multiple occurrences of c in string to 1 occurrence */
{
  char* cp = string;
  int cnt = 0;
  while (*string != '\0')
  {
    if (*string != c)
    {
      *cp++ = *string; cnt = 0;
    }
    else if (cnt++ == 0) *cp++ = *string;
    string++;
  }
  *cp = '\0';
}

char* supp_tb(char* string) /* suppress trailing blanks in string */
{
  int l = strlen(string), j;
  for (j = l-1; j >= 0; j--)
  {
    if (string[j] != ' ') break;
    string[j] = '\0';
  }
  return string;
}

void termination_handler(int signum)
{
  if (strcmp(myfree_caller, "none") == 0)
    puts("+++ memory access outside program range, fatal +++");
  else
    printf("+++ illegal call to free memory from routine: %s +++\n",
           myfree_caller);
  puts("good bye");
  exit(1);
}

double tgrndm(double cut)
  /* returns random variable from normal distribution cat at 'cut' sigmas */
{
  double ret = zero;
  if (cut > zero)
  {
    ret = grndm();
    while (fabs(ret) > fabs(cut))  ret = grndm();
  }
  return ret;
}

double vdot(int* n, double* v1, double* v2)
  /* returns dot product of vectors v1 and v2 */
{
  int i;
  double dot = 0;
  for (i = 0; i < *n; i++)  dot += v1[i] * v2[i];
  return dot;
}

double vmod(int* n, double* v)
{
  int i;
  double mod = 0;
  for (i = 0; i < *n; i++)  mod += v[i] * v[i];
  return sqrt(mod);
}

int v_length(char* form)
{
  int ret = 0;
  if      (form[1] == 'I') sscanf(int_format, "%d", &ret);
  else if (form[1] == 'F') sscanf(float_format, "%d", &ret);
  return ret;
}

char* v_format(char* string)
  /* copies string to gloval variable var_form
     replacing  %S, %I, and %F by the user defined formats;
     %NF and %NI are replaced by the field lengths (!) of the defined formats */
{
  char *p, *q = string, *s = string, *t;
  char c;
  *var_form = '\0';
  while ((p = strpbrk(s, "NIFS")))
  {
    if ((uintptr_t)p > (uintptr_t)q)
    {
      t = p; t--;
      if (*t == '%')
      {
        c = *p;
        strncat(var_form, q, (uintptr_t)p - (uintptr_t)q);
        if (c == 'N')
        {
          sprintf(&var_form[strlen(var_form)], "%d", v_length(p));
          p++;
        }
        else if (c == 'F')  strcat(var_form, float_format);
        else if (c == 'S')  strcat(var_form, string_format);
        else if (c == 'I')  strcat(var_form, int_format);
        q = p; q++;
      }
    }
    s = ++p;
  }
  strcat(var_form, q);
  return var_form;
}

void write_elems(struct el_list* ell, struct command_list* cl, FILE* file)
{
  int i;
  for (i = 0; i < ell->curr; i++)
  {
    if (pass_select_list(ell->elem[i]->name, cl))
      export_element(ell->elem[i], ell, file);
  }
}

void write_elems_8(struct el_list* ell, struct command_list* cl, FILE* file)
{
  int i;
  for (i = 0; i < ell->curr; i++)
  {
    if (pass_select_list(ell->elem[i]->name, cl))
      export_elem_8(ell->elem[i], ell, file);
  }
}

void write_nice(char* string, FILE* file)
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

void write_nice_8(char* string, FILE* file)
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

void write_sequs(struct sequence_list* sql,struct command_list* cl, FILE* file)
{
  /* exports sequences in order of their nest level, flat first etc. */
  int i, j, max_nest = 0;
  for (i = 0; i < sql->curr; i++)
    if(sql->sequs[i]->nested > max_nest) max_nest = sql->sequs[i]->nested;
  for (j = 0; j <= max_nest; j++)
  {
    for (i = 0; i < sql->curr; i++)
      if(sql->sequs[i]->nested == j)
      {
        if (pass_select_list(sql->sequs[i]->name, cl))
          export_sequence(sql->sequs[i], file);
      }
  }
}

void write_vars(struct var_list* varl, struct command_list* cl, FILE* file)
{
  int i;
  for (i = 0; i < varl->curr; i++)
  {
    if (predef_var(varl->vars[i]) == 0
        && pass_select_list(varl->vars[i]->name, cl))
      export_variable(varl->vars[i], file);
  }
}

void write_vars_8(struct var_list* varl, struct command_list* cl, FILE* file)
{
  int i;
  for (i = 0; i < varl->curr; i++)
  {
    if (predef_var(varl->vars[i]) == 0
        && pass_select_list(varl->vars[i]->name, cl))
      export_var_8(varl->vars[i], file);
  }
}

void write_table(struct table* t, char* filename)
  /* writes rows with columns listed in row and col */
{
  char l_name[NAME_L];
  char sys_name[200];
  char* pc = 0x0;
  struct int_array* col = t->col_out;
  struct int_array* row = t->row_out;
  int i, j, k, tmp;
  time_t now;
  struct tm* tm;
#ifndef _WIN32
  struct utsname u;
  i = uname(&u); /* get system name */
  strcpy(sys_name, u.sysname);
#endif
#ifdef _WIN32
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
    fprintf(out_file,
            "@ NAME             %%%02ds \"%s\"\n", strlen(t->name),
            stoupper(l_name));

    strcpy(l_name, t->type);
    fprintf(out_file,
            "@ TYPE             %%%02ds \"%s\"\n", strlen(t->type),
            stoupper(l_name));

    if (t->header != NULL)
    {
      for (j = 0; j < t->header->curr; j++)
        fprintf(out_file, "%s\n", t->header->p[j]);
    }
    if (title != NULL)
      fprintf(out_file,
              "@ TITLE            %%%02ds \"%s\"\n", strlen(title), title);

    fprintf(out_file,
            "@ ORIGIN           %%%02ds \"%s %s\"\n",
            strlen(myversion)+strlen(sys_name)+1, myversion, sys_name);

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
              pc = mymalloc("write_table",1);
              *pc =       '\"';
              *(c_dum->c) = '\"';
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

void zero_double(double* a, int n)
  /* sets first n values in double array a to zero */
{
  int j;
  for (j = 0; j < n; j++)  a[j] = zero;
}

int zero_string(char* string) /* returns 1 if string defaults to '0', else 0 */
{
  int i, l = strlen(string);
  char c;
  for (i = 0; i < l; i++)
    if ((c = string[i]) != '0' && c != ' ' && c != '.') return 0;
  return 1;
}

void cf77flush()
#include <stdio.h>
{
  fflush(stdout);
  fflush(stderr);
}

/************************************************************************/
/* The functions below belongs to matchc2, but were moved here to make
   mpars compile */

/************************************************************************/
int match2_augmentnmacros()
{
  /*makes place in the working arrays for a new macro
    Piotr Skowronski Mar 2007
  */
  int i,j;
  char fn[]={"match2_augmentnmacros"};
  char**   new_match2_macro_name;
  char* ** new_match2_cons_name;
  double** new_match2_cons_value;
  double** new_match2_cons_value_rhs;
  double** new_match2_cons_value_lhs;
  double** new_match2_cons_weight;
  char**   new_match2_cons_sign;
  struct expression* ** new_match2_cons_rhs;
  struct expression* ** new_match2_cons_lhs;

  if(MAX_MATCH_MACRO == 0)
  {
    error("match2_augmentnconstraints","match with use_maco was not initialized");
    return 1;
  }

  new_match2_macro_name     = (char**)  mycalloc(fn,MAX_MATCH_MACRO+1,sizeof(char*));
  new_match2_cons_name      = (char* **)mycalloc(fn,MAX_MATCH_MACRO+1,sizeof(char**));
  new_match2_cons_value     = (double**)mycalloc(fn,MAX_MATCH_MACRO+1,sizeof(double*));
  new_match2_cons_value_rhs = (double**)mycalloc(fn,MAX_MATCH_MACRO+1,sizeof(double*));
  new_match2_cons_value_lhs = (double**)mycalloc(fn,MAX_MATCH_MACRO+1,sizeof(double*));
  new_match2_cons_weight    = (double**)mycalloc(fn,MAX_MATCH_MACRO+1,sizeof(double*));
  new_match2_cons_sign      = (char**)  mycalloc(fn,MAX_MATCH_MACRO+1,sizeof(char*));
  new_match2_cons_rhs       = (struct expression* **) mycalloc(fn,MAX_MATCH_MACRO+1,sizeof(struct expression**));
  new_match2_cons_lhs       = (struct expression* **) mycalloc(fn,MAX_MATCH_MACRO+1,sizeof(struct expression**));

  /*copy old pointers to arrays*/
  for(i=0;i<MAX_MATCH_MACRO;i++)
  {
    new_match2_cons_name[i]      = match2_cons_name[i];

    new_match2_cons_value[i]     = match2_cons_value[i];

    new_match2_cons_value_rhs[i] = match2_cons_value_rhs[i];
    new_match2_cons_value_lhs[i] = match2_cons_value_lhs[i];
    new_match2_cons_weight[i]    = match2_cons_weight[i];
    new_match2_cons_sign[i]      = match2_cons_sign[i];

    new_match2_cons_rhs[i]       = match2_cons_rhs[i];
    new_match2_cons_lhs[i]       = match2_cons_lhs[i];

    new_match2_macro_name[i]     = match2_macro_name[i];
  }

  /*free the old arrays*/
  myfree(fn,match2_cons_name);
  myfree(fn,match2_cons_value);
  myfree(fn,match2_cons_value_rhs);
  myfree(fn,match2_cons_value_lhs);
  myfree(fn,match2_cons_weight);
  myfree(fn,match2_cons_sign);
  myfree(fn,match2_cons_rhs);
  myfree(fn,match2_cons_lhs);
  myfree(fn,match2_macro_name);


  /*assign freed pointers to the new arrays*/

  match2_cons_name = new_match2_cons_name  ;
  match2_cons_value = new_match2_cons_value  ;
  match2_cons_value_rhs = new_match2_cons_value_rhs  ;
  match2_cons_value_lhs = new_match2_cons_value_lhs  ;
  match2_cons_weight = new_match2_cons_weight  ;
  match2_cons_sign = new_match2_cons_sign  ;
  match2_cons_rhs = new_match2_cons_rhs  ;
  match2_cons_lhs = new_match2_cons_lhs  ;
  match2_macro_name = new_match2_macro_name  ;



  /*make arrays in the new row*/

  match2_cons_name[MAX_MATCH_MACRO]      = (char**)mycalloc(fn,MAX_MATCH_CONS,sizeof(char*));
  match2_cons_value[MAX_MATCH_MACRO]     = (double*)mycalloc(fn,MAX_MATCH_CONS,sizeof(double));

  match2_cons_value_rhs[MAX_MATCH_MACRO] = (double*)mycalloc(fn,MAX_MATCH_CONS,sizeof(double));
  match2_cons_value_lhs[MAX_MATCH_MACRO] = (double*)mycalloc(fn,MAX_MATCH_CONS,sizeof(double));
  match2_cons_weight[MAX_MATCH_MACRO]    = (double*)mycalloc(fn,MAX_MATCH_CONS,sizeof(double));
  match2_cons_sign[MAX_MATCH_MACRO]      = (char*)mycalloc(fn,MAX_MATCH_CONS,sizeof(char));

  match2_cons_rhs[MAX_MATCH_MACRO]       = (struct expression**) mycalloc(fn,MAX_MATCH_CONS,sizeof(struct expression*));
  match2_cons_lhs[MAX_MATCH_MACRO]       = (struct expression**) mycalloc(fn,MAX_MATCH_CONS,sizeof(struct expression*));


  /*initializes arrays in the last row*/

  match2_macro_name[MAX_MATCH_MACRO]=NULL;

  for(j=0;j<MAX_MATCH_CONS;j++)
  {
    match2_cons_name     [MAX_MATCH_MACRO][j]=0x0;

    match2_cons_value    [MAX_MATCH_MACRO][j]=0.0;
    match2_cons_value_lhs[MAX_MATCH_MACRO][j]=0.0;
    match2_cons_value_rhs[MAX_MATCH_MACRO][j]=0.0;
    match2_cons_weight   [MAX_MATCH_MACRO][j]=0.0;
    match2_cons_sign     [MAX_MATCH_MACRO][j]='n';

    match2_cons_rhs      [MAX_MATCH_MACRO][j]=0x0;
    match2_cons_lhs      [MAX_MATCH_MACRO][j]=0x0;
  }


  return ++MAX_MATCH_MACRO;

}
/************************************************************************/

void match2_delete_arrays()
{
  /*clean the stuff;*/
  int i;
  char fn[]={"match2_delete_arrays"};

  if(MAX_MATCH_MACRO <= 0) return;

  for(i=0;i<MAX_MATCH_MACRO;i++)
  {
    if(match2_cons_name[i] == 0x0) break;
    myfree(fn,match2_cons_name     [i]);
    myfree(fn,match2_cons_value    [i]);
    myfree(fn,match2_cons_value_lhs[i]);
    myfree(fn,match2_cons_value_rhs[i]);
    myfree(fn,match2_cons_weight   [i]);
    myfree(fn,match2_cons_sign     [i]);
    myfree(fn,match2_cons_rhs      [i]);
    myfree(fn,match2_cons_lhs      [i]);
  }

  myfree(fn,match2_cons_name);
  myfree(fn,match2_cons_value);
  myfree(fn,match2_cons_value_rhs);
  myfree(fn,match2_cons_value_lhs);
  myfree(fn,match2_cons_weight);
  myfree(fn,match2_cons_sign);
  myfree(fn,match2_cons_rhs);
  myfree(fn,match2_cons_lhs);
  myfree(fn,match2_macro_name);


  match2_cons_name = 0x0;
  match2_cons_value = 0x0;
  match2_cons_value_rhs = 0x0;
  match2_cons_value_lhs = 0x0;
  match2_cons_weight = 0x0;
  match2_cons_sign = 0x0;
  match2_cons_rhs = 0x0;
  match2_cons_lhs = 0x0;
  match2_macro_name = 0x0;

  /*for security so we cannot add more constraints if the module is not initialized*/
  MAX_MATCH_CONS =  0;
  MAX_MATCH_MACRO = 0;

}
/************************************************************************/


void match2_delete_expressions()
{
  char rout_name[] = "match2_delete_expressions";

  int i,j;

  for(i=0;i<MAX_MATCH_MACRO;i++)
  {
    if ( match2_cons_name[i][0] == 0x0) break;
    for(j=0;j<MAX_MATCH_CONS;j++)
    {
      if ( match2_cons_name[i][j] == 0x0) break;
      myfree(rout_name,match2_cons_name[i][j]);
      delete_expression(match2_cons_rhs[i][j]);
      delete_expression(match2_cons_lhs[i][j]);
      match2_cons_rhs[i][j] = 0x0;
      match2_cons_lhs[i][j] = 0x0;
    }
  }

}
