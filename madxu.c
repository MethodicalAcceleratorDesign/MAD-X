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
  struct expr_list* clone = new_expr_list(p->curr);
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

struct name_list* clone_name_list(struct name_list* p)
{
  int i, l = p->curr > 0 ? p->curr : 1;
  struct name_list* clone = new_name_list(l);
  for (i = 0; i < p->curr; i++) clone->index[i] = p->index[i];
  for (i = 0; i < p->curr; i++) clone->inform[i] = p->inform[i];
  for (i = 0; i < p->curr; i++) clone->names[i] = p->names[i];
  clone->curr = p->curr;
  return clone;
}

struct node* clone_node(struct node* p)
{
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
  return clone;
}

void copy_name_list(struct name_list* out, struct name_list* in)
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
  if (pa == NULL)  return NULL;
  if (stamp_flag && pa->stamp != 123456) 
     fprintf(stamp_file, "d_c_a double delete --> %s\n", pa->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", pa->name);
  if (pa->c != NULL)  free(pa->c);
  free(pa);
  return NULL;
}

struct char_p_array* delete_char_p_array(struct char_p_array* pa, int flag)
     /* flag = 0: delete only pointer array, = 1: delete char. buffers, too */
{
  int i;
  if (pa == NULL)  return NULL;
  if (stamp_flag && pa->stamp != 123456) 
     fprintf(stamp_file, "d_c_p_a double delete --> %s\n", pa->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", pa->name);
  if (flag)
    {
     for (i = 0; i < pa->curr; i++)  if (pa->p[i]) free(pa->p[i]);
    }
  if (pa->p != NULL)  free(pa->p);
  free(pa);
  return NULL;
}

struct command* delete_command(struct command* cmd)
{
  if (cmd == NULL) return NULL;
  if (stamp_flag && cmd->stamp != 123456) 
     fprintf(stamp_file, "d_c double delete --> %s\n", cmd->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", cmd->name);
  if (cmd->par != NULL)  delete_command_parameter_list(cmd->par);
  if (cmd->par_names != NULL) delete_name_list(cmd->par_names);
  free(cmd);
  return NULL; 
}

struct command_list* delete_command_list(struct command_list* cl)
{
  int i;
  if (cl == NULL) return NULL;
  if (stamp_flag && cl->stamp != 123456) 
     fprintf(stamp_file, "d_c_l double delete --> %s\n", cl->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", cl->name);
  if (cl->list != NULL) delete_name_list(cl->list);
  for (i = 0; i < cl->curr; i++) delete_command(cl->commands[i]);
  if (cl->commands) free(cl->commands);
  free(cl);
  return NULL;
}

struct command_parameter* 
       delete_command_parameter(struct command_parameter* par)
{
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
  free(par);
  return NULL; 
}

struct command_parameter_list* 
       delete_command_parameter_list(struct command_parameter_list* parl)
{
  int i;
  if (parl == NULL) return NULL;
  if (stamp_flag && parl->stamp != 123456) 
     fprintf(stamp_file, "d_c_p_l double delete --> %s\n", parl->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", parl->name);
  if (parl->parameters != NULL) 
    {
     for (i = 0; i < parl->curr; i++)
	if (parl->parameters[i] != NULL) 
        parl->parameters[i] = delete_command_parameter(parl->parameters[i]);
     if (parl->parameters)  free(parl->parameters);
    }
  free(parl);
  return NULL;
}

struct constraint* delete_constraint(struct constraint* cst)
{
  if (cst == NULL)  return NULL;
  if (stamp_flag && cst->stamp != 123456) 
     fprintf(stamp_file, "d_c double delete --> %s\n", cst->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", "constraint");
  free(cst);
  return NULL;
}

struct constraint_list* delete_constraint_list(struct constraint_list* cl)
{
  if (cl == NULL)  return NULL;
  if (stamp_flag && cl->stamp != 123456) 
     fprintf(stamp_file, "d_c_l double delete --> %s\n", cl->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", "constraint_list");
  free(cl);
  return NULL;
}

struct double_array* delete_double_array(struct double_array* a)
{
  if (a == NULL)  return NULL;
  if (stamp_flag && a->stamp != 123456) 
     fprintf(stamp_file, "d_d_a double delete --> %s\n", a->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", a->name);
  if (a->a != NULL) free(a->a);
  free(a);
  return NULL;
}

struct element* delete_element(struct element* el)
{
  if (el == NULL)  return NULL;
  if (stamp_flag && el->stamp != 123456) 
     fprintf(stamp_file, "d_e double delete --> %s\n", el->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", el->name);
  free(el);
  return NULL;
}

struct el_list* delete_el_list(struct el_list* ell)
{
  if (ell->list == NULL) return NULL;
  if (stamp_flag && ell->stamp != 123456) 
     fprintf(stamp_file, "d_e_l double delete --> %s\n", ell->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", ell->name);
  delete_name_list(ell->list);
  if (ell->elem != NULL) free(ell->elem);
  free(ell);
  return NULL;
}

struct expression* delete_expression(struct expression* expr)
{
  if (expr == NULL) return NULL;
  if (stamp_flag && expr->stamp != 123456) 
     fprintf(stamp_file, "d_ex double delete --> %s\n", expr->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", expr->name);
  if (expr->polish != NULL) expr->polish = delete_int_array(expr->polish);
  if (expr->string != NULL) free(expr->string);
  free(expr);
  return NULL;
}

struct expr_list* delete_expr_list(struct expr_list* exprl)
{
  int i;
  if (exprl == NULL) return NULL;
  if (stamp_flag && exprl->stamp != 123456) 
     fprintf(stamp_file, "d_ex_l double delete --> %s\n", exprl->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", exprl->name);
  if (exprl->list != NULL)
    {
     for (i = 0; i < exprl->curr; i++)
	if (exprl->list[i] != NULL)  delete_expression(exprl->list[i]);
     free(exprl->list);
    }
  free(exprl);
  return NULL;
}

struct in_cmd* delete_in_cmd(struct in_cmd* cmd)
{
  if (cmd == NULL) return NULL;
  if (stamp_flag && cmd->stamp != 123456) 
     fprintf(stamp_file, "d_i_c double delete --> %s\n", cmd->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", cmd->name);
  if (cmd->tok_list != NULL) 
        cmd->tok_list = delete_char_p_array(cmd->tok_list, 0);
  free(cmd);
  return NULL;
}

struct int_array* delete_int_array(struct int_array* i)
{
  if (i == NULL)  return NULL;
  if (stamp_flag && i->stamp != 123456) 
     fprintf(stamp_file, "d_i_a double delete --> %s\n", i->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", i->name);
  if (i->i != NULL) free(i->i);
  free(i);
  return NULL;
}

struct macro* delete_macro(struct macro* macro)
{
  if (macro == NULL)  return NULL;
  if (stamp_flag && macro->stamp != 123456) 
     fprintf(stamp_file, "d_m double delete --> %s\n", macro->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", macro->name);
  if (macro->formal != NULL) delete_char_p_array(macro->formal, 0);
  if (macro->tokens != NULL) delete_char_p_array(macro->tokens, 0);
  if (macro->body != NULL) delete_char_array(macro->body);
  free(macro);
  return NULL;
}

struct name_list* delete_name_list(struct name_list* l)
{
  if (l == NULL) return NULL;
  if (stamp_flag && l->stamp != 123456) 
     fprintf(stamp_file, "d_n_l double delete --> %s\n", l->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", l->name);
  if (l->index != NULL)  free(l->index);
  if (l->inform != NULL)  free(l->inform);
  if (l->names != NULL)  free(l->names);
  free(l);
  return NULL;
}

struct node* delete_node(struct node* p)
{
  if (p == NULL) return NULL;
  if (stamp_flag && p->stamp != 123456) 
     fprintf(stamp_file, "d_n double delete --> %s\n", p->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", p->name);
  if (p->p_al_err) p->p_al_err = delete_double_array(p->p_al_err);
  if (p->p_fd_err) p->p_fd_err = delete_double_array(p->p_fd_err);
  free(p);
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
  if (l == NULL)  return NULL;
  if (stamp_flag && l->stamp != 123456) 
     fprintf(stamp_file, "d_no_l double delete --> %s\n", l->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", l->name);
  if (l->nodes != NULL)  free(l->nodes);
  if (l->list != NULL)  delete_name_list(l->list);
  free(l);
  return NULL;
}

struct sequence_list* delete_sequence_list(struct sequence_list* sql)
{
  if (sql == NULL) return NULL;
  if (stamp_flag && sql->stamp != 123456) 
     fprintf(stamp_file, "d_s_l double delete --> %s\n", sql->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", sql->name);
  if (sql->list != NULL) delete_name_list(sql->list);
  if (sql->sequs != NULL) free(sql->sequs);
  free(sql);
  return NULL;
}

struct table* delete_table(struct table* t)
{
  int i, j;
  if (t == NULL) return NULL;
  if (stamp_flag && t->stamp != 123456) 
     fprintf(stamp_file, "d_t double delete --> %s\n", t->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", "table");
  if (t->header != NULL) t->header = delete_char_p_array(t->header, 1);
  if (t->col_out != NULL) t->col_out = delete_int_array(t->col_out);
  if (t->row_out != NULL) t->row_out = delete_int_array(t->row_out);
  t->node_nm = delete_char_p_array(t->node_nm, 0);
  for (i = 0; i < t->curr; i++)
    {
     if (t->l_head[i] != NULL) 
        t->l_head[i] = delete_char_p_array(t->l_head[i], 1);
    }
  if (t->l_head)  free(t->l_head);
  if (t->p_nodes) free(t->p_nodes);
  if (t->d_cols)
    {
     for (i = 0; i < t->num_cols; i++)
        if (t->columns->inform[i] < 3 && t->d_cols[i]) free(t->d_cols[i]);
     free(t->d_cols);
    }
  if (t->s_cols)
    {
     for (i = 0; i < t->num_cols; i++)
       {
        if (t->columns->inform[i] == 3 && t->s_cols[i])
	  {
	   for (j = 0; j < t->curr; j++) 
              if (t->s_cols[i][j]) free(t->s_cols[i][j]);
           free(t->s_cols[i]);
	  }
       }
     free (t->s_cols);
    }
  t->columns = delete_name_list(t->columns);
  free(t);
  return NULL; 
}

struct variable* delete_variable(struct variable* var)
{
  if (var == NULL)  return NULL;
  if (stamp_flag && var->stamp != 123456) 
     fprintf(stamp_file, "d_v double delete --> %s\n", var->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", var->name);
  if (var->expr != NULL) delete_expression(var->expr);
  if (var->string != NULL) free(var->string);
  free(var);
  return NULL;
}

struct var_list* delete_var_list(struct var_list* varl)
{
  if (varl == NULL) return NULL;
  if (stamp_flag && varl->stamp != 123456) 
     fprintf(stamp_file, "d_v_l double delete --> %s\n", varl->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", varl->name);
  if (varl->list != NULL) delete_name_list(varl->list);
  if (varl->vars != NULL) free(varl->vars);
  free(varl);
  return NULL;
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
  free(p_loc);
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
  free(p_loc);
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
  free(c_loc);
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
  free(c_loc);
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
  free(c_loc);
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
  free(c_loc);
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
  free(a_loc);
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
  free(e_loc);
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
  free(e_loc);
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
  free(e_loc);
  p->input_files = mycalloc(rout_name, new, sizeof(FILE*));
  for (j = 0; j < p->curr; j++) p->input_files[j] = f_loc[j];
  free(f_loc);
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
  free(c_loc);
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
  free(i_loc);
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
  free(n_loc);
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
  free(n_loc);
  free(l_ind);
  free(l_inf);
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
  free(n_loc);
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
  free(sloc);
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
  free(pa_loc);
  t->node_nm->curr = t->curr; free(t_loc);
  for (j = 0; j < t->num_cols; j++)
    {
     if ((s_loc = t->s_cols[j]) != NULL)
       {
        t->s_cols[j] = (char**) mycalloc(rout_name,new, sizeof(char*));
        for (i = 0; i < t->curr; i++) t->s_cols[j][i] = s_loc[i];
        free(s_loc);
       }
    }
  for (j = 0; j < t->num_cols; j++)
    {
     if ((d_loc = t->d_cols[j]) != NULL)
       {
        t->d_cols[j] = (double*) mycalloc(rout_name,new, sizeof(double));
        for (i = 0; i < t->curr; i++) t->d_cols[j][i] = d_loc[i];
        free(d_loc);
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
  free(t_loc);
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
  free(v_loc);
}

void* mycalloc(char* caller, size_t nelem, size_t size)
{
  /* calls calloc, checks for memory granted */
  void* p;
  if ((p = calloc(nelem, size)) == NULL)
     fatal_error("memory overflow, called from routine:", caller);
  return p;
}

void mycpy(char* sout, char* sin)
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

void* mymalloc(char* caller, size_t size)
{
  /* calls malloc, checks for memory granted */
  void* p;
  if ((p = malloc(size)) == NULL) 
    fatal_error("memory overflow, called from routine:", caller);
  return p;
}

char* mystrchr(char* string, char c)
     /* returns strchr, but only outside strings */
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

char* mystrstr(char* string, char* s)
     /* returns strstr, but only outside strings */
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

struct char_array* new_char_array(int length)
{
  char rout_name[] = "new_char_array";
  struct char_array* il = 
       (struct char_array*) mycalloc(rout_name,1, sizeof(struct char_array));
  strcpy(il->name, "char_array");
  il->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", il->name);
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
  return il;
}

struct command* new_command(char* name, int nl_length, int pl_length,
                            char* module, char* group, int link, int mad_8)
{
  char rout_name[] = "new_command";
  struct command* new 
   = (struct command*) mycalloc(rout_name,1, sizeof(struct command));
  new->stamp = 123456;
  strcpy(new->name, name); 
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", new->name);
  strcpy(new->module, module);
  strcpy(new->group, group);
  new->link_type = link;
  new->mad8_type = mad_8;
  if (nl_length == 0) nl_length = 1;
  new->par_names = new_name_list(nl_length);
  new->par = new_command_parameter_list(pl_length);
  return new;
}

struct command_list* new_command_list(int length)
{
  char rout_name[] = "new_command_list";
  struct command_list* il = 
    (struct command_list*) mycalloc(rout_name,1, sizeof(struct command_list));
  strcpy(il->name, "command_list");
  il->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", il->name);
  il->curr = 0;
  il->max = length;
  il->list = new_name_list(length);
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
  il->list = new_name_list(length);
  il->command_lists 
   = (struct command_list**) 
     mycalloc(rout_name,length, sizeof(struct command_list*));
  return il;
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
  strcpy(il->name, "double_array");
  il->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", il->name);
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
  ell->list = new_name_list(length);
  ell->elem  
     = (struct element**) mycalloc(rout_name,length, sizeof(struct element*));
  ell->max = length;
  return ell;
}

struct node* new_elem_node(struct element* el, int occ_cnt)
{
  struct node* p;
  p = new_node(compound(el->name, occ_cnt));
  p->p_elem = el;
  p->length = el->length;
  p->base_name = el->base_type->name;
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
  il->labels = new_name_list(length);
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
  m->body = new_char_array(++length);
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
  nll->list = new_name_list(length);
  nll->macros  
     = (struct macro**) mycalloc(rout_name,length, sizeof(struct macro*));
  nll->max = length;
  return nll;
}

struct name_list* new_name_list(int length)
{
  char rout_name[] = "new_name_list";
  struct name_list* il = 
       (struct name_list*) mycalloc(rout_name,1, sizeof(struct name_list));
  strcpy(il->name, "name_list");
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
  nll->list = new_name_list(length);
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
  s->list = new_name_list(length);
  s->sequs 
    = (struct sequence**) mycalloc(rout_name,length, sizeof(struct sequence*));
  return s;
}

struct node* new_sequ_node(struct sequence* sequ, int occ_cnt)
{
  struct node* p;
  p = new_node(compound(sequ->name, occ_cnt));
  p->p_sequ = sequ;
  p->length = sequ->length;
  p->base_name = permbuff("sequence");
  return p;
}

struct table* new_table(char* name, int rows, struct name_list* cols)
{
  char rout_name[] = "new_table";
  int i, n = cols->curr;
  struct table* t 
     = (struct table*) mycalloc(rout_name,1, sizeof(struct table));
  strcpy(t->name, name);
  t->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", "table");
  t->columns = cols;
  t->num_cols = n;
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
  tl->names = new_name_list(size);
  tl->tables 
    = (struct table**) mycalloc(rout_name,size, sizeof(struct table*));
  return tl;
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
  var->list = new_name_list(length);
  var->vars  
    = (struct variable**) mycalloc(rout_name,length, sizeof(struct variable*));
  var->max = length;
  return var;
}
