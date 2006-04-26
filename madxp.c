/* Production version of MAD-X, version number: see madxd.h */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>
#include <sys/types.h>
#ifndef _WIN32
#include <sys/utsname.h>
#include <unistd.h>
#endif
#include <sys/timeb.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#ifdef _CATCH_MEM
#include <signal.h>
#endif
#include "madxl.h"
#include "madx.h"
#include "madxreg.h"
#include "madxd.h"
#include "madxdict.h"

/* JMJ 7/11/2002 moved this here from c6t.c */
/* FS & TdA 15.03.2004 plot upgrade, bugs correction, ptc_twiss upgrade, touschek preparation */
#include "c6t.h"

static const int kSkowronDebug = 0;

void madx()
{
#ifdef _CATCH_MEM
  /* provide a termination routine for access to memory outside scope */
  if (signal(SIGSEGV, termination_handler) == SIG_IGN)
    signal(SIGSEGV, SIG_IGN);
#endif
  madx_start();
  madx_init();
  main_input(0);
  madx_finish();
}

#ifdef _FULL
#include "madxn.c"
#endif

#include "madxu.c"

#include "madxreg.c"

int act_special(int type, char* statement)
  /* acts on special commands (IF{..} etc.) */
{
  struct char_array* loc_buff = NULL;
  struct char_array* loc_w = NULL;
  int cnt_1, start_2, rs, re, level = pro->curr, ls = strlen(statement);
  int ret_val = 0;
  struct char_p_array* logic;
  int logex = 0;
  char *cp = statement;
  if (ls < IN_BUFF_SIZE) ls = IN_BUFF_SIZE;
  if (level == pro->max) grow_in_buff_list(pro);
  if (pro->buffers[level] == NULL)
    pro->buffers[level] = new_in_buffer(ls);
  else
  {
    while(pro->buffers[level]->c_a->max < ls)
      grow_char_array(pro->buffers[level]->c_a);
  }
  if (type == 5) /* macro */ return make_macro(statement);
  else if (type == 6) /* line */ return make_line(statement);
  logic = new_char_p_array(1000);
  loc_buff = new_char_array(ls);
  loc_w = new_char_array(ls);
  get_bracket_range(statement, '{', '}', &rs, &re);
  if (re < 0) fatal_error("missing '{' or '}' in statement:",statement);
  cnt_1 = rs; start_2 = rs + 1;
  mystrcpy(loc_buff, statement); loc_buff->c[re] =  '\0';
  while(aux_buff->max < cnt_1) grow_char_array(aux_buff);
  strncpy(aux_buff->c, statement, cnt_1); aux_buff->c[cnt_1] = '\0';
  switch (type)
  {
    case 1:  /* if */
      pro->buffers[level]->flag = 0;
    case 3:  /* else if */
      if (pro->buffers[level]->flag < 0)
      {
        ret_val = -1;
        break;
      }
      if (pro->buffers[level]->flag == 0)
      {
        pre_split(aux_buff->c, loc_w, 0);
        mysplit(loc_w->c, tmp_l_array);
        get_bracket_t_range(tmp_l_array->p, '(', ')', 0, tmp_l_array->curr,
                            &rs, &re);
        rs++;
        if ((logex = logic_expr(re-rs, &tmp_l_array->p[rs])) > 0)
        {
          pro->buffers[level]->flag = 1;
          pro->curr++;
          /* now loop over statements inside {...} */
          pro_input(&loc_buff->c[start_2]);
          pro->curr--;
        }
        else if (logex < 0) warning("illegal if construct set false:", cp);
      }
      break;
    case 2: /* else */
      if (pro->buffers[level]->flag < 0)
      {
        ret_val = -1;
        break;
      }
      if (pro->buffers[level]->flag == 0)
      {
        pro->curr++;
        /* now loop over statements inside {...} */
        pro_input(&loc_buff->c[start_2]);
        pro->curr--;
        pro->buffers[level]->flag = -1;
      }
      break;
    case 4: /* while */
      pre_split(aux_buff->c, loc_w, 0);
      mysplit(loc_w->c, logic);
      get_bracket_t_range(logic->p, '(', ')', 0, logic->curr,
                          &rs, &re);
      pro->curr++; rs++;
      while ((logex = logic_expr(re-rs, &logic->p[rs])) > 0)
      {
        /* now loop over statements inside {...} */
        pro_input(&loc_buff->c[start_2]);
      }
      pro->curr--;
      break;
    default:
      ret_val = -1;
  }
  if (loc_buff != NULL) delete_char_array(loc_buff);
  if (loc_w != NULL) delete_char_array(loc_w);
  delete_char_p_array(logic, 0);
  return ret_val;
}

double act_value(int pos, struct name_list* chunks)
  /* returns the actual value of a variable, element, or command parameter */
{
  char* name = chunks->names[pos];
  char comm[NAME_L];
  char par[NAME_L];
  char *p, *n = name, *q = comm;
  double val = zero;
  struct element* el;
  struct command* cmd = NULL;
  if ((p = strstr(name, "->")) == NULL) /* variable */
  {
    if ((current_variable = find_variable(name, variable_list)) == NULL)
    {
      if (get_option("verify"))
        warning("undefined variable set to zero:", name);
      current_variable = new_variable(name, zero, 1, 1, NULL, NULL);
      val = zero;
      add_to_var_list(current_variable, variable_list, 0);
    }
    else val = variable_value(current_variable);
  }
  else /* element or command parameter */
  {
    while (n < p)  *(q++) = *(n++);
    *q = '\0';
    q = par; n++; n++;
    while (*n != '\0')  *(q++) = *(n++);
    *q = '\0';
    if (strncmp(comm, "beam", 4) == 0)
    {
      cmd = current_beam = find_command("default_beam", beam_list);
      if ((p = strchr(comm, '%')) != NULL)
      {
        if ((current_beam = find_command(++p, beam_list)) == NULL)
          current_beam = cmd;
      }
      val = command_par_value(par, current_beam);
    }
    else if ((el = find_element(comm, element_list)) != NULL)
      val = el_par_value(par, el);
    else if ((cmd = find_command(comm, stored_commands)) != NULL)
      val = command_par_value(par, cmd);
    else if ((cmd = find_command(comm, beta0_list)) != NULL)
      val = command_par_value(par, cmd);
    else if ((cmd = find_command(comm, defined_commands)) != NULL)
      val = command_par_value(par, cmd);
  }
  return val;
}

char* alias(char* par_string) /* returns main parameter for alias */
{
  if (strcmp(par_string, "filename") == 0)  return file_string;
  else return par_string;
}

void all_node_pos(struct sequence* sequ)
  /* calculates all node positions in an expanded sequence */
{
  struct node* node = sequ->start;
  while (node != NULL)
  {
    if (node->p_elem != NULL)
      node->length = node->p_elem->length
        = element_value(node, "l");
    else if (node->p_sequ != NULL)
      node->length = node->p_sequ->length;
    else fatal_error("node is neither element nor sequence:",
                     node->name);
    if ((node->position = get_node_pos(node, sequ)) < zero)
      node->position += sequ->length;
    if (node == sequ->end) break;
    node = node->next;
  }
}

int aperture_count(struct sequence* sequ)
  /* returns max. number of aperture parameters needed for sequence */
{
  int i = 0, n = 0;
  char* p;
  struct node* node = sequ->start;
  while (node != NULL)
  {
    if ((node->p_elem)&&(p = command_par_string("apertype", node->p_elem->def)))
    {
      while(aperture_types[i][0] != ' ')
      {
        if (strcmp(p, aperture_types[i]) == 0)
        {
          if (n < aperture_npar[i]) n = aperture_npar[i];
          break;
        }
        i++;
      }
    }
    if (node == sequ->end) break;
    node = node->next;
  }
  return n;
}

int belongs_to_class(struct element* el, char* class)
  /* returns 1 if an element belongs to a class, else 0 */
{
  int in = 0;
  if (strcmp(el->name, class) == 0) in = 1;
  else
  {
    while (el->parent != el)
    {
      if (strcmp(el->parent->name, class) == 0)
      {
        in = 1; break;
      }
      el = el->parent;
    }
  }
  return in;
}

struct in_cmd* buffered_cmd(struct in_cmd* cmd)
  /* returns a buffered command if found */
{
  int k;
  if ((k = name_list_pos(cmd->tok_list->p[cmd->decl_start],
                         buffered_cmds->labels)) > -1)
    return buffered_cmds->in_cmds[k];
  else return cmd;
}

void buffer_in_cmd(struct in_cmd* cmd)
  /* stores an input command in a buffer */
{
  int i;
  if (buffered_cmds->curr == buffered_cmds->max)
    grow_in_cmd_list(buffered_cmds);
  cmd->label = permbuff(cmd->label);
  add_to_name_list(cmd->label, 0, buffered_cmds->labels);
  buffered_cmds->in_cmds[buffered_cmds->curr++] = cmd;
  for (i = 0; i < cmd->tok_list->curr; i++)
    cmd->tok_list->p[i] = permbuff(cmd->tok_list->p[i]);
}

void check_table(char* string)
  /* replaces argument of "table" if any by a string variable */
{

  char *pa, *pb, *pt, *pl, *pr, *sv, *quote;
  char* qpos[1000];
  int npos = 0;
  pa = string;
  while (*pa) /* find all strings in quotes, keep positions */
  {
    if (*pa == '\'' || *pa == '\"')
    {
      quote = pa;
      qpos[npos++] = pa++;
      if (!(pa = strchr(pa, *quote))) return; /* quote string not closed */
      qpos[npos++] = pa;
    }
    pa++;
  }
  pa = string;
  while ((pb = strstr(pa, "table")) != NULL)
  {
    if (!inbounds(pb, npos, qpos))
    {
      mystrcpy(c_join, pa);
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
    pa = ++pb;
  }
}

int cmd_match(int cnt, char** toks, int* cmd_pos, int* decl_start)
  /* matches input (user) command with the symbolic description
     from the list in madxl.h and returns is type */
{
  int i, i2, j, k, lp;
  for (i = 0; i < n_match; i++)
  {
    k = 0; lp = -1;
    i2 = t_match[i];
    for (j = s_match[i2]; j < s_match[i2+1]; j++)
    {
      if (k == cnt)  break;
      if (*cmd_match_base[j] == '@')
      {
        if (strcmp(cmd_match_base[j], "@cmd") == 0)
        {
          if ((lp = name_list_pos(toks[k],
                                  defined_commands->list)) < 0) break;
        }
        else if (isalpha(*toks[k]) == 0) break;
      }
      else if (strcmp(cmd_match_base[j], toks[k]) != 0)  break;
      k++;
    }
    
    if (j == s_match[i2+1]) goto found;
  }
  return -3;
  found:
  *cmd_pos = lp; *decl_start = s_match[i2+1] - s_match[i2];
  
  
  return i2;
}

struct expression* command_par_expr(char* parameter, struct command* cmd)
  /* returns a command parameter expression if found, else NULL */
{
  struct expression* expr = NULL;
  int i;
  if ((i = name_list_pos(parameter, cmd->par_names)) > -1)
    expr = cmd->par->parameters[i]->expr;
  return expr;
}

char* command_par_string(char* parameter, struct command* cmd)
  /* returns a command parameter string if found, else NULL */
{
  struct command_parameter* cp;
  char* p = NULL;
  int i;
  if ((i = name_list_pos(parameter, cmd->par_names)) > -1)
  {
    cp = cmd->par->parameters[i];
    if (cp->type == 3) p = cp->string;
  }
  return p;
}

int command_par_value2(char* parameter, struct command* cmd, double* val)
  /* returns a command parameter value val 
     if found returns 1, else 0 */
{
  struct command_parameter* cp;
  int i;
  int ret = 0;
  
  *val = zero;
  if ((i = name_list_pos(parameter, cmd->par_names)) > -1)
  {
    cp = cmd->par->parameters[i];
    if (cp->type < 3)
    {
      if (cp->expr == NULL)  *val = cp->double_value;
      else *val = expression_value(cp->expr, 2);
      ret = 1;
    }
  }
  
  return ret;
}


double command_par_value(char* parameter, struct command* cmd)
  /* returns a command parameter value if found, else zero */
{
  struct command_parameter* cp;
  double val = zero;
  int i;
  if ((i = name_list_pos(parameter, cmd->par_names)) > -1)
  {
    cp = cmd->par->parameters[i];
    if (cp->type < 3)
    {
      if (cp->expr == NULL)  val = cp->double_value;
      else val = expression_value(cp->expr, 2);
    }
  }
  return val;
}

char* compound(char* e_name, int occ)
  /* makes node name from element name and occurrence count */
{
  sprintf(c_dum->c,"%s:%d", e_name, occ);
  return c_dum->c;
}

void control(struct in_cmd* cmd)
  /* executes so-called "control" commands */
{
  char** toks = cmd->tok_list->p;
  int k = cmd->decl_start - 1;
  if      (strcmp(toks[k], "assign")      == 0) exec_assign(cmd);
  else if (strcmp(toks[k], "beam")        == 0) exec_beam(cmd, 0);
  else if (strcmp(toks[k], "call")        == 0) exec_call(cmd);
  else if (strcmp(toks[k], "option")      == 0) exec_option();
  else if (strcmp(toks[k], "resbeam")     == 0) exec_beam(cmd, 1);
  else if (strcmp(toks[k], "save")        == 0) exec_save(cmd);
#ifdef _FULL
  else if (strcmp(toks[k], "dumpsequ")    == 0) exec_dumpsequ(cmd);
  else if (strcmp(toks[k], "set")         == 0) store_set(cmd->clone, 1);
  else if (strcmp(toks[k], "sodd")        == 0) exec_sodd(cmd);
  else if (strcmp(toks[k], "threader")    == 0) store_threader(cmd);
  else if (strcmp(toks[k], "use")         == 0) use_sequ(cmd);
  else if (strcmp(toks[k], "write")       == 0) exec_dump(cmd);
  else if (strcmp(toks[k], "beta0")       == 0) store_beta0(cmd);
  else if (strcmp(toks[k], "coguess")     == 0) exec_store_coguess(cmd);
  else if (strcmp(toks[k], "create")      == 0) exec_create_table(cmd);
  else if (strcmp(toks[k], "fill")        == 0) exec_fill_table(cmd);
  else if (strcmp(toks[k], "plot")        == 0) exec_plot(cmd);
  else if (strcmp(toks[k], "print")       == 0) exec_print(cmd);
  else if (strcmp(toks[k], "readtable")   == 0) read_table(cmd);
  else if (strcmp(toks[k], "savebeta")    == 0) store_savebeta(cmd);
  else if (strcmp(toks[k], "select")      == 0) store_select(cmd);
  else if (strcmp(toks[k], "deselect")    == 0) store_deselect(cmd);
#endif
#ifndef _FULL
  puts("++++++++++++++ command skipped in parser version");
  /* insert your proper command action here */
#endif
}

int decode_command () /* compares command with templates, fills this_cmd
                         return: 0 command from list (except below);
                         1 element definition (as well inside sequence);
                         2 variable definition;
                         3 sequence or endsequence;
                         4 element reference inside sequence (no def.);
                         -1 illegal in this context or form;
                         -2 label is protected keyword;
                         -3 statement not recognised */
{

  
  int i, aux_pos, cmd_pos, decl_start, type;
  int n = this_cmd->tok_list->curr;
  char** toks = this_cmd->tok_list->p;
  this_cmd->type = -3;
  if ((i = cmd_match(n, toks, &cmd_pos, &decl_start)) < 0)  return i;
  

  this_cmd->sub_type = i;
  this_cmd->decl_start = decl_start;
  switch (i)
  {
    case 0:
      if (n > 1 && *toks[1] == ':') /* label is (protected) key */ return -2;
      this_cmd->cmd_def = defined_commands->commands[cmd_pos];
      if (strcmp(toks[0], "endsequence") == 0)
      {
        if (sequ_is_on == 0) return -1;
        this_cmd->type = 3;
        sequ_is_on = 0;
      }
      else if (strcmp(this_cmd->cmd_def->module, "element") == 0
               || strcmp(this_cmd->cmd_def->module, "sequence") == 0) return -1;
      else this_cmd->type = 0;
      break;
    case 1:
    case 16:
      if (i == 1) aux_pos = 0;
      else          aux_pos = 1;
      this_cmd->cmd_def = defined_commands->commands[cmd_pos];
      this_cmd->type = 0;
      this_cmd->label = permbuff(toks[aux_pos]);
      if (strcmp(this_cmd->cmd_def->module, "element") == 0)
        this_cmd->type = 1; /* element definition */
      else if (strcmp(this_cmd->cmd_def->module, "sequence") == 0)
      {
        if (strcmp(toks[aux_pos+2], "sequence") == 0)
        {
          this_cmd->type = 3;
          sequ_is_on = 1;
        }
        else return -1;
      }
      type = this_cmd->cmd_def->link_type;
      if (type == 1)
      {
        if (group_is_on || sequ_is_on) return -1; /* group in group illegal */
        group_is_on = 1;
        current_link_group = this_cmd->cmd_def->group;
      }
      else if (group_is_on)
      {
        if (strcmp(none, this_cmd->cmd_def->group) != 0
            && strcmp(current_link_group, this_cmd->cmd_def->group) != 0)
          return -1; /* command does not belong to this group */
        if (type == 2)
        {
          current_link_group = none; /* end of group */
          group_is_on = 0;
        }
      }
      break;
    case 14:
      this_cmd->type = 1;
      this_cmd->decl_start++; /* skip (as yet unknown) element class pos. */
      break;
    case 15:
      this_cmd->type = 4; /* element or sequence reference inside sequence
                             or element parameter definition */
      this_cmd->decl_start--; /* second name is already declaration start */
      break;
    default:
      type = 0;
      this_cmd->type = 2;
  }
  return this_cmd->type;
}

int decode_par(struct in_cmd* cmd, int start, int number, int pos, int log)
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
          if     (strcmp(toks[i+2], "true") == 0)  ival = 1;
          else if(strcmp(toks[i+2], "false") == 0) ival = 0;
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
          if (name_list_pos(alias(toks[j]),
                            cmd->cmd_def->par_names) >= 0) break;
        if (*toks[j-1] == ',') j--;
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
    else if (lp->type == 4)  /* one constraint */
    {
      if (i+1 < number && *toks[i+1] == ':')
      {
        con_flag = 1; i++;
      }
      else con_flag = 0; /* if != zero, := or :< or :> */
      if (i+2 < number)
      {
        if (*toks[i+1] == '=') c_type = 4;
        else if (*toks[i+1] == '>') c_type = 1;
        else if (*toks[i+1] == '<') c_type = 2;
        if (c_type)
        {
          start = i + 2;
          for (t_num = start; t_num < number; t_num++) if(*toks[t_num] == ',')
            break;
          if ((e_type = loc_expr(toks, t_num, start, &end)) == 0)
            return -i;
          tot_end = end;
          if (e_type == 1) /* simple number */
          {
            val = simple_double(toks, start, end);
            expr = NULL;
          }
          else /* expression */
          {
            if ((expr =
                 make_expression(end + 1 - start, &toks[start])) == NULL)
              return -i;
            val = expression_value(expr, 2);
          }
        }
      }
      else
      {
        c_type = 4;
        val = zero;
        expr = NULL;
        tot_end = i;
      }
    }
    if (clp->c_type == 1 && c_type == 2) clp->c_type = 3; /* min present */
    else if (clp->c_type == 2 && c_type == 1) clp->c_type = 3; /* max */
    else clp->c_type = c_type;
    if (con_flag == 0)  expr = NULL;
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

void deco_init()
  /* initializes Polish decoding */
{
  expr_chunks = new_name_list("expr_chunks", 2000);
  cat = new_int_array(MAX_ITEM);
  deco = new_int_array(MAX_ITEM);
  d_var = new_int_array(MAX_ITEM);
  oper = new_int_array(MAX_ITEM);
  func = new_int_array(MAX_ITEM);
  cat_doubles = new_double_array(MAX_ITEM);
  doubles = new_double_array(MAX_D_ITEM);
  twiss_deltas = new_double_array(MAX_ITEM);
}

double element_value(struct node* node, char* par)
  /* all element parameter values except vectors are modified here
     resp. in el_par_value if any */
{
  double e_val;
  struct element* el = node->p_elem;
  if (strcmp(par, "mad8_type") == 0) e_val = el->def->mad8_type;
  else  e_val = el_par_value(par, el);
  return e_val;
}

int down_unit(char* file_name)
  /* makes a called file the current input unit */
{
  FILE* new;
  if ((new = fopen(file_name, "r")) == NULL)
  {
    if (interactive) warning("cannot open input file:", file_name);
    else             fatal_error("cannot open input file:", file_name);
    return 0;
  }
  if (in->curr+1 == in->max) grow_in_buff_list(in);
  in->input_files[++in->curr] = new;
  return 1;
}

int element_vector(struct element* el, char* par, double* vector)
  /* returns length + vector of parameter par for element el */
{
  int i, l = 0;
  struct double_array* da;
  struct expr_list* ell;
  if ((i = name_list_pos(par, el->def->par_names)) > -1)
  {
    if ((da = el->def->par->parameters[i]->double_array) != NULL)
    {
      if ((ell = el->def->par->parameters[i]->expr_list) != NULL)
        update_vector(ell, da);
      l = da->curr;
      copy_double(da->a, vector, l);
    }
  }
  return l;
}

double el_par_value(char* par, struct element* el)
  /* returns an element parameter value */
{
  int k = 0, n;
  char tmp[8];
  double val = zero, angle = zero, l, vec[100];
  double fact = strcmp(el->base_type->name, "rbend") == 0 ? one : zero;
  int mult = strcmp(el->base_type->name, "multipole") == 0 ? 1 : 0;
  if (fact != zero || strcmp(el->base_type->name, "sbend") == 0) /* bend */
  {
    if ((l = command_par_value("l", el->def)) == zero)
      fatal_error("bend with zero length:",el->name);
    angle = command_par_value("angle", el->def);
    if (strcmp(par, "angle") == 0)  val = angle;
    else if (strcmp(par, "tilt") == 0)
      val = command_par_value("tilt", el->def);
    else if (strcmp(par, "k0") == 0) val = command_par_value("k0", el->def);
    else if (strcmp(par, "k0s") == 0) val = command_par_value("k0s", el->def);
    else if (strcmp(par, "l") == 0)
    {
      if (fact != zero && get_option("rbarc") && angle != zero)
        val = l * angle / (two * sin(angle/two));
      else val = l;
    }
    else if (strcmp(par, "e1") == 0)
      val = command_par_value("e1", el->def) + fact * angle / two;
    else if (strcmp(par, "e2") == 0)
      val = command_par_value("e2", el->def) + fact * angle / two;
    else if (strcmp(par, "rhoinv") == 0) val = angle / l;
    else if (strcmp(par, "blen") == 0) val = l;
    else val = command_par_value(par, el->def);
  }
  /* all elements except bends */
  else if (strcmp(par, "rhoinv") == 0) val = zero;
  else if (strcmp(par, "blen") == 0) val = zero;
  else if (mult)  /* multipole */
  {
    if (strcmp(par, "l") == 0) val = zero;
    else if (par[0] == 'k' && isdigit(par[1]) && par[strlen(par)-1] == 'l')
      /* single component requested for multipole */
    {
      if (strchr(par, 's')) strcpy(tmp, "ksl");
      else                  strcpy(tmp, "knl");
      sscanf(&par[1], "%d", &k);
      if ((n = element_vector(el, tmp, vec)) > k)  val = vec[k];
    }
    else val = command_par_value(par, el->def);
  }
  else val = command_par_value(par, el->def);
  /* extra code for kickers */
  if (val == zero && strcmp(el->base_type->name, "hkicker") == 0)
  {
    if (strcmp(par,"hkick") == 0)
      val = command_par_value("kick", el->def);
    else if (strcmp(par,"kick") == 0)
      val = command_par_value("hkick", el->def);
  }
  else if (val == zero && strcmp(el->base_type->name, "vkicker") == 0)
  {
    if (strcmp(par,"vkick") == 0)
      val = command_par_value("kick", el->def);
    else if (strcmp(par,"kick") == 0)
      val = command_par_value("vkick", el->def);
  }
  return val;
}

void enter_element(struct in_cmd* cmd)
  /* enters an element in the list (and the sequence if applicable) */
{
  struct name_list* nl;
  struct command_parameter_list* pl;
  char** toks = cmd->tok_list->p;
  struct element *el, *parent;
  struct command* comm;
  int pos, flag = 0, k = cmd->type == 1 ? 2 : 0;
  if ((parent = find_element(toks[k], element_list)) == NULL)
  {
    fatal_error("unknown class type:", toks[k]);
  }
  else
  {
    cmd->cmd_def = parent->def;
    cmd->clone = clone_command(cmd->cmd_def);
    strcpy(cmd->clone->name, toks[0]);
    scan_in_cmd(cmd);
    if (k == 0 || strcmp(toks[0], toks[2]) == 0) el = parent;
    else
    {
      if ((el = make_element(toks[0], parent->name,
                             cmd->clone, 1+sequ_is_on)) == NULL) return;
      el->def_type = sequ_is_on;
      flag = 1; /* new element - definition only once in sequence allowed */
    }
    nl = cmd->clone->par_names;
    pl = cmd->clone->par;
    if (el != parent && (pos = name_list_pos("bv", nl)) > -1)
    {
      if (nl->inform[pos]) el->bv = command_par_value("bv", cmd->clone);
      else if ((comm = find_command(el->parent->name, defined_commands))
               != NULL && (pos = name_list_pos("bv", comm->par_names)) > -1)
        el->bv = command_par_value("bv", comm);
      else el->bv = parent->bv;
    }
    if (sequ_is_on) enter_elm_reference(cmd, el, flag);
  }
}

void enter_elm_reference(struct in_cmd* cmd, struct element* el, int flag)
  /* enters an element in a sequence */
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int i, pos, k = 1;
  double at;
  if (strcmp(el->base_type->name, "rfcavity") == 0 &&
      find_element(el->name, current_sequ->cavities) == NULL)
    add_to_el_list(&el, 0, current_sequ->cavities, 0);
  if (nl->inform[name_list_pos("at", nl)] == 0)
    fatal_error("element reference without 'at':",
                join(cmd->tok_list->p, cmd->tok_list->curr));
  at = command_par_value("at", cmd->clone);
  if ((i = name_list_pos(el->name, occ_list)) < 0)
    i = add_to_name_list(el->name, k, occ_list);
  else if (flag)
    fatal_error("multiple element definition inside sequence:", el->name);
  else k = ++occ_list->inform[i];
  make_elem_node(el, k);
  current_node->at_value = at;
  current_node->at_expr = command_par_expr("at", cmd->clone);
  pos = name_list_pos("from", nl);
  if (nl->inform[pos])
    current_node->from_name = permbuff(pl->parameters[pos]->string);
}

void enter_sequence(struct in_cmd* cmd)
  /* handles sequence start and end on input */
{
  struct name_list* nl;
  struct command_parameter_list* pl;
  int i, k = 0, pos, aux_pos;
  char** toks = cmd->tok_list->p;
  struct element* el;
  struct command* clone;
  aux_pos = strcmp(toks[0], "shared") == 0 ? 1 : 0;
  if (strcmp(toks[0], "endsequence") == 0)
  {
    pos = name_list_pos("marker", defined_commands->list);
    clone = clone_command(defined_commands->commands[pos]);
    sprintf(c_dum->c, "%s$end", current_sequ->name);
    el = make_element(c_dum->c, "marker", clone, 0);
    make_elem_node(el, 1);
    current_node->at_value = current_sequ->length;
    current_sequ->end = current_node;
    current_sequ->start->previous = current_sequ->end;
    current_sequ->end->next = current_sequ->start;
  }
  else if (strcmp(toks[aux_pos+2], "sequence") == 0)
  {
    for (i = aux_pos+3; i < cmd->tok_list->curr; i++)
    {
      if (strcmp(toks[i], "refer") == 0)
      {
        if (i+2 < cmd->tok_list->curr)
        {
          if (strcmp(toks[i+2], "entry") == 0)  k = 1;
          else if (strcmp(toks[i+2], "exit") == 0)  k = -1;
        }
        break;
      }
    }
    if ((pos = name_list_pos(toks[aux_pos], sequences->list)) >= 0)
    {
      remove_from_sequ_list(sequences->sequs[pos], sequences);
      delete_sequence(sequences->sequs[pos]);
    }
    current_sequ = new_sequence(toks[aux_pos], k);
    add_to_sequ_list(current_sequ, sequences);
    cmd->clone = clone_command(cmd->cmd_def);
/* prevent a line with this name from expansion */
    disable_line(current_sequ->name, line_list);
    scan_in_cmd(cmd);
    nl = cmd->clone->par_names;
    pl = cmd->clone->par;
    if ((current_sequ->length
         = command_par_value("l", cmd->clone)) == zero)
      fatal_error("missing length for sequence:", toks[aux_pos]);
    current_sequ->l_expr = command_par_expr("l", cmd->clone);
    pos = name_list_pos("refpos", nl);
    if (nl->inform[pos])
      current_sequ->refpos = permbuff(pl->parameters[pos]->string);
    current_node = NULL;
    if (occ_list == NULL)
      occ_list = new_name_list("occ_list", 10000);  /* for occurrence count */
    else occ_list->curr = 0;
    if (current_sequ->cavities != NULL)  current_sequ->cavities->curr = 0;
    else current_sequ->cavities = new_el_list(100);
    pos = name_list_pos("marker", defined_commands->list);
    clone = clone_command(defined_commands->commands[pos]);
    sprintf(c_dum->c, "%s$start", current_sequ->name);
    el = make_element(c_dum->c, "marker", clone, 0);
    make_elem_node(el, 1);
    current_sequ->start = current_node;
    current_sequ->share = aux_pos;
  }
}

void enter_variable(struct in_cmd* cmd) /* stores variable contained in cmd */
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
          expr = new_expression(join(&cmd->tok_list->p[start],
                                     end + 1 - start), deco);
          val = expression_value(expr, type);
        }
        else
        {
          expr = NULL;
          val = polish_value(deco);
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

void enter_sequ_reference(struct in_cmd* cmd, struct sequence* sequ)
  /* enters a sequence in a sequence */
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int i, pos, k = 1;
  double at;
  if (nl->inform[name_list_pos("at", nl)] == 0)
    fatal_error("sequence reference without 'at':",
                join(cmd->tok_list->p, cmd->tok_list->curr));
  at = command_par_value("at", cmd->clone);
  if ((i = name_list_pos(cmd->tok_list->p[0], occ_list)) < 0)
    i = add_to_name_list(permbuff(cmd->tok_list->p[0]), k, occ_list);
  else k = ++occ_list->inform[i];
  make_sequ_node(sequ, k);
  current_node->at_value = at;
  current_node->at_expr = command_par_expr("at", cmd->clone);
  pos = name_list_pos("from", nl);
  if (nl->inform[pos])
    current_node->from_name = permbuff(pl->parameters[pos]->string);
  if (current_sequ->nested <= sequ->nested)
    current_sequ->nested = sequ->nested + 1;
}

void exec_assign(struct in_cmd* cmd)
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

void exec_beam(struct in_cmd* cmd, int flag)
  /* chooses correct beam for beam definitions, upgrades, and resets */
{
  char* name;
  char name_def[] = "default_beam";
  struct command* keep_beam = current_beam;
  struct command_parameter_list* pl = cmd->clone->par;
  struct name_list* nl = cmd->clone->par_names;
  int pos = name_list_pos("sequence", nl);
  int bpos = name_list_pos("sequence", current_beam->par_names);
  if (nl->inform[pos])
  {
    name = pl->parameters[pos]->string;
    if ((current_beam = find_command(name, beam_list)) == NULL)
    {
      set_defaults("beam");
      add_to_command_list(name, current_beam, beam_list, 0);
    }
  }
  else
  {
    name = name_def;
    current_beam = find_command(name, beam_list);
  }
  current_beam->par->parameters[bpos]->string = permbuff(name);
  current_beam->beam_def = 1;
  if (flag == 0) update_beam(cmd->clone);
  else if (flag == 1)  set_defaults("beam");
  current_beam = keep_beam;
}

void exec_call(struct in_cmd* cmd)
  /* handles calling external files */
{
  struct command_parameter_list* pl = cmd->clone->par;
  struct name_list* nl = cmd->clone->par_names;
  int pos = name_list_pos("file", nl);
  int top = in->curr;
  if (nl->inform[pos])
  {
    if (down_unit(pl->parameters[pos]->string))  main_input(top);
  }
  else warning("call without filename:", "ignored");
}

void exec_command()
  /* executes one command */
{
  char** toks;
  char* cmd_name;
  struct in_cmd* p = this_cmd;
  struct in_cmd* pp;
  int ret, izero = 0, pos;
  if (p->cmd_def != NULL)
  {
    while (strcmp(p->cmd_def->name, "exec") == 0)
    {
      if ((pos = name_list_pos(p->tok_list->p[p->decl_start],
                               macro_list->list)) > -1)
      {
        exec_macro(p, pos);
        return;
      }
      pp = p;
      if ((p = buffered_cmd(pp)) == pp) break;
    }
    this_cmd = p;
    toks = p->tok_list->p;
    cmd_name = p->cmd_def->name;
    if (strcmp(cmd_name, "stop") == 0 || strcmp(cmd_name, "quit") == 0
        || strcmp(cmd_name, "exit") == 0)
    {
      madx_finish(); stop_flag = 1; return;
    }
    else if (strcmp(cmd_name, "help") == 0) exec_help(p);
    else if (strcmp(cmd_name, "show") == 0) exec_show(p);
    else if (strcmp(cmd_name, "return") == 0)  return_flag = 1;
    else if (strcmp(cmd_name, "value") == 0)
    {
      print_value(p);
    }
    else if (strcmp(cmd_name, "system") == 0)
      ret = system(noquote(toks[p->decl_start]));
    else if (strcmp(cmd_name, "title") == 0)
      title = permbuff(noquote(toks[p->decl_start]));
    else if (strcmp(cmd_name, "resplot") == 0)
    {
      plot_options = delete_command(plot_options);
      set_defaults("setplot");
    }
    else
    {
      if (get_option("trace")) time_stamp(cmd_name);
      /* clones with defaults for most commands */
      if (strcmp(cmd_name, "option") == 0
          && options != NULL)
      {
        set_option("tell", &izero); /* reset every time */
        p->clone = options; p->clone_flag = 1;
      }
#ifdef _FULL
      else if (strcmp(cmd_name, "setplot") == 0
               && plot_options != NULL)
      {
        p->clone = plot_options; p->clone_flag = 1;
      }
#endif
      else p->clone = clone_command(p->cmd_def);
      scan_in_cmd(p); /* match input command with clone + fill */
      current_command = p->clone;
      if (strcmp(p->cmd_def->module, "control") == 0) control(p);
#ifdef _FULL
      else if (strcmp(p->cmd_def->module, "c6t") == 0) conv_sixtrack(p);
      else if (strcmp(p->cmd_def->module, "edit") == 0) seq_edit_main(p);
      else if (strcmp(p->cmd_def->module, "ibs") == 0)
      {
        current_ibs = p->clone;
        pro_ibs(p);
      }
      else if (strcmp(p->cmd_def->module, "aperture") == 0)
      {
        pro_aperture(p);
      }
      else if (strcmp(p->cmd_def->module, "touschek") == 0)
      {
        current_touschek = p->clone;
        pro_touschek(p);
      }
      else if (strcmp(p->cmd_def->module, "makethin") == 0) makethin(p);
      else if (strcmp(p->cmd_def->module, "match") == 0)
      {
        current_match = p->clone; /* OB 23.1.2002 */
        pro_match(p);
      }
      else if (strcmp(p->cmd_def->module, "correct") == 0)
      {
        pro_correct(p);
      }
      else if (strcmp(p->cmd_def->module, "emit") == 0)
      {
        pro_emit(p);
      }
      else if (strcmp(p->cmd_def->module, "error") == 0)
      {
        current_error = p->clone;
        pro_error(p);
      }
      else if (strcmp(p->cmd_def->module, "ptc_create_universe") == 0)
      {
        w_ptc_create_universe_();
        curr_obs_points = 1;  /* default: always observe at machine end */
      }
      else if (strcmp(p->cmd_def->module, "ptc_create_layout") == 0)
      {
        w_ptc_create_layout_();
      }
      else if (strcmp(p->cmd_def->module, "ptc_move_to_layout") == 0)
      {
        w_ptc_move_to_layout_();
      }
      else if (strcmp(p->cmd_def->module, "ptc_align") == 0)
      {
        w_ptc_align_();
      }
      else if (strcmp(p->cmd_def->module, "ptc_twiss") == 0)
      {
        current_twiss = p->clone;
        pro_ptc_twiss();
      }
      else if (strcmp(p->cmd_def->module, "ptc_normal") == 0)
      {
        w_ptc_normal_();
      }
      else if (strcmp(p->cmd_def->module, "select_ptc_normal") == 0)
      {
        select_ptc_normal(p);
      }
      else if (strcmp(p->cmd_def->module, "ptc_trackline") == 0)
      {
        if (kSkowronDebug) printf("madxp.c: Command is ptc_trackline, calling pro_ptc_trackline\n");
        pro_ptc_trackline(p);
      }
      else if (strcmp(p->cmd_def->module, "ptc_dumpmaps") == 0)
      {
        if (kSkowronDebug) printf("madxp.c: Command is ptc_trackline, calling pro_ptc_dumpmaps\n");
        ptc_dumpmaps(p);
      }
      else if (strcmp(p->cmd_def->module, "ptc_twiss_linac") == 0)
      {
        if (kSkowronDebug) printf("madxp.c: Command is ptc_twiss_linac, calling pro_ptc_trackline\n");
        current_twiss = p->clone;
        pro_ptc_twiss_linac(p);
      }
      else if (strcmp(p->cmd_def->module, "ptc_track") == 0)
      {
        if (kSkowronDebug) printf("madxp.c: Command is ptc_track, calling pro_ptc_track\n");
        pro_ptc_track(p);
      }
      else if (strcmp(p->cmd_def->module, "ptc_setswitch") == 0)
      {
        if (kSkowronDebug) printf("madxp.c: Command is ptc_setswitch, calling pro_ptc_setswitch\n");
        pro_ptc_setswitch(p);
      }
      else if (strcmp(p->cmd_def->module, "ptc_select") == 0)
      {
        if (kSkowronDebug) printf("madxp.c: Command is ptc_select, calling pro_ptc_select\n");
        pro_ptc_select(p);
      }
      else if (strcmp(p->cmd_def->module, "ptc_script") == 0)
      {
        if (kSkowronDebug) printf("madxp.c: Command is ptc_script, calling pro_ptc_script\n");
        pro_ptc_script(p);
      }
      else if (strcmp(p->cmd_def->module, "ptc_observe") == 0)
      {
        ptc_track_observe(p);
      }
      else if (strcmp(p->cmd_def->module, "ptc_start") == 0)
      {
        track_is_on = 1;
        track_start(p->clone);
        p->clone_flag = 1;
        /* w_ptc_start_(); */
      }
      else if (strcmp(p->cmd_def->module, "ptc_track_end") == 0)
      {
        ptc_track_end();
      }
      else if (strcmp(p->cmd_def->module, "ptc_end") == 0)
      {
        if(track_is_on) ptc_track_end();
        w_ptc_end_();
      }
      else if (strcmp(p->cmd_def->module, "sxf") == 0)
      {
        pro_sxf(p);
      }
      else if (strcmp(p->cmd_def->module, "survey") == 0)
      {
        current_survey = p->clone;
        pro_survey(p);
      }
      else if (strcmp(p->cmd_def->module, "track") == 0)
      {
        pro_track(p);
      }
      else if (strcmp(p->cmd_def->module, "twiss") == 0)
      {
        current_twiss = p->clone;
        pro_twiss();
      }
#endif
    }
  }
}

void exec_dumpsequ(struct in_cmd* cmd)
  /* writes a sequence out */
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int level, spos, pos = name_list_pos("sequence", nl);
  struct sequence* sequ = NULL;
  char* name = NULL;
  if (nl->inform[pos] == 0)  sequ = current_sequ;
  else
  {
    name = pl->parameters[pos]->string;
    if ((spos = name_list_pos(name, sequences->list)) >= 0)
      sequ = sequences->sequs[spos];
  }
  pos = name_list_pos("level", nl);
  if (nl->inform[pos] > 0) level = pl->parameters[pos]->double_value;
  else level = 0;
  if (sequ != NULL) dump_exp_sequ(sequ, level);
}

void exec_help(struct in_cmd* cmd)
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

void exec_macro(struct in_cmd* cmd, int pos)
  /* executes a macro */
{
  int i, rs, re, sum = 0, any = 0, level = pro->curr;
  int n = macro_list->macros[pos]->n_formal;
  char** toks = cmd->tok_list->p;
  if (level == pro->max) grow_in_buff_list(pro);
  if (pro->buffers[level] == NULL)
    pro->buffers[level] = new_in_buffer(IN_BUFF_SIZE);
  pro->curr++;
  strcpy(pro->buffers[level]->c_a->c, macro_list->macros[pos]->body->c);
  if (n)
  {
    get_bracket_t_range(toks, '(', ')', cmd->decl_start+1,
                        cmd->tok_list->curr-1, &rs, &re);
    any = re - rs - 1; rs++;
    if (any < 0) any = 0;
    else if (any > n) any = n;
    for (i = 0; i < any; i++)  sum += strlen(toks[rs+i]);
    while (l_wrk->max < (strlen(pro->buffers[level]->c_a->c)+sum))
      grow_char_array(l_wrk);
    for (i = 0; i < any; i++)
    {
      myrepl(macro_list->macros[pos]->formal->p[i], toks[rs+i],
             pro->buffers[level]->c_a->c, l_wrk->c);
      mystrcpy(pro->buffers[level]->c_a, l_wrk->c);
    }
  }
  pro_input(pro->buffers[level]->c_a->c);
  pro->curr--;
}

void exec_option()
{
  if (get_option("reset")) set_defaults("option");
  if (get_option("tell")) print_command(options);

}

void exec_save(struct in_cmd* cmd)
  /* save a sequence with all necessary parameters and sub-sequences */
{
  int i, n = 0, pos, prev = 0, beam_save = log_val("beam", cmd->clone),
    mad8 = log_val("mad8", cmd->clone),
    bare = log_val("bare", cmd->clone), all_sequ = 0;
  char *name, *filename;
  struct element* el;
  struct el_list* ell;
  struct node* c_node;
  struct sequence* sequ;
  struct sequence_list *sql, *sqo;
  struct var_list* varl;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  struct command_parameter* clp;
  default_beam_saved = 0;
  i = name_list_pos("file", nl);
  if (nl->inform[i] == 0)
  {
    warning("save without file:", "ignored");
    return;
  }
  filename = pl->parameters[i]->string;
  if ((out_file = fopen(filename, "w")) == NULL)
  {
    warning("cannot open output file:", filename);
    return;
  }
  warning("SAVE makes all previous USE invalid !", " ");
  pos = name_list_pos("sequence", nl);
  clp = cmd->clone->par->parameters[pos];
  if (nl->inform[pos] == 0)  /* no sequence given, all sequences saved */
  {
    sqo = sequences; all_sequ = 1;
  }
  else
  {
    sqo = new_sequence_list(20);
    for (i = 0; i < clp->m_string->curr; i++)
    {
      name = clp->m_string->p[i];
      if ((pos = name_list_pos(name, sequences->list)) < 0)
        warning("unknown sequence ignored:", name);
      else add_to_sequ_list(sequences->sequs[pos], sqo);
    }
  }
  /* now do it */
  sql = new_sequence_list(20);
  ell = new_el_list(10000);
  if (all_sequ == 0)  varl = new_var_list(2000);
  else varl = clone_var_list(variable_list); /* write all variables */
  for (pos = 0; pos < sqo->curr; pos++)
  {
    sequ = sqo->sequs[pos];
    add_to_sequ_list(sequ, sql);
    /* check for inserted sequences, flatten if necessary  - HG 23.3.04 */
    c_node = sequ->start;
    while(c_node != NULL)
    {
      if (c_node->p_sequ != NULL)
      {
        warning("structured sequence flattened:", sequ->name);
        seq_edit_ex(sequ);
        seq_flatten(edit_sequ);
        seq_end_ex();
        break;
      }
      if (c_node == sequ->end) break;
      c_node = c_node->next;
    }
    /* end mod - HG 23.3.04 */
    if (beam_save && bare == 0)
    {
      if (mad8 == 0) save_beam(sequ, out_file); /* only mad-X */
      else warning("when mad-8 format requested,","beam not saved");
    }
  }
  for (i = sql->curr-1; i >= 0; i--) /* loop over sequences, get elements */
  {
    c_node = sql->sequs[i]->start;
    while (c_node != NULL)
    {
      if ((el = c_node->p_elem) != NULL && strchr(el->name, '$') == NULL
          && strcmp(el->base_type->name, "drift") != 0)
      {
        if (el->def_type != 0) el = el->parent;
        while (el->base_type != el)
        {
          add_to_el_list(&el, 0, ell, 0);
          el = el->parent;
        }
      }
      if (c_node == sql->sequs[i]->end) break;
      c_node = c_node->next;
    }
  }
  if (all_sequ == 0)
  {
    while (prev < ell->curr) /* loop over elements, get variables -
                                recursive, since elements may be added */
    {
      prev = ell->curr;
      for (i = n; i < ell->curr; i++)
        fill_elem_var_list(ell->elem[i], ell, varl);
      n = prev;
    }
    fill_sequ_var_list(sql, ell, varl); /* variables for positions */
  }
  if (mad8)
  {
    if (bare == 0)
    {
      write_vars_8(varl, save_select, out_file);
      write_elems_8(ell, save_select, out_file);
    }
    for (pos = 0; pos < sql->curr; pos++)
    {
      sequ = sql->sequs[pos];
      all_node_pos(sequ);
      sequ->ex_nodes = new_node_list(2*sequ->nodes->curr);
      expand_sequence(sequ, 0);
      export_sequ_8(sequ, save_select, out_file);
      sequ->ex_nodes = delete_node_list(sequ->ex_nodes);
    }
  }
  else
  {
    if (bare == 0)
    {
      write_vars(varl, save_select, out_file);
      write_elems(ell, save_select, out_file);
    }
    write_sequs(sql, save_select, out_file);
  }
  fclose(out_file);
  if (sqo != sequences) sqo = delete_sequence_list(sqo);
  sql = delete_sequence_list(sql);
  ell = delete_el_list(ell);
  varl = delete_var_list(varl);
  current_sequ = NULL;
}

void exec_show(struct in_cmd* cmd)
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
      else fprintf(prt_file, "%s not found\n", toks[i]);
    }
  }

}

struct node* expand_node(struct node* node, struct sequence* top_sequ,
                         struct sequence* sequ, double position)
  /* replaces a (sequence) node by a sequence of nodes - recursive */
{
  struct sequence* nodesequ = node->p_sequ;
  struct node *p, *q = nodesequ->start;
  int i;
  p = clone_node(q, 0);
  if ((i = name_list_pos(p->p_elem->name, occ_list)) < 0)
    i = add_to_name_list(p->p_elem->name, 1, occ_list);
  else strcpy(p->name, compound(p->p_elem->name, ++occ_list->inform[i]));
  add_to_node_list(p, 0, top_sequ->ex_nodes);
  p->previous = node->previous;
  p->previous->next = p;
  p->position = position + get_node_pos(p, nodesequ)
    - get_refpos(nodesequ);
  while (p != NULL)
  {
    if (q == nodesequ->end) break;
    p->next = clone_node(q->next, 0);
    p->next->previous = p;
    p = p->next;
    q = q->next;
    p->position = position + get_node_pos(p, nodesequ)
      - get_refpos(nodesequ);
    if (p->p_sequ != NULL) p = expand_node(p, top_sequ,
                                           p->p_sequ, p->position);
    else
    {
      p->share = top_sequ->share;
      add_to_node_list(p, 0, top_sequ->ex_nodes);
    }
  }
  delete_node(node);
  return p;
}

void expand_sequence(struct sequence* sequ, int flag)
  /* expands a sequence into nodes, expands sequence nodes */
{
  /* Transfers errors from original nodes if flag != 0;
     this is needed for SXF input  */
  struct node *p, *q = sequ->start;
  p = sequ->ex_start = clone_node(sequ->start, 0);
  add_to_node_list(p, 0, sequ->ex_nodes);
  while (p != NULL)
  {
    if (q == sequ->end) break;
    p->next = clone_node(q->next, flag);
    p->next->previous = p;
    p = p->next;
    q = q->next;
    if (p->p_sequ != NULL) p = expand_node(p, sequ, sequ, p->position);
    else add_to_node_list(p, 0, sequ->ex_nodes);
  }
  sequ->ex_end = p;
  sequ->ex_end->next = sequ->ex_start;
  sequ->ex_start->previous = sequ->ex_end;
  p = sequ->ex_start;
  while (p)
  {
    if (strstr(p->base_name, "kicker") ||strstr(p->base_name, "monitor"))
      p->enable = 1; /* flag for orbit correction module */
    if (p == sequ->ex_end) break;
    p = p->next;
  }
}

double expression_value(struct expression* expr, int flag) /* recursive */
  /* returns the value of an expression if valid, else zero */
{
  double val = zero;
  if (expr->status == 0 || flag == 2)
  {
    if (expr->polish != NULL)
    {
      val = expr->value = polish_value(expr->polish);
      expr->status = 1;
    }
  }
  else val = expr->value;
  return val;
}

void fatal_error(char* t1, char* t2)
  /*prints fatal error message, halts program */
{
  printf("+=+=+= fatal: %s %s\n",t1,t2); exit(1);
}

struct command* find_command(char* name, struct command_list* cl)
{
  int pos;
  if ((pos = name_list_pos(name, cl->list)) < 0)
    return NULL;
  return cl->commands[pos];
}

void fill_elem_var_list(struct element* el, struct el_list* ell,
                        struct var_list* varl)
  /* puts all variables an element depends on, in a list */
{
  struct command* cmd = el->def;
  int i;
  for (i = 0; i < cmd->par->curr; i++)
    fill_par_var_list(ell, cmd->par->parameters[i], varl);
}

void fill_expr_list(char** toks, int s_start, int s_end,
                    struct expr_list* p)
{
  int start = s_start, nitem = s_end + 1, end, cnt = 0, nc;
  while (start < nitem)
  {
    if ((nc = next_char(',', toks, start, nitem)) < 0) nc = nitem;
    if (loc_expr(toks, nc, start, &end))
    {
      if (cnt == p->max)  grow_expr_list(p);
      p->list[cnt++] = make_expression(end + 1 - start, &toks[start]);
      start = nc + 1;
    }
    else break;
  }
  p->curr = cnt;
}

void fill_expr_var_list(struct el_list* ell,
                        struct expression* expr, struct var_list* varl)
  /* puts all variables an expression depends on, in a list */
{
  struct variable* var;
  struct element* el;
  char name[2*NAME_L];
  char* p;
  int i, k, kc;
  struct int_array* deco = expr->polish;
  for (i = 0; i < deco->curr; i++)   /* decoding loop */
  {
    if ((k = deco->i[i]) > 4 && (kc = k / 100000000) == 1)
    {
      k -= 100000000 * kc;  strcpy(name, expr_chunks->names[k]);
      if ((p = strstr(name, "->")) != NULL)
      {
        *p = '\0';
        if ((el = find_element(name, element_list)) != NULL)
          add_to_el_list(&el, 0, ell, 0);
      }
      else if ((var = find_variable(name, variable_list)) != NULL)
      {
        add_to_var_list(var, varl, 2);
        if (var->type == 2 && var->expr != NULL)
          fill_expr_var_list(ell, var->expr, varl);
      }
    }
  }
}

void fill_par_var_list(struct el_list* ell,
                       struct command_parameter* par, struct var_list* varl)
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

void fill_sequ_var_list(struct sequence_list* sql, struct el_list* ell,
                        struct var_list* varl)
  /* puts all variables a sequence depends on, in a list */
{
  int i;
  struct sequence* sequ;
  struct node* c_node;
  for (i = 0; i < sql->curr; i++)
  {
    sequ = sql->sequs[i];
    if (sequ->l_expr != NULL) fill_expr_var_list(ell, sequ->l_expr, varl);
    c_node = sequ->start;
    while(c_node != NULL)
    {
      if (c_node->at_expr != NULL)
        fill_expr_var_list(ell, c_node->at_expr, varl);
      if (c_node == sequ->end)  break;
      c_node = c_node->next;
    }
  }
}

struct element* find_element(char* name, struct el_list* ell)
{
  int pos;
  if ((pos = name_list_pos(name, ell->list)) < 0)
    return NULL;
  return ell->elem[pos];
}

struct variable* find_variable(char* name, struct var_list* varl)
{
  int pos;
  if ((pos = name_list_pos(name, varl->list)) < 0)
    return NULL;
  return varl->vars[pos];
}

double get_aperture(struct node* node, char* par)
  /* returns aperture parameter 'i' where i is integer at the end of par;
     e.g. aper_1 gives i = 1 etc. (count starts at 1) */
{
  int i, k, n = strlen(par);
  double val = zero, vec[100];
  for (i = 0; i < n; i++)  if(isdigit(par[i])) break;
  if (i == n) return val;
  sscanf(&par[i], "%d", &k); k--;
  if ((n = element_vector(node->p_elem, "aperture", vec)) > k)  val = vec[k];
  return val;
}

void get_bracket_range(char* string, char lb, char rb, int* rs, int* re)
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

void get_bracket_t_range(char* toks[], char lb, char rb,
                         int start, int end, int* rs, int* re)
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

void get_defined_commands()
  /* reads + stores the commands defined in madxdict.h */
{
  char rout_name[] = "get_defined_commands";
  int i;
  char** p;
  int n = char_cnt(';', command_def);
  p = (char**) mymalloc(rout_name,n * sizeof(char*));
  p[0] = strtok(command_def, ";");
  for (i = 1; i < n; i++) /* make temporary list - strtok is called again */
    p[i] = strtok(NULL, ";");
  for (i = 0; i < n; i++) store_command_def(p[i]);
  myfree(rout_name, p);
}

void get_defined_constants()
{
  /* reads + stores the constants defined in madxdict.h */
  supp_char('\n', predef_constants);
  pro_input(predef_constants);
  start_var = variable_list->curr;
}

struct element* get_drift(double length)
  /* makes a drift space with the required length */
{
  struct element *p, *bt;
  struct command* clone;
  char key[NAME_L];
  int i;
  for (i = 0; i < drift_list->curr; i++)
  {
    p = drift_list->elem[i];
    if (fabs(p->length - length) < ten_m_12) return p;
  }
  sprintf(key, "drift_%d", drift_list->curr);
  bt = find_element("drift", base_type_list);
  clone = clone_command(bt->def);
  store_comm_par_value("l", length, clone);
  p = make_element(key, "drift", clone, 0);
  add_to_el_list(&p, 1, drift_list, 0);
  return p;
}

char* get_new_name()
  /* makes a new internal element or variable name */
{
  char name[NAME_L] = "__";
  sprintf(&name[2], "%d", new_name_count++);
  strcat(name, "__");
  return permbuff(name);
}

double get_node_pos(struct node* node, struct sequence* sequ) /*recursive */
  /* returns node position from declaration for expansion */
{
  double fact = 0.5 * sequ->ref_flag; /* element half-length offset */
  double pos, from = 0;
  if (loop_cnt++ == MAX_LOOP)
  {
    sprintf(c_dum->c, "%s   occurrence: %d", node->p_elem->name,
            node->occ_cnt);
    fatal_error("circular call in position of", c_dum->c);
  }
  if (node->at_expr == NULL) pos = node->at_value;
  else                       pos = expression_value(node->at_expr, 2);
  pos += fact * node->length; /* make centre position */
  if (node->from_name != NULL)
  {
    if ((from = hidden_node_pos(node->from_name, sequ)) == INVALID)
      fatal_error("'from' reference to unknown element:", node->from_name);
  }
  loop_cnt--;
  pos += from;
  return pos;
}

int get_option(char* str)
{
/* This function is called by fortran to get option of a command */
  int i, k;
  mycpy(c_dum->c, str);
  if (options != NULL
      && (i = name_list_pos(c_dum->c, options->par_names)) > -1)
    return (k = options->par->parameters[i]->double_value);
  else if (strcmp(c_dum->c, "warn") == 0) return init_warn;
  else return 0;
}

double get_refpos(struct sequence* sequ)
  /* returns the position of a refpos element, or zero */
{
  int i;
  if (sequ != NULL && sequ->refpos != NULL)
  {
    sprintf(c_dum->c, "%s:1", sequ->refpos);
    if ((i = name_list_pos(c_dum->c, sequ->nodes->list)) < 0)
      fatal_error("'refpos' reference to unknown element:", sequ->refpos);
    return get_node_pos(sequ->nodes->nodes[i], sequ);
  }
  else return zero;
}

int get_stmt(FILE* file, int supp_flag)
  /* Returns complete command(s) in in->buffers[inbuf_level] */
  /* return value:  0   no more input */
  /*                >0  OK */

{
  struct char_array* ca;
  char *c_st, *c_cc, *c_ex;
  char end = ';';
  int in_comment = 0, spec_test = 0, curl_level = 0, go_cond;
  if (in->buffers[in->curr] == NULL)
    in->buffers[in->curr] = new_in_buffer(IN_BUFF_SIZE);
  ca = in->buffers[in->curr]->c_a;
  ca->curr = 0;
  do /* read lines until complete statement(s) */
  {
    next:
    if (ca->max - ca->curr < MAX_LINE) grow_char_array(ca);
    if (!fgets(&ca->c[ca->curr], MAX_LINE, file)) return 0;
    if (get_option("echo")) puts(&ca->c[ca->curr]);
    c_cc = mystrstr(&ca->c[ca->curr], "//");
    c_ex = mystrchr(&ca->c[ca->curr], '!');
    if (c_cc != NULL && c_ex != NULL)
    {
      c_cc = (int)c_cc < (int)c_ex ? c_cc : c_ex; *c_cc = '\0';
    }
    else if(c_cc != NULL)
    {
      if (c_cc == &ca->c[ca->curr]) goto next;
      else *c_cc = '\0';
    }
    else if(c_ex != NULL)
    {
      if (c_ex == &ca->c[ca->curr]) goto next;
      else *c_ex = '\0';
    }
    if (in_comment)
    {
      if (mystrstr(&ca->c[ca->curr], "*/") == NULL)  goto next;
      else
      {
        remove_upto(&ca->c[ca->curr], "*/"); in_comment = 0;
      }
    }
    else if((c_st = mystrstr(&ca->c[ca->curr], "/*")) != NULL)
    {
      if (mystrstr(&ca->c[ca->curr], "*/") == NULL)
      {
        *c_st = '\0'; ca->curr += supp_lt(&ca->c[ca->curr], supp_flag);
        in_comment = 1; goto next;
      }
      else remove_range(&ca->c[ca->curr], "/*", "*/");
    }
    if (spec_test == 0 && mystrchr(&ca->c[ca->curr], '{') != NULL)
    {
      spec_test = 1; if (in_spec_list(ca->c)) end = '}';
    }
    curl_level += char_cnt('{', &ca->c[ca->curr])
      - char_cnt('}', &ca->c[ca->curr]);
    if ((ca->curr += supp_lt(&ca->c[ca->curr], supp_flag)) == 0) go_cond = 1;
    else if ((go_cond = curl_level) == 0)
      go_cond = ca->c[ca->curr-1] == ';' || ca->c[ca->curr-1] == end ? 0 : 1;
  }
  while (go_cond);
  return 1;
}

int get_string(char* name, char* par, char* string)
  /* returns string for  value "par" of command or store "name" if present,
     length = string length, else length = 0 if not present */
{
  struct name_list* nl = NULL;
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
    if (current_survey != NULL) nl = current_survey->par_names;
    if (nl != NULL && nl->inform[name_list_pos(c_dum->c, nl)])
    {
      if ((p = command_par_string(c_dum->c, current_survey)) != NULL)
      {
        strcpy(string, p); length = strlen(p);
      }
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
    if (current_twiss != NULL) nl = current_twiss->par_names;
    if (nl != NULL && nl->inform[name_list_pos(c_dum->c, nl)])
    {
      if ((p = command_par_string(c_dum->c, current_twiss)) != NULL)
      {
        strcpy(string, p); length = strlen(p);
      }
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
    if ((cmd = find_command(c_dum->c, stored_commands)) != NULL)
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

void get_sxf_names()
  /* reads and stores names for SXF I/O from madxl.h */
{
  int i = 0;
  while (sxf_table_names[i][0] != ' ')
  {
    add_to_name_list(sxf_table_names[i++], 0, sxf_list);
  }
}

int get_val_num(char* in_string, int start, int end)
{
  int j, dot = 0, exp = 0, sign = 0;
  char c;
  for (j = start; j < end; j++)
  {
    c = in_string[j];
    if(!isdigit(c))
    {
      if ((c = in_string[j]) == '.')
      {
        if (dot || exp) return (j - 1);
        dot = 1;
      }
      else if (c == 'e' || c == 'd')
      {
        if (c == 'd') in_string[j] = 'e';
        if (exp) return (j - 1);
        else exp = j+1;
      }
      else if(strchr("+-", c))
      {
        if (exp != j || sign) return (j - 1);
        sign = 1;
      }
      else return (j - 1);
    }
  }
  return (j - 1);
}

double get_value(char* name, char* par)
  /* this function is used by fortran to get the parameters values
     returns parameter value "par" for command or store "name" if present,
     else INVALID */
{
/* WHY IT IS NOT A SWITCH?????????*/
  struct name_list* nl = NULL;
  mycpy(c_dum->c, name);
  mycpy(aux_buff->c, par);
  if (strcmp(c_dum->c, "beam") == 0)
    return command_par_value(aux_buff->c, current_beam);
  else if (strcmp(c_dum->c, "probe") == 0)
    return command_par_value(aux_buff->c, probe_beam);
  else if (strcmp(c_dum->c, "survey") == 0)
  {
    if (current_survey != NULL) nl = current_survey->par_names;
    if (nl != NULL && nl->inform[name_list_pos(aux_buff->c, nl)])
      return command_par_value(aux_buff->c, current_survey);
    else return zero;
  }
  else if (strcmp(c_dum->c, "twiss") == 0)
  {
    if (current_twiss != NULL) nl = current_twiss->par_names;
    if (nl != NULL && nl->inform[name_list_pos(aux_buff->c, nl)])
      return command_par_value(aux_buff->c, current_twiss);
    else return zero;
  }
  else if (strcmp(c_dum->c, "sequence") == 0)
  {
    if (strcmp(aux_buff->c, "l") == 0) return current_sequ->length;
    else if (strcmp(aux_buff->c, "range_start") == 0)
      return (current_sequ->range_start->position
              - 0.5 * current_sequ->range_start->length);
    else return INVALID;
  }
  else if (current_command != NULL
           && strcmp(c_dum->c, current_command->name) == 0)
    return command_par_value(aux_buff->c, current_command);
  else return INVALID;
}

double get_variable(char* name)
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
    if ((el = find_element(comm, element_list)) != NULL)
      val = command_par_value(par, el->def);
    else if ((cmd = find_command(comm, stored_commands)) != NULL)
      val = command_par_value(par, cmd);
    else if ((cmd = find_command(comm, beta0_list)) != NULL)
      val = command_par_value(par, cmd);
    else if ((cmd = find_command(comm, defined_commands)) != NULL)
      val = command_par_value(par, cmd);
  }
  return val;
}

double hidden_node_pos(char* name, struct sequence* sequ) /*recursive */
  /* (for 'from' calculation:) returns the position of a node
     in the current or a sub-sequence (hence hidden) */
{
  double pos;
  struct node* c_node;
  char tmp[NAME_L];
  int i;
  strcpy(tmp, name);
  square_to_colon(tmp);
  if ((i = name_list_pos(tmp, sequ->nodes->list)) > -1)
    return get_node_pos(sequ->nodes->nodes[i], sequ);
  else
  {
    c_node = sequ->start;
    while (c_node != NULL)
    {
      if (c_node->p_sequ != NULL)
      {
        if ((pos = hidden_node_pos(name, c_node->p_sequ)) != INVALID)
        {
          pos += get_node_pos(c_node, sequ) - get_refpos(c_node->p_sequ);
          return pos;
        }
      }
      if (c_node == sequ->end) break;
      c_node = c_node->next;
    }
    return INVALID;
  }
}

int in_spec_list(char* string)
  /* checks for presence of special commands IF() etc. */
{
  char* cp;
  int i = 0, n = mymin((int)strlen(string), 100);
  strncpy(c_dum->c, string, n); c_dum->c[n] = '\0'; stolower(c_dum->c);
  supp_char(' ', c_dum->c);
  while (special_comm_cnt[i])
  {
    if (special_comm_desc[i][0] == '>')
    {
      if ((cp = strchr(c_dum->c, special_comm_desc[i][1])) != NULL)
      {
        if (strncmp(++cp, &special_comm_desc[i][2], special_comm_cnt[i])
            == 0)  return i+1;
      }
    }
    else if (strncmp(c_dum->c, &special_comm_desc[i][0],special_comm_cnt[i])
             == 0)  return i+1;
    i++;
  }
  return 0;
}

void link_in_front(struct node* new, struct node* el)
  /* links a node 'new' in front of a node 'el' */
{
  el->previous->next = new;
  new->previous = el->previous; new->next = el;
  el->previous = new;
}

int loc_expr(char** items, int nit, int start, int* end)
  /* Returns the type and end of an expression, or 0 if illegal */
{
  char c;
  int i, e_type = 1, par_level = 0, ltog = -1;
  *end = start - 1;
  if (nit > start && is_expr_start(*items[start]))
  {
    c = *items[start];
    for (i = start; i < nit; i++)
    {
      c = *items[i];
      if (c == '(')  {par_level++; e_type = 2;}
      else if (c == ')')
      {
        if (par_level == 0) return 0;
        par_level--; ltog = 0;
      }
      else if (par_level == 0)
      {
        if (ltog < 0)  ltog = is_operator(c) ? 1 : 0;
        else if ((ltog == 0 && is_operator(c))
                 || (ltog != 0 && is_operand(c)))  ltog = 1 - ltog;
        else return 0;
      }
      *end = i;
      if ((*end > start && ltog > 0) || (isalpha(c) || c == '_'))  e_type = 2;
    }
    return e_type;
  }
  else return 0;
}

int logic_expr(int nit, char* toks[])
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

int log_val(char* name, struct command* cmd)
  /* returns 0 = flase, 1 = true for a logical command parameter */
{
  struct name_list* nl = cmd->par_names;
  struct command_parameter_list* pl = cmd->par;
  int pos = name_list_pos(name, nl);
  if (pos > -1 && nl->inform[pos]) /* "name" has beem read */
    return pl->parameters[pos]->double_value == zero ? 0 : 1;
  else return 0;
}

void madx_finish()
  /* write the termination message */
{
  int nwarnings = warn_numb+warn_numbf;
  /* should work with Lahey on windows 24.03.2004 */

  if (final_message == 0)
  {
    final_message = 1;
#ifdef _FULL
    if (plots_made)
    {
      gxterm_();
    }
#endif
    printf("\n  Number of warnings: %d\n",nwarnings);
    if (nwarnings > 0)
     {
        printf("%d in C and %d in Fortran\n",warn_numb,warn_numbf);
     }
    if (get_option("trace")) time_stamp("end");
    printf("\n  +++++++++++++++++++++++++++++++++++\n");
    printf("  + %s finished normally +\n", myversion);
    printf("  +++++++++++++++++++++++++++++++++++\n");
  }
}
void madx_init()
  /* initializes program */
{
  struct variable* var;
  int ione = 1;
#ifdef _WIN32
  interactive = 1;
#endif
#ifndef _WIN32
  interactive = intrac();
#endif
  init55(123456789);          /* random generator */
  if (watch_flag == 1)  debug_file = fopen("madx.debug", "w");
  else if (watch_flag == 2)  debug_file = stdout;
  if (stamp_flag == 1)  stamp_file = fopen("madx.stamp", "w");
  else if (stamp_flag == 2)  stamp_file = stdout;
  in = new_in_buff_list(100); /* list of input buffers, dynamic */
  in->input_files[0] = stdin;
  prt_file = stdout;
  pro = new_in_buff_list(100); /* list of process buffers, dynamic */
  pro->buffers[0] = new_in_buffer(IN_BUFF_SIZE);
  pro->curr = 1;
  c_dum = new_char_array(AUX_LG);
  c_join = new_char_array(AUX_LG);
  work = new_char_array(AUX_LG);
  l_wrk = new_char_array(AUX_LG);
  char_buff = new_char_array_list(100); /* list of character arrays, dynamic */
  char_buff->ca[char_buff->curr++] = new_char_array(CHAR_BUFF_SIZE);
  aux_buff = new_char_array(AUX_LG);  /* dynamic temporary buffer */
  drift_list = new_el_list(1000); /* dynamic list for internal drifts */
  variable_list = new_var_list(2000); /* dynamic list of variables */
  comm_constraints = new_constraint_list(10); /* dynamic constraint list */
  beam_list = new_command_list("beam_list", 10); /* dynamic beam list */
  stored_track_start = new_command_list("track_start", 100); /* dynamic */
  table_deselect = new_command_list_list(10); /* dynamic */
  table_select = new_command_list_list(10); /* dynamic */
  defined_commands = new_command_list("defined_commands", 100); /* dynamic */
  stored_commands = new_command_list("stored_commands", 500); /* dynamic */
  line_list = new_macro_list(100); /* dynamic */
  macro_list = new_macro_list(100); /* dynamic */
  base_type_list = new_el_list(60); /* dynamic */
  element_list = new_el_list(20000); /* dynamic */
  buffered_cmds = new_in_cmd_list(10000); /* dynamic */
  sequences = new_sequence_list(20); /* dynamic */
  match_sequs = new_sequence_list(2);
  selected_ranges = new_node_list(10000); /* dynamic */
  selected_elements = new_el_list(10000); /* dynamic */
  tmp_p_array = new_char_p_array(1000); /* dynamic */
  tmp_l_array = new_char_p_array(1000); /* dynamic */
  sxf_list = new_name_list("sxf_list", 50); /* dynamic */
  deco_init();
  get_defined_constants();
  get_defined_commands();
  get_sxf_names();
  pi = get_variable("pi");
  twopi = two * pi;
  degrad = 180 / pi;
  raddeg = pi / 180;
  e = get_variable("e");
  clight = get_variable("clight");
  hbar = get_variable("hbar");
  var = new_variable("twiss_tol", 1.e-6, 1, 1, NULL, NULL);
  add_to_var_list(var, variable_list, 1);
  title = permbuff("no-title");
  set_defaults("option");
  set_defaults("beam");
  add_to_command_list("default_beam", current_beam, beam_list, 0);
  set_defaults("set");
  set_defaults("setplot");
  set_defaults("threader");
  table_register = new_table_list(10); /* dynamic */
  beta0_list = new_command_list("beta0_list", 10); /* dynamic */
  savebeta_list = new_command_list("savebeta_list", 10); /* dynamic */
  seqedit_select = /* dynamic - for "select seqedit" commands */
    new_command_list("seqedit_select", 10);
  error_select = /* dynamic - for "select error" commands */
    new_command_list("error-select", 10);
  save_select = /* dynamic - for "select save" commands */
    new_command_list("save_select", 10);
  slice_select = /* dynamic - for "select makethin" commands */
    new_command_list("slice_select", 10);
  sector_select = /* dynamic - for "select sectormap" commands */
    new_command_list("sector_select", 10);
  s_range = new_int_array(10); /* dynamic */
  e_range = new_int_array(10); /* dynamic */
  sd_range = new_int_array(10); /* dynamic */
  ed_range = new_int_array(10); /* dynamic */
  zero_double(orbit0, 6);
  zero_double(disp0, 6);
  zero_double(guess_orbit,6);
  zero_double(oneturnmat, 36);
  set_option("twiss_print", &ione);
}

void madx_start()
  /* prints start message after having read madxdict.h */
{
  struct tm* tm;

  /*  setbuf(stdout,(char *)0); */ /* no buffering - for debugging */
  time(&start_time); /* initialize timing */
  tm = localtime(&start_time); /* split system time */
  last_time = start_time;
  printf("\n  +++++++++++++++++++++++++++++++++++++++++++\n");
  printf("  +              %s              +\n", myversion);
  printf("  + %s      +\n",code_mod_date);
  printf("  + Execution Time Stamp: %02d.%02d.%02d %02d.%02d.%02d +\n",
         tm->tm_mday, tm->tm_mon+1, tm->tm_year%100,
         tm->tm_hour, tm->tm_min, tm->tm_sec);
  printf("  +++++++++++++++++++++++++++++++++++++++++++\n");
}

void main_input(int top)
  /* loops over input until end of execution */
{
  while (in_stop == 0)
  {
    if (interactive && in->curr == 0) puts("X: ==>");
    if (return_flag || get_stmt(in->input_files[in->curr], 0) == 0)
    {
      if (in->curr == 0) return;
      fclose(in->input_files[in->curr--]);
      return_flag = 0;
      if (in->curr == top) return;
    }
    else
    {
      stolower_nq(in->buffers[in->curr]->c_a->c);
      pro_input(in->buffers[in->curr]->c_a->c);
      if (stop_flag)  return;
    }
  }
}

struct element* make_element(char* name, char* parent,
                             struct command* def, int flag)
  /* makes a new element from declaration, stores in list */
{
/*  double length; */
  struct element* el = new_element(name);
  el->def = def;
  if (strcmp(name, parent) == 0)  /* basic element type like drift etc. */
  {
    add_to_el_list(&el, def->mad8_type, base_type_list, 1);
    el->parent = el->base_type = el;
  }
  else
  {
    if((el->parent = find_element(parent, element_list)) == NULL)
      fatal_error("unknown class type:", parent);
    el->base_type = el->parent->base_type;
    el->length = el_par_value("l", el);
  }
  add_to_el_list(&el, def->mad8_type, element_list, flag);
  return el;
}

void make_elem_node(struct element* el, int occ_cnt)
  /* makes + links a new node at the end of the current sequence */
{
  prev_node = current_node;
  current_node = new_elem_node(el, occ_cnt);
  current_node->occ_cnt = occ_cnt;
  add_to_node_list(current_node, 0, current_sequ->nodes);
  if (prev_node != NULL) prev_node->next = current_node;
  current_node->previous = prev_node;
  current_node->next = NULL;
}

struct expression* make_expression(int n, char** toks)
  /* makes an expression from a list of tokens */
{
  struct expression* expr = NULL;

  if (polish_expr(n, toks) == 0)
    expr = new_expression(join(toks, n), deco);
  else warning("Invalid expression starting at:", join_b(toks, n));
  return expr;
}

int make_line(char* statement)
  /* makes a new line from input command, stores name in name list */
{
  struct macro* m;
  char** toks = tmp_l_array->p;
  char *prs, *psem;
  int i, n, rs, re;
  while(strlen(statement) >= aux_buff->max) grow_char_array(aux_buff);
  strcpy(aux_buff->c, statement);
  if ((prs = strchr(aux_buff->c, '=')) == NULL) return -3;
  *prs = '\0'; prs++;
  pre_split(aux_buff->c, l_wrk, 0);
  mysplit(l_wrk->c, tmp_l_array);
  get_bracket_t_range(toks, '(', ')', 0, tmp_l_array->curr-1, &rs, &re);
  if ((n = re - rs - 1) < 0) n = 0; /* number of formal arguments if any */
  m = new_macro(n, 2*strlen(prs), 50);
  strcpy(m->name, toks[0]); rs++;
  for (i = 0; i < n; i++) m->formal->p[i] = permbuff(toks[rs+i]);
  if (n > 0) m->formal->curr = n;
  if ((psem = strchr(prs, ';')) != NULL) *psem = '\0';
  mystrcpy(l_wrk, prs);
  pre_split(l_wrk->c, m->body, 0);
  m->body->curr = strlen(m->body->c);
  mysplit(m->body->c, m->tokens);
  n = 0;
  for (i = 0; i < m->tokens->curr; i++)
  {
    if (*m->tokens->p[i] == '(')  n++;
    else if (*m->tokens->p[i] == ')')  n--;
    else if(i > 0 && *m->tokens->p[i] == '*'
            && isdigit(*m->tokens->p[i-1]) == 0) return -3;
  }
  if (n) return -3; /* bracket closing error */
  add_to_macro_list(m, line_list);
  return 0;
}

int make_macro(char* statement)
  /* makes a new macro from input command, stores name in name list */
{
  struct macro* m;
  char** toks = tmp_l_array->p;
  int i, n, rs, re, start_2;
  while(strlen(statement) >= aux_buff->max) grow_char_array(aux_buff);
  strcpy(aux_buff->c, statement);
  get_bracket_range(aux_buff->c, '{', '}', &rs, &re);
  start_2 = rs + 1;
  aux_buff->c[rs] = '\0'; aux_buff->c[re] = '\0'; /* drop '{' and '}' */
  pre_split(aux_buff->c, l_wrk, 0);
  mysplit(l_wrk->c, tmp_l_array);
  get_bracket_t_range(toks, '(', ')', 0, tmp_l_array->curr-1, &rs, &re);
  if ((n = re - rs - 1) < 0) n = 0; /* number of formal arguments if any */
  m = new_macro(n, strlen(&aux_buff->c[start_2]), 0);
  strcpy(m->name, toks[0]); rs++;
  for (i = 0; i < n; i++) m->formal->p[i] = permbuff(toks[rs+i]);
  if (n > 0) m->formal->curr = n;
  strcpy(m->body->c, &aux_buff->c[start_2]); m->body->curr = strlen(m->body->c);
  add_to_macro_list(m, macro_list);
  return 0;
}

void make_occ_list(struct sequence* sequ)
  /* makes the node occurrence list */
{
  struct node* c_node = sequ->start;
  int i;
  while (c_node != NULL)
  {
    if (c_node->p_elem != NULL)
    {
      if ((i = name_list_pos(c_node->p_elem->name, occ_list)) < 0)
        i = add_to_name_list(c_node->p_elem->name, 1, occ_list);
      else ++occ_list->inform[i];
    }
    if (c_node == sequ->end) break;
    c_node = c_node->next;
  }
}

void make_sequ_node(struct sequence* sequ, int occ_cnt)
  /* makes + links a node pointing to a sub-sequence */
{
  prev_node = current_node;
  current_node = new_sequ_node(sequ, occ_cnt);
  current_node->occ_cnt = occ_cnt;
  add_to_node_list(current_node, 0, current_sequ->nodes);
  prev_node->next = current_node;
  current_node->previous = prev_node;
  current_node->next = NULL;
}

char* make_string_variable(char* string)
  /* creates + stores a variable containing a character string */
{
  char* name = get_new_name();
  struct variable* var = new_variable(name, zero, 3, 0, NULL, string);
  add_to_var_list(var, variable_list, 0);
  return var->name;
}

int mysplit(char* buf, struct char_p_array* list)
{
  /* splits into tokens */
  int j;
  char* p;
  if ((p =strtok(buf, " \n")) == NULL) return 0;
  list->curr = 0;
  list->p[list->curr++] = p;
  while ((p = strtok(NULL, " \n")) != NULL)
  {
    if (list->curr == list->max) grow_char_p_array(list);
    list->p[list->curr++] = p;
  }
  /* remove '@' in strings */
  for (j = 0; j < list->curr; j++)
    if(*list->p[j] == '\"' || *list->p[j] == '\'') /* quote */
      replace(list->p[j], '@', ' ');
  return list->curr;
}

int next_char(char c, char** toks, int start, int nitem)
  /* returns the number of the token starting with c after token start */
{
  int i;
  for (i = start; i < nitem; i++) if(*toks[i] == c)  return i;
  return -1;
}

double node_value(char* par)
  /* returns value for parameter par of current element */
{
  double value;
  char lpar[NAME_L];
  mycpy(lpar, par);
  if (strcmp(lpar, "l") == 0) value = current_node->length;
  else if (strcmp(lpar, "dipole_bv") == 0) value = current_node->dipole_bv;
  else if (strcmp(lpar, "other_bv") == 0) value = current_node->other_bv;
  else if (strcmp(lpar, "chkick") == 0) value = current_node->chkick;
  else if (strcmp(lpar, "cvkick") == 0) value = current_node->cvkick;
  else if (strcmp(lpar, "obs_point") == 0) value = current_node->obs_point;
  else if (strcmp(lpar, "sel_sector") == 0) value = current_node->sel_sector;
  else if (strcmp(lpar, "enable") == 0) value = current_node->enable;
  else if (strcmp(lpar, "occ_cnt") == 0) value = current_node->occ_cnt;
  else value =  element_value(current_node, lpar);
  return value;
}

int pass_select(char* name, struct command* sc)
  /* checks name against class (if element) and pattern that may
     (but need not) be contained in command sc;
     0: does not pass, 1: passes */
{
  struct name_list* nl = sc->par_names;
  struct command_parameter_list* pl = sc->par;
  struct element* el = find_element(strip(name), element_list);
  int pos, in = 0, any = 0;
  char *class, *pattern;
  pos = name_list_pos("class", nl);
  if (pos > -1 && nl->inform[pos])  /* parameter has been read */
  {
    el = find_element(strip(name), element_list);
    if (el != NULL)
    {
      class = pl->parameters[pos]->string;
      in = belongs_to_class(el, class);
      if (in == 0) return 0;
    }
  }
  any = in = 0;
  pos = name_list_pos("pattern", nl);
  if (pos > -1 && nl->inform[pos])  /* parameter has been read */
  {
    any = 1;
    pattern = stolower(pl->parameters[pos]->string);
    if(myregex(pattern, strip(name)) == 0)  in = 1;
  }
  if (any == 0) return 1;
  else return in;
}

int pass_select_list(char* name, struct command_list* cl)
  /* returns 0 (does not pass) or 1 (passes) for a list of selects */
{
  int i, ret = 0;
  if (cl->curr == 0)  return 1;
  for (i = 0; i < cl->curr; i++)
  {
    if ((ret = pass_select(name, cl->commands[i]))) break;
  }
  return ret;
}

char* permbuff(char* string)  /* string -> general buffer, returns address */
{
  int n, k = char_buff->ca[char_buff->curr-1]->curr;
  if (string == NULL)  return NULL;
  n = strlen(string)+1;
  if (k + n >= char_buff->ca[char_buff->curr-1]->max)
  {
    if (char_buff->curr == char_buff->max) grow_char_array_list(char_buff);
    char_buff->ca[char_buff->curr++] = new_char_array(CHAR_BUFF_SIZE);
    k = 0;
  }
  strcpy(&char_buff->ca[char_buff->curr-1]->c[k], string);
  char_buff->ca[char_buff->curr-1]->curr += n;
  return &char_buff->ca[char_buff->curr-1]->c[k];
}

int polish_expr(int c_item, char** item)   /* split input */
  /* deco output array containing:
     expression in Polish notation of length deco->curr,
     coded as 0-, 1+, 2*, 3/, 4^ (power),
     6 evaluate function
     100000000 + n = variable n (refers to vars),
     200000000 + n = function n (refers to functs),
     400000000 + n = real n (refers to doubles)
     -- Example: suppose a, b are variables 0 and 4, exp is function 3:
     then     3 * a * q[l] * q[k1] / exp((b - 1.57)^2) + 1.57
     would result in
     400000000 100000000 2 100000001 2 100000002 2
     100000003 400000001 0 400000002 3 200000003 3 400000001 1
     where 3 = real 0, 1.57 = real 1, 2 = real 2
     a = vars 0, q[l] vars 1, q[k1] vars 2, exp functs 3
  */
{
  int i, j, error, op, id, stack = 0, l_deco, l_double;
  int up[100][3] = {{-1, -1, -1}};

  l_deco = deco->curr = 0;
  l_double = doubles->curr;
  error = scan_expr(c_item, item);
  if (error) return error;
  for (i = 0; i < cat->curr; i++)
  {

    /* categories: 1: variable, 3: floating constant, 4: operator
       6: left par., 7: right par.     */
    switch (cat->i[i])
    {
      case 1:                              /* variable */
        if (l_deco == deco->max) grow_int_array(deco);
        deco->i[l_deco++] = 100000000 + d_var->i[i];
        break;
      case 3:                              /* constant */
        if (l_deco == deco->max) grow_int_array(deco);
        if (l_double == doubles->max) grow_double_array(doubles);
        doubles->a[l_double] = cat_doubles->a[i];
        deco->i[l_deco++] = 400000000 + l_double++;
        doubles->curr = l_double;
        break;
      case 4:
        if ((op = oper->i[i]) < 5)           /* operator */
        {
          id = op / 2;
          for (j = 2; j >= id; j--)
          {
            if (up[stack][j] > -1)
            {
              if (l_deco == deco->max) grow_int_array(deco);
              deco->i[l_deco++] = up[stack][j];
              up[stack][j] = -1;
            }
          }
          up[stack][id] = op;
        }
        else
        {
          if (l_deco == deco->max) grow_int_array(deco);
          deco->i[l_deco++] = 200000000 + func->i[i];  /* function */
        }
        break;
      case 6:      /*  '(' */
        stack++;
        for (j = 0; j < 3; j++)  up[stack][j] = -1;
        break;
      case 7:      /*  ')' */
        for (j = 2; j >= 0; j--)
        {
          if (up[stack][j] > -1)
          {
            if (l_deco == deco->max) grow_int_array(deco);
            deco->i[l_deco++] = up[stack][j];
          }
        }
        stack--;
        break;
      default:
        return 9;
    }   /* end switch */
  }     /* end loop over categories */
  for (j = 2; j >= 0; j--)   /* clear stack */
  {
    if (up[stack][j] > -1)
    {
      if (l_deco == deco->max) grow_int_array(deco);
      deco->i[l_deco++] = up[stack][j];
    }
  }
  deco->curr = l_deco;
  return 0;
}

double polish_value(struct int_array* deco)  /* coded input (see below) */
  /* description see polish_expression */
{
  int i, k, kc, c_stack = -1;
  double stack[MAX_ITEM];
  char tmp[20];

  if (++polish_cnt > MAX_LOOP)
    fatal_error("circular call in expression","evaluation");
  stack[0] = 0;
  for (i = 0; i < deco->curr; i++)   /* decoding loop */
  {
    k = deco->i[i];
    if ( k < 5)     /* operator */
    {
      if (c_stack < 0)
      {
        fatal_error("polish_value","stack underflow in Polish decoding");
      }
      else if (c_stack == 0)
      {
        stack[1] = stack[0]; stack[0] = 0;
      }
      else  c_stack--;
      
      
      switch(k)
      {
        case 0:
          stack[c_stack] -= stack[c_stack+1];
          break;
        case 1:
          stack[c_stack] += stack[c_stack+1];
          break;
        case 2:
          stack[c_stack] *= stack[c_stack+1];
          break;
        case 3:
          if (stack[c_stack+1] == 0.0)
           {
             warning("polish_value: division by zero","k=3, stack[c_stack+1]=%f. Putting result to 0!",stack[c_stack+1]);
             /*fatal_error("Division by zero is not defined","Aborting!");*/
             stack[c_stack] = 0.0;
             /*stack[c_stack] = nan("");*/
             break;
           }  
          stack[c_stack] /= stack[c_stack+1];
          break;
        case 4:
          stack[c_stack] = pow(stack[c_stack],stack[c_stack+1]);
          break;
        default:
          fatal_error("polish_value","illegal operator in Polish decoding");
      }
    }
    else
    {
      kc = k / 100000000;  k -= 100000000 * kc;
      switch(kc)
      {
        case 1:            /* variable */
          stack[++c_stack] = act_value(k, expr_chunks);
          break;
        case 4:            /* real constant */
          stack[++c_stack] = doubles->a[k];
          break;
        case 2:            /* function */
          switch(k-1)      /* the offset is due to dummyfunction */
          {
            case 0:
              stack[c_stack] = fabs(stack[c_stack]);
              break;
            case 1:
              stack[c_stack] = sqrt(stack[c_stack]);
              break;
            case 2:
              stack[c_stack] = exp(stack[c_stack]);
              break;
            case 3:
              stack[c_stack] = log(stack[c_stack]);
              break;
            case 4:
              stack[c_stack] = log10(stack[c_stack]);
              break;
            case 5:
              stack[c_stack] = sin(stack[c_stack]);
              break;
            case 6:
              stack[c_stack] = cos(stack[c_stack]);
              break;
            case 7:
              stack[c_stack] = tan(stack[c_stack]);
              break;
            case 8:
              stack[c_stack] = asin(stack[c_stack]);
              break;
            case 9:
              stack[c_stack] = acos(stack[c_stack]);
              break;
            case 10:
              stack[c_stack] = atan(stack[c_stack]);
              break;
            case 11:
              stack[c_stack] = sinh(stack[c_stack]);
              break;
            case 12:
              stack[c_stack] = cosh(stack[c_stack]);
              break;
            case 13:
              stack[c_stack] = tanh(stack[c_stack]);
              break;
            case 14:
              stack[c_stack] = frndm();
              break;
            case 15:
              stack[c_stack] = grndm();
              break;
            case 16:
              stack[c_stack] = tgrndm(stack[c_stack]);
              break;
            case 17:
              stack[c_stack] = table_value();
              break;
            default:
              fatal_error("polish_value",
                          "illegal function in Polish decoding");
          }
          break;
        default:
          /* if you get here, then most likely someone has created
             more than 100000000 double precision numbers */
          sprintf(tmp, "%d", k-1);
          fatal_error("illegal type in Polish decoding: ", tmp);
          exit(1);
      }
    }
  }       /* end of decoding loop */
  polish_cnt--;
  return stack[0];
}

void pre_split(char* inbuf, struct char_array* outbuf, int fill_flag)
  /* inserts blanks between tokens */
  /* fill_flag != 0 makes a 0 to be inserted into an empty "()" */
{
  char c, cp = ' ', cpnb = ' ', quote = ' ';
  int j, k, kn, sl = strlen(inbuf), cout = 0, quote_level = 0, rb_level = 0;
  int left_b = 0, new_string = 1, c_digit = 0, f_equal = 0, comm_cnt = 0;
  while (2*strlen(inbuf) > outbuf->max) grow_char_array(outbuf);
  for (k = 0; k < sl; k++)
  {
    c = inbuf[k];
    if (quote_level > 0)
    {
      if (c == quote)
      {
        quote_level--; outbuf->c[cout++] = c; outbuf->c[cout++] = ' ';
      }
      else outbuf->c[cout++] = c == ' ' ? '@' : c;
    }
    else
    {
      c = inbuf[k];
      switch (c)
      {
        case '\"':
        case '\'':
          quote = c;
          quote_level++; outbuf->c[cout++] = ' '; outbuf->c[cout++] = c;
          break;
        case '-':
          if (inbuf[k+1] == '>')
          {
            outbuf->c[cout++] = c; break;
          }
        case '+':
          if (left_b > 0)
          {
            outbuf->c[cout++] = ' ';
            outbuf->c[cout++] = '0';
            outbuf->c[cout++] = ' ';
            left_b = 0;
          }
          if (!(new_string > 0 && c_digit > 0 && strchr("ed",cp)) && cout > 0)
            outbuf->c[cout++] = ' ';
          outbuf->c[cout++] = c;
          if (!(new_string > 0 && c_digit > 0 && strchr("ed",cp)))
          {
            outbuf->c[cout++] = ' ';
            new_string = 1;
          }
          break;
        case '(':
          rb_level++;
          left_b = 1;
          new_string = 1;
          outbuf->c[cout++] = ' ';
          outbuf->c[cout++] = c;
          outbuf->c[cout++] = ' ';
          break;
        case '>':
          if (cout > 0 && outbuf->c[cout-1] == '-')
          {
            outbuf->c[cout++] = c; break;
          }
        case ')':
          rb_level--;
          if (fill_flag && cpnb == '(') outbuf->c[cout++] = '0';
        case '<':
        case ':':
        case '*':
        case '/':
        case '^':
        case '{':
        case '}':
        case '[':
        case ']':
        case '|':
        case '&':
          left_b = 0;
          new_string = 1;
          outbuf->c[cout++] = ' ';
          outbuf->c[cout++] = c;
          outbuf->c[cout++] = ' ';
          break;
        case '=':
          f_equal = 1;
          left_b = 0;
          new_string = 1;
          outbuf->c[cout++] = ' ';
          outbuf->c[cout++] = c;
          outbuf->c[cout++] = ' ';
          break;
        case ',': /* kept behind first "=", or if not first "," */
          /* not kept inside round brackets before '=' */
          left_b = 0;
          new_string = 1;
          outbuf->c[cout++] = ' ';
          if (f_equal || (comm_cnt && rb_level == 0))
          {
            outbuf->c[cout++] = c;
            outbuf->c[cout++] = ' ';
          }
          comm_cnt++;
          break;
        case ';':
          left_b = 0;
          new_string = 1;
          outbuf->c[cout++] = ' ';
          break;
        default:
          if (c == ' ') outbuf->c[cout++] = c;
          else
          {
            left_b = 0;
            if (new_string && (isdigit(c) || c == '.'))
            {
              kn = get_val_num(inbuf, k, sl);
              for (j = k; j <= kn; j++) outbuf->c[cout++] = inbuf[j];
              outbuf->c[cout++] = ' ';
              k = kn;
            }
            else
            {
              new_string = 0;
              if (cout > 0 || c != ' ') outbuf->c[cout++] = c;
            }
          }
      }
      cp = c; if (c != ' ') cpnb = c;
    }
  }
  outbuf->c[cout] = '\0';
}

void process()  /* steering routine: processes one command */
{
  int pos;
  char* name;
  struct element* el;
  if (this_cmd != NULL)
  {
    switch (this_cmd->type)
    {
      case 0: /* executable commands */
        exec_command();
        if (stop_flag)
        {
          if (this_cmd)
          {
            if (this_cmd->clone != NULL)
              this_cmd->clone = delete_command(this_cmd->clone);
            this_cmd = delete_in_cmd(this_cmd);
          }
          return;
        }
        break;
      case 1: /* element definition */
        enter_element(this_cmd);
        buffer_in_cmd(this_cmd);
        break;
      case 2: /* variable definition */
        enter_variable(this_cmd);
        break;
      case 3: /* sequence start or end */
        enter_sequence(this_cmd);
        break;
      case 4:
        name = this_cmd->tok_list->p[0];
        if (sequ_is_on)
          /* element or sequence reference in sequence */
        {
          if ((pos =
               name_list_pos(name, sequences->list)) < 0)
            enter_element(this_cmd);
          else
          {
            this_cmd->cmd_def  = find_command("sequence", defined_commands);
            this_cmd->clone = clone_command(this_cmd->cmd_def);
            strcpy(this_cmd->clone->name, name);
            scan_in_cmd(this_cmd);
            enter_sequ_reference(this_cmd, sequences->sequs[pos]);
          }
        }
        else
          /* element parameter definition */
        {
          if ((el = find_element(name, element_list)) == NULL)
            warning("skipped, command or element unknown:", name);
          else
          {
            this_cmd->cmd_def = el->def;
            this_cmd->clone = clone_command(this_cmd->cmd_def);
            strcpy(this_cmd->clone->name, name);
            scan_in_cmd(this_cmd);
            update_element(el, this_cmd->clone);
          }
        }
        break;
      default:
        warning("unknown command type:",
                join_b(this_cmd->tok_list->p,
                       this_cmd->tok_list->curr));
    }
    if (this_cmd != NULL && (this_cmd->type == 0 || this_cmd->type == 2))
    {
      if (this_cmd->clone != NULL)
      {
        if (this_cmd->clone_flag == 0)
          this_cmd->clone = delete_command(this_cmd->clone);
        else add_to_command_list(this_cmd->clone->name,
                                 this_cmd->clone, stored_commands, 0);
      }
      if (this_cmd->label != NULL) buffer_in_cmd(this_cmd);
      else this_cmd = delete_in_cmd(this_cmd);
    }
  }
}

void pro_input(char* statement)
  /* processes one special (IF() etc.), or one normal statement after input */
{
  int type, code, nnb, ktmp;
  char* sem;
  int rs, re, start = 0, l = strlen(statement);
  while (start < l)
  {
    if ((type = in_spec_list(&statement[start])))
    {
      if (type == 6)
      {
        get_bracket_range(&statement[start], '(', ')', &rs, &re);
        ktmp = re+1;
        if (re > rs && strchr(&statement[ktmp], ':')) /* formal arg.s */
        {
          get_bracket_range(&statement[ktmp], '(', ')', &rs, &re);
          rs += ktmp; re += ktmp;
        }
      }
      else get_bracket_range(&statement[start], '{', '}', &rs, &re);
      if (re > rs)
      {
        re += start + 1;
        if (re < l && next_non_blank(&statement[re]) == ';')
          re += next_non_blank_pos(&statement[re]) + 1;
      }
      if((code = act_special(type, &statement[start])) < 0)
      {
        if (get_option("warn"))
        {
          switch (code)
          {
            case -1:
              warning("statement illegal in this context,", "skipped");
              break;
            case -2:
              warning("statement label is protected keyword,", "skipped");
            case -3:
              warning("statement not recognised:", statement);
              break;
            default:
              fatal_error("illegal return code","from act_special");
          }
        }
      }
      if (re > rs && re < l)
      {
        if (next_non_blank_pos(&statement[re]) < 0) start = l;
        else start = re;
      }
      else start = l;
    }
    else
    {
      if ((sem = mystrchr(&statement[start], ';')) == NULL) return;
      if ((int) sem > (int) &statement[start]) /* skip empty ';' */
      {
        *sem = '\0';
        this_cmd = new_in_cmd(400);
        pre_split(&statement[start], work, 1);
        check_table(work->c);
        this_cmd->tok_list->curr = mysplit(work->c, this_cmd->tok_list);
        if ((type = decode_command()) < 0) /* error */
        {
          if (get_option("warn"))
          {
            switch (type)
            {
              case -1:
                warning("statement illegal in this context,", "skipped");
                break;
              case -2:
                warning("statement label is protected keyword,","skipped");
              case -3:
                warning("statement not recognised:", work->c);
                break;
              default:
                fatal_error("illegal return code","from decode_command");
            }
          }
        }
        else process();
        if (stop_flag)  return;
        *sem = ';';
      }
      sem++;
      start = (int)sem - (int)statement;
      if (start < l)
      {
        if ((nnb = next_non_blank_pos(sem)) < 0)  start = l;
        else start += nnb;
      }
    }
  }
}

void put_info(char* t1, char* t2)
{
  if (get_option("info") && get_option("warn"))
    printf("++++++ info: %s %s\n",t1,t2);
}

void remove_range(char* string, char* s1, char* s2)
  /* remove portion s1...s2 (included) in string */
{
  char* ps1 = strstr(string, s1);
  char* ps2 = strstr(string, s2);
  if (ps1 != NULL && ps2 != NULL && (int) ps1 < (int) ps2)
  {
    ps2++; ps2++;
    while (*ps2 != '\0')  *ps1++ = *ps2++;
    *ps1 = '\0';
  }
}

void remove_upto(char* string, char* s1)
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

void resequence_nodes(struct sequence* sequ)
  /* resequences occurrence list */
{
  struct node* c_node = sequ->start;
  int i, cnt;
  while (c_node != NULL)
  {
    if (c_node->p_elem != NULL)
    {
      if ((i = name_list_pos(c_node->p_elem->name, occ_list)) < 0)
      {
        i = add_to_name_list(c_node->p_elem->name, 1, occ_list);
        cnt = 1;
      }
      else cnt = ++occ_list->inform[i];
      strcpy(c_node->name, compound(c_node->p_elem->name, cnt));
    }
    if (c_node == sequ->end) break;
    c_node = c_node->next;
  }
}

void save_beam(struct sequence* sequ, FILE* file)
{
  struct command* comm;
  char beam_buff[AUX_LG];
  int i, def = 0;
  if ((comm = find_command(sequ->name, beam_list)) == NULL)
  {
    if (default_beam_saved == 0)
    {
      def = default_beam_saved = 1;
      comm = find_command("default_beam", beam_list);
    }
  }
  if (comm != NULL)
  {
    beam_buff[0] = '\0';
    strcat(beam_buff, "beam");
    for (i = 0; i < comm->par->curr; i++)
    {
      if (comm->par_names->inform[i])
      {
        if (strcmp(comm->par_names->names[i], "sequence") != 0
            || def == 0)
          export_comm_par(comm->par->parameters[i], beam_buff);
      }
    }
    write_nice(beam_buff, file);
  }
}

int scan_expr(int c_item, char** item)   /* split input */

  /* scans expressions for parameters, elements, numbers, and functions

  categories: 1: variable, 3: floating constant, 4: operator
  operator types:
  1 = +, 2 = -, 3 = *, 4 = /, 5 = ^ (power), 6 = function (from functs)
  */

{
  int i, lp, lx = -1, l_cat = 0, level = 0, pos, f_level[MAX_ITEM];
  char c;
  char* bf;

  for (i = 0; i < c_item; i++)  /* loop over input items */
  {
    c = item[i][0];  /* first character of item i */
    if (c == '(')
    {
      f_level[level++] = 0;
      if (l_cat > 0)
      {
        if (cat->i[l_cat-1] < 4)  return 2;  /* error: missing operator */
        if (cat->i[l_cat-1] == 5)   /* function */
        {
          f_level[level-1] = func->i[--l_cat];
          if (l_cat == func->max)  grow_int_array(func);
          func->i[l_cat] = 0;
        }
      }
      if (l_cat == cat->max)  grow_int_array(cat);
      cat->i[l_cat++] = 6;
    }
    else if (c == ')')
    {
      if (level == 0)  return 1;  /* error: too many right brackets */
      if (l_cat == cat->max)  grow_int_array(cat);
      cat->i[l_cat++] = 7;
      level--;
      if (f_level[level] != 0)
      {
        if (l_cat == oper->max)  grow_int_array(oper);
        if (l_cat == func->max)  grow_int_array(func);
        if (l_cat == cat->max)  grow_int_array(cat);
        oper->i[l_cat] = 6;
        func->i[l_cat] = f_level[level];
        cat->i[l_cat++] = 4;
      }
    }
    else if (isalpha(c) || c == '_')  /* start of variable or function */
    {
      lp = 0;
      while (strlen(functs[lp]))
      {
        lx = lp;
        if (strcmp(item[i], functs[lp]) == 0)  break;
        lp++;
      }
      if (lx == lp)    /* function found */
      {
        if (l_cat == cat->max)  grow_int_array(cat);
        if (l_cat == func->max)  grow_int_array(func);
        cat->i[l_cat] = 5;
        func->i[l_cat++] = lp;
      }
      else
      {
        if (l_cat == cat->max)  grow_int_array(cat);
        if (l_cat == d_var->max)  grow_int_array(d_var);
        cat->i[l_cat] = 1;
        if ((pos = name_list_pos(item[i], expr_chunks)) < 0)
        {
          bf = permbuff(item[i]);
          d_var->i[l_cat++] = add_to_name_list(bf, 0, expr_chunks);
        }
        else d_var->i[l_cat++] = pos;
      }
    }
    else if (isdigit(c) || c == '.')  /* number */
    {
      if (l_cat == cat->max)  grow_int_array(cat);
      if (l_cat == cat_doubles->max) grow_double_array(cat_doubles);
      cat->i[l_cat] = 3;
      cat_doubles->a[l_cat++] = atof(item[i]);
    }
    else if (is_operator(c))
    {
      if (l_cat == cat->max) grow_int_array(cat);
      if (l_cat == oper->max) grow_int_array(oper);
      cat->i[l_cat] = 4;
      oper->i[l_cat++] = str_pos(op_string, c);
      /* oper->i[l_cat++] = (int)strchr(op_string, c) -(int)op_string;*/
    }
    else   return 2;  /* illegal character */
  }
  if (level != 0)  return 1;  /* unclosed parentheses */
  cat->curr = l_cat;
  return 0;
}

void scan_in_cmd(struct in_cmd* cmd)
  /* reads a command into a clone of the original */
{
  int cnt = 0, /* gives position in command (from 1) */
    i, k, log, n;
  struct name_list* nl = cmd->clone->par_names;
  for (i = 0; i < nl->curr; i++) nl->inform[i] = 0; /* set when read */
  n = cmd->tok_list->curr;
  i = cmd->decl_start;
  cmd->tok_list->p[n] = blank;
  while (i < n)
  {
    log = 0;
    if (i+1 < n && *cmd->tok_list->p[i] == '-')
    {
      log = 1; i++;
    }
    if (*cmd->tok_list->p[i] != ',')
    {
      if ((k = name_list_pos(cmd->tok_list->p[i],
                             cmd->cmd_def->par_names)) < 0)  /* try alias */
      {
        if ((k = name_list_pos(alias(cmd->tok_list->p[i]),
                               cmd->cmd_def->par_names)) < 0)
          warning("illegal keyword:", cmd->tok_list->p[i]);
        break;
      }
      else if ((i = decode_par(cmd, i, n, k, log)) < 0)
      {
        warning("illegal format near:", cmd->tok_list->p[-i]);
        break;
      }
      cmd->clone->par_names->inform[k] = ++cnt; /* mark parameter as read */
    }
    i++;
  }
}

void seq_edit_ex(struct sequence* seq)
{
  edit_sequ = seq;
  edit_is_on = 1;
  seqedit_install = seqedit_move = seqedit_remove = 0;
  if (edit_sequ->ex_start != NULL)
  {
    edit_sequ->ex_nodes = delete_node_list(edit_sequ->ex_nodes);
    edit_sequ->ex_start = delete_node_ring(edit_sequ->ex_start);
  }
  if (occ_list == NULL)
    occ_list = new_name_list("occ_list", 10000);  /* for occurrence count */
  else occ_list->curr = 0;
  resequence_nodes(edit_sequ);
  all_node_pos(edit_sequ);
}

void seq_end_ex()
{
  occ_list->curr = 0;
  resequence_nodes(edit_sequ);
  selected_ranges->curr = 0;
  selected_ranges->list->curr = 0;
  edit_is_on = 0;
}

void seq_flatten(struct sequence* sequ)
  /* executes flatten command */
{
  struct node* c_node;
  struct node_list* nl;
  if (occ_list == NULL)
    occ_list = new_name_list("occ_list", 10000);  /* for occurrence count */
  else occ_list->curr = 0;
  make_occ_list(sequ);
  all_node_pos(sequ);
  sequ->ex_nodes = new_node_list(2*sequ->nodes->curr);
  expand_sequence(sequ, 0);
  sequ->nested = 0;
  nl = sequ->nodes;
  sequ->nodes = sequ->ex_nodes;
  sequ->ex_nodes = delete_node_list(nl);
  sequ->start = sequ->ex_start; sequ->ex_start = NULL;
  sequ->end = sequ->ex_end; sequ->ex_end = NULL;
  c_node = sequ->start;
  while (c_node != NULL)
  {
    c_node->at_value = c_node->position;
    c_node->at_expr = NULL;
    c_node->from_name = NULL;
    if (c_node == sequ->end) break;
    c_node = c_node->next;
  }
}

void set_defaults(char* string) /* reset options, beam etc. to defaults */
{
  int i, pos;
  struct command* beam_clone;

  if ((pos = name_list_pos(string, defined_commands->list)) > -1)
  {
    if (strcmp(string, "option") == 0)
    {
      if (options != NULL) delete_command(options);
      options = clone_command(defined_commands->commands[pos]);
    }
    else if (strcmp(string, "set") == 0)
      store_set(defined_commands->commands[pos], 0);
    else if (strcmp(string, "setplot") == 0)
    {
      if (plot_options != NULL) delete_command(plot_options);
      plot_options = clone_command(defined_commands->commands[pos]);
    }
    else if (strcmp(string, "threader") == 0)
    {
      if (threader_par != NULL)  delete_command(threader_par);
      threader_par = clone_command(defined_commands->commands[pos]);
    }
    else if (strcmp(string, "beam") == 0)
    {
      if (current_beam == NULL)
        current_beam = clone_command(defined_commands->commands[pos]);
      beam_clone = clone_command(defined_commands->commands[pos]);
      for (i = 0; i < beam_clone->par_names->curr; i++)
        beam_clone->par_names->inform[i] = 1; /* mark as "read" */
      update_beam(beam_clone);
      delete_command(beam_clone);
    }
  }
}

void set_option(char* str, int* opt)
  /* sets an (old or new) option with name "str",
     value *opt (0 flase, 1 true) */
{
  int i, j, k;
  char* bc;
  mycpy(c_dum->c, str); bc = permbuff(c_dum->c);
  if ((i = name_list_pos(bc, options->par_names)) < 0)
  {
    j = add_to_name_list(bc, 0, options->par_names);
    if ((k = options->par->curr) == options->par->max)
      grow_command_parameter_list(options->par);
    options->par->parameters[options->par->curr++]
      = new_command_parameter(bc, 0);
    options->par->parameters[k]->double_value = *opt;
  }
  else options->par->parameters[i]->double_value = *opt;
}

void set_sub_variable(char* comm, char* par, struct in_cmd* cmd)
{
  char* p;
  struct element* el;
  struct command *command, *keep_beam = current_beam;
  int end, start = cmd->decl_start, t_num, exp_type;
  double val = 0;
  for (t_num = start; t_num < cmd->tok_list->curr; t_num++)
    if (*(cmd->tok_list->p[t_num]) == ',') break;
  exp_type = loc_expr(cmd->tok_list->p, t_num,
                      start, &end);
  if (exp_type == 1) /* literal constant */
    val = simple_double(cmd->tok_list->p, start, end);
  else if (polish_expr(end + 1 - start, &cmd->tok_list->p[start]) == 0)
    val = polish_value(deco);
  if (strncmp(comm, "beam", 4) == 0)
  {
    command = current_beam = find_command("default_beam", beam_list);
    if ((p = strchr(comm, '%')) != NULL)
    {
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

void set_command_par_value(char* parameter, struct command* cmd, double val)
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

void show_beam(char* tok)
{
  struct command* comm;
  if (strlen(tok) > 5 && tok[4] == '%')
    comm = find_command(&tok[5], beam_list);
  else comm = find_command("default_beam", beam_list);
  if (comm != NULL) dump_command(comm);
}

double simple_double(char** toks, int start, int end)
{
  if (start > end && start + 1 != end)  return INVALID;
  else if (start == end)  return atof(toks[start]);
  else
  {
    if (*toks[start] == '-') return -atof(toks[end]);
    else if (*toks[start] == '+') return atof(toks[end]);
    else return INVALID;
  }
}

int simple_logic_expr(int nit, char* toks[])
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
      val1 = polish_value(deco);
    else return -1;
  }
  else  val1 = simple_double(toks, l1_start, l1_end);
  if ((t2 = loc_expr(toks, nit, l2_start, &l2_end)) == 0) return -1;
  if (t2 == 2)
  {
    if (polish_expr(l2_end + 1 - l2_start, &toks[l2_start]) == 0)
      val2 = polish_value(deco);
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

char* spec_join(char** it_list, int n)
  /* replaces variable in table(...) by original string */
{
  char rout_name[] = "spec_join";
  int j;
  char** p;
  struct variable* var;
  *c_join->c = '\0';
  if (n > 0)
  {
    p = (char**) mymalloc(rout_name,n*sizeof(char*));
    for (j = 0; j < n; j++) p[j] = it_list[j];
    for (j = 0; j < n; j++)
      if (strcmp(p[j], "table") == 0 && j+3 < n
          && (var = find_variable(p[j+2], variable_list)) != NULL)
        p[j+2] = var->string;
    for (j = 0; j < n; j++) strcat(c_join->c, p[j]);
    myfree(rout_name, p);
  }
  return c_join->c;
}

void store_command_def(char* cmd_string)  /* processes command definition */
{
  int i, n, j, b_s = 0, r_start, r_end, b_cnt;
  struct element* el;
  struct command* cmd;
  struct command_parameter* p;
  struct in_cmd* tmp_cmd = new_in_cmd(1000);
  struct char_p_array* toks = tmp_cmd->tok_list;

  pre_split(cmd_string, work, 0);
  n = mysplit(work->c, toks);
  if (n < 6 || *toks->p[1] != ':') fatal_error("illegal command:", cmd_string);
  if (defined_commands->curr == defined_commands->max)
    grow_command_list(defined_commands);
  cmd = defined_commands->commands[defined_commands->curr++] =
    new_command(toks->p[0], 40, 40, toks->p[2], toks->p[3],
                atoi(toks->p[4]), atoi(toks->p[5]));
  
/*  
    printf("skowrondebug: store_command_def.c command name %s\n",cmd->name);
    if (strcmp(cmd->name,"twcavity") == 0)
    {
    printf("skowrondebug: store_command_def.c I have got TWCAVITY\n");
    }
*/  
  i = add_to_name_list(cmd->name, 0, defined_commands->list);
  if (n > 6)
  {
    b_cnt = string_cnt('[', n, toks->p);
    for (i = 0; i < b_cnt; i++)  /* loop over parameter definitions */
    {
      get_bracket_t_range(toks->p, '[', ']', b_s, n, &r_start, &r_end);
      if (r_start < b_s || r_end - r_start < 2)
        fatal_error("empty or illegal cmd parameter definition:",
                    cmd->name);
      if (cmd->par->curr == cmd->par->max)
        grow_command_parameter_list(cmd->par);
      p = store_comm_par_def(toks->p, r_start+1, r_end-1);
      if (p == NULL) fatal_error("illegal cmd parameter definition:",
                                 cmd->name);
      cmd->par->parameters[cmd->par->curr++] = p;
      j = add_to_name_list(p->name, 1, cmd->par_names);
      b_s = r_end + 1;
    }
  }
  if (strcmp(toks->p[2], "element") == 0)
    el = make_element(toks->p[0], toks->p[0], cmd, 0);
  delete_in_cmd(tmp_cmd);
}

struct command_parameter* store_comm_par_def(char* toks[], int start, int end)
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
      for (j = 0; j <= mymin((end - start),2); j++)
      {
        if (strcmp(toks[start+j], "true") == 0)  pl[jj]->double_value = 1;
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
            if (*toks[s_start+j] != ',' &&
                strcmp(toks[s_start+j], none) != 0)
              pl[0]->m_string->p[i++] = permbuff(toks[s_start+j]);
          }
          pl[0]->m_string->curr = i;
          pl[0]->type += 10;
        }
      }
  }
  return pl[0];
}

void store_comm_par_value(char* parameter, double val, struct command* cmd)
{
  struct command_parameter* cp;
  int i;
  if ((i = name_list_pos(parameter, cmd->par_names)) > -1)
  {
    cp = cmd->par->parameters[i];
    cp->type = 2;
    if(cp->expr != NULL) cp->expr = delete_expression(cp->expr);
    cp->double_value = val;
  }
}

void store_node_value(char* par, double* value)
  /* stores value for parameter par at current node */
{
  char lpar[NAME_L];
  struct element* el = current_node->p_elem;

  mycpy(lpar, par);
  if (strcmp(lpar, "chkick") == 0) current_node->chkick = *value;
  else if (strcmp(lpar, "cvkick") == 0) current_node->cvkick = *value;
  else if (strcmp(lpar, "dipole_bv") == 0) current_node->dipole_bv = *value;
  else if (strcmp(lpar, "other_bv") == 0) current_node->other_bv = *value;
  else if (strcmp(lpar, "obs_point") == 0) current_node->obs_point = *value;
  else if (strcmp(lpar, "sel_sector") == 0) current_node->sel_sector = *value;
  else if (strcmp(lpar, "enable") == 0) current_node->enable = *value;

  /* added by E. T. d'Amico 27 feb 2004 */

  else if (strcmp(lpar, "e1") == 0) store_comm_par_value("e1",*value,el->def);
  else if (strcmp(lpar, "e2") == 0) store_comm_par_value("e2",*value,el->def);
  else if (strcmp(lpar, "angle") == 0) store_comm_par_value("angle",*value,el->def);

  /* added by E. T. d'Amico 12 may 2004 */

  else if (strcmp(lpar, "h1") == 0) store_comm_par_value("h1",*value,el->def);
  else if (strcmp(lpar, "h2") == 0) store_comm_par_value("h2",*value,el->def);
  else if (strcmp(lpar, "fint") == 0) store_comm_par_value("fint",*value,el->def);
  else if (strcmp(lpar, "fintx") == 0) store_comm_par_value("fintx",*value,el->def);
  else if (strcmp(lpar, "hgap") == 0) store_comm_par_value("hgap",*value,el->def);

  /* end of additions */
}

void store_set(struct command* comm, int flag)
{
  char* p;
  char* name;
  struct command_parameter* cp;
  struct name_list* nl = comm->par_names;
  int i, lp, n = 0, posf = name_list_pos("format", nl),
    poss = name_list_pos("sequence", nl);
  if (flag == 0 || (posf > -1 && nl->inform[posf]))
  {
    n++;
    cp = comm->par->parameters[posf];
    for (i = 0; i < cp->m_string->curr; i++)
    {
      p = noquote(cp->m_string->p[i]);
      if (strchr(p, 's')) strcpy(string_format, p);
      else if (strpbrk(p, "id")) strcpy(int_format, p);
      else if (strpbrk(p, "feEgG")) strcpy(float_format, p);
    }
  }
  if (flag != 0 && poss > -1 && nl->inform[poss])
  {
    n++;
    name = comm->par->parameters[poss]->string;
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

int table_row(struct table* table, char* name)
{
  int i, j, ret = -1;
  for (i = 0; i < table->num_cols; i++)
   {
     if(table->columns->inform[i] == 3) 
      {
        if (debuglevel > 2) 
          printf("table_row: Column %d named <<%s>> is of strings. We use it to find the name.\n",
                 i,table->columns->names[i]);
        break;
      } 
   }  
  
  if (i < table->num_cols)
   {
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
  if(ret==-1) warning("table_row","Name of row not found: %s,", name);
  return ret;
}

double table_value()
{
  double val = zero;
  int ntok, pos, col, row;
  char** toks;
  struct table* table;
  if (current_variable != NULL && current_variable->string != NULL)
  {
    strcpy(c_dum->c, current_variable->string);
    supp_char(',', c_dum->c);
    mysplit(c_dum->c, tmp_p_array);
    toks = tmp_p_array->p; ntok = tmp_p_array->curr;
    if (ntok > 1)
    {
      if ((pos = name_list_pos(toks[0], table_register->names)) > -1)
      {
        table = table_register->tables[pos];
        if ((col = name_list_pos(toks[ntok-1], table->columns)) > -1)
        {
          if (ntok > 2)  /* find row - else current (dynamic), or 0 */
          {
            row = table_row(table, toks[1]);
          }
          else if (table->dynamic)  row = table->curr;
          else row = 0;
          val = table->d_cols[col][row];
        }
        else if ((ntok == 3) && ((col = name_list_pos(toks[1], table->columns)) > -1))
        {
          row = atoi(toks[2])-1;
          if(row < table->curr) val = table->d_cols[col][row];
        }
      }
    }
  }
  return val;
}

int tab_name_code(char* name, char* t_name)
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

void time_stamp(char* place)
{
  time_t now;
  int k, l;
  time(&now);    /* get system time */
  k = (int)now - (int)start_time;
  l = (int)now - (int)last_time;
  last_time = now;
  fprintf(prt_file, "sec.s since start: %d   since last call: %d\n", k, l);
}

char* tmpbuff(char* string)
  /* buffers string in a temporary (i.e. allocated) buffer */
{
  char* p;
  p = (char*) mymalloc("tmpbuff",strlen(string)+1);
  strcpy(p, string);
  return p;
}

void update_beam(struct command* comm)
  /* calculates consistent values for modified beam data set.
     beam command values are evaluated in the order:
     particle->(mass+charge)
     energy->pc->gamma->beta
     ex->exn
     ey->eyn
     current->npart
     et->sigt->sige
     where any item to the left takes precendence over the others;
     for ions, the input energy is multiplied by the charge, and the
  */
{
  struct name_list* nlc = comm->par_names;
  struct command_parameter_list* plc = comm->par;
  struct command_parameter_list* pl = current_beam->par;
  int pos, lp;
  char* name = blank;
  double energy = 0, beta = 0, gamma = 0, charge = 0, freq0 = 0, bcurrent = 0,
    npart = 0, mass = 0, pc = 0, ex, exn, ey, eyn, alfa, circ = 0,
    arad = 0;
  pos = name_list_pos("particle", nlc);
  if (nlc->inform[pos])  /* parameter has been read */
  {
    pl->parameters[pos]->string = name
      = plc->parameters[pos]->string;
    if ((lp = name_list_pos(name, defined_commands->list)) > -1)
    {
      mass = command_par_value("mass", defined_commands->commands[lp]);
      charge = command_par_value("charge", defined_commands->commands[lp]);
    }
    else /* unknown particle, then mass and charge must be given as well */
    {
      pos = name_list_pos("mass", nlc);
      if (nlc->inform[pos]) mass = command_par_value("mass", comm);
      else
      {
        warning("emass given to unknown particle:", name);
        mass = get_variable("emass");
      }
      pos = name_list_pos("charge", nlc);
      if (nlc->inform[pos]) charge = command_par_value("charge", comm);
      else
      {
        warning("charge +1 given to unknown particle:", name);
        charge = 1;
      }
    }
  }
  else name = pl->parameters[pos]->string;
  if (strcmp(name, "ion") == 0)
  {
    pos = name_list_pos("mass", nlc);
    if (nlc->inform[pos]) mass = command_par_value("mass", comm);
    pos = name_list_pos("charge", nlc);
    if (nlc->inform[pos]) charge = command_par_value("charge", comm);
    else charge = command_par_value("charge", current_beam);
  }
  if (mass == zero) mass = command_par_value("mass", current_beam);
  if (charge == zero) charge = command_par_value("charge", current_beam);
  arad = ten_m_16 * charge * charge * get_variable("qelect")
    * clight * clight / mass;
  if ((pos = name_list_pos("energy", nlc)) > -1 && nlc->inform[pos])
  {
    energy = command_par_value("energy", comm);
    if (energy <= mass) fatal_error("energy must be","> mass");
    pc = sqrt(energy*energy - mass*mass);
    gamma = energy / mass;
    beta = pc / energy;
  }
  else if((pos = name_list_pos("pc", nlc)) > -1 && nlc->inform[pos])
  {
    pc = command_par_value("pc", comm);
    energy = sqrt(pc*pc + mass*mass);
    gamma = energy / mass;
    beta = pc / energy;
  }
  else if((pos = name_list_pos("gamma", nlc)) > -1 && nlc->inform[pos])
  {
    if ((gamma = command_par_value("gamma", comm)) <= one)
      fatal_error("gamma must be","> 1");
    energy = gamma * mass;
    pc = sqrt(energy*energy - mass*mass);
    beta = pc / energy;
  }
  else if((pos = name_list_pos("beta", nlc)) > -1 && nlc->inform[pos])
  {
    if ((beta = command_par_value("beta", comm)) >= one)
      fatal_error("beta must be","< 1");
    gamma = one / sqrt(one - beta*beta);
    energy = gamma * mass;
    pc = sqrt(energy*energy - mass*mass);
  }
  else
  {
    energy = command_par_value("energy", current_beam);
    if (energy <= mass) fatal_error("energy must be","> mass");
    pc = sqrt(energy*energy - mass*mass);
    gamma = energy / mass;
    beta = pc / energy;
  }
  if (nlc->inform[name_list_pos("ex", nlc)])
  {
    ex = command_par_value("ex", comm);
    exn = ex * 4 * beta * gamma;
  }
  else if (nlc->inform[name_list_pos("exn", nlc)])
  {
    exn = command_par_value("exn", comm);
    ex = exn / (4 * beta * gamma);
  }
  else
  {
    ex = command_par_value("ex", current_beam);
    exn = command_par_value("exn", current_beam);
  }
  if (nlc->inform[name_list_pos("ey", nlc)])
  {
    ey = command_par_value("ey", comm);
    eyn = ey * 4 * beta * gamma;
  }
  else if (nlc->inform[name_list_pos("eyn", nlc)])
  {
    eyn = command_par_value("eyn", comm);
    ey = eyn / (4 * beta * gamma);
  }
  else
  {
    ey = command_par_value("ey", current_beam);
    eyn = command_par_value("eyn", current_beam);
  }
  alfa = one / (gamma * gamma);
  if (nlc->inform[name_list_pos("circ", nlc)])
  {
    circ = command_par_value("circ", comm);
    if (circ > zero) freq0 = (beta * clight) / (ten_p_6 * circ);
  }
  else if (nlc->inform[name_list_pos("freq0", nlc)])
  {
    freq0 = command_par_value("eyn", comm);
    if (freq0 > zero) circ = (beta * clight) / (ten_p_6 * freq0);
  }
  else if ((pos = name_list_pos(name, sequences->list)) >= 0)
  {
    circ = sequences->sequs[pos]->length;
    freq0 = (beta * clight) / (ten_p_6 * circ);
  }
  if (nlc->inform[name_list_pos("bcurrent", nlc)])
  {
    bcurrent = command_par_value("bcurrent", comm);
    if (freq0 > zero)
      npart = bcurrent / (beta * freq0 * ten_p_6 * get_variable("qelect"));
  }
  else if (nlc->inform[name_list_pos("npart", nlc)])
  {
    npart = command_par_value("npart", comm);
    bcurrent = npart * beta * freq0 * ten_p_6 * get_variable("qelect");
  }
  pos = name_list_pos("bunched", nlc);
  if (nlc->inform[pos])
    pl->parameters[pos]->double_value = plc->parameters[pos]->double_value;
  pos = name_list_pos("radiate", nlc);
  if (nlc->inform[pos])
    pl->parameters[pos]->double_value = plc->parameters[pos]->double_value;
  pos = name_list_pos("et", nlc);
  if (nlc->inform[pos])
    pl->parameters[pos]->double_value = plc->parameters[pos]->double_value;
  pos = name_list_pos("sigt", nlc);
  if (nlc->inform[pos])
    pl->parameters[pos]->double_value = plc->parameters[pos]->double_value;
  pos = name_list_pos("sige", nlc);
  if (nlc->inform[pos])
    pl->parameters[pos]->double_value = plc->parameters[pos]->double_value;
  pos = name_list_pos("kbunch", nlc);
  if (nlc->inform[pos])
    pl->parameters[pos]->double_value = plc->parameters[pos]->double_value;
  pos = name_list_pos("bv", nlc);
  if (nlc->inform[pos])
    pl->parameters[pos]->double_value = plc->parameters[pos]->double_value;
  pos = name_list_pos("pdamp", nlc);
  if (nlc->inform[pos]) copy_double(plc->parameters[pos]->double_array->a,
                                    pl->parameters[pos]->double_array->a, 3);
  store_comm_par_value("mass", mass, current_beam);
  store_comm_par_value("charge", charge, current_beam);
  store_comm_par_value("energy", energy, current_beam);
  store_comm_par_value("pc", pc, current_beam);
  store_comm_par_value("gamma", gamma, current_beam);
  store_comm_par_value("ex", ex, current_beam);
  store_comm_par_value("exn", exn, current_beam);
  store_comm_par_value("ey", ey, current_beam);
  store_comm_par_value("eyn", eyn, current_beam);
  store_comm_par_value("npart", npart, current_beam);
  store_comm_par_value("bcurrent", bcurrent, current_beam);
  store_comm_par_value("freq0", freq0, current_beam);
  store_comm_par_value("circ", circ, current_beam);
  store_comm_par_value("beta", beta, current_beam);
  store_comm_par_value("alfa", alfa, current_beam);
  store_comm_par_value("arad", arad, current_beam);
}

void update_element(struct element* el, struct command* update)
  /* updates the parameters of el from those read into update */
{
  struct command_parameter_list* e_pl = el->def->par;
  struct command_parameter_list* pl = update->par;
  struct command_parameter *e_par, *par;
  int pos;
  for (pos = 0; pos < e_pl->curr; pos++)
  {
    if (update->par_names->inform[pos])  /* parameter has been read */
    {
      el->def->par_names->inform[pos]=update->par_names->inform[pos]; /*hbu activate this parameter in the element */
      e_par = e_pl->parameters[pos];
      par = pl->parameters[pos];
      switch (par->type)
      {
        case 0:
        case 1:
        case 2:
          e_par->double_value = par->double_value;
          e_par->expr = clone_expression(par->expr);
          /* fix for bv flag start */
          if (strcmp(e_par->name, "bv") == 0)
            el->bv = e_par->double_value;
          /* fix for bv flag end */
          break;
        case 3:
          e_par->string = permbuff(par->string);
          break;
        case 11:
        case 12:
          e_par->double_array = clone_double_array(par->double_array);
          e_par->expr_list = clone_expr_list(par->expr_list);
      }
    }
  }
}

void update_vector(struct expr_list* ell, struct double_array* da)
{
  int i;
  for (i = 0; i < ell->curr; i++)
  {
    if (ell->list[i] != NULL)
    {
      while (da->max < i) grow_double_array(da);
      da->a[i] = expression_value(ell->list[i], 2);
    }
  }
  if (da->curr < ell->curr)  da->curr = ell->curr;
}

double variable_value(struct variable* var)
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

void warning(char* t1, register char* fmt, ...)
{
/*prints warning on the standard error and accepts parameters printout with std C formatting*/
/*Piotr Skowronski CERN*/
  va_list         list;

  warn_numb++; /*I think that warnings should be counted even if the user does not want to see them*/
  fflush(0); /*flushes all the buffers -> so the warning appears in a correct place*/

  if (get_option("warn") == 0)
   {
     return;
   }

  va_start( list, fmt );

  fprintf(stderr,"++++++ warning: %s : ",t1); /*prints first part to the STDERR and +++....*/
  vfprintf(stderr, fmt, list); /*prints the second part and variables*/
  fprintf(stderr,"\n"); /*prints end of line*/
  fflush(stderr); /*flushes STDERR*/
  va_end(list);
}

void error(char* t1, register char* fmt, ...)
{
/*prints warning on the standard error and accepts parameters printout with std C formatting*/
/*Piotr Skowronski CERN*/
  va_list         list;

  warn_numb++; /*I think that warnings should be counted even if the user does not want to see them*/
  fflush(0); /*flushes all the buffers -> so the warning appears in a correct place*/

  va_start( list, fmt );

  fprintf(stderr,"++++++ Error: %s : ",t1); /*prints first part to the STDERR and +++....*/
  vfprintf(stderr, fmt, list); /*prints the second part and variables*/
  fprintf(stderr,"\n"); /*prints end of line*/
  fflush(stderr); /*flushes STDERR*/
  va_end(list);
}



void warningOld(char* t1, char* t2)
{
  if (get_option("warn")) 
  {
    printf("++++++ warning: %s %s\n",t1,t2);
    warn_numb++;
  }
}

void augmentfwarn() 
{
/*increases counter of the fortran warnings*/
   warn_numbf++;
}

