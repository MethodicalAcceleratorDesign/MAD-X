/* Production version of MAD-X, version number: see madxd.h */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#ifndef _WIN32
#include <sys/utsname.h>
#include <unistd.h>
#endif
#include <sys/timeb.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include "madxl.h"
#include "madx.h"
#include "madxreg.h"
#include "madxd.h"
#include "madxdict.h"

/* JMJ 7/11/2002 moved this here from c6t.c */ 
#include "c6t.h"

void madx()
{
  madx_start();
  madx_init();
  main_input(0);
  madx_finish();
}

#include "c6t.c"

#include "madxe.c"

#include "madxc.c"

#include "madxreg.c"

#include "sxf.c"

#include "makethin.c"

#include "matchc.c"

#include "madxu.c"

int act_special(int type, char* statement)
     /* acts on special commands (IF{..} etc.) */
{
  char* loc_buff = NULL;
  char* loc_w = NULL;
  int cnt_1, start_2, rs, re, level = pro->curr, ls = strlen(statement);
  int ret_val = 0;
  struct char_p_array* logic = new_char_p_array(1000);
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
  loc_buff = (char*) mymalloc("act_special", ls);  
  loc_w = (char*) mymalloc("act_special", ls);
  get_bracket_range(statement, '{', '}', &rs, &re);
  if (re < 0) fatal_error("missing '{' or '}' in statement:",statement); 
  cnt_1 = rs; start_2 = rs + 1;
  strcpy(loc_buff, statement); loc_buff[re] =  '\0';
  strncpy(l_dummy, statement, cnt_1); l_dummy[cnt_1] = '\0';
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
        pre_split(l_dummy, loc_w, 0);
        mysplit(loc_w, tmp_l_array);
        get_bracket_t_range(tmp_l_array->p, '(', ')', 0, tmp_l_array->curr,
                          &rs, &re);
        rs++;
        if ((logex = logic_expr(re-rs, &tmp_l_array->p[rs])) > 0)
          {
           pro->buffers[level]->flag = 1;
	   pro->curr++;
	   /* now loop over statements inside {...} */
           pro_input(&loc_buff[start_2]);
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
        pro_input(&loc_buff[start_2]);
        pro->curr--;
        pro->buffers[level]->flag = -1;
       }
     break;
    case 4: /* while */
     pre_split(l_dummy, loc_w, 0);
     mysplit(loc_w, logic);
     get_bracket_t_range(logic->p, '(', ')', 0, logic->curr,
                          &rs, &re);
     pro->curr++; rs++;
     while ((logex = logic_expr(re-rs, &logic->p[rs])) > 0)
       {
	/* now loop over statements inside {...} */
        pro_input(&loc_buff[start_2]);
       }
     pro->curr--;
     break;
    default:
      ret_val = -1;
    }
  if (loc_buff != NULL) free(loc_buff);
  if (loc_w != NULL) free(loc_w);
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

void adjust_beam()
     /* adjusts beam parameters to current beta, gamma, bcurrent, npart */
{
  struct name_list* nl = current_beam->par_names;
  double circ = one, freq0, alfa, beta, gamma, bcurrent, npart = 0;
  if (current_sequ != NULL && current_sequ->length != zero) 
      circ = current_sequ->length;
  beta = command_par_value("beta", current_beam);
  gamma = command_par_value("gamma", current_beam);
  alfa = one / (gamma * gamma);
  freq0 = (beta * clight) / (ten_p_6 * circ);
  if (nl->inform[name_list_pos("bcurrent", nl)] &&
    (bcurrent = command_par_value("bcurrent", current_beam)) > zero)
     npart = bcurrent / (freq0 * ten_p_6 * get_variable("qelect"));
  else if (nl->inform[name_list_pos("npart", nl)] &&
           (npart = command_par_value("npart", current_beam)) > zero)
     bcurrent = npart * freq0 * ten_p_6 * get_variable("qelect");
  store_comm_par_value("alfa", alfa, current_beam);
  store_comm_par_value("freq0", freq0, current_beam);
  store_comm_par_value("circ", circ, current_beam);
  store_comm_par_value("npart", npart, current_beam);
  store_comm_par_value("bcurrent", bcurrent, current_beam);
}

void adjust_probe(double delta_p)
     /* adjusts beam parameters to the current deltap */
{
  int j;
  double etas, slope, qs, fact, tmp, ds = oneturnmat[34];
  double alfa, beta, gamma, dtbyds, circ, deltat, freq0;
  double betas, gammas, et, sigt, sige;
  et = command_par_value("et", current_beam);
  sigt = command_par_value("sigt", current_beam);
  sige = command_par_value("sige", current_beam);
  beta = command_par_value("beta", current_beam);
  gamma = command_par_value("gamma", current_beam);
  circ = command_par_value("circ", current_beam);
  for (j = 0; j < 4; j++) ds += oneturnmat[4 + 6*j] * disp0[j];
  tmp = - beta * beta * ds / circ;
  freq0 = (clight * ten_m_6 * beta) / (circ * (one + tmp * delta_p));
  etas = beta * gamma * (one + delta_p);
  gammas = sqrt(one + etas * etas);
  betas = etas / gammas;
  tmp = - betas * betas * ds / circ;
  alfa = one / (gammas * gammas) + tmp;
  dtbyds = delta_p * tmp / betas;
  deltat = circ * dtbyds;
  store_comm_par_value("freq0", freq0, probe_beam);
  store_comm_par_value("alfa", alfa, probe_beam);
  store_comm_par_value("beta", betas, probe_beam);
  store_comm_par_value("gamma", gammas, probe_beam);
  store_comm_par_value("dtbyds", dtbyds, probe_beam);
  store_comm_par_value("deltap", delta_p, probe_beam);
  slope = -rfc_slope();
  qs = sqrt(fabs((tmp * slope) / (twopi * betas)));
  if (qs != zero)
    {
     fact = (tmp * circ) / (twopi * qs);
     if (et > zero)
        {
         sigt = sqrt(fabs(et * fact));
         sige = sqrt(fabs(et / fact));
  	}
     else if (sigt > zero)
  	{
         sige = sigt / fact;
         et = sige * sigt;
  	}
     else if (sige > zero)
       {
        sigt = sige * fact;
        et = sige * sigt;
       }
    }
  if (sigt < ten_m_15)
    {
     put_info("Zero value of SIGT", "replaced by 1.");
     sigt = one;
    }
  if (sige < ten_m_15)
    {
     put_info("Zero value of SIGE", "replaced by 1/1000.");
     sigt = ten_m_3;
    }
  store_comm_par_value("qs", qs, probe_beam);
  store_comm_par_value("et", et, probe_beam);
  store_comm_par_value("sigt", sigt, probe_beam);
  store_comm_par_value("sige", sige, probe_beam);
}

void adjust_rfc()
{
  /* adjusts rfc frequency to given harmon number */
  double freq0, harmon, freq;
  int i;
  struct element* el;
  freq0 = command_par_value("freq0", probe_beam);
  for (i = 0; i < current_sequ->cavities->curr; i++)
    {
     el = current_sequ->cavities->elem[i];
     if ((harmon = command_par_value("harmon", el->def)) > zero)
       {
	freq = freq0 * harmon;
        store_comm_par_value("freq", freq, el->def);
       }
    }
}

int advance_node()
     /* advances to next node in expanded sequence;
        returns 0 if end of range, else 1 */
{
  if (current_node == current_sequ->range_end)  return 0;
  current_node = current_node->next;
  return 1;
}

int advance_to_pos(char* table, int* t_pos)
     /* advances current_node to node at t_pos in table */
{
  struct table* t;
  int pos, cnt = 0, ret = 0;
  mycpy(c_dummy, table);
  if ((pos = name_list_pos(c_dummy, table_register->names)) > -1)
    {
     ret = 1;
     t = table_register->tables[pos];
     if (t->origin == 1)  return 0; /* table is read, has no node pointers */
     while (current_node)
       {
        if (current_node == t->p_nodes[*t_pos-1]) break;
        if ((current_node = current_node->next) 
             == current_sequ->ex_start) cnt++;
        if (cnt > 1) return 0;
       }
    }
  return ret;
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
     if ((p = command_par_string("apertype", node->p_elem->def)))
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

int attach_beam(struct sequence* sequ)
     /* attaches the beam belonging to the current sequence */
{
  if ((current_beam = find_command(sequ->name, beam_list)) == NULL)
    current_beam = find_command("default_beam", beam_list);
  return current_beam->beam_def;
}

void augment_count(char* table) /* increase table occ. by 1, fill missing */
{
  int pos;
  struct table* t;
  mycpy(c_dummy, table);
  if ((pos = name_list_pos(c_dummy, table_register->names)) > -1)
    t = table_register->tables[pos];
  else return;
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

char* buffer(char* string)  /* replaced by permbuff */
{
  return permbuff(string);
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

int char_from_table(char* table, char* name, int* row, char* val)
     /* OB 2.4.2002 */
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
  mycpy(c_dummy, table);
  if ((pos = name_list_pos(c_dummy, table_register->names)) > -1)
    t = table_register->tables[pos];
  else return -1;
  mycpy(c_dummy, name);
  if ((pos = name_list_pos(c_dummy, t->columns)) < 0) return -2;
  if (*row > t->curr)  return -3;
  strncpy(val,t->node_nm->p[*row-1],NAME_L);
  while (strlen(val)<=NAME_L) val[strlen(val)]=' ';
  return 0;
}

void check_table(char* string)
     /* replaces argument of "table" if any by a string variable */
{
  char *pt, *pl, *pr, *sv;
  int start = 0;
  while (strstr(&string[start], "table") != NULL)
    {
     strcpy(c_join, &string[start]);
     pt = strstr(c_join, "table");
     if ((pl = strchr(pt, '(')) == NULL) return;
     if ((pr = strchr(pl, ')')) == NULL) return;
     *pl = '\0';
     *pr = '\0';
     sv = make_string_variable(++pl);
     string[start] ='\0';
     strcat(string, c_join);
     strcat(string, " ( ");
     strcat(string, sv);
     strcat(string, " ) ");
     start = strlen(string);
     strcat(string, ++pr);
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

struct double_array* command_par_array(char* parameter, struct command* cmd)
     /* returns an updated command parameter array if found, else NULL */
{
  struct command_parameter* cp;
  struct double_array* arr = NULL;
  int i;
  if ((i = name_list_pos(parameter, cmd->par_names)) > -1)
    {
     cp = cmd->par->parameters[i];
     if (cp->type == 11 || cp->type == 12) 
       {
        arr = cp->double_array;
        if (cp->expr_list != NULL) update_vector(cp->expr_list, arr);
       }
    }
  return arr;
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

int command_par_vector(char* parameter, struct command* cmd, double* vector)
     /* returns the length of, and an updated command parameter vector 
        if found, else 0 */
{
  struct command_parameter* cp;
  int i;
  if ((i = name_list_pos(parameter, cmd->par_names)) > -1)
    {
     cp = cmd->par->parameters[i];
     if (cp->double_array != NULL)
       {
        if (cp->expr_list != NULL) 
            update_vector(cp->expr_list, cp->double_array);
	copy_double(cp->double_array->a, vector, cp->double_array->curr);
        return cp->double_array->curr;
       }
    }
  return 0;
}

void comment_to_table(char* table, char* comment, int* length)
     /* Saves the comment string at the current line.
        This comment is then printed in front of this line.
        Several calls to the same current line are possible. */
{
  int pos;
  struct table* t;
  mycpy(c_dummy, table);
  if ((pos = name_list_pos(c_dummy, table_register->names)) > -1)
    t = table_register->tables[pos];
  else return;
  strncpy(c_dummy, comment, *length); c_dummy[*length] = '\0';
  if (t->l_head[t->curr] == NULL)  
      t->l_head[t->curr] = new_char_p_array(2);
  else if (t->l_head[t->curr]->curr == t->l_head[t->curr]->max)
    grow_char_p_array(t->l_head[t->curr]);
  t->l_head[t->curr]->p[t->l_head[t->curr]->curr++] = tmpbuff(c_dummy);
}

void comm_para(char* name, int* n_int, int* n_double, int* n_string,
               int* int_array, double* double_array, char* strings, 
               int* string_lengths)
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
  int i, l, pos;
  struct command_parameter* cp;
  struct double_array* arr = NULL;
  *n_int = *n_double = *n_string = 0;
  mycpy(c_dummy, name);
  if (this_cmd != NULL && this_cmd->clone != NULL)
    {
     if ((pos = name_list_pos(c_dummy, this_cmd->clone->par_names)) > -1)
       {
        cp = this_cmd->clone->par->parameters[pos];
        switch (cp->type)
          {
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
            if (cp->string != NULL)
	      {
               *n_string = 1;
               l = *string_lengths = strlen(cp->string);
               strncpy(strings, cp->string, l);
	      }
            break;
          case 11:
          case 12:
            arr = cp->double_array;
            if (cp->expr_list != NULL) update_vector(cp->expr_list, arr);
            if (cp->type == 11)
	      {
	       for (i = 0; i < arr->curr; i++) int_array[i] = arr->a[i];
               *n_int = arr->curr;
	      }
            else
	      {
	       for (i = 0; i < arr->curr; i++) double_array[i] = arr->a[i];
               *n_double = arr->curr;
	      }
            break;
          case 13:
            for (i = 0; i < cp->m_string->curr; i++)
	      {
	       string_lengths[i] = l = strlen(cp->m_string->p[i]);
               strncpy(strings, cp->m_string->p[i], l);
               strings += l;
	      }
            *n_string = cp->m_string->curr;
          }
       }
    }
}

void complete_twiss_table(struct table* t)
     /* fills all items missing after "twiss" into twiss table */
{
  int i, j, mult, n;
  double el, val;
  struct node* c_node;
  char tmp[16];

  if (t == NULL) return;
  i = t->curr;
  c_node = current_node;
  mult = strcmp(c_node->base_name, "multipole") == 0 ? 1 : 0;
  t->s_cols[0][i] = tmpbuff(c_node->name);
  t->s_cols[1][i] = tmpbuff(c_node->base_name);
  t->s_cols[twiss_fill_end+1][i] = tmpbuff(c_node->p_elem->parent->name);
  for (j = twiss_opt_end+1; j<= twiss_fill_end; j++)
    {
     el = c_node->length;
     if (strcmp(twiss_table_cols[j], "l") == 0) val = el;
     else if(mult) 
       {
        val = mult_par(twiss_table_cols[j], c_node->p_elem);
        if (strstr(twiss_table_cols[j], "k0")) val *= c_node->dipole_bv;
        else val *= c_node->other_bv; 
       }
     else
       {
	strcpy(tmp, twiss_table_cols[j]);
        n = strlen(tmp) - 1;
        if (n > 1 && tmp[0] == 'k' && isdigit(tmp[1]) && tmp[n] == 'l') 
            tmp[n] = '\0'; /* suppress trailing l in k0l etc. */
        val = el_par_value(tmp, c_node->p_elem);
        if ((strstr(tmp, "k0") || strstr(tmp, "kick"))&& c_node->dipole_bv) 
            val *= c_node->dipole_bv;
        else if(c_node->other_bv) val *= c_node->other_bv; 
        if (strstr(tmp,"kick") == NULL && el != zero) val *= el;
       }
     t->d_cols[j][i] = val;
    }
}

char* compound(char* e_name, int occ)
     /* makes node name from element name and occurrence count */
{
  sprintf(c_dummy,"%s:%d", e_name, occ);
  return c_dummy;
}

struct expression* compound_expr(struct expression* e1, double v1,
                      char* oper, struct expression* e2, double v2)
     /* make one out of two expressions, using oper to connect them */
{
  char** toks = tmp_l_array->p;
  struct expression* expr = NULL;
  char tmp[30];
  int n;
  char lb[] = "(", rb[] = ")";
  if (e1 != NULL || e2 != NULL)
    {
     if (e1 != NULL)
       {
        if (e2 != NULL) 
          {
           toks[0] = lb; toks[1] = e1->string; toks[2] = rb;
           toks[3] = oper;
           toks[4] = lb; toks[5] = e2->string; toks[6] = rb;
          }
        else
          {
	   sprintf(tmp, "%e", v2);
           toks[0] = lb; toks[1] = e1->string; toks[2] = rb;
           toks[3] = oper;
           toks[4] = lb; toks[5] = tmp; toks[6] = rb;
          }
       }
     else
       {
	sprintf(tmp, "%e", v1);
        toks[0] = lb; toks[1] = tmp; toks[2] = rb;
        toks[3] = oper;
        toks[4] = lb; toks[5] = e2->string; toks[6] = rb;
       }
     join(toks, 7);
     pre_split(c_join, l_work, 0);
     n = mysplit(l_work, tmp_l_array);
     expr = make_expression(n, toks);
    }
  return expr;
}

void control(struct in_cmd* cmd)
     /* executes so-called "control" commands */
{
  char** toks = cmd->tok_list->p;
  int k = cmd->decl_start - 1;
  if      (strcmp(toks[k], "assign")      == 0) exec_assign(cmd);
  else if (strcmp(toks[k], "beam")        == 0) exec_beam(cmd, 0);
  else if (strcmp(toks[k], "beta0")       == 0) store_beta0(cmd);
  else if (strcmp(toks[k], "call")        == 0) exec_call(cmd);
  else if (strcmp(toks[k], "coguess")     == 0) exec_store_coguess(cmd);
  else if (strcmp(toks[k], "create")      == 0) exec_create_table(cmd);
  else if (strcmp(toks[k], "dumpsequ")    == 0) exec_dumpsequ(cmd);
  else if (strcmp(toks[k], "fill")        == 0) exec_fill_table(cmd);
  else if (strcmp(toks[k], "option")      == 0) exec_option();
  else if (strcmp(toks[k], "plot")        == 0) exec_plot(cmd);
  else if (strcmp(toks[k], "print")       == 0) exec_print(cmd);
  else if (strcmp(toks[k], "readtable")   == 0) read_table(cmd);
  else if (strcmp(toks[k], "resbeam")     == 0) exec_beam(cmd, 1);
  else if (strcmp(toks[k], "save")        == 0) exec_save(cmd);
  else if (strcmp(toks[k], "savebeta")    == 0) store_savebeta(cmd);
  else if (strcmp(toks[k], "select")      == 0) store_select(cmd);
  else if (strcmp(toks[k], "threader")    == 0) store_threader(cmd);
  else if (strcmp(toks[k], "use")         == 0) use_sequ(cmd);
  else if (strcmp(toks[k], "write")       == 0) exec_dump(cmd);
}

void deco_init()
     /* initializes Polish decoding */
{
  expr_chunks = new_name_list(2000);
  cat = new_int_array(MAX_ITEM);
  deco = new_int_array(MAX_ITEM);
  d_var = new_int_array(MAX_ITEM);
  oper = new_int_array(MAX_ITEM);
  func = new_int_array(MAX_ITEM);
  cat_doubles = new_double_array(MAX_ITEM);
  doubles = new_double_array(MAX_D_ITEM);
  twiss_deltas = new_double_array(MAX_ITEM);
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
  int j, k, ks, i = start, e_type, ival, end, e_end, tot_end, c_type = 0,
      val_type = 0, cnt = 0, con_flag;
  double val;
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
           if ((e_type = loc_expr(toks, number, start, &end)) == 0) return -i;
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
              if ((e_type = loc_expr(toks, number, start, &end)) == 0) 
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

double double_from_expr(char** toks, int s_start, int s_end)
     /* returns the value of an expression if valid, else INVALID */
{
  int end, nitem = s_end + 1;
  int type = loc_expr(toks, nitem, s_start, &end);
  if (type == 1) /* simple number */ 
     return simple_double(toks, s_start, end);
  else if (polish_expr(end + 1 - s_start, &toks[s_start]) == 0)
    return polish_value(deco);
  else return INVALID;
}

void double_to_table(char* table, char* name, double* val)
     /* puts val at current position in column with name "name".
        The table count is increased separately with "augment_count" */
{
  int pos;
  struct table* t;

  mycpy(c_dummy, table);
  if ((pos = name_list_pos(c_dummy, table_register->names)) > -1)
    t = table_register->tables[pos];
  else return;
  mycpy(c_dummy, name);
  if ((pos = name_list_pos(c_dummy, t->columns)) >= 0
      && t->columns->inform[pos] < 3) t->d_cols[pos][t->curr] = *val;
}

int double_from_table(char* table, char* name, int* row, double* val)
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
  mycpy(c_dummy, table);
  if ((pos = name_list_pos(c_dummy, table_register->names)) > -1)
    t = table_register->tables[pos];
  else return -1;
  mycpy(c_dummy, name);
  if ((pos = name_list_pos(c_dummy, t->columns)) < 0) return -2;
  if (*row > t->curr)  return -3;
  *val = t->d_cols[pos][*row-1];
  return 0;
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

void dynap_tables_create(struct in_cmd* cmd)
     /* creates the dynamic tables for DYNAP execution */
{
  int npart = stored_track_start->curr;

  struct table* t;
  t = make_table("tracksumm", "tracksumm", tracksumm_table_cols, 
  		 tracksumm_table_types, 2*stored_track_start->curr);
  add_to_table_list(t, table_register);
  t = make_table("dynap", "dynap", dynap_table_cols, dynap_table_types, 10);
  add_to_table_list(t, table_register);
  t = make_table("dynaptune", "dynaptune", dynaptune_table_cols, 
                 dynaptune_table_types, npart);
  add_to_table_list(t, table_register);
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

double el_par_value(char* par, struct element* el)
     /* returns an element parameter value */
{
  int k = 0, n;
  char tmp[8];
  double val = zero, angle = zero, tilt = zero, l, k0, k0s, vec[100];
  double fact = strcmp(el->base_type->name, "rbend") == 0 ? one : zero;
  int noang = 0, mult = strcmp(el->base_type->name, "multipole") == 0 ? 1 : 0;
  if (fact != zero || strcmp(el->base_type->name, "sbend") == 0)
    {
     if ((l = command_par_value("l", el->def)) == zero) 
        fatal_error("bend with zero length:",el->name);
     if ((angle = command_par_value("angle", el->def)) == zero)
       {
	noang = 1;
        k0 = command_par_value("k0", el->def);
        k0s = command_par_value("k0s", el->def);
        angle = l * sqrt(k0*k0+k0s*k0s);
        if (k0 < zero) 
	  {
	   angle = -angle; tilt = -atan2(k0s, fabs(k0));
	  }
        else tilt = atan2(k0s,k0);
       }
     else tilt = command_par_value("tilt", el->def);
     if (strcmp(par, "angle") == 0)  val = angle;
     else if (strcmp(par, "tilt") == 0)  val = tilt;
     else if (strcmp(par, "k0") == 0)
       {
        if ((k0 = command_par_value("k0", el->def)) == zero && noang == 0)
           k0 = cos(tilt) * angle / l;
        val = k0;
       }
     else if (strcmp(par, "k0s") == 0)
       {
        if ((k0s = command_par_value("k0s", el->def)) == zero && noang == 0)
           k0s = sin(tilt) * angle / l;
        val = k0s;
       }
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
  else if (strcmp(par, "rhoinv") == 0) val = zero;
  else if (strcmp(par, "blen") == 0) val = zero;
  else if (mult)
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

void element_name(char* name, int* l)
     /* returns current node element name in Fortran format */
     /* l is max. allowed length in name */
{
  int ename_l = strlen(current_node->p_elem->name);
  int i, ncp = ename_l < *l ? ename_l : *l;
  int nbl = *l - ncp;
  for (i = 0; i < ncp; i++) name[i] = current_node->p_elem->name[i];
  for (i = 0; i < nbl; i++) name[ncp+i] = ' ';
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

void enter_elm_reference(struct in_cmd* cmd, struct element* el, int flag)
     /* enters an element in a sequence */
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int i, pos, k = 1;
  double at;
  if (strcmp(el->base_type->name, "rfcavity") == 0 &&
      find_element(el->name, current_sequ->cavities) == NULL)
    add_to_el_list(el, 0, current_sequ->cavities, 0);
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
     sprintf(c_dummy, "%s$end", current_sequ->name);
     el = make_element(c_dummy, "marker", clone, 0);
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
       current_sequ = sequences->sequs[pos];
     else
       {
        current_sequ = new_sequence(toks[aux_pos], k);
        add_to_sequ_list(current_sequ, sequences);
       }
     cmd->clone = clone_command(cmd->cmd_def);
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
           occ_list = new_name_list(10000);  /* for occurrence count */
     else occ_list->curr = 0;
     if (current_sequ->cavities != NULL)  current_sequ->cavities->curr = 0;
     else current_sequ->cavities = new_el_list(100);
     pos = name_list_pos("marker", defined_commands->list);
     clone = clone_command(defined_commands->commands[pos]);
     sprintf(c_dummy, "%s$start", current_sequ->name);
     el = make_element(c_dummy, "marker", clone, 0);
     make_elem_node(el, 1);
     current_sequ->start = current_node;
     current_sequ->share = aux_pos;
    }
}

void enter_variable(struct in_cmd* cmd) /* stores variable contained in cmd */
{
  struct variable* var;
  struct expression* expr = NULL;
  int k, end, type, name_pos, start = cmd->decl_start, val_type;
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
  current_beam->par->parameters[pos]->string = permbuff(name);
  current_beam->beam_def = 1;
  if (flag == 0) update_beam(cmd->clone);
  else if (flag == 1)
    {
     set_defaults("beam");
     add_to_command_list(name, current_beam, beam_list, 0);
    }
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
     else if (strcmp(cmd_name, "resplot") == 0)
       { 
        plot_options = delete_command(plot_options);
	set_defaults("setplot");
       }
     else if (strcmp(cmd_name, "value") == 0)
       {
	print_value(p);
       }
     else if (strcmp(cmd_name, "system") == 0)
	ret = system(noquote(toks[p->decl_start]));
     else if (strcmp(cmd_name, "title") == 0)
	title = permbuff(noquote(toks[p->decl_start]));
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
	else if (strcmp(cmd_name, "setplot") == 0
            && plot_options != NULL) 
	  {
           p->clone = plot_options; p->clone_flag = 1;
	  }
        else p->clone = clone_command(p->cmd_def);
	scan_in_cmd(p); /* match input command with clone + fill */
        current_command = p->clone;
        if (strcmp(p->cmd_def->module, "control") == 0) control(p);
        else if (strcmp(p->cmd_def->module, "c6t") == 0) conv_sixtrack(p);
        else if (strcmp(p->cmd_def->module, "edit") == 0) seq_edit_main(p);
        else if (strcmp(p->cmd_def->module, "ibs") == 0) 
	  {
	    current_ibs = p->clone;
	    pro_ibs(p);
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
        /* else if (strcmp(p->cmd_def->module, "ptc") == 0) ttwm_(); */
      }
    }
}

void exec_create_table(struct in_cmd* cmd)
     /* makes a user defined table */
{
  struct table* t;
  int* t_types;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  struct char_p_array* m;
  char** t_c;
  int j, pos = name_list_pos("table", nl);
  char* name = NULL;
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
  t_types = malloc(m->curr*sizeof(int));
  t_c = malloc((m->curr+1)*sizeof(char*));
  for (j = 0; j < m->curr; j++) 
    {
     t_types[j] = 2; /* type double */
     t_c[j] = permbuff(m->p[j]);
    }
  t_c[m->curr] = blank;
  t = make_table(name, "user", t_c, t_types, USER_TABLE_LENGTH);
  t->org_cols = 0;  /* all entries are "added" */
  add_to_table_list(t, table_register);
  free(t_c); free(t_types);
}

void exec_store_coguess(struct in_cmd* cmd)
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

void exec_dump(struct in_cmd* cmd)
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
            || *f == '0') strcpy(filename, name);
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

void exec_fill_table(struct in_cmd* cmd)
     /* adds variables to a table */
{
  struct table* t;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int pos = name_list_pos("table", nl);
  char* name = NULL;
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
     t = table_register->tables[pos];
     add_vars_to_table(t);
     if (++t->curr == t->max) grow_table(t);
    }
  else warning("table not found: ", "ignored");
     return;
}

void exec_macro(struct in_cmd* cmd, int pos)
     /* executes a macro */
{
  int i, rs, re, any = 0, level = pro->curr;
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
     for (i = 0; i < any; i++) 
       {
        my_repl(macro_list->macros[pos]->formal->p[i], toks[rs+i],
        pro->buffers[level]->c_a->c, l_work);
        strcpy(pro->buffers[level]->c_a->c, l_work);
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

void exec_plot(struct in_cmd* cmd)
{
  int ierr, nt = strcmp(title,"no-title") == 0 ? 1 : 0;
  char* pt = title;

 /* <JMJ 7/11/2002> The following ifndef exclusion is a quick fix so that 
     the WIN32 version 
     does not try to do X11 graphics. However this has the consequence that 
	 the program will not make Postscript files.  HG needs to separate these things.
  </JMJ 7/11/2002> */

#ifndef _WIN32
  if (nt && current_sequ != NULL) title = current_sequ->name;
  pesopt_(&ierr);
  if (ierr == 0) pefill_(&ierr);
  if (ierr == 0) 
    {
     pemima_();
     plotit_(&plots_made);
     plots_made = 1;
    }
  if (nt) title = pt;
#endif
}

void exec_print(struct in_cmd* cmd)
     /* prints text from "print" command to current output unit */
{
  struct command_parameter_list* pl = cmd->clone->par;
  struct name_list* nl = cmd->clone->par_names;
  int pos = name_list_pos("text", nl);
  if (nl->inform[pos]) fprintf(prt_file,"%s\n", pl->parameters[pos]->string);
}

void exec_save(struct in_cmd* cmd)
     /* save a sequence with all necessary parameters and sub-sequences */
{
  int i, n = 0, pos, prev = 0, beam_save = log_val("beam", cmd->clone), 
      mad8 = log_val("mad8", cmd->clone), all_sequ = 0;
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
  else                varl = variable_list; /* write all variables */
  for (pos = 0; pos < sqo->curr; pos++)
    {
     sequ = sqo->sequs[pos];
     fill_sequ_list(sequ, sql);
     if (beam_save) 
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
                 add_to_el_list(el, 0, ell, 0);
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
     write_vars_8(varl, save_select, out_file);
     write_elems_8(ell, save_select, out_file);
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
     write_vars(varl, save_select, out_file);
     write_elems(ell, save_select, out_file);
     write_sequs(sql, save_select, out_file);
    }
  fclose(out_file);
  if (sqo != sequences) sqo = delete_sequence_list(sqo);
  sql = delete_sequence_list(sql);
  ell = delete_el_list(ell);
  varl = delete_var_list(varl);
  current_sequ = NULL;
}

void exec_savebeta()
     /* stores twiss values in a beta0 structure */
{
  struct name_list* nl;
  struct command_parameter_list* pl;
  struct node* nodes[2];
  struct command* beta0;
  char* label;
  int i, pos;
  for (i = 0; i < savebeta_list->curr; i++)
    {
     nl = savebeta_list->commands[i]->par_names;
     pl = savebeta_list->commands[i]->par;
     pos = name_list_pos("label", nl);
     label = pl->parameters[pos]->string;
     if (find_command(label, beta0_list) == NULL) /* fill only once */
       {
        pos = name_list_pos("sequence", nl);
        if (nl->inform[pos] == 0 
         || strcmp(pl->parameters[pos]->string, current_sequ->name) == 0)
          {
           pos = name_list_pos("place", nl);
           if (get_ex_range(pl->parameters[pos]->string, current_sequ, nodes))
	     {
	      pos = name_list_pos("beta0", defined_commands->list);
              beta0 = clone_command(defined_commands->commands[pos]);
              fill_beta0(beta0, nodes[0]);
              add_to_command_list(label, beta0, beta0_list, 0);
	     }
	  }
       }
    }
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
           else fprintf(prt_file, "%s = %-18.10g ;\n", toks[i], var->value);
          }
        else fprintf(prt_file, "%s not found\n", toks[i]);
       }
    }

}

void expand_line(struct char_p_array* l_buff)
     /* expands a beam line, applies rep. count and inversion */
{
  /* first get all bracket pairs with their level; keep max. level */
  int add, i, j, k, n, number, dummy, rep, pos;
  int level = 0, l_max = 0, b_cnt = 0;
  char* p;
  struct int_array* lbpos = new_int_array(l_buff->curr);
  struct int_array* rbpos = new_int_array(l_buff->curr);
  struct int_array* b_level = new_int_array(l_buff->curr);

  for (i = 0; i < l_buff->curr; i++) 
    {
     if (*l_buff->p[i] == '(')
       {
        lbpos->i[b_cnt] = i;
        b_level->i[b_cnt++] = level++;
        if (level > l_max) l_max = level;
       }
     else if (*l_buff->p[i] == ')')  level--;
    }
  l_max--;
  for (i = 0; i < b_cnt; i++)
     get_bracket_t_range(l_buff->p, '(', ')', lbpos->i[i], 
                         l_buff->curr-1, &dummy, &rbpos->i[i]);
  lbpos->curr = rbpos->curr = b_level->curr = b_cnt;
  /* now loop over level from highest down to zero, expand '*' in each pair */
  for (level = l_max; level >=0; level--)
    {
     for (i = 0; i < b_cnt; i++)
       {
	if (b_level->i[i] == level && (pos = lbpos->i[i]) > 1)
	  {
	   if (*l_buff->p[pos-1] == '*')
	     {
	      sscanf(l_buff->p[pos-2], "%d", &rep);
              add = rep - 1;
              number = rbpos->i[i] - pos - 1; /* inside bracket */
              n = number * add; /* extra tokens */
              while (l_buff->curr + n >= l_buff->max) 
                    grow_char_p_array(l_buff);
              for (j = l_buff->curr; j > pos + number; j--) /* shift upwards */
		l_buff->p[j+n] = l_buff->p[j];
              l_buff->curr += n;
              for (k = 1; k <= add; k++)
		{
		 for (j = pos+1; j <= pos+number; j++) 
		   l_buff->p[j+k*number] = l_buff->p[j];
		}
              for (j = 0; j < b_cnt; j++)  /* reset bracket pointers */
		{
		 if (lbpos->i[j] > pos + number) lbpos->i[j] += n;
		 if (rbpos->i[j] > pos + number) rbpos->i[j] += n;
		}
              l_buff->p[pos-1] = l_buff->p[pos-2] = blank;
	     }
	  }
       }
    }
  /* loop over buffer, expand simple element repetition */
  for (pos = 2; pos < l_buff->curr; pos++)
    {
     if (*l_buff->p[pos] == '*')
       {
	sscanf(l_buff->p[pos-2], "%d", &rep);
        n = add = rep - 1;
        while (l_buff->curr + n >= l_buff->max) grow_char_p_array(l_buff);
        for (j = l_buff->curr; j > pos + 1; j--) /* shift upwards */
	  l_buff->p[j+n] = l_buff->p[j];
        l_buff->curr += n;
        for (k = 1; k <= add; k++)
	  {
	   j = pos+1; 
	   l_buff->p[j+k] = l_buff->p[j];
	  }
        for (j = 0; j < b_cnt; j++)  /* reset bracket pointers */
	  {
	   if (lbpos->i[j] > pos + 1) lbpos->i[j] += n;
	   if (rbpos->i[j] > pos + 1) rbpos->i[j] += n;
	  }
        l_buff->p[pos-1] = l_buff->p[pos-2] = blank;
       }
    }
  /* now loop over level from highest down to zero, invert if '-' */
  for (level = l_max; level >= 0; level--)
    {
     for (i = 0; i < b_cnt; i++)
       {
	pos = lbpos->i[i];
	if (b_level->i[i] == level)
	  {
	   p = blank;
	   for (j = pos - 1; j > 0; j--)
	     {
	      p = l_buff->p[j];
              if (*p != ' ')  break;
	     }
           if (*p == '-')
	     {
              number = rbpos->i[i] - pos - 1; 
              n = number / 2;
              for (j = 0; j < n; j++)
		{
		 p = l_buff->p[pos+j+1];
                 l_buff->p[pos+j+1] = l_buff->p[pos+number-j];
                 l_buff->p[pos+number-j] = p;
		}
	     }
	  }
       }
    }
  /* finally remove all non-alpha tokens */
  n = 0;
  for (i = 0; i < l_buff->curr; i++) 
    if (isalpha(*l_buff->p[i]))  l_buff->p[n++] = l_buff->p[i];
  l_buff->curr = n;
  lbpos = delete_int_array(lbpos);
  rbpos = delete_int_array(rbpos);
  b_level = delete_int_array(b_level);
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

void expand_curr_sequ(int flag)
     /* expands the current sequence, i.e. flattens it, inserts drifts etc. */
{
  struct node* c_node;
  int j;
  if (current_sequ->ex_start != NULL)
    {
     current_sequ->ex_nodes = delete_node_list(current_sequ->ex_nodes);
     current_sequ->ex_start = delete_node_ring(current_sequ->ex_start);
     current_sequ->orbits = delete_vector_list(current_sequ->orbits);
    }
  if (current_sequ->ex_start == NULL)
    {
     use_count++;
     if (occ_list == NULL) 
        occ_list = new_name_list(10000);  /* for occurrence count */
     else occ_list->curr = 0;
     make_occ_list(current_sequ);
     all_node_pos(current_sequ);
     current_sequ->ex_nodes = new_node_list(2*current_sequ->nodes->curr);
     expand_sequence(current_sequ, flag);
     current_sequ->n_nodes = 
       add_drifts(current_sequ->ex_start, current_sequ->ex_end);
     if (current_sequ->all_nodes != NULL) free(current_sequ->all_nodes);
     current_sequ->all_nodes 
        = (struct node**) malloc(current_sequ->n_nodes * sizeof(struct node*));
     c_node = current_sequ->ex_start;
     for (j = 0; j < current_sequ->n_nodes; j++)
       {
	current_sequ->all_nodes[j] = c_node;
        c_node = c_node->next;
       }
    }
  set_node_bv(current_sequ); /* set bv factors for all nodes */
  if (current_range) set_range(current_range, current_sequ);
  else
    {
     current_sequ->range_start = current_sequ->ex_start;
     current_sequ->range_end = current_sequ->ex_end;
    }
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

void fill_beta0(struct command* beta0, struct node* node)
{
  /* makes uses of fact that beta0 and twiss_table have the same variables
     at the same place (+2)  up to energy inclusive */
  struct command_parameter_list* pl = beta0->par;
  struct name_list* nl = beta0->par_names;
  int i = -1, pos;
  if (twiss_table == NULL) return;
  for (pos = 0; pos < twiss_table->curr; pos++)
    {
     if (twiss_table->p_nodes[pos] == node)  break;
    }
  if (pos < twiss_table->curr)
    {
     do
       {
	i++;
	pl->parameters[i]->double_value = twiss_table->d_cols[i+3][pos];
       }
     while (strcmp(nl->names[i], "energy") != 0);
    }
}

void fill_constraint_list(int type /* 1 node, 2 global */, 
                          struct command* cd, struct constraint_list* cl)
{
  struct command_parameter_list* pl = cd->par;
  struct name_list* nl = cd->par_names;
  struct constraint* l_cons;
  int j;
  for (j = 0; j < pl->curr; j++)
    {
     if (nl->inform[j] && pl->parameters[j]->type == 4) 
	{
	 l_cons = make_constraint(type, pl->parameters[j]);
         add_to_constraint_list(l_cons, cl);
        }
    }
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
              add_to_el_list(el, 0, ell, 0);
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

void fill_orbit_table(struct table* t_out, struct table* t_in)
     /* fills a table with orbit values at monitor positions */
{
  int i, j, pos;
  t_out->curr = 0;
  for (i = 0; i < t_in->curr; i++)
    {
     if (strstr(t_in->s_cols[1][i], "monitor"))
       {
	for (j = 0; j < t_out->num_cols; j++)
	  {
	   if ((pos = name_list_pos(t_out->columns->names[j], 
                t_in->columns)) > -1)
	     {
	      if (t_out->columns->inform[j] < 3) 
		t_out->d_cols[j][t_out->curr] = t_in->d_cols[pos][i];
              else t_out->s_cols[j][t_out->curr] 
                   = tmpbuff(t_in->s_cols[pos][i]);
	     }
           else
	     {
	      if (t_out->columns->inform[j] < 3) 
		t_out->d_cols[j][t_out->curr] = zero;
              else t_out->s_cols[j][t_out->curr] = tmpbuff(blank);
	     }
	  }
        t_out->curr++;
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

void fill_sequ_list(struct sequence* sequ, struct sequence_list* sql)
     /* puts all sequences depending on sequ into a sequ_list, recursively */
{
  struct node* c_node;
  add_to_sequ_list(sequ, sql);
  c_node = sequ->start;
  while(c_node != NULL)
    {
     if (c_node->p_sequ != NULL) fill_sequ_list(c_node->p_sequ, sql);
     if (c_node == sequ->end) break;
     c_node = c_node->next;
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

void fill_twiss_header(struct table* t)
     /* puts beam parameters etc. at start of twiss table */
{
  int i, pos, h_length = 33; /* change adding header lines ! */
  double dtmp;
  struct table* s;
  char tmp[16];

  if (t == NULL) return;
  /* ATTENTION: if you add header lines, augment h_length accordingly */
  if (t->header == NULL)  t->header = new_char_p_array(h_length);
  strcpy(tmp, t->org_sequ->name);
  sprintf(c_dummy, "@ SEQUENCE         %%%02ds \"%s\"", strlen(tmp), 
          stoupper(tmp));
  t->header->p[t->header->curr++] = tmpbuff(c_dummy);
  i = get_string("beam", "particle", tmp);
  sprintf(c_dummy, "@ PARTICLE         %%%02ds \"%s\"", i, stoupper(tmp));
  t->header->p[t->header->curr++] = tmpbuff(c_dummy);
  dtmp = get_value("beam", "mass");
  sprintf(c_dummy, "@ MASS             %%le  %22.12g", dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dummy);
  dtmp = get_value("beam", "charge");
  sprintf(c_dummy, "@ CHARGE           %%le  %22.12g", dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dummy);
  dtmp = get_value("beam", "energy");
  sprintf(c_dummy, "@ ENERGY           %%le  %22.12g", dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dummy);
  dtmp = get_value("beam", "pc");
  sprintf(c_dummy, "@ PC               %%le  %22.12g", dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dummy);
  dtmp = get_value("beam", "gamma");
  sprintf(c_dummy, "@ GAMMA            %%le  %22.12g", dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dummy);
  dtmp = get_value("beam", "kbunch");
  sprintf(c_dummy, "@ KBUNCH           %%le  %22.12g", dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dummy);
  dtmp = get_value("beam", "bcurrent");
  sprintf(c_dummy, "@ BCURRENT         %%le  %22.12g", dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dummy);
  dtmp = get_value("beam", "sige");
  sprintf(c_dummy, "@ SIGE             %%le  %22.12g", dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dummy);
  dtmp = get_value("beam", "sigt");
  sprintf(c_dummy, "@ SIGT             %%le  %22.12g", dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dummy);
  dtmp = get_value("beam", "npart");
  sprintf(c_dummy, "@ NPART            %%le  %22.12g", dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dummy);
  dtmp = get_value("beam", "ex");
  sprintf(c_dummy, "@ EX               %%le  %22.12g", dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dummy);
  dtmp = get_value("beam", "ey");
  sprintf(c_dummy, "@ EY               %%le  %22.12g", dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dummy);
  dtmp = get_value("beam", "et");
  sprintf(c_dummy, "@ ET               %%le  %22.12g", dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dummy);
  if ((pos = name_list_pos("summ", table_register->names)) > -1)
    {
     s = table_register->tables[pos];
     pos = name_list_pos("length", s->columns);
     dtmp = s->d_cols[pos][0];
     sprintf(c_dummy, "@ LENGTH           %%le  %22.12g", dtmp);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     pos = name_list_pos("alfa", s->columns);
     dtmp = s->d_cols[pos][0];
     sprintf(c_dummy, "@ ALFA             %%le  %22.12g", dtmp);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     pos = name_list_pos("gammatr", s->columns);
     dtmp = s->d_cols[pos][0];
     sprintf(c_dummy, "@ GAMMATR          %%le  %22.12g", dtmp);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     pos = name_list_pos("q1", s->columns);
     dtmp = s->d_cols[pos][0];
     sprintf(c_dummy, "@ Q1               %%le  %22.12g", dtmp);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     pos = name_list_pos("q2", s->columns);
     dtmp = s->d_cols[pos][0];
     sprintf(c_dummy, "@ Q2               %%le  %22.12g", dtmp);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     pos = name_list_pos("dq1", s->columns);
     dtmp = s->d_cols[pos][0];
     sprintf(c_dummy, "@ DQ1              %%le  %22.12g", dtmp);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     pos = name_list_pos("dq2", s->columns);
     dtmp = s->d_cols[pos][0];
     sprintf(c_dummy, "@ DQ2              %%le  %22.12g", dtmp);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     pos = name_list_pos("dxmax", s->columns);
     dtmp = s->d_cols[pos][0];
     sprintf(c_dummy, "@ DXMAX            %%le  %22.12g", dtmp);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     pos = name_list_pos("dymax", s->columns);
     dtmp = s->d_cols[pos][0];
     sprintf(c_dummy, "@ DYMAX            %%le  %22.12g", dtmp);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     pos = name_list_pos("xcomax", s->columns);
     dtmp = s->d_cols[pos][0];
     sprintf(c_dummy, "@ XCMAX            %%le  %22.12g", dtmp);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     pos = name_list_pos("ycomax", s->columns);
     dtmp = s->d_cols[pos][0];
     sprintf(c_dummy, "@ YCMAX            %%le  %22.12g", dtmp);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     pos = name_list_pos("betxmax", s->columns);
     dtmp = s->d_cols[pos][0];
     sprintf(c_dummy, "@ BETXMAX          %%le  %22.12g", dtmp);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     pos = name_list_pos("betymax", s->columns);
     dtmp = s->d_cols[pos][0];
     sprintf(c_dummy, "@ BETYMAX          %%le  %22.12g", dtmp);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     pos = name_list_pos("xcorms", s->columns);
     dtmp = s->d_cols[pos][0];
     sprintf(c_dummy, "@ XCRMS            %%le  %22.12g", dtmp);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     pos = name_list_pos("ycorms", s->columns);
     dtmp = s->d_cols[pos][0];
     sprintf(c_dummy, "@ YCRMS            %%le  %22.12g", dtmp);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     pos = name_list_pos("dxrms", s->columns);
     dtmp = s->d_cols[pos][0];
     sprintf(c_dummy, "@ DXRMS            %%le  %22.12g", dtmp);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     pos = name_list_pos("dyrms", s->columns);
     dtmp = s->d_cols[pos][0];
     sprintf(c_dummy, "@ DYRMS            %%le  %22.12g", dtmp);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     pos = name_list_pos("deltap", s->columns);
     dtmp = s->d_cols[pos][0];
     sprintf(c_dummy, "@ DELTAP           %%le  %22.12g", dtmp);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
    }
}

struct command* find_command(char* name, struct command_list* cl)
{
  int pos;
  if ((pos = name_list_pos(name, cl->list)) < 0) 
      return NULL;
  return cl->commands[pos];
}

struct command_list* find_command_list(char* name, 
                                       struct command_list_list* sl)
{
  int pos;
  if ((pos = name_list_pos(name, sl->list)) < 0) 
      return NULL;
  return sl->command_lists[pos];
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

void madx_finish()
     /* write the termination message */
{
  if (final_message == 0)
    {
     final_message = 1;
     if (plots_made) gxterm_();
     if (get_option("trace")) time_stamp("end");
     printf("\n  ++++++++++++++++++++++++++++++++\n");
     printf("  + %s finished normally +\n", myversion);
     printf("  ++++++++++++++++++++++++++++++++\n");
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
  char_buff = new_char_array_list(100);
  char_buff->ca[char_buff->curr++] = new_char_array(CHAR_BUFF_SIZE);
  drift_list = new_el_list(1000);
  variable_list = new_var_list(2000);
  comm_constraints = new_constraint_list(10);
  beam_list = new_command_list(10);
  stored_track_start = new_command_list(100);
  table_select = new_command_list_list(10);
  defined_commands = new_command_list(100);
  stored_commands = new_command_list(500);
  line_list = new_macro_list(100);
  macro_list = new_macro_list(100);
  base_type_list = new_el_list(60);
  element_list = new_el_list(20000);
  buffered_cmds = new_in_cmd_list(10000);
  sequences = new_sequence_list(20);
  match_sequs = new_sequence_list(2);
  selected_ranges = new_node_list(10000);
  selected_elements = new_el_list(10000);
  tmp_p_array = new_char_p_array(1000);
  tmp_l_array = new_char_p_array(1000);
  line_buffer = new_char_p_array(1000);
  sxf_list = new_name_list(50);
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
  set_defaults("setplot");
  set_defaults("threader");
  table_register = new_table_list(10);
  beta0_list = new_command_list(10);
  savebeta_list = new_command_list(10);
  seqedit_select = new_command_list(10); /* for "select seqedit" commands */
  error_select = new_command_list(10); /* for "select error" commands */
  save_select = new_command_list(10); /* for "select save" commands */
  slice_select = new_command_list(10); /* for "select makethin" commands */
  sector_select = new_command_list(10); /* for "select sectormap" commands */
  s_range = new_int_array(10);
  e_range = new_int_array(10);
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
  printf("\n  ++++++++++++++++++++++++++++++++\n");
  printf("  + %s %02d/%02d/%02d %02d.%02d.%02d +\n", myversion,
         tm->tm_mday, tm->tm_mon+1, tm->tm_year%100,
         tm->tm_hour, tm->tm_min, tm->tm_sec);
  printf("  ++++++++++++++++++++++++++++++++\n");
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

void get_disp0(double* disp)
{
  copy_double(disp0, disp, 6);
}

char* get_new_name()
     /* makes a new internal element or variable name */
{
  char name[NAME_L] = "__";
  sprintf(&name[2], "%d", new_name_count++);
  strcat(name, "__");
  return permbuff(name);
}

int get_select_ex_ranges(struct sequence* sequ, struct command_list* select,
                         struct node_list* s_ranges)
     /* makes a list of nodes of an expanded sequence that pass the range
        selection */
{
  /*returns 0 if invalid sequence pointer
            1 if nodes in s_ranges (including 0) */
  struct name_list* nl;
  struct command* cd;
  struct command_parameter_list* pl;
  char* name;
  int full = 0, i, k, pos;
  struct node* c_node;
  struct node* nodes[2];
  if (sequ == NULL) return 0; 
  for (i = 0; i < select->curr; i++)
    {
     cd = select->commands[i];
     nl = cd->par_names;
     pl = cd->par;
     pos = name_list_pos("full", nl);
     if ((pos = name_list_pos("full", nl)) > -1 && nl->inform[pos] 
	 && command_par_value("full", cd) != zero) full = 1;
     if (full == 0 && (pos = name_list_pos("range", nl)) > -1 
                   && nl->inform[pos])
       {
        name = pl->parameters[pos]->string;
        if ((k = get_ex_range(name, sequ, nodes)) == 0) return 0;
       }
     else
       {
	if ((nodes[0] = sequ->ex_start) == NULL ||
	     (nodes[1] = sequ->ex_end) == NULL) return 0;
       }
     c_node = nodes[0];
     while (c_node != NULL)
       {
	if (full != 0 || pass_select(c_node->p_elem->name, cd) != 0)
              add_to_node_list(c_node, 0, s_ranges);
        if (c_node == nodes[1]) break;
        c_node = c_node->next;
       }
     if (full != 0) break;
    }
  return 1;
}

int get_select_ranges(struct sequence* sequ, struct command_list* select,
                      struct node_list* s_ranges)
     /* makes a list of nodes of a sequence that pass the range selection */
{
  struct name_list* nl;
  struct command_parameter_list* pl;
  char* name;
  int i, k, pos;
  struct node* c_node;
  struct node* nodes[2]; 
  for (i = 0; i < select->curr; i++)
    {
     nl = select->commands[i]->par_names;
     pl = select->commands[i]->par;
     pos = name_list_pos("range", nl);
     if (pos > -1 && nl->inform[pos])  /* parameter has been read */
       {
        name = pl->parameters[pos]->string;
        if ((k = get_range(name, sequ, nodes)) > 0)
	  {
	   c_node = nodes[0];
           while (c_node != NULL)
	     {
	      add_to_node_list(c_node, 0, s_ranges);
              if (c_node == nodes[1]) break;
              c_node = c_node->next;
	     }
	  }
       }
    }
  return s_ranges->curr;
}

void get_select_t_ranges(struct command_list* select, struct table* t)
     /* makes a list of table rows that pass the range selection */
{
  int rows[2];
  struct name_list* nl;
  struct command_parameter_list* pl;
  int i, pos, any = 0;
  s_range->curr = 0; e_range->curr = 0;
  while (s_range->max < select->curr) grow_int_array(s_range);
  while (e_range->max < select->curr) grow_int_array(e_range);
  for (i = 0; i < select->curr; i++)
    {
     nl = select->commands[i]->par_names;
     pl = select->commands[i]->par;
     pos = name_list_pos("range", nl);
     if (pos > -1 && nl->inform[pos])  /* parameter has been read */
       {
	if (get_table_range(pl->parameters[pos]->string, t, rows))
	  {
	   any = 1;
	   if (rows[0] <= rows[1])
	     {
              s_range->i[s_range->curr++] = rows[0];
              e_range->i[e_range->curr++] = rows[1];
	     }
	   else
	     {
              s_range->i[s_range->curr++] = 0;
              e_range->i[e_range->curr++] = t->curr - 1;
	     }
	  }
       }
    }
  if (any == 0)
    {
     for (i = 0; i < select->curr; i++)
       {
        s_range->i[s_range->curr++] = 0;  
        e_range->i[e_range->curr++] = t->curr - 1;
       }
    } 
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
  free(p);
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
  add_to_el_list(p, 1, drift_list, 0);
  return p;
}

int get_node_count(struct node* node)
     /* finds the count of a node in the current expanded sequence */
{
  int cnt = 0;
  current_node = current_sequ->ex_start;
  while (current_node != NULL)
    {
     if (current_node == node) return cnt;
     cnt++;
     if (current_node == current_sequ->ex_end) break;
     current_node = current_node->next;
    }
  return -1;
}

double get_node_pos(struct node* node, struct sequence* sequ) /*recursive */
     /* returns node position from declaration for expansion */
{
  double fact = 0.5 * sequ->ref_flag; /* element half-length offset */
  double pos, from = 0;
  if (loop_cnt++ == MAX_LOOP) 
    {
     sprintf(c_dummy, "%s   occurrence: %d", node->p_elem->name, 
             node->occ_cnt);
     fatal_error("circular call in position of", c_dummy);
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
  int i, k;
  mycpy(c_dummy, str);
  if (options != NULL 
      && (i = name_list_pos(c_dummy, options->par_names)) > -1) 
     return (k = options->par->parameters[i]->double_value);
  else if (strcmp(c_dummy, "warn") == 0) return init_warn;
  else return 0;
}

void get_node_vector(char* par, int* length, double* vector)  
/* returns vector for parameter par of current element */
{
  char lpar[NAME_L];
  mycpy(lpar, par);
  if (strcmp(lpar, "orbit0") == 0) copy_double(orbit0, vector, 6);
  else if (strcmp(lpar, "obs_orbit") == 0)
    {
     if (current_node->obs_orbit)
       {
	*length = current_node->obs_orbit->curr;
        copy_double(current_node->obs_orbit->a, vector, *length);
       }
     else *length = 0;
    }
  else if (strcmp(lpar, "orbit_ref") == 0)
    {
     if (current_node->orbit_ref)
       {
        *length = current_node->orbit_ref->curr;
        copy_double(current_node->orbit_ref->a, vector, *length);
       }
    }
  else *length = element_vector(current_node->p_elem, lpar, vector);
}

int get_ex_range(char* range, struct sequence* sequ, struct node** nodes)
     /* returns start and end node (nodes[0] and nodes[1])
        of a range in the full expanded sequence */
{
  int i, n, pos;
  char* c[2];
  char tmp[NAME_L];
  if (sequ == NULL) return 0;
  strcpy(c_dummy, range); stolower(c_dummy);
  c[0] = strtok(c_dummy, "/");
  if ((c[1] = strtok(NULL,"/")) == NULL) /* only one element given */
    n = 1;
  else n = 2;
  for (i = 0; i < n; i++)
    {
     if (*c[i] == '#')
       {
	if (strncmp(c[i], "#s", 2) == 0) nodes[i] = sequ->ex_start;
	else if (strncmp(c[i], "#e", 2) == 0) nodes[i] = sequ->ex_end;
        else
	  {
	   warning("illegal expand range ignored:", range);
           return 0;
	  }
       }
     else
       {
	strcpy(tmp, c[i]);
        if (square_to_colon(tmp) == 0)
	  {
	   warning("illegal expand range ignored:", range);
           return 0;
	  }
        if ((pos = 
            name_list_pos(tmp, sequ->ex_nodes->list)) > -1)
         nodes[i] = sequ->ex_nodes->nodes[pos];
        else
          {
	   warning("illegal expand range ignored:", range);
           return 0;
          }
       }
    }
  if (n == 1) nodes[1] = nodes[0];
  return n;
}

int get_sub_range(char* range, struct sequence* sequ, struct node** nodes)
{
     /* returns start and end node (nodes[0] and nodes[1])
        of a range between range_start and range_end of an expanded sequence */
  int i, n;
  char* c[2];
  struct node* c_node;
  char tmp[NAME_L];
  if (sequ == NULL) return 0;
  strcpy(c_dummy, range); stolower(c_dummy);
  c[0] = strtok(c_dummy, "/");
  if ((c[1] = strtok(NULL,"/")) == NULL) /* only one element given */
    n = 1;
  else n = 2;
  for (i = 0; i < n; i++)
    {
     if (*c[i] == '#')
       {
	if (strncmp(c[i], "#s", 2) == 0) nodes[i] = sequ->range_start;
	else if (strncmp(c[i], "#e", 2) == 0) nodes[i] = sequ->range_end;
        else
	  {
	   warning("illegal expand range ignored:", range);
           return 0;
	  }
       }
     else
       {
	strcpy(tmp, c[i]);
        if (square_to_colon(tmp) == 0)
	  {
	   warning("illegal expand range ignored:", range);
           return 0;
	  }
        c_node = sequ->range_start;
        while(c_node)
	  {
	   if (strcmp(c_node->name, tmp) == 0) break;
           if ((c_node = c_node->next) == sequ->range_end)
	     {
	      warning("illegal expand range ignored:", range);
              return 0;
	     }
	  }
	nodes[i] = c_node;
       }
    }
  if (n == 1) nodes[1] = nodes[0];
  return n;
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

double plot_option(char* name)
  /* returns the value of setplot parameters */
{
  double val = zero;
  int i;
  mycpy(c_dummy, name);
  if (plot_options != NULL 
      && (i = name_list_pos(c_dummy, plot_options->par_names)) > -1) 
     val = plot_options->par->parameters[i]->double_value;
  return val;
}

int get_range(char* range, struct sequence* sequ, struct node** nodes)
     /* returns start and end node (nodes[0] and nodes[1])
        of a range in the non-expanded sequence */
{
  int i, n, pos;
  char* c[2];
  char tmp[NAME_L];
  if (sequ == NULL) return 0;
  strcpy(c_dummy, range); stolower(c_dummy);
  c[0] = strtok(c_dummy, "/");
  if ((c[1] = strtok(NULL,"/")) == NULL) /* only one element given */
    n = 1;
  else n = 2;
  for (i = 0; i < n; i++)
    {
     if (*c[i] == '#')
       {
	if (strncmp(c[i], "#s", 2) == 0) nodes[i] = sequ->start;
	else if (strncmp(c[i], "#e", 2) == 0) nodes[i] = sequ->end;
        else
	  {
	   warning("illegal range ignored:", range);
           return 0;
	  }
       }
     else
       {
	strcpy(tmp, c[i]);
        if (square_to_colon(tmp) == 0)
	  {
	   warning("illegal range ignored:", range);
           return 0;
	  }
        if ((pos = 
            name_list_pos(tmp, sequ->nodes->list)) > -1)
         nodes[i] = sequ->nodes->nodes[pos];
        else
          {
	   warning("illegal range ignored:", range);
           return 0;
          }
       }
    }
  if (n == 1) nodes[1] = nodes[0];
  return n;
}

double get_refpos(struct sequence* sequ)
     /* returns the position of a refpos element, or zero */
{
  int i;
  if (sequ != NULL && sequ->refpos != NULL)
    {
     sprintf(c_dummy, "%s:1", sequ->refpos);
     if ((i = name_list_pos(c_dummy, sequ->nodes->list)) < 0)
	fatal_error("'refpos' reference to unknown element:", sequ->refpos);
     return get_node_pos(sequ->nodes->nodes[i], sequ);
    }
  else return zero;
}

int get_table_range(char* range, struct table* table, int* rows)
     /* returns start and end row (rows[0] and rows[1])
        of a range in a table; 0 if not found, 1 (1 row) or 2 ( > 1) */
{
  int i, n;
  char* c[2];
  char tmp[NAME_L], dumtex[3*NAME_L];;
  rows[0] = rows[1] = 0;
  mycpy(c_dummy, range); stolower(c_dummy); strcpy(dumtex, c_dummy); 
  c[0] = strtok(c_dummy, "/");
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
     if (fgets(&ca->c[ca->curr], ca->max - ca->curr, file) == NULL) return 0;
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
  mycpy(c_dummy, name);
  if (strcmp(c_dummy, "beam") == 0)
    {
     mycpy(c_dummy, par);
     if ((p = command_par_string(c_dummy, current_beam)) != NULL)
       {
	strcpy(string, p); length = strlen(p);
       }
    }
  else if (strcmp(c_dummy, "probe") == 0)
    {
     mycpy(c_dummy, par);
     if ((p = command_par_string(c_dummy, probe_beam)) != NULL)
       {
	strcpy(string, p); length = strlen(p);
       }
    }
  else if (strcmp(c_dummy, "survey") == 0)
    {
     mycpy(c_dummy, par);
     if (current_survey != NULL) nl = current_survey->par_names;
     if (nl != NULL && nl->inform[name_list_pos(c_dummy, nl)])
       {
        if ((p = command_par_string(c_dummy, current_survey)) != NULL)
          {
	   strcpy(string, p); length = strlen(p);
          }
       }
    }
  else if (strcmp(c_dummy, "twiss") == 0)
    {
     mycpy(c_dummy, par);
     if (current_twiss != NULL) nl = current_twiss->par_names;
     if (nl != NULL && nl->inform[name_list_pos(c_dummy, nl)])
       {
        if ((p = command_par_string(c_dummy, current_twiss)) != NULL)
          {
	   strcpy(string, p); length = strlen(p);
          }
       }
    }
  else if (strcmp(c_dummy, "sequence") == 0)
    {
     mycpy(c_dummy, par);
     if (current_sequ != NULL && strcmp(c_dummy, "name") == 0)
       {
	p = current_sequ->name;
	strcpy(string, p); length = strlen(p);
       }
    }
  else if (strcmp(c_dummy, "element") == 0)
    {
     mycpy(c_dummy, par);
     if (current_sequ != NULL && strcmp(c_dummy, "name") == 0)
       {
	p = current_node->p_elem->name;
	strcpy(string, p); length = strlen(p);
       }
    }
  else if ((cmd = find_command(c_dummy, stored_commands)) != NULL)
    {
     mycpy(c_dummy, par);
     if ((p = command_par_string(c_dummy, cmd)) != NULL)
       {
	strcpy(string, p); length = strlen(p);
       }
    }
  return length;
}

void get_title(char* tlt, int* l)
     /* copies title from buffer into tl without trailing '\0' */
{
  *l = 0;
  if (title != NULL)
    {
     *l = strlen(title);
     strncpy(tlt, title, *l);
    }
}

double get_value(char* name, char* par)
     /* returns parameter value "par" for command or store "name" if present,
        else INVALID */
{
  struct name_list* nl = NULL;
  mycpy(c_dummy, name);
  mycpy(l_dummy, par);
  if (strcmp(c_dummy, "beam") == 0)
     return command_par_value(l_dummy, current_beam);
  else if (strcmp(c_dummy, "probe") == 0)
     return command_par_value(l_dummy, probe_beam);
  else if (strcmp(c_dummy, "survey") == 0)
    {
     if (current_survey != NULL) nl = current_survey->par_names;
     if (nl != NULL && nl->inform[name_list_pos(l_dummy, nl)])
        return command_par_value(l_dummy, current_survey);
     else return zero;
    }
  else if (strcmp(c_dummy, "twiss") == 0)
    {
     if (current_twiss != NULL) nl = current_twiss->par_names;
     if (nl != NULL && nl->inform[name_list_pos(l_dummy, nl)])
        return command_par_value(l_dummy, current_twiss);
     else return zero;
    }
  else if (strcmp(c_dummy, "sequence") == 0)
    {
     if (strcmp(l_dummy, "l") == 0) return current_sequ->length;
     else if (strcmp(l_dummy, "range_start") == 0)
        return (current_sequ->range_start->position 
                - 0.5 * current_sequ->range_start->length);
     else return INVALID;
    }
  else if (current_command != NULL 
            && strcmp(c_dummy, current_command->name) == 0)
     return command_par_value(l_dummy, current_command);
  else return INVALID;
}

int get_vector(char* name, char* par, double* vector)
     /* returns double "vector" for "par" of command or store "name";
        length is returned as function value (0 if not found) */
{
  mycpy(c_dummy, name);
  mycpy(l_dummy, par);
  if (strcmp(c_dummy, "threader") == 0)
     return command_par_vector(l_dummy, threader_par, vector);
  else return 0;
}

void get_version(char* tlt, int* l)
     /* returns version number */
{
  time_t tmp;
  struct tm* tm;
  int n = strlen(myversion);
  time(&tmp);
  tm = localtime(&tmp); /* split system time */
  strncpy(tlt, myversion, n);
  tlt += n;
  sprintf(tlt, "  %02d/%02d/%02d %02d.%02d.%02d\n",
         tm->tm_mday, tm->tm_mon+1, tm->tm_year%100,
         tm->tm_hour, tm->tm_min, tm->tm_sec);
  *l = n + 19;
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

int int_in_array(int k, int n, int* array)
     /* returns 1 if k in first n elements of array, else 0 */
{
  int j;
  for (j = 0; j < n; j++)  if (k == array[j])  return 1;
  return 0;
}

int in_spec_list(char* string)
     /* checks for presence of special commands IF() etc. */
{
  char* cp;
  int i = 0, n = mymin((int)strlen(string), 100);
  strncpy(c_dummy, string, n); c_dummy[n] = '\0'; stolower(c_dummy);
  supp_char(' ', c_dummy);
  while (special_comm_cnt[i])
    {
     if (special_comm_desc[i][0] == '>')
       {
	if ((cp = strchr(c_dummy, special_comm_desc[i][1])) != NULL)
	  {
	   if (strncmp(++cp, &special_comm_desc[i][2], special_comm_cnt[i])
               == 0)  return i+1;
	  }
       }
     else if (strncmp(c_dummy, &special_comm_desc[i][0],special_comm_cnt[i]) 
               == 0)  return i+1;
     i++;
    }
  return 0;
}

void insert_elem(struct sequence* sequ, struct node* node)
     /* inserts an element in a sequence as function of its position */
{
  struct node* c_node = sequ->start;
  while (c_node != NULL)
    {
     if (node->position <= c_node->position || c_node == sequ->end) break;
     c_node = c_node->next;
    }
  link_in_front(node, c_node);
}

void install_one(struct element* el, char* from_name, double at_value,
                 struct expression* at_expr, double position)
     /* adds an element to a sequence */
{
  struct node* node;
  int i, occ = 1;
  if (strcmp(el->base_type->name, "rfcavity") == 0 &&
      find_element(el->name, edit_sequ->cavities) == NULL)
    add_to_el_list(el, 0, edit_sequ->cavities, 0);
  if ((i = name_list_pos(el->name, occ_list)) < 0)
    i = add_to_name_list(el->name, occ, occ_list);
  else occ = ++occ_list->inform[i];
  node = new_elem_node(el, occ);
  add_to_node_list(node, 0, edit_sequ->nodes);
  node->position = position;
  node->at_value = at_value;
  node->at_expr = at_expr;
  node->from_name = from_name;
  insert_elem(edit_sequ, node);
}

double line_nodes(struct char_p_array* flat)
     /* creates a linked node list from a flat element list of a line */
{
  int i, j, k;
  double pos = zero;
  struct element* el;
  for (j = 0; j < flat->curr; j++)
    {
     if ((el = find_element(flat->p[j], element_list)) == NULL)
       fatal_error("line contains unknown element:", flat->p[j]);
     if (strcmp(el->base_type->name, "rfcavity") == 0 &&
         find_element(el->name, current_sequ->cavities) == NULL)
       add_to_el_list(el, 0, current_sequ->cavities, 0);
     pos += el->length / 2;
     k = 1;
     if ((i = name_list_pos(el->name, occ_list)) < 0)
         i = add_to_name_list(el->name, k, occ_list);
     else k = ++occ_list->inform[i];
     make_elem_node(el, k);
     current_node->at_value = current_node->position = pos;
     pos += el->length / 2;
    }
  return pos;
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
     else break;
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
  char c;

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

struct constraint* make_constraint(int type, struct command_parameter* par)
     /* makes + stores a constraint from command parameter */
{
  struct constraint* new = new_constraint(par->c_type);
  strcpy(new->name, par->name);
  switch(par->c_type)
    {
    case 1: /* minimum */
    case 3: /* both */
      if (par->min_expr == NULL) new->c_min = par->c_min;
      else 
	{
         new->c_min = expression_value(par->min_expr, 2);
         new->ex_c_min = par->min_expr;
	}
      if (par->c_type == 1) break;
    case 2: /* maximum */
      if (par->max_expr == NULL) new->c_max = par->c_max;
      else 
	{
         new->c_max = expression_value(par->max_expr, 2);
         new->ex_c_max = par->max_expr;
	}
      break;
    case 4: /* value */
      if (par->expr == NULL) new->value = par->double_value;
      else 
	{
         new->value = expression_value(par->expr, 2);
         new->ex_value = par->expr;
	}
    }
  if (type == 1) new->weight = command_par_value(new->name, current_weight);
  else           new->weight = command_par_value(new->name, current_gweight);
  return new;
}

struct element* make_element(char* name, char* parent, 
                             struct command* def, int flag)
     /* makes a new element from declaration, stores in list */
{
  double length;
  struct element* el = new_element(name);
  if ((length = command_par_value("l", def)) != INVALID) el->length = length;
  el->def = def;
  if (strcmp(name, parent) == 0)  /* basic element type like drift etc. */
    {
     add_to_el_list(el, def->mad8_type, base_type_list, 1);
     el->parent = el->base_type = el;
    }
  else 
    {
     if((el->parent = find_element(parent, element_list)) == NULL)
       fatal_error("unknown class type:", parent);
     el->base_type = el->parent->base_type;
    }
  add_to_el_list(el, def->mad8_type, element_list, flag);
  return el;
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

int make_line(char* statement)
     /* makes a new line from input command, stores name in name list */
{
  struct macro* m;
  char** toks = tmp_l_array->p;
  char *prs, *psem;
  int i, n, rs, re;
  strcpy(l_dummy, statement);
  if ((prs = strchr(l_dummy, '=')) == NULL) return -3;
  *prs = '\0'; prs++;
  pre_split(l_dummy, l_work, 0);
  mysplit(l_work, tmp_l_array);
  get_bracket_t_range(toks, '(', ')', 0, tmp_l_array->curr-1, &rs, &re);
  if ((n = re - rs - 1) < 0) n = 0; /* number of formal arguments if any */
  m = new_macro(n, 2*strlen(prs), 50);
  strcpy(m->name, toks[0]); rs++;
  for (i = 0; i < n; i++) m->formal->p[i] = permbuff(toks[rs+i]);
  if (n > 0) m->formal->curr = n;
  if ((psem = strchr(prs, ';')) != NULL) *psem = '\0';
  strcpy(l_work, prs);
  pre_split(l_work, m->body->c, 0);
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
  strcpy(l_dummy, statement);
  get_bracket_range(l_dummy, '{', '}', &rs, &re);
  start_2 = rs + 1;
  l_dummy[rs] = '\0'; l_dummy[re] = '\0'; /* drop '{' and '}' */
  pre_split(l_dummy, l_work, 0);
  mysplit(l_work, tmp_l_array);
  get_bracket_t_range(toks, '(', ')', 0, tmp_l_array->curr-1, &rs, &re);
  if ((n = re - rs - 1) < 0) n = 0; /* number of formal arguments if any */
  m = new_macro(n, strlen(&l_dummy[start_2]), 0);
  strcpy(m->name, toks[0]); rs++;
  for (i = 0; i < n; i++) m->formal->p[i] = permbuff(toks[rs+i]);
  if (n > 0) m->formal->curr = n;
  strcpy(m->body->c, &l_dummy[start_2]); m->body->curr = strlen(m->body->c);
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

void make_sequ_from_line(char* name)
     /* converts a line into a sequence from actual line definition */
{
  char** tmp = NULL;
  int pos = name_list_pos(name, line_list->list);
  int spos;
  struct sequence* old_sequ = NULL;
  struct macro* line;
  int mpos = name_list_pos("marker", defined_commands->list);
  struct command* clone = clone_command(defined_commands->commands[mpos]);
  struct element* el;
  if (pos < 0) fatal_error("unknown line: ", name);
  line = line_list->macros[pos];
  line_buffer->curr = 0;
  replace_lines(line, 0, tmp); /* replaces all referenced lines */
  expand_line(line_buffer); /* act on '-' and rep. count */
  current_sequ = new_sequence(name, 0); /* node positions = centre */
  if ((spos = name_list_pos(name, sequences->list)) >= 0) 
    old_sequ = sequences->sequs[spos];
  add_to_sequ_list(current_sequ, sequences);
  if (old_sequ) old_sequ = delete_sequence(old_sequ);
  if (current_sequ->cavities != NULL)  current_sequ->cavities->curr = 0;
  else current_sequ->cavities = new_el_list(100);
  if (occ_list == NULL) 
        occ_list = new_name_list(10000);  /* for occurrence count */
  else occ_list->curr = 0;
  sprintf(c_dummy, "%s$start", current_sequ->name);
  el = make_element(c_dummy, "marker", clone, 0);
  current_node = NULL;
  make_elem_node(el, 1);
  current_sequ->start = current_node;
  current_sequ->length = line_nodes(line_buffer);
  sprintf(c_dummy, "%s$end", current_sequ->name);
  el = make_element(c_dummy, "marker", clone, 0);
  make_elem_node(el, 1);
  current_node->at_value = current_node->position = current_sequ->length;
  current_sequ->end = current_node;
  current_sequ->start->previous = current_sequ->end;
  current_sequ->end->next = current_sequ->start;
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

struct table* make_table(char* name, char* type, char** table_cols, 
                         int* table_types, int rows)
{
  struct table* t;
  struct name_list *cols;
  struct command_list* scl;
  int i, n = 0;
  while (*table_cols[n] != ' ') {n++;}
  cols = new_name_list(n);
  for (i = 0; i < n; i++) 
       add_to_name_list(table_cols[i], table_types[i], cols);
  if ((scl = find_command_list(name, table_select)) != NULL && scl->curr > 0) 
          add_table_vars(cols, scl);
  t = new_table(name, type, rows, cols);
  t->org_cols = n;
  return t;
}

double mult_par(char* par, struct element* el)
     /* returns multipole parameter for par = "k0l" or "k0sl" etc. */
{
  char tmp[12];
  char* p;
  double val = zero, vect[FIELD_MAX];
  int k, l, skew = 0;
  strcpy(tmp, par);
  if (*tmp == 'k' && (p = strchr(tmp, 'l')) != NULL)
    {
     *p = '\0';  /* suppress trailing l */
     if ((p = strchr(tmp, 's')) != NULL)
       {
	skew = 1; *p = '\0';
       }
     sscanf(&tmp[1], "%d", &k);
     if (skew) l = element_vector(el, "ksl", vect);
     else      l = element_vector(el, "knl", vect);
     if (k < l) val = vect[k];
    }
  return val;
}

int next_char(char c, char** toks, int start, int nitem)
     /* returns the number of the token starting with c after token start */
{
  int i;
  for (i = start; i < nitem; i++) if(*toks[i] == c)  return i;
  return -1;
}

int next_constraint(char* name, int* name_l, int* type, double* value, 
                    double* c_min, double* c_max, double* weight)
     /* returns the parameters of the next constraint; 0 = none, else count */
{
  int i, ncp, nbl;
  struct constraint* c_c;
  if (current_node->cl == NULL) return 0;
  if (current_node->con_cnt == current_node->cl->curr)
    {
     current_node->con_cnt = 0; return 0;
    }
  c_c = current_node->cl->constraints[current_node->con_cnt];
  ncp = strlen(c_c->name) < *name_l ? strlen(c_c->name) : *name_l;
  nbl = *name_l - ncp;
  strncpy(name, c_c->name, ncp);
  for (i = 0; i < nbl; i++) name[ncp+i] = ' ';
  *type = c_c->type;
  if (c_c->ex_value == NULL) *value = c_c->value;
  else                       *value = expression_value(c_c->ex_value,2);
  if (c_c->ex_c_min == NULL) *c_min = c_c->c_min;
  else                       *c_min = expression_value(c_c->ex_c_min,2);
  if (c_c->ex_c_max == NULL) *c_max = c_c->c_max;
  else                       *c_max = expression_value(c_c->ex_c_max,2);
  *weight = c_c->weight;
  return ++current_node->con_cnt;
}

int next_global(char* name, int* name_l, int* type, double* value, 
                    double* c_min, double* c_max, double* weight)
     /* returns the parameters of the next global constraint; 
        0 = none, else count */
{
  int i, ncp, nbl;
  struct constraint* c_c;
  if (current_sequ->cl == NULL) return 0;
  if (current_sequ->con_cnt == current_sequ->cl->curr)
    {
     current_sequ->con_cnt = 0; return 0;
    }
  c_c = current_sequ->cl->constraints[current_sequ->con_cnt];
  ncp = strlen(c_c->name) < *name_l ? strlen(c_c->name) : *name_l;
  nbl = *name_l - ncp;
  strncpy(name, c_c->name, ncp);
  for (i = 0; i < nbl; i++) name[ncp+i] = ' ';
  *type = c_c->type;
  if (c_c->ex_value == NULL) *value = c_c->value;
  else                       *value = expression_value(c_c->ex_value,2);
  if (c_c->ex_c_min == NULL) *c_min = c_c->c_min;
  else                       *c_min = expression_value(c_c->ex_c_min,2);
  if (c_c->ex_c_max == NULL) *c_max = c_c->c_max;
  else                       *c_max = expression_value(c_c->ex_c_max,2);
  *weight = c_c->weight;
  return ++current_sequ->con_cnt;
}

int next_start(double* x,double* px,double* y,double* py,double* t,
               double* deltae,double* fx,double* phix,double* fy,double* phiy,
               double* ft,double* phit)
     /* returns the parameters of the next particle to track;
        0 = none, else count */
{
  struct command* comm;
  if (start_cnt == stored_track_start->curr)
    {
     start_cnt = 0; return 0;
    }
  comm = stored_track_start->commands[start_cnt];
  *x = command_par_value("x", comm);
  *px = command_par_value("px", comm);
  *y = command_par_value("y", comm);
  *py = command_par_value("py", comm);
  *t = command_par_value("t", comm);
  *deltae = command_par_value("pt", comm);
  *fx = command_par_value("fx", comm);
  *phix = command_par_value("phix", comm);
  *fy = command_par_value("fy", comm);
  *phiy = command_par_value("phiy", comm);
  *ft = command_par_value("ft", comm);
  *phit = command_par_value("phit", comm);
  return ++start_cnt;
}

int next_vary(char* name, int* name_l, 
                    double* lower, double* upper, double* step)
     /* returns the next variable to be varied during match;
        0 = none, else count */
{
  int i, pos, ncp, nbl;
  double l_step;
  char* v_name;
  struct name_list* nl;
  struct command* comm;
  struct command_parameter_list* pl;
  if (vary_cnt == stored_match_var->curr)
    {
     vary_cnt = 0; return 0;
    }
  comm = stored_match_var->commands[vary_cnt];
  nl = comm->par_names;
  pl = comm->par;
  pos = name_list_pos("name", nl);
  v_name = pl->parameters[pos]->string;
  ncp = strlen(v_name) < *name_l ? strlen(v_name) : *name_l;
  nbl = *name_l - ncp;
  strncpy(name, v_name, ncp);
  for (i = 0; i < nbl; i++) name[ncp+i] = ' ';
  *lower = command_par_value("lower", comm);
  *upper = command_par_value("upper", comm);
  if ((l_step = command_par_value("step", comm)) < ten_m_12) l_step = ten_m_12;
  *step = l_step;
  return ++vary_cnt;
}

int node_al_errors(double* errors)
     /* returns the alignment errors of a node */
{
  if (current_node->p_al_err == NULL) return 0;
  else
    {
     copy_double(current_node->p_al_err->a, errors, 
                 current_node->p_al_err->curr);
     return current_node->p_al_err->curr;
    }
}

int node_fd_errors(double* errors)
     /* returns the field errors of a node */
{
  if (current_node->p_fd_err == NULL) return 0;
  else
    {
     copy_double(current_node->p_fd_err->a, errors, 
                 current_node->p_fd_err->curr);
     return current_node->p_fd_err->curr;
    }
}

void node_string(char* key, char* string, int* l)
     /* returns current node string value for "key" in Fortran format */
     /* l is max. allowed length in string */
{
  char tmp[2*NAME_L];
  char* p;
  int i, l_p, nbl, ncp = 0;
  mycpy(tmp, key);
  if ((p = command_par_string(tmp, current_node->p_elem->def)))
    {
     l_p = strlen(p);
     ncp = l_p < *l ? l_p : *l;
    }
  nbl = *l - ncp;
  for (i = 0; i < ncp; i++) string[i] = p[i];
  for (i = 0; i < nbl; i++) string[ncp+i] = ' ';
}

double spec_node_value(char* par, int* number)  
/* returns value for parameter par of specified node (start = 1 !!) */
{
  double value = zero;
  struct node* node = current_node;
  int n = *number + current_sequ->start_node - 1;
  if (0 <= n && n < current_sequ->n_nodes)
    {
     current_node = current_sequ->all_nodes[n];
     value = node_value(par);
     current_node = node;
    }
  return value;
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

void out_table(char* tname, struct table* t, char* filename)
     /* output of a table */
{
  int j;
  struct command_list* scl;
  while (t->num_cols > t->col_out->max) 
         grow_int_array(t->col_out);
  while (t->curr > t->row_out->max) 
         grow_int_array(t->row_out);
  t->row_out->curr = t->curr;
  if ((scl = find_command_list(tname, table_select)) != NULL
       && scl->curr > 0 && par_present("full", NULL, scl) == 0) 
          prepare_table_file(t, scl);
  else
    {
     for (j = 0; j < t->curr; j++) t->row_out->i[j] = 1;
     for (j = 0; j < t->num_cols; j++) 
          t->col_out->i[j] = j;
     t->col_out->curr = t->num_cols;
    }
  write_table(t, filename);
}

int par_present(char* par, struct command* cmd, struct command_list* c_list)
     /* returns 1 if in cmd or in c_list par is read, else returns 0 */
{
  struct name_list* nl;
  int i, pos;
  if (cmd != NULL)
    {
     nl = cmd->par_names;
     pos = name_list_pos(par, nl);
     if (pos > -1 && nl->inform[pos] > 0)  return 1;
    }
  if (c_list != NULL)
    {
     for (i = 0; i < c_list->curr; i++)
       {
        nl = c_list->commands[i]->par_names;
        pos = name_list_pos(par, nl);
        if (pos > -1 && nl->inform[pos] > 0)  return 1;
       }
    }
  return 0;
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

void prepare_table_file(struct table* t, struct command_list* scl)
     /* prepares a table file for output */
{
  int j;
  if (scl->curr != 0) set_selected_columns(t, scl);
  if (scl->curr == 0 || par_present("full", NULL, scl))
    {
      for (j = 0; j < t->row_out->curr; j++) /* select all rows */
          t->row_out->i[j] = 1;
    }
  else set_selected_rows(t, scl);
}

void pre_split(char* inbuf, char* outbuf, int fill_flag)
     /* inserts blanks between tokens */
     /* fill_flag != 0 makes a 0 to be inserted into an empty "()" */
{
  char c, cp, cpnb = ' ', quote;
  int k, sl = strlen(inbuf), cout = 0, quote_level = 0, rb_level = 0;
  int left_b = 0, in_num = 0, c_digit = 0, f_equal = 0, comm_cnt = 0;
  for (k = 0; k < sl; k++)
   {
    c = inbuf[k];
    if (quote_level > 0)
      {
       if (c == quote) 
	 {
          quote_level--; outbuf[cout++] = c; outbuf[cout++] = ' ';
	 }
       else outbuf[cout++] = c == ' ' ? '@' : c;
      }
    else
      {
       c = inbuf[k];
       switch (c)
	 {
	 case '\"':
	 case '\'':
            quote = c;
            quote_level++; outbuf[cout++] = ' '; outbuf[cout++] = c;
            break;
          case '-':
            if (inbuf[k+1] == '>')
	      {
               outbuf[cout++] = c; break;
	      }
          case '+':
            if (left_b > 0)
	      {
               outbuf[cout++] = ' ';
               outbuf[cout++] = '0';
               outbuf[cout++] = ' ';
               left_b = 0;
	      }
            if (!(in_num > 0 && c_digit > 0 && strchr("ed",cp)) && cout > 0)
               outbuf[cout++] = ' ';
            outbuf[cout++] = c;
            if (!(in_num > 0 && c_digit > 0 && strchr("ed",cp)))
	      {
               outbuf[cout++] = ' ';
               in_num = 1;
	      }
            break;
          case '(':
            rb_level++;
            left_b = 1;
            in_num = 1;
            outbuf[cout++] = ' ';
            outbuf[cout++] = c;
            outbuf[cout++] = ' ';
            break;
          case '>':
            if (cout > 0 && outbuf[cout-1] == '-')
	      {
               outbuf[cout++] = c; break;
	      }
          case ')':
            rb_level--;
            if (fill_flag && cpnb == '(') outbuf[cout++] = '0';
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
            in_num = 1;
            outbuf[cout++] = ' ';
            outbuf[cout++] = c;
            outbuf[cout++] = ' ';
            break;
          case '=':
            f_equal = 1;
            left_b = 0;
            in_num = 1;
            outbuf[cout++] = ' ';
            outbuf[cout++] = c;
            outbuf[cout++] = ' ';
            break;
	  case ',': /* kept behind first "=", or if not first "," */
	    /* not kept inside round brackets before '=' */
            left_b = 0;
            in_num = 1;
            outbuf[cout++] = ' ';
            if (f_equal || (comm_cnt && rb_level == 0))
	      {
               outbuf[cout++] = c;
               outbuf[cout++] = ' ';
	      }
            comm_cnt++;
            break;
	  case ';':
            left_b = 0;
            in_num = 1;
            outbuf[cout++] = ' ';
            break;
          default:
	    if (c != ' ')  left_b = 0;
          if (cout > 0 || c != ' ') outbuf[cout++] = c;
          c_digit += isdigit(c);
          if (strchr(" ,=",c) || is_operator(c))
	    { in_num = 1; c_digit = 0; }
          else  in_num = 
          (isdigit(c) || c == '.' 
          || (strchr("ed",c) && c_digit > 0)) ? in_num : 0;
	 }
       cp = c; if (c != ' ') cpnb = c;
      }
   }
  outbuf[cout] = '\0';
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
         exec_command(); if (stop_flag)  return;
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

void pro_emit(struct in_cmd* cmd)
     /* calls the emit module */
{
  struct command* emit = cmd->clone;
  double e_deltap, e_tol, u0;
  int j, error, keep;
  double* tt;
  double emit_v[3], nemit_v[3], bmax[9], gmax[9], dismax[4], tunes[3],
    sig_v[4], pdamp[3];
  char tmp[16];

  if (current_sequ == NULL || current_sequ->ex_start == NULL)
    {
     warning("sequence not active,", "EMIT ignored");
     return;
    }
  fprintf(prt_file, "enter EMIT module\n");
  if (attach_beam(current_sequ) == 0)
    fatal_error("EMIT - sequence without beam:", current_sequ->name);
  e_deltap = command_par_value("deltap", emit);
  e_tol = command_par_value("tol", emit);
  keep = get_option("twiss_print");
  j = 0;
  set_option("twiss_print", &j);
  zero_double(orbit0, 6);
  zero_double(disp0, 6);
  zero_double(oneturnmat, 36);
  tt = (double*) mycalloc("pro_emit", 216, sizeof(double));
  adjust_beam();
  probe_beam = clone_command(current_beam);
  tmrefe_(oneturnmat); /* one-turn linear transfer map */
  adjust_probe(e_deltap); /* sets correct gamma, beta, etc. */
  print_global(e_deltap);
  adjust_rfc(); /* sets freq in rf-cavities from probe */
  printf("guess: %d %f %f\n",guess_flag, guess_orbit[0],guess_orbit[1]);
  if (guess_flag) copy_double(guess_orbit, orbit0, 6);
  getclor_(orbit0, oneturnmat, tt, &error); /* closed orbit */
  free(tt);
  if (error == 0)
    {
     current_node = current_sequ->ex_start;
     emit_(&e_deltap, &e_tol, orbit0, disp0, oneturnmat, &u0, emit_v, nemit_v,
           bmax, gmax, dismax, tunes, sig_v, pdamp);
     if (e_deltap == zero)
       {
        store_comm_par_value("ex", emit_v[0], current_beam);
        store_comm_par_value("exn", nemit_v[0], current_beam);
        store_comm_par_value("ey", emit_v[1], current_beam);
        store_comm_par_value("eyn", nemit_v[1], current_beam);
        store_comm_par_value("et", emit_v[2], current_beam);
        store_comm_par_value("sigt", sig_v[2], current_beam);
        store_comm_par_value("sige", sig_v[3], current_beam);
        store_comm_par_value("u0", u0, current_beam);
        store_comm_par_value("qs", tunes[2], current_beam);
        store_comm_par_vector("pdamp", pdamp, current_beam);
       }
     else 
       {
	sprintf(tmp, "%14.6f", e_deltap);
        warning("EMIT: beam not updated, non-zero deltap: ", tmp);
       }     
     print_rfc();
    }
  set_option("twiss_print", &keep);
}

void pro_ibs(struct in_cmd* cmd)
  /* control for IBS module */
{
  struct command* keep_beam = current_beam;
  struct name_list* nl = current_ibs->par_names;
  struct command_parameter_list* pl = current_ibs->par;
  char *filename, *table_name;
  int pos, w_file;

  if (twiss_table == NULL)
      warning("no TWISS table present","IBS command ignored");
  else 
    {
     if ((current_beam
          = find_command(twiss_table->org_sequ->name, beam_list)) == NULL) 
       current_beam = find_command("default_beam", beam_list);
     if (probe_beam != NULL) delete_command(probe_beam);
      probe_beam = clone_command(current_beam);
      pos = name_list_pos("file", nl);
      if (nl->inform[pos])
	{
	  if ((filename = pl->parameters[pos]->string) == NULL)
	    {
	      if (pl->parameters[pos]->call_def != NULL)
		filename = pl->parameters[pos]->call_def->string;
	    }
	  if (filename == NULL) filename = permbuff("dummy");
	  w_file = 1;
	}
      else w_file = 0;
      set_option("ibs_table", &w_file); /* fill only if output */
      if (w_file)
	{
         table_name = permbuff("ibs");
         ibs_table = make_table(table_name, "ibs", ibs_table_cols, 
			     ibs_table_types, current_sequ->n_nodes);
         add_to_table_list(ibs_table, table_register);
	}
      adjust_probe(zero); /* sets correct gamma, beta, etc. */
      ibs_();
      if (w_file) out_table(table_name, ibs_table, filename);
    }
  if (probe_beam) delete_command(probe_beam);
  current_beam = keep_beam;
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
           check_table(work);
           this_cmd->tok_list->curr = mysplit(work, this_cmd->tok_list);
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
                    warning("statement not recognised:", work);
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

void pro_match(struct in_cmd* cmd)
     /* controls the matching module */
{
  /* OB 12.2.2002: changed the sequence of if statements so that MAD
                   can go through the whole matching sequence */
  if (strcmp(cmd->tok_list->p[0], "match") == 0)
    {
     match_match(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "cell") == 0)
    {
     warning("CELL command no longer valid, ","use MATCH");
     return;
    }
  else if (match_is_on == 0)
    {
     warning("no MATCH command seen,","ignored");
     return;
    }
  else if (strcmp(cmd->tok_list->p[0], "endmatch") == 0)
    {
     match_end(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "migrad") == 0 ||
           strcmp(cmd->tok_list->p[0], "lmdif") == 0 ||
           strcmp(cmd->tok_list->p[0], "simplex") == 0)
    {
     match_action(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "constraint") == 0)
    {
     match_constraint(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "couple") == 0)
    {
     match_couple(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "fix") == 0)
    {
     match_fix(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "global") == 0)
    {
     match_global(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "level") == 0)
    {
     match_level(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "vary") == 0)
    {
     match_vary(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "weight") == 0)
    {
     match_weight(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "gweight") == 0)
    {
     match_gweight(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "rmatrix") == 0)
    {
     match_rmatrix(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "tmatrix") == 0)
    {
     match_tmatrix(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "global") == 0)
    {
     match_global(cmd);
    }
}

void pro_survey(struct in_cmd* cmd)
     /* calls survey module */
{
  struct name_list* nl = current_survey->par_names;
  struct command_parameter_list* pl = current_survey->par;
  char *filename, *table_name;
  int pos, w_file;
  int iarc = 1, keep;
  if (current_sequ == NULL)
    {
     warning("SURVEY, but no active sequence:", "ignored");
     return;
    }
  fprintf(prt_file, "enter Survey module\n");
  keep = get_option("rbarc");
  set_option("rbarc", &iarc);
  pos = name_list_pos("file", nl);
  if (nl->inform[pos])
    {
     if ((filename = pl->parameters[pos]->string) == NULL)
       {
        if (pl->parameters[pos]->call_def != NULL)
        filename = pl->parameters[pos]->call_def->string;
       }
     if (filename == NULL) filename = permbuff("dummy");
     w_file = 1;
    }
  else w_file = 0;
  pos = name_list_pos("table", nl);
  if(nl->inform[pos]) /* table name specified - overrides save */
    {
     if ((table_name = pl->parameters[pos]->string) == NULL)
       table_name = pl->parameters[pos]->call_def->string;
    }
  else table_name = permbuff("survey");
  survey_table = make_table(table_name, "survey", survey_table_cols, 
                            survey_table_types, current_sequ->n_nodes);
  add_to_table_list(survey_table, table_register);
  survey_();
  if (w_file) out_table(table_name, survey_table, filename);
  set_option("rbarc", &keep);
}

void pro_track(struct in_cmd* cmd)
     /* controls track module */
{
  if (current_sequ == NULL || current_sequ->ex_start == NULL)
    {
     warning("TRACK, but no active sequence:", "ignored");
     return;
    }
  if (strcmp(cmd->tok_list->p[0], "track") == 0)
    {
     track_track(cmd);
    }
  if (strcmp(cmd->tok_list->p[0], "dynap") == 0)
    {
     track_dynap(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "endtrack") == 0)
    {
     track_end(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "observe") == 0)
    {
     track_observe(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "run") == 0)
    {
     track_run(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "ripple") == 0)
    {
     track_ripple(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "start") == 0)
    {
     track_start(cmd->clone);
     cmd->clone_flag = 1;
    }
}

void pro_twiss()
     /* controls twiss module */
{
  struct command* keep_beam = current_beam;
  struct name_list* nl = current_twiss->par_names;
  struct command_parameter_list* pl = current_twiss->par;
  struct int_array* tarr;
  struct node *nodes[2], *use_range[2];
  char *filename, *name, *table_name, *sector_name;
  double tol, tol_keep;
  int i, j, l, lp, k_orb, u_orb, pos, k = 1, ks, w_file, beta_def;
  int keep_info = get_option("info");
  i = keep_info * get_option("twiss_print");
  set_option("info", &i);
  /*
         start command decoding
  */
  pos = name_list_pos("sequence", nl);
  if(nl->inform[pos]) /* sequence specified */
    {
     name = pl->parameters[pos]->string;
     if ((lp = name_list_pos(name, sequences->list)) > -1)
       {
	current_sequ = sequences->sequs[lp];
       }
     else
       {
        warning("unknown sequence ignored:", name);
        return;
       }
    }
  if (current_sequ == NULL || current_sequ->ex_start == NULL)
    {
     warning("sequence not active,", "Twiss ignored");
     return;
    }
  if(get_option("twiss_print")) fprintf(prt_file, "enter Twiss module\n");
  if (attach_beam(current_sequ) == 0)
    fatal_error("TWISS - sequence without beam:", current_sequ->name);
  pos = name_list_pos("table", nl);
  if(nl->inform[pos]) /* table name specified - overrides save */
    {
     if ((table_name = pl->parameters[pos]->string) == NULL)
      table_name = pl->parameters[pos]->call_def->string;
    }
  else if((pos = name_list_pos("save", nl)) > -1 &&
          nl->inform[pos]) /* save name specified */
    {
     if ((table_name = pl->parameters[pos]->string) == NULL)
      table_name = pl->parameters[pos]->call_def->string;
    }
  else table_name = "twiss";
  if ((ks = get_value(current_command->name,"sectormap")) != 0)
    {
     set_option("twiss_sector", &k);
     pos = name_list_pos("sectorfile", nl);
     if(nl->inform[pos]) 
       {
        if ((sector_name = pl->parameters[pos]->string) == NULL)
         sector_name = pl->parameters[pos]->call_def->string;
       }
     else  sector_name = pl->parameters[pos]->call_def->string;
     if ((sec_file = fopen(sector_name, "w")) == NULL)
          fatal_error("cannot open output file:", sector_name);
    }
  use_range[0] = current_sequ->range_start;
  use_range[1] = current_sequ->range_end;
  if ((pos = name_list_pos("range", nl)) > -1 && nl->inform[pos])
    {
     if (get_sub_range(pl->parameters[pos]->string, current_sequ, nodes))
       {
	current_sequ->range_start = nodes[0];
	current_sequ->range_end = nodes[1];
       }
     else warning("illegal range ignored:", pl->parameters[pos]->string);
    }
  for (j = 0; j < current_sequ->n_nodes; j++)
    {
     if (current_sequ->all_nodes[j] == current_sequ->range_start) break;
    }
  if((pos = name_list_pos("useorbit", nl)) > -1 &&nl->inform[pos]) 
    /* orbit specified */
    {
     if (current_sequ->orbits == NULL)  
        warning("orbit not found, ignored: ", pl->parameters[pos]->string);
     else
       {
        name = pl->parameters[pos]->string;
        if ((u_orb = name_list_pos(name, current_sequ->orbits->names)) < 0)
            warning("orbit not found, ignored: ", name);
        else set_option("useorbit", &k);
       }
    }
  pos = name_list_pos("keeporbit", nl);
  if(nl->inform[pos]) /* orbit specified */
    {
     name = pl->parameters[pos]->string;
     if (current_sequ->orbits == NULL)  
       current_sequ->orbits = new_vector_list(10);
     else if (current_sequ->orbits->curr == current_sequ->orbits->max)
	      grow_vector_list(current_sequ->orbits);
     if ((k_orb = name_list_pos(name, current_sequ->orbits->names)) < 0)
       {
        k_orb = add_to_name_list(permbuff(name), 0, 
                                 current_sequ->orbits->names);
        current_sequ->orbits->vectors[k_orb] = new_double_array(6);
       }
     set_option("keeporbit", &k);
    }
  pos = name_list_pos("file", nl);
  if (nl->inform[pos])
    {
     if ((filename = pl->parameters[pos]->string) == NULL)
       {
        if (pl->parameters[pos]->call_def != NULL)
        filename = pl->parameters[pos]->call_def->string;
       }
     if (filename == NULL) filename = permbuff("dummy");
     w_file = 1;
    }
  else w_file = 0;
  tol_keep = get_variable("twiss_tol");
  pos = name_list_pos("tolerance", nl);
  if (nl->inform[pos])
    {
     tol = command_par_value("tolerance", current_twiss);
    }
  /*
             end of command decoding
  */
  current_sequ->start_node = j;
  zero_double(orbit0, 6);
  zero_double(disp0, 6);
  zero_double(oneturnmat, 36);
  if ((beta_def = twiss_input(current_twiss)) < 0)
    {
     if (beta_def == -1) warning("unknown beta0,", "Twiss ignored");
     else if (beta_def == -2) 
         warning("betx or bety missing,", "Twiss ignored");
     return;
    }
  set_option("twiss_inval", &beta_def);
  set_option("twiss_summ", &k);
  pos = name_list_pos("chrom", nl);
  set_option("twiss_chrom", &nl->inform[pos]);
  set_option("twiss_save", &k);
  set_twiss_deltas(current_twiss);
  adjust_beam();
  probe_beam = clone_command(current_beam);
  tmrefe_(oneturnmat); /* one-turn linear transfer map */
  summ_table = make_table("summ", "summ", summ_table_cols, summ_table_types,
                          twiss_deltas->curr+1);
  add_to_table_list(summ_table, table_register);
  l = strlen(table_name);
  tarr = new_int_array(l+1);
  conv_char(table_name, tarr);
  if (get_option("twiss_sector"))
    {
     reset_sector(current_sequ, 0);
     set_sector();
    }
  if (get_option("useorbit"))
      copy_double(current_sequ->orbits->vectors[u_orb]->a, orbit0, 6);
  else if (guess_flag)
    {
     for (i = 0; i < 6; i++)
       {
        if (guess_orbit[i] != zero) orbit0[i] = guess_orbit[i];
       }
    }
  for (i = 0; i < twiss_deltas->curr; i++)
    {
     twiss_table = make_table(table_name, "twiss", twiss_table_cols, 
                            twiss_table_types, current_sequ->n_nodes);
     twiss_table->dynamic = 1; /* flag for table row access to current row */
     add_to_table_list(twiss_table, table_register);
     current_sequ->tw_table = twiss_table;
     twiss_table->org_sequ = current_sequ;
     adjust_probe(twiss_deltas->a[i]); /* sets correct gamma, beta, etc. */
     adjust_rfc(); /* sets freq in rf-cavities from probe */
     current_node = current_sequ->ex_start;
     twiss_(oneturnmat, disp0, tarr->i);
     if ((twiss_success = get_option("twiss_success")))
       {
        if (get_option("keeporbit"))  copy_double(orbit0, 
                        current_sequ->orbits->vectors[k_orb]->a, 6);
        fill_twiss_header(twiss_table);
        if (i == 0) exec_savebeta(); /* fill beta0 at first delta_p only */
        if (w_file) out_table(table_name, twiss_table, filename);
       }
     else puts("Twiss failed");
    }
  if (sec_file)  
    {
     fclose(sec_file); sec_file = NULL;
    }
  tarr = delete_int_array(tarr);
  if (twiss_success && get_option("twiss_print")) print_table(summ_table);
  /* cleanup */
  current_beam = keep_beam;
  probe_beam = delete_command(probe_beam);
  k = 0;
  set_option("couple", &k);
  set_option("chrom", &k);
  set_option("rmatrix", &k);
  set_option("centre", &k);
  set_option("twiss_sector", &k);
  set_option("keeporbit", &k);
  set_option("useorbit", &k);
  set_option("info", &keep_info);
  set_variable("twiss_tol", &tol_keep);
  current_sequ->range_start = use_range[0];
  current_sequ->range_end = use_range[1];
}

void put_info(char* t1, char* t2) 
{
  if (get_option("info") && get_option("warn")) 
      printf("++++++ info: %s %s\n",t1,t2);
}

struct table* read_table(struct in_cmd* cmd)
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
  while (fgets(l_dummy, AUX_LG, tab_file))
    {
     cc = strtok(l_dummy, " \"\n");
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
	tnl = new_name_list(20);
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
           if (strcmp(tmp,"%s") == 0) t->s_cols[i][t->curr] = tmpbuff(cc);
           else if (strcmp(tmp,"%d") == 0 || strcmp(tmp,"%hd") == 0)
	     {
	      sscanf(cc, tmp, &k); t->d_cols[i][t->curr] = k;
	     }
           else sscanf(cc, tmp, &t->d_cols[i][t->curr]);
           if (i+1 < tnl->curr)
	     {
              if ((cc =strtok(NULL, " \"\n")) == NULL)
	        {
	         warning("incomplete table line starting with:", l_dummy);
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

void remove_from_command_list(char* label, struct command_list* list)
{
  int i;
  if ((i = remove_from_name_list(label, list->list)) > -1)
    {
     if (i < --list->curr)
       {
	delete_command(list->commands[i]);
        list->commands[i] = list->commands[list->curr];
       }
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

void remove_from_node_list(struct node* node, struct node_list* nodes)
{
  int i;
  if ((i = remove_from_name_list(node->name, nodes->list)) > -1)
     nodes->nodes[i] = nodes->nodes[--nodes->curr];
}

int remove_one(struct node* node)
{
  int pos;
  /* removes a node from a sequence being edited */
  if ((pos = name_list_pos(node->p_elem->name, occ_list)) < 0)  return 0;
  if (node->previous != NULL) node->previous->next = node->next;
  if (node->next != NULL) node->next->previous = node->previous;
  if (occ_list->inform[pos] == 1)
    {
     remove_from_node_list(node, edit_sequ->nodes);
     remove_from_name_list(node->p_elem->name, occ_list);
    }
  else --occ_list->inform[pos];
  free(node);
  return 1;
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

void replace_one(struct node* node, struct element* el)
     /* replaces an existing node by a new one made from el */
{
  int i, k = 1;
  remove_from_node_list(node, edit_sequ->nodes);
  if ((i = name_list_pos(el->name, occ_list)) < 0)
  i = add_to_name_list(el->name, 1, occ_list);
  else k = ++occ_list->inform[i];
  strcpy(node->name, compound(el->name, k));
  add_to_node_list(node, 0, edit_sequ->nodes);
  node->p_elem = el;
  node->base_name = el->base_type->name;
  node->length = el->length;
  if (strcmp(el->base_type->name, "rfcavity") == 0 &&
      find_element(el->name, edit_sequ->cavities) == NULL)
    add_to_el_list(el, 0, edit_sequ->cavities, 0);
}

void replace_lines(struct macro* org, int replace, char** reps)
     /* replaces lines in line by elements - recursive */
{
  int i, j, k, l, n, pos; 
  int mf = replace < org->n_formal ? replace : org->n_formal;
  char* p;
  struct macro* line;
  if (org->tokens == NULL) fatal_error("line not split:", org->name);
  line = clone_macro(org);
  for (j = 0; j < mf; j++)
    {
     for (i = 0; i < line->tokens->curr; i++)
       {
        p = line->tokens->p[i];
        if (isalpha(*p) && strcmp(line->formal->p[j], p) == 0)  
          line->tokens->p[i] = reps[j];
       }
    }
  for (i = 0; i < line->tokens->curr; i++)
    {
     p = line->tokens->p[i];
     if (isalpha(*p) && (pos = name_list_pos(p, line_list->list)) > -1)
       {
	if (*line->tokens->p[i+1] == '(') /* formal arguments */
	  {
           for (k = i+2; k < line->tokens->curr; k++)
	    if (*line->tokens->p[k] == ')') break;
	   n = k - i - 2;
           l = k;
	  }
        else 
	  {
	   n = 0; l = i;
	  }
        replace_lines(line_list->macros[pos], n, &line->tokens->p[i+2]);
        i = l;
       }
     else
       {
        if (line_buffer->curr == line_buffer->max) 
            grow_char_p_array(line_buffer);
        line_buffer->p[line_buffer->curr++] = p;
       }
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

void reset_count(char* table) /* resets table counter to zero */
{
  int pos;
  struct table* t;
  mycpy(c_dummy, table);
  if ((pos = name_list_pos(c_dummy, table_register->names)) > -1)
    t = table_register->tables[pos];
  else return;
  t->curr = 0;
}

void reset_errors(struct sequence* sequ)
     /* zeros the sel_err node flag for all nodes of an expanded sequence */
{
  struct node* c_node;
  if (sequ != NULL && sequ->ex_start != NULL && sequ->ex_end != NULL)
    {
     c_node = sequ->ex_start;
     while (c_node != NULL)
       {
	c_node->sel_err = 0;
        if (c_node == sequ->ex_end) break;
        c_node = c_node->next;
       }
    }
}

void reset_sector(struct sequence* sequ, int val)
     /* sets node->sel_sector = val for all nodes of an expanded sequence */
{
  struct node* c_node;
  if (sequ != NULL && sequ->ex_start != NULL && sequ->ex_end != NULL)
    {
     c_node = sequ->ex_start;
     while (c_node != NULL)
       {
	c_node->sel_sector = val;
        if (c_node == sequ->ex_end) break;
        c_node = c_node->next;
       }
    }
}

int restart_sequ()
{
  current_node = current_sequ->range_start;
  return 1;
}

int retreat_node()
     /* replaces current node by previous node; 0 = already at start, else 1 */
{
  if (current_node == current_sequ->range_start)  return 0;
  current_node = current_node->previous;
  return 1;
}

double rfc_slope()
     /* calculates the accumulated "slope" of all cavities */
{
  double slope = zero, lag, volt, harmon, charge, pc;
  struct node* c_node = current_sequ->range_start;
  struct element* el;
  charge = command_par_value("charge", current_beam);
  pc = command_par_value("pc", current_beam);
  do
    {
     el = c_node->p_elem;
     if (strcmp(el->base_type->name, "rfcavity") == 0 &&
         (harmon = command_par_value("harmon", el->def)) > zero)
       {
	volt = command_par_value("volt", el->def);
	lag = command_par_value("lag", el->def);
        slope += ten_m_3 * charge * volt * harmon * cos(twopi * lag) / pc;
       }
     if (c_node == current_sequ->range_end) break;
     c_node = c_node->next;
    }
  while (c_node != NULL);
  return slope;
}

int scan_expr(int c_item, char** item)   /* split input */

     /* scans expressions for parameters, elements, numbers, and functions

        categories: 1: variable, 3: floating constant, 4: operator
        operator types:
        1 = +, 2 = -, 3 = *, 4 = /, 5 = ^ (power), 6 = function (from functs)
     */
    
{
  int i, lp, lx, l_cat = 0, level = 0, pos, f_level[MAX_ITEM];
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

void sector_out(double* pos, double* kick, double* rmatrix, double* tmatrix)
     /* writes a sector map to sec_file */
{
  int i;
  fprintf(sec_file, " %-20.6g   %s\n", *pos, current_node->p_elem->name);
  for (i = 0; i < 6; i++) fprintf(sec_file, "%15.8e ", kick[i]);
  fprintf(sec_file,"\n");
  for (i = 0; i < 36; i++) 
    {
     fprintf(sec_file, "%15.8e ", rmatrix[i]);
     if ((i+1)%6 == 0)  fprintf(sec_file,"\n");
    }
  for (i = 0; i < 216; i++) 
    {
     fprintf(sec_file, "%15.8e ", tmatrix[i]);
     if ((i+1)%6 == 0)  fprintf(sec_file,"\n");
    }
}

void seq_cycle(struct in_cmd* cmd)
     /* cycles a sequence */
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  struct node *node, *clone;
  char* name = NULL;
  int pos = name_list_pos("start", nl);
  if (nl->inform[pos] && (name = pl->parameters[pos]->string) != NULL)
    {
     sprintf(c_dummy, "%s:1", name);
     if ((pos = name_list_pos(c_dummy, edit_sequ->nodes->list)) > -1)
       {
	node = edit_sequ->nodes->nodes[pos];
        sprintf(c_dummy, "%s$_p_", node->name);
	if (strchr(node->previous->name, '$') == NULL)
	  {
	   clone = clone_node(node, 0);
           strcpy(clone->name, c_dummy);
           link_in_front(clone, node);
	  }
        edit_sequ->start = node;
        edit_sequ->end = node->previous;
        set_new_position(edit_sequ);
        all_node_pos(edit_sequ);
       }
     else warning("cycle: unknown element ignored:", name);
    }
  else warning("cycle: no start given,","ignored");
}

void seq_edit_main(struct in_cmd* cmd)
     /* controls sequence editing */
{
  int k = cmd->decl_start - 1;
  char** toks = cmd->tok_list->p;
  if (strcmp(toks[k], "seqedit") == 0)  seq_edit(cmd);
  else if(edit_is_on)
    {
     if (strcmp(toks[k], "install") == 0)  seq_install(cmd);
     else if (strcmp(toks[k], "move") == 0)  seq_move(cmd);
     else if (strcmp(toks[k], "remove") == 0)  seq_remove(cmd);
     else if (strcmp(toks[k], "cycle") == 0)  seq_cycle(cmd);
     else if (strcmp(toks[k], "flatten") == 0)  seq_flatten(edit_sequ);
     else if (strcmp(toks[k], "reflect") == 0)  seq_reflect(cmd);
     else if (strcmp(toks[k], "replace") == 0)  seq_replace(cmd);
     else if (strcmp(toks[k], "endedit") == 0)  seq_end(cmd);
    }
  else warning("seqedit command outside edit", "ignored");
}

void seq_edit(struct in_cmd* cmd)
     /* executes seqedit command */
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  char* name = NULL;
  int pos;
  pos = name_list_pos("sequence", nl);
  if (nl->inform[pos] && (name = pl->parameters[pos]->string) != NULL)
    {
     if ((pos = name_list_pos(name, sequences->list)) >= 0)
  	  {
	   edit_is_on = 1;
           seqedit_install = seqedit_move = seqedit_remove = 0;
  	   edit_sequ = sequences->sequs[pos];
           if (edit_sequ->ex_start != NULL)
  	     {
              edit_sequ->ex_nodes = delete_node_list(edit_sequ->ex_nodes);
              edit_sequ->ex_start = delete_node_ring(edit_sequ->ex_start);
  	     }
           if (occ_list == NULL) 
               occ_list = new_name_list(10000);  /* for occurrence count */
           else occ_list->curr = 0;
           resequence_nodes(edit_sequ);
           all_node_pos(edit_sequ);
  	  }
     else warning("unknown sequence:", "ignored");
    }
  else warning("seqedit without sequence:", "ignored");
}

void seq_end(struct in_cmd* cmd)
     /* executes endedit command */
{
  char tmp[8];
  sprintf(tmp, "%d", seqedit_install);
  put_info("seqedit - number of elements installed: ", tmp);
  sprintf(tmp, "%d", seqedit_move);
  put_info("seqedit - number of elements moved:     ", tmp);
  sprintf(tmp, "%d", seqedit_remove);
  put_info("seqedit - number of elements removed:   ", tmp);
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
        occ_list = new_name_list(10000);  /* for occurrence count */
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

void seq_install(struct in_cmd* cmd)
     /* executes install command */
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  struct element *cl, *el;
  struct node* c_node;
  struct expression* expr = NULL;
  double at, from = zero;
  char name[NAME_L], *pname, *name_e = NULL, *name_c = NULL, *from_name = NULL;
  int k, pos, any = 0;
  int pos_e = name_list_pos("element", nl);
  int pos_c = name_list_pos("class", nl);
  if (nl->inform[pos_e] && (name_e = pl->parameters[pos_e]->string) != NULL)
    {
     if (nl->inform[pos_c] && (name_c = pl->parameters[pos_c]->string) != NULL)
       {
        if ((cl = find_element(name_c, element_list)) == NULL)
	  {
	   warning("ignored because of unknown class:", name_c);
           return;
	  }
        else
	  {
	   el = clone_element(cl);
           strcpy(el->name, name_e);
           add_to_el_list(el, cl->def->mad8_type, element_list, 1);
	  }
       }
     else if ((el = find_element(name_e, element_list)) == NULL)
       { 
	warning("ignored, unknown command or element:", name_c); return;
       }
    }
  else
    {
     warning("no element specified,","ignored"); return;
    }
  if (nl->inform[name_list_pos("at", nl)] == 0)
    {
     warning("no 'at':", "ignored"); return;
    }
  at = command_par_value("at", cmd->clone);
  expr = clone_expression(command_par_expr("at", cmd->clone));
  pos = name_list_pos("from", nl);
  if (nl->inform[pos]) 
    {
     from_name = pl->parameters[pos]->string;
     if (strcmp(from_name, "selected") == 0)
       {
        if (seqedit_select->curr == 0)
          {
           warning("no active select commands:", "ignored"); return;
          }
        else
          {
	   if (get_select_ranges(edit_sequ, seqedit_select, selected_ranges) 
               == 0) any = 1;
           c_node = edit_sequ->start;
           while (c_node != NULL)
	     {
	      if (any 
                  || name_list_pos(c_node->name, selected_ranges->list) > -1)
		{
		 for (k = 0; k < seqedit_select->curr; k++)
		   {
	            my_repl(":", "[", c_node->name, name);
                    strcat(name, "]");
                    if (strchr(name, '$') == NULL &&
		        pass_select(c_node->name, 
                          seqedit_select->commands[k])) break;
		   }
                 if (k < seqedit_select->curr)
		   {
                    from = get_node_pos(c_node, edit_sequ);
                    pname = permbuff(name);
                    install_one(el, pname, at, expr, at+from);
                    seqedit_install++;
		   }
		}
               if (c_node == edit_sequ->end) break;
               c_node = c_node->next;
	     }
	  }
       }
     else
       {
        from_name = permbuff(pl->parameters[pos]->string);
        if ((from = hidden_node_pos(from_name, edit_sequ)) == INVALID)
          {
	   warning("ignoring 'from' reference to unknown element:", from_name);
           return;
          }
        install_one(el, from_name, at, expr, at+from);
        seqedit_install++;
       }
    }
  else 
    {
     install_one(el, from_name, at, expr, at);
     seqedit_install++;
    }
}

void seq_move(struct in_cmd* cmd)
     /* executes move command */
{
  char *name, *from_name;
  double at, by, to, from;
  int any = 0, k;
  struct node *node;
  struct element* el;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int pos = name_list_pos("element", nl);
  if (nl->inform[pos] && (name = pl->parameters[pos]->string) != NULL)
    {
     if (strcmp(name, "selected") == 0)
       {
        if (seqedit_select->curr == 0)
          {
           warning("no active select commands:", "ignored"); return;
          }
        else
          {
           if (nl->inform[name_list_pos("by", nl)] == 0)
	     {
              warning("no 'by' given,", "ignored"); return;
	     }
	   by = command_par_value("by", cmd->clone);
	   if (get_select_ranges(edit_sequ, seqedit_select, selected_ranges) 
               == 0) any = 1;
           node = edit_sequ->start;
           while (node != NULL)
	     {
	      if (any 
                  || name_list_pos(node->name, selected_ranges->list) > -1)
		{
		 name = NULL;
		 for (k = 0; k < seqedit_select->curr; k++)
		   {
		    if (node->p_elem != NULL) name = node->p_elem->name;
                    if (name != NULL && strchr(name, '$') == NULL &&
		        pass_select(name, 
                          seqedit_select->commands[k])) break;
		   }
                 if (k < seqedit_select->curr)
		   {
		    at = node->position + by;
		    el = node->p_elem;
                    if (remove_one(node) > 0)
		      {
                       install_one(el, NULL, at, NULL, at);
                       seqedit_move++;
		      }
		   }
		}
               if (node == edit_sequ->end) break;
               node = node->next;
	     }
	  }
       }
     else
       {
        strcpy(c_dummy, name);
        square_to_colon(c_dummy);
        if ((pos = name_list_pos(c_dummy, edit_sequ->nodes->list)) > -1)
          {
	   node = edit_sequ->nodes->nodes[pos];
           if (nl->inform[name_list_pos("by", nl)] == 0)
	     {
              if (nl->inform[name_list_pos("to", nl)] == 0)
                {
                 warning("no position given,", "ignored"); return;
                }
              to = command_par_value("to", cmd->clone);
              pos = name_list_pos("from", nl);
              if (nl->inform[pos])
		{
                 from_name = pl->parameters[pos]->string;
                 if ((from = hidden_node_pos(from_name, edit_sequ)) == INVALID)
                    {
	             warning("ignoring 'from' reference to unknown element:", 
                     from_name);
                     return;
                    }
		}
              at = to + from;
	     }
           else
	     {
	      by = command_par_value("by", cmd->clone);
	      at = node->position + by;
	     }
	   el = node->p_elem;
           if (remove_one(node) > 0)
	     {
              install_one(el, NULL, at, NULL, at);
              seqedit_move++;
	     }
	  }
       }
    }
}

void seq_reflect(struct in_cmd* cmd)
     /* executes reflect command */
{
  struct node *tmp, *c_node;
  c_node = edit_sequ->start;
  while (c_node != NULL)
    {
     tmp = c_node->next;
     c_node->next = c_node->previous;
     c_node->previous = tmp;
     if (c_node == edit_sequ->end) break;
     c_node = tmp;
    }
  tmp = edit_sequ->start;
  edit_sequ->start = edit_sequ->end;
  edit_sequ->end = tmp;
  c_node = edit_sequ->start;
  edit_sequ->range_start = edit_sequ->start;
  edit_sequ->range_end = edit_sequ->end;
  while (c_node != NULL)
    {
     c_node->at_expr = NULL;
     c_node->from_name = NULL;
     c_node->position = c_node->at_value 
             = edit_sequ->length - c_node->position;
     if (c_node == edit_sequ->end) break;
     c_node = c_node->next;
    }
}

void seq_remove(struct in_cmd* cmd)
     /* executes remove command */
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  struct node *c_node;
  char *name;
  int k, any = 0;
  int pose = name_list_pos("element", nl);
  if (nl->inform[pose] && (name = pl->parameters[pose]->string) != NULL)
    {
     if (strcmp(name, "selected") == 0)
       {
        if (seqedit_select->curr == 0)
          {
           warning("no active select commands:", "ignored"); return;
          }
        else
          {
	   if (get_select_ranges(edit_sequ, seqedit_select, selected_ranges) 
               == 0) any = 1;
           c_node = edit_sequ->start;
           while (c_node != NULL)
	     {
	      if (any 
                  || name_list_pos(c_node->name, selected_ranges->list) > -1)
		{
		 name = NULL;
		 for (k = 0; k < seqedit_select->curr; k++)
		   {
		    if (c_node->p_elem != NULL) name = c_node->p_elem->name;
                    if (name != NULL && strchr(name, '$') == NULL &&
		        pass_select(name, 
                          seqedit_select->commands[k])) break;
		   }
                 if (k < seqedit_select->curr)
		   {
                    seqedit_remove += remove_one(c_node);
		   }
		}
               if (c_node == edit_sequ->end) break;
               c_node = c_node->next;
	     }
	  }
       }
     else
       {
	strcpy(c_dummy, name);
        square_to_colon(c_dummy);
        if ((pose = name_list_pos(c_dummy, edit_sequ->nodes->list)) > -1)
	  {
           seqedit_remove += remove_one(edit_sequ->nodes->nodes[pose]);
	  }
        else warning("ignored because of unknown element:", name);
       }
    }
  else  warning("no element specified,","ignored");
}

void seq_replace(struct in_cmd* cmd)
     /* executes replace command */
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  struct node *node, *c_node;
  char* name;
  struct element* el;
  int any = 0, k, pos = name_list_pos("element", nl);
  if (nl->inform[pos] && (name = pl->parameters[pos]->string) != NULL)
    {
     if (strcmp(name, "selected") == 0)
       {
        if (seqedit_select->curr == 0)
          {
           warning("no active select commands:", "ignored"); return;
          }
        else
          {
           pos = name_list_pos("by", nl);
           if (nl->inform[pos] && (name = pl->parameters[pos]->string) != NULL)
              {
               if ((el = find_element(name, element_list)) == NULL)
		 {
                  warning("ignoring unknown 'by' element:",name);
                  return;
		 }
	      }
           else 
	     {
              warning("'by' missing, ","ignored");
              return;
	     }
	   if (get_select_ranges(edit_sequ, seqedit_select, selected_ranges) 
               == 0) any = 1;
           c_node = edit_sequ->start;
           while (c_node != NULL)
	     {
	      if (any 
                  || name_list_pos(c_node->name, selected_ranges->list) > -1)
		{
		 name = NULL;
		 for (k = 0; k < seqedit_select->curr; k++)
		   {
		    if (c_node->p_elem != NULL) name = c_node->p_elem->name;
                    if (name != NULL && strchr(name, '$') == NULL &&
		        pass_select(name, 
                          seqedit_select->commands[k])) break;
		   }
                 if (k < seqedit_select->curr) replace_one(node, el);
		}
               if (c_node == edit_sequ->end) break;
               c_node = c_node->next;
	     }
	  }
       }
     else
       {
        strcpy(c_dummy, name);
        square_to_colon(c_dummy);
        if ((pos = name_list_pos(c_dummy, edit_sequ->nodes->list)) > -1)
          {
	   node = edit_sequ->nodes->nodes[pos];
           pos = name_list_pos("by", nl);
           if (nl->inform[pos] && (name = pl->parameters[pos]->string) != NULL)
              {
               if ((el = find_element(name, element_list)) != NULL)
	         replace_one(node, el);
               else warning("ignoring unknown 'by' element: ",name);
	      }
           else warning("'by' missing, ","ignored");
          }
        else warning("ignored because of unknown element: ", name);
       }
    }
  else  warning("no element specified, ","ignored");
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
     else if (strcmp(string, "setplot") == 0)
       {
	if (plot_options != NULL) delete_command(plot_options);
        plot_options = clone_command(defined_commands->commands[pos]);
       }
     else if (strcmp(string, "threader") == 0)
       {
        threader_par = clone_command(defined_commands->commands[pos]);
       }
     else if (strcmp(string, "beam") == 0)
       {
        current_beam = clone_command(defined_commands->commands[pos]);
        beam_clone = clone_command(defined_commands->commands[pos]);
        for (i = 0; i < beam_clone->par_names->curr; i++)
          beam_clone->par_names->inform[i] = 1; /* mark as "read" */
        update_beam(beam_clone);
        delete_command(beam_clone);
       }
    }
}

void sequence_name(char* name, int* l)
     /* returns current sequence name in Fortran format */
{
  int sname_l = strlen(current_sequ->name);
  int i, ncp = sname_l < *l ? sname_l : *l;
  int nbl = *l - ncp;
  for (i = 0; i < ncp; i++) name[i] = current_sequ->name[i];
  for (i = 0; i < nbl; i++) name[ncp+i] = ' ';
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
       }
    }
}

void set_new_position(struct sequence* sequ)
     /* sets a new node position for all nodes */
{
  struct node* c_node = sequ->start;
  double zero_pos = c_node->position;
  int flag = 0;
  while (c_node != NULL)
    {
     if (c_node->from_name == NULL)
       {
        c_node->position -= zero_pos;
        if (c_node->position < zero || (flag && c_node->position == zero))
	   c_node->position += sequ->length;
        if (c_node->position > zero) flag = 1;
        c_node->at_value = c_node->position;
        c_node->at_expr = NULL;
       }
     if (c_node == sequ->end) break;
     c_node = c_node->next;
    }
  c_node->position = c_node->at_value = sequ->length;
}

void set_node_bv(struct sequence* sequ)
     /* sets bv flag for all nodes */
{
  struct node* c_node = sequ->ex_start;
  double beam_bv;
  beam_bv = command_par_value("bv", current_beam);
  while (c_node != NULL)
    {
     if(command_par_value("magnet", c_node->p_elem->def))
       {
        c_node->other_bv = beam_bv;
        if (c_node->p_elem->bv) c_node->dipole_bv = beam_bv;
        else                    c_node->dipole_bv = 1;
       }
     if (c_node == sequ->ex_end) break;
     c_node = c_node->next;
    }
}

void set_option(char* str, int* opt)
     /* sets an (old or new) option with name "str", 
        value *opt (0 flase, 1 true) */
{
  int i, j, k;
  char* bc;
  mycpy(c_dummy, str); bc = permbuff(c_dummy);
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

void set_value(char* name, char* par, double* value)
     /* sets parameter value "par" for command or store "name" if present */
{
  mycpy(c_dummy, name);
  mycpy(l_dummy, par);
  if (strcmp(c_dummy, "beam") == 0)
     set_command_par_value(l_dummy, current_beam, *value);
  else if (strcmp(c_dummy, "probe") == 0)
     set_command_par_value(l_dummy, probe_beam, *value);
  else if (strcmp(c_dummy, "survey") == 0)
     set_command_par_value(l_dummy, current_survey, *value);
  else if (strcmp(c_dummy, "twiss") == 0)
     set_command_par_value(l_dummy, current_twiss, *value);
  else if (current_command != NULL 
            && strcmp(c_dummy, current_command->name) == 0)
     set_command_par_value(l_dummy, current_command, *value);
}

double sss_variable(char* name)
{
  char comm[NAME_L];
  char par[NAME_L];
  double val = zero;
  struct variable* var;
  struct element* el;
  struct command* cmd;
  char *p, *n = c_dummy, *q = comm;
  mycpy(c_dummy, name);
  if ((p = strstr(c_dummy, "->")) == NULL) /* variable */
    {
     if ((var = find_variable(c_dummy, variable_list)) != NULL)
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

void set_variable(char* name, double* value)
{
  /* sets variable name to value */
  char comm[NAME_L];
  char par[NAME_L];
  struct variable* var;
  double val = *value;
  struct element* el;
  struct command* cmd;
  char *p, *n = c_dummy, *q = comm;
  mycpy(c_dummy, name);
  if ((p = strstr(c_dummy, "->")) == NULL) /* variable */
    {
     if ((var = find_variable(c_dummy, variable_list)) != NULL)
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
        var = new_variable(c_dummy, val, 1, 1, NULL, NULL);
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

char* spec_join(char** it_list, int n)
     /* replaces variable in table(...) by original string */
{
  char rout_name[] = "spec_join";
  int j;
  char** p;
  struct variable* var;
  *c_join = '\0';
  if (n > 0)
    {
     p = (char**) mymalloc(rout_name,n*sizeof(char*));
     for (j = 0; j < n; j++) p[j] = it_list[j];
     for (j = 0; j < n; j++) 
        if (strcmp(p[j], "table") == 0 && j+3 < n
	    && (var = find_variable(p[j+2], variable_list)) != NULL)
	  p[j+2] = var->string;
     for (j = 0; j < n; j++) strcat(c_join, p[j]);
     free(p);
    }
  return c_join;
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

int set_enable(char* type, struct in_cmd* cmd)
{
  char* name;
  struct command_parameter* cp;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  struct sequence* sequ;
  struct node* nodes[2];
  struct node* c_node;
  int pos, k, n, status, count = 0;
  pos = name_list_pos("sequence", nl);
  if(nl->inform[pos]) /* sequence specified */
    {
     cp = cmd->clone->par->parameters[pos];
     if ((n = name_list_pos(cp->string, sequences->list)) >= 0)
       sequ = sequences->sequs[n];
     else
       {
        warning(cp->string," :sequence not found, skipped");
        return 0;
       }
    }
  else sequ = current_sequ;
  if (sequ->ex_start == NULL)
    {
     warning(sequ->name," :sequence not USEed, skipped");
     return 0;
    }
  pos = name_list_pos("status", nl);
  if (pos > -1 && nl->inform[pos])  /* parameter has been read */
    {
     name = pl->parameters[pos]->string;
     status = strcmp(name, "on") == 0 ? 1 : 0;
    }
  else status = 1;
  pos = name_list_pos("range", nl);
  if (pos > -1 && nl->inform[pos])  /* parameter has been read */
    {
     name = pl->parameters[pos]->string;
     if ((k = get_ex_range(name, sequ, nodes)) == 0)
  	  {
  	   nodes[0] = NULL; nodes[1] = NULL;
  	  }
    }
  else
    {
     nodes[0] = sequ->ex_start; nodes[1] = sequ->ex_end;
    }
  c_node = nodes[0];
  while (c_node)
    {
     if (strstr(c_node->base_name, type) &&
         pass_select(c_node->p_elem->name, cmd->clone) != 0)
       {
	c_node->enable = status; count++;
       }
     if (c_node == nodes[1]) break;
     c_node = c_node->next;
    }
  return count;
}

void set_range(char* range, struct sequence* sequ)
{
  struct node* nodes[2];
  current_sequ->range_start = current_sequ->ex_start;
  current_sequ->range_end = current_sequ->ex_end;
  if (get_ex_range(range, sequ, nodes) == 0) return;
  current_sequ->range_start = nodes[0];
  current_sequ->range_end = nodes[1];
}

void set_selected_columns(struct table* t, struct command_list* select)
{
  int i, j, pos, k, n = 0;
  char* p;
  struct name_list* nl;
  struct command_parameter_list* pl;
  if (par_present("column", NULL, select) == 0)
    {
     for (j = 0; j < t->num_cols; j++)  /* select all columns */
          t->col_out->i[j] = j;
     t->col_out->curr = t->num_cols;
    }
  else
    {
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
    }
  t->col_out->curr = n;
}

void set_selected_elements()
{
  struct name_list* nl;
  struct command* cd;
  struct command_parameter_list* pl;
  struct element* el;
  int i, j, pos, slice;
  selected_elements->curr = 0;
  selected_elements->list->curr = 0;
  for (j = 0; j < element_list->curr; j++)
    {
     el = element_list->elem[j];
     for (i = 0; i < slice_select->curr; i++)
       {
        cd = slice_select->commands[i];
        nl = cd->par_names;
        pos = name_list_pos("slice", nl);
        pl = cd->par;
        if (pos > -1 && nl->inform[pos])  /* parameter has been read */
          slice = pl->parameters[pos]->double_value;
        else slice = 1;
        if (pass_select(el->name, cd) != 0)
          {
           if ((pos = name_list_pos(el->name, selected_elements->list)) > -1)
	     {
	      if (selected_elements->list->inform[pos] < slice)
		selected_elements->list->inform[pos] = slice;
	     }
           else add_to_el_list(el, slice, selected_elements, 0);
           break;
          }
       }
    }
}

void set_selected_errors()
{
  int i, flag;
  if ((flag = 
      get_select_ex_ranges(current_sequ, error_select, selected_ranges)) != 0)
    {
     for (i = 0; i < selected_ranges->curr; i++)
          selected_ranges->nodes[i]->sel_err = 1;
    }
}

void set_selected_rows(struct table* t, struct command_list* select)
{
  int i, j;
  c_range_start = get_node_count(current_sequ->range_start);
  c_range_end = get_node_count(current_sequ->range_end);
  get_select_t_ranges(select, t);
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

void set_twiss_deltas(struct command* comm)
{
  char* string;
  int i, k = 0, n = 0, pos;
  double s, sign = one, ar[3];
  struct name_list* nl = comm->par_names;
  pos = name_list_pos("deltap", nl);
  twiss_deltas->curr = 1;
  twiss_deltas->a[0] = zero;
  if ((pos = name_list_pos("deltap", nl)) >= 0 && nl->inform[pos] 
           && (string = comm->par->parameters[pos]->string) != NULL)
    {
     pre_split(string, c_dummy, 0);
     mysplit(c_dummy, tmp_p_array);
     while (k < tmp_p_array->curr)
       {
	for (i = k; i < tmp_p_array->curr; i++)
	   if (*tmp_p_array->p[i] == ':') break;
        ar[n++] = double_from_expr(tmp_p_array->p, k, i-1);
	k = i + 1;
       }
     if (n == 1) twiss_deltas->a[0] = ar[0];
     else  /* there is a range given - fill array */
       {
	if (n == 2) ar[n++] = ar[1] - ar[0];
        if (ar[2] == zero) twiss_deltas->a[0] = ar[0];
        else if (ar[2] * (ar[1] - ar[0]) < zero)
	  warning("illegal deltap range ignored:", string);
        else
	  {
	   twiss_deltas->a[0] = ar[0];
	   if (ar[2] < zero) sign = -sign;
	   for (s = sign * (ar[0] + ar[2]); 
                s <= sign * ar[1]; s+= sign * ar[2])
	     {
	      if (twiss_deltas->curr == twiss_deltas->max)
		{
		 sprintf(c_dummy, "%d values", twiss_deltas->max);
		 warning("deltap loop cut at", c_dummy);
                 break;
		}
              twiss_deltas->a[twiss_deltas->curr] 
                = twiss_deltas->a[twiss_deltas->curr-1] + ar[2];
              twiss_deltas->curr++;
	     }
	  }
       }
    }
}

void set_sub_variable(char* comm, char* par, struct in_cmd* cmd)
{
  char* p;
  struct element* el;
  struct command *command, *keep_beam = current_beam;
  int end, start = cmd->decl_start;
  int exp_type = loc_expr(cmd->tok_list->p, cmd->tok_list->curr,
                          start, &end);
  double val = 0;
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

void set_sector()
{
  int i, flag;
  if (sector_select->curr == 0) reset_sector(current_sequ, 1);
  else 
    {
     sector_ranges->curr = 0; sector_ranges->list->curr = 0; 
     if ((flag = 
      get_select_ex_ranges(current_sequ, sector_select, sector_ranges)) != 0)
       {
        for (i = 0; i < sector_ranges->curr; i++)
              sector_ranges->nodes[i]->sel_sector = 1;
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

void store_beta0(struct in_cmd* cmd)
{
  int k = cmd->decl_start - 1;
  if (k == 0) warning("beta0 without label:", "ignored");
  else
    {
     cmd->clone_flag = 1; /* do not delete */
     add_to_command_list(cmd->tok_list->p[0], cmd->clone, beta0_list, 0);
    }
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
  n = mysplit(work, toks);
  if (n < 6 || *toks->p[1] != ':') fatal_error("illegal command:", cmd_string);
  if (defined_commands->curr == defined_commands->max)
      grow_command_list(defined_commands);
  cmd = defined_commands->commands[defined_commands->curr++] =
      new_command(toks->p[0], 40, 40, toks->p[2], toks->p[3], 
                  atoi(toks->p[4]), atoi(toks->p[5]));
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
  int i, j, jj, k, dummy, type, s_start, s_end, ss_end;
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
               if ((dummy = loc_expr(toks, end+1, start, &k)) > 1)
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
     case 3:
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
     cp->expr = NULL;
     cp->double_value = val;
    }
}

void store_comm_par_vector(char* parameter, double* val, struct command* cmd)
{
  struct command_parameter* cp;
  int i;
  if ((i = name_list_pos(parameter, cmd->par_names)) > -1)
    {
     cp = cmd->par->parameters[i];
     if (cp->double_array != NULL)
       {
	copy_double(val, cp->double_array->a, cp->double_array->curr);
        if (cp->expr_list != NULL) 
            cp->expr_list = delete_expr_list(cp->expr_list);
       }
    }
}

void store_savebeta(struct in_cmd* cmd)
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int pos;
  char* name = NULL;
  struct command* comm;
  if (log_val("clear", cmd->clone)) 
    {
     delete_command_list(savebeta_list);
     savebeta_list = new_command_list(10);
     delete_command_list(beta0_list);
     beta0_list = new_command_list(10);
    }
  else if (nl->inform[name_list_pos("place", nl)] == 0)
  warning("savebeta without place:", "ignored");
  else
    {
     pos = name_list_pos("label", nl);
     if (nl->inform[pos])  name = pl->parameters[pos]->string;
     else warning("savebeta without label:", "ignored");
     if (name != NULL)
       {
        cmd->clone_flag = 1; /* do not delete */
        if ((comm = find_command(name, beta0_list)))
           remove_from_command_list(name, beta0_list);
        add_to_command_list(permbuff(name), cmd->clone, savebeta_list, 0);
       }
    }
}

void store_select(struct in_cmd* cmd)
{
  char* flag_name;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  struct command_list* scl;
  int pos = name_list_pos("flag", nl);
  if (nl->inform[pos] == 0 ||
      (flag_name = pl->parameters[pos]->string) == NULL)
    {
     warning("no FLAG specified", "ignored");
     return;
    }
  if (strcmp(flag_name, "seqedit") == 0) 
    {
     if (log_val("clear", cmd->clone)) 
       {
        delete_command_list(seqedit_select);
        seqedit_select = new_command_list(10);
       }
     else
       {
	if (seqedit_select->curr == seqedit_select->max) 
            grow_command_list(seqedit_select);
        seqedit_select->commands[seqedit_select->curr++] = cmd->clone;
        cmd->clone_flag = 1; /* do not drop */
       }
    }
  else if (strcmp(flag_name, "error") == 0) 
    {
     if (log_val("clear", cmd->clone)) 
       {
	delete_command_list(error_select);  
        error_select = new_command_list(10);
        selected_ranges->curr = 0;
        selected_ranges->list->curr = 0;
        reset_errors(current_sequ);
       }
     else
       {
	if (error_select->curr == error_select->max) 
            grow_command_list(error_select);
        error_select->commands[error_select->curr++] = cmd->clone;
        cmd->clone_flag = 1; /* do not drop */
       }
    }
  else if (strcmp(flag_name, "makethin") == 0) 
    {
     if (log_val("clear", cmd->clone)) 
       {
        slice_select->curr = 0;
       }
     else
       {
	if (slice_select->curr == slice_select->max) 
            grow_command_list(slice_select);
        slice_select->commands[slice_select->curr++] = cmd->clone;
        cmd->clone_flag = 1; /* do not drop */
       }
    }
  else if (strcmp(flag_name, "save") == 0) 
    {
     if (log_val("clear", cmd->clone)) 
       {
        save_select->curr = 0;
       }
     else
       {
	if (save_select->curr == save_select->max) 
            grow_command_list(save_select);
        save_select->commands[save_select->curr++] = cmd->clone;
        cmd->clone_flag = 1; /* do not drop */
       }
    }
  else if (strcmp(flag_name, "sectormap") == 0) 
    {
     if (sector_ranges == NULL)   sector_ranges = new_node_list(10000);
     if (log_val("clear", cmd->clone)) 
       {
	delete_command_list(sector_select);  
        sector_select = new_command_list(10);
        sector_ranges->curr = 0;
        sector_ranges->list->curr = 0;
       }
     else
       {
	if (sector_select->curr == sector_select->max) 
            grow_command_list(sector_select);
        sector_select->commands[sector_select->curr++] = cmd->clone;
        cmd->clone_flag = 1; /* do not drop */
       }
    }
  else /* store select for all tables */
    {
     if ((scl = find_command_list(flag_name, table_select)) == NULL)
       {
	scl = new_command_list(10);
        add_to_command_list_list(flag_name, scl, table_select);
       }
     if (log_val("clear", cmd->clone)) 
       {
	scl = new_command_list(10);
        add_to_command_list_list(flag_name, scl, table_select);
       }
     else
       {
	if (scl->curr == scl->max) grow_command_list(scl);
        scl->commands[scl->curr++] = cmd->clone;
        cmd->clone_flag = 1; /* do not drop */
       }
    }
}

void store_node_value(char* par, double* value)  
/* stores value for parameter par at current node */
{
  char lpar[NAME_L];
  mycpy(lpar, par);
  if (strcmp(lpar, "chkick") == 0) current_node->chkick = *value;
  else if (strcmp(lpar, "cvkick") == 0) current_node->cvkick = *value;
  else if (strcmp(lpar, "dipole_bv") == 0) current_node->dipole_bv = *value;
  else if (strcmp(lpar, "other_bv") == 0) current_node->other_bv = *value;
  else if (strcmp(lpar, "obs_point") == 0) current_node->obs_point = *value;
  else if (strcmp(lpar, "sel_sector") == 0) current_node->sel_sector = *value;
  else if (strcmp(lpar, "enable") == 0) current_node->enable = *value;
}

void store_node_vector(char* par, int* length, double* vector)  
/* stores vector at node */
{
  char lpar[NAME_L];
  mycpy(lpar, par);
  if (strcmp(lpar, "orbit0") == 0)  copy_double(vector, orbit0, 6);
  else if (strcmp(lpar, "orbit_ref") == 0)
    {
     if (current_node->orbit_ref)
       {
	while (*length > current_node->orbit_ref->max) 
             grow_double_array(current_node->orbit_ref);
       }
     else current_node->orbit_ref = new_double_array(*length);
     copy_double(vector, current_node->orbit_ref->a, *length);
     current_node->orbit_ref->curr = *length;
    }
}

void store_orbit(struct command* comm, double* orbit)
{
  struct name_list* nl = comm->par_names;
  if (nl->inform[name_list_pos("x", nl)]) 
      orbit[0] = command_par_value("x",comm);
  if (nl->inform[name_list_pos("px", nl)]) 
      orbit[1] = command_par_value("px",comm);
  if (nl->inform[name_list_pos("y", nl)]) 
      orbit[2] = command_par_value("y",comm);
  if (nl->inform[name_list_pos("py", nl)]) 
      orbit[3] = command_par_value("py",comm);
  if (nl->inform[name_list_pos("t", nl)]) 
      orbit[4] = command_par_value("t",comm);
  if (nl->inform[name_list_pos("pt", nl)]) 
      orbit[5] = command_par_value("pt",comm);
}

void store_threader(struct in_cmd* cmd)
{
  threader_par = cmd->clone;
  cmd->clone_flag = 1;
  dump_command(threader_par);
}

void string_to_table(char* table, char* name, char* string)
     /* buffers + puts "string"
        at current position in column with name "name".
        The table count is increased separately with "augment_count" */
{
  int pos;
  struct table* t;

  mycpy(c_dummy, table);
  if ((pos = name_list_pos(c_dummy, table_register->names)) > -1)
    t = table_register->tables[pos];
  else return;
  mycpy(c_dummy, name);
  if ((pos = name_list_pos(c_dummy, t->columns)) >= 0
      && t->columns->inform[pos] == 3) 
    {
     mycpy(c_dummy, string);
     if (strcmp(c_dummy, "name") == 0) 
       t->s_cols[pos][t->curr] = tmpbuff(current_node->name);
     else t->s_cols[pos][t->curr] = tmpbuff(c_dummy);
    }
}

int table_length(char* table)
     /* returns no. of rows in table */
{
  int pos;
  int length = 0;
  mycpy(c_dummy, table);
  if ((pos = name_list_pos(c_dummy, table_register->names)) > -1)
     length = table_register->tables[pos]->curr;
  return length;
}

int table_org(char* table)
     /* returns origin: 0  this job, 1 read or unknown */
{
  int pos;
  int org = 1;
  mycpy(c_dummy, table);
  if ((pos = name_list_pos(c_dummy, table_register->names)) > -1)
     org = table_register->tables[pos]->origin;
  return org;
}

void table_range(char* table, char* range, int* rows)
     /* returns first and last row numbers (start=1) in rows
        or 0 if table or range invalid */
{
  int pos;
  struct table* t;

  rows[0] = rows[1] = 0;
  mycpy(c_dummy, table);
  if ((pos = name_list_pos(c_dummy, table_register->names)) > -1)
    {
     t = table_register->tables[pos];
     get_table_range(range, t, rows);
     rows[0]++; rows[1]++;
    }
}


int table_row(struct table* table, char* name)
{
  int i, j, ret = 0;
  for (i = 0; i < table->num_cols; i++) 
     if(table->columns->inform[i] == 3) break;
  if (i < table->num_cols)
    {
     for (j = 0; j < table->curr; j++)
       if (tab_name_code(name, table->s_cols[i][j])) break;
     if (j < table->curr) ret = j;
    }
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
     strcpy(c_dummy, current_variable->string);
     supp_char(',', c_dummy);
     mysplit(c_dummy, tmp_p_array);
     toks = tmp_p_array->p; ntok = tmp_p_array->curr;
     if (ntok > 1)
       {
	if ((pos = name_list_pos(toks[0], table_register->names)) > -1)
	  {
	   table = table_register->tables[pos];
           if ((col = name_list_pos(toks[ntok-1], table->columns)) > -1)
	     {
              if (ntok > 2)  /* find row - else current (dynamic), or 0 */
	        row = table_row(table, toks[1]);
              else if (table->dynamic)  row = table->curr;
              else row = 0;
              val = table->d_cols[col][row];
	     }
	  }
       } 
    }
  return val;
}

int tab_name_code(char* name, char* t_name)
     /* returns 1 if name corresponds to t_name, else 0 */
{
  char tmp[NAME_L];
  char *p, *n = one_string;
  strcpy(tmp, name); 
  if ((p = strstr(tmp, "->")) != NULL)
    {
     *p = '\0'; p = strstr(name, "->"); p++; n = ++p;
    }
  strcat(tmp, ":"); strcat(tmp, n);
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

void track_dynap(struct in_cmd* cmd)
{
  char rout_name[] = "track_dynap";
  int e_flag, flag = 2, izero = 0,
      turns = command_par_value("turns", cmd->clone),
      npart = stored_track_start->curr;
  int *ibuf1, *ibuf2, *ibuf3;
  double orbit[6];
  double *buf1, *buf2, *buf_dxt, *buf_dyt, *buf3, *buf4, *buf5, *buf6, 
    *buf7, *buf8, *buf9, *buf10, *buf11;
  struct table* t;
  if (track_is_on == 0)
    {
     warning("track_dynap: no TRACK command seen yet", "ignored");
     return;
    }
  if (npart == 0)
    {
     warning("track_dynap: no START command seen yet", "ignored");
     return;
    }
  if (turns < 64)
    {
     warning("track_dynap: turns cannot be < 64", "set to 64");
     turns = 64;
    }
  adjust_beam();
  if (probe_beam) probe_beam = delete_command(probe_beam);
  probe_beam = clone_command(current_beam);
  adjust_probe(track_deltap); /* sets correct gamma, beta, etc. */
  adjust_rfc(); /* sets freq in rf-cavities from probe */
  zero_double(orbit0, 6);
  zero_double(oneturnmat, 36);
  if (get_option("onepass") == 0) 
    {
      tmrefo_(&izero,orbit0,orbit,oneturnmat); 
       /* closed orbit and one-turn linear transfer map */
    }
  dynap_tables_create(cmd);
  /* allocate buffers */
  ibuf1 = (int*) mymalloc(rout_name,npart*sizeof(int));
  ibuf2 = (int*) mymalloc(rout_name,npart*sizeof(int));
  ibuf3 = (int*) mymalloc(rout_name, current_sequ->n_nodes*sizeof(int));
  buf1 = (double*) mymalloc(rout_name,npart*sizeof(double));
  buf2 = (double*) mymalloc(rout_name,6*npart*sizeof(double));
  buf_dxt = (double*) mymalloc(rout_name,npart*sizeof(double));
  buf_dyt = (double*) mymalloc(rout_name,npart*sizeof(double));
  buf3 = (double*) mymalloc(rout_name,6*npart*sizeof(double));
  buf4 = (double*) mymalloc(rout_name,36*sizeof(double)); /* eigenvectors */
  buf5 = (double*) mymalloc(rout_name,6*npart*(turns+1)*sizeof(double));
  buf6 = (double*) mymalloc(rout_name, current_sequ->n_nodes*sizeof(double));
  buf7 = (double*) mymalloc(rout_name, turns*sizeof(double));
  buf8 = (double*) mymalloc(rout_name, 6*turns*sizeof(double));
  buf9 = (double*) mymalloc(rout_name, 2*turns*sizeof(double));
  buf10 = (double*) mymalloc(rout_name, turns*sizeof(double));
  buf11 = (double*) mymalloc(rout_name, turns*sizeof(double));
  trrun_(&flag, &turns,orbit0, oneturnmat, ibuf1, ibuf2, buf1, buf2,
	buf_dxt, buf_dyt, buf3, buf4, buf5, &e_flag, ibuf3, buf6);
  t = 
  table_register->tables[name_list_pos("tracksumm", table_register->names)];
  print_table(t);
  if (e_flag)
    {
     warning("track_dynap: particle lost before last turn,", "ignored");
     return;
    }
  dynap_(buf4, buf5, &turns, &npart, buf7, buf8, buf9, buf10, buf11);
  /*
  table_register->tables[name_list_pos("dynapsumm", table_register->names)];
  print_table(t);
  if (get_option("dynap_dump")) dynap_tables_dump();
  */
  /* free buffers */
  free(ibuf1); free(ibuf2); free(ibuf3); free(buf1); free(buf2);
  free(buf_dxt); free(buf_dyt); free(buf3); free(buf4); free(buf5); free(buf6);
  free(buf7); free(buf8); free(buf9); free(buf10); free(buf11);
}

void track_end(struct in_cmd* cmd)
{
  int i;
  struct node* c_node;
  if (track_is_on == 0)
    {
     warning("track_end: no TRACK command seen yet", "ignored");
     return;
    }
  for (i = 0; i < stored_track_start->curr; i++)
    stored_track_start->commands[i] = 
         delete_command(stored_track_start->commands[i]);
  stored_track_start->curr = 0;
  c_node = current_sequ->ex_start;
  while(c_node != NULL) /* clean observation points */
    {
     c_node->obs_point = 0;
     c_node->obs_orbit = delete_double_array(c_node->obs_orbit);
     if (c_node == current_sequ->ex_end)  break;
     c_node = c_node->next;
    }
  track_is_on = 0;
  fprintf(prt_file, "exit TRACK module\n\n");
}

void track_observe(struct in_cmd* cmd)
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  struct node* nodes[2];
  int pos;
  if (track_is_on == 0)
    {
     warning("track_observe: no TRACK command seen yet,", "ignored");
     return;
    }
  pos = name_list_pos("place", nl);
  if (get_ex_range(pl->parameters[pos]->string, current_sequ, nodes))
    {
     nodes[0]->obs_point = ++curr_obs_points;
     nodes[0]->obs_orbit = new_double_array(6);
     nodes[0]->obs_orbit->curr = 6;
     if (get_option("onepass") == 0) 
       {
        tmrefo_(&curr_obs_points,orbit0,nodes[0]->obs_orbit->a,oneturnmat);
        /* closed orbit and one-turn linear transfer map */
       }
    }
  else
    {
     warning("track_observe: unknown place,", "ignored");
     return;
    }
}

void track_pteigen(double* eigen)
{
  int i, j, pos;
  struct table* t;
  double tmp;
  if ((pos = name_list_pos("trackone", table_register->names)) > -1)
    {
     t = table_register->tables[pos];
     if (t->header == NULL)  t->header = new_char_p_array(45);
     sprintf(c_dummy, "@ XC               %%le  %22.12g", orbit0[0]);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     sprintf(c_dummy, "@ PXC              %%le  %22.12g", orbit0[1]);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     sprintf(c_dummy, "@ YC               %%le  %22.12g", orbit0[2]);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     sprintf(c_dummy, "@ PYC              %%le  %22.12g", orbit0[3]);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     sprintf(c_dummy, "@ TC               %%le  %22.12g", orbit0[4]);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     sprintf(c_dummy, "@ PTC              %%le  %22.12g", orbit0[5]);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     tmp = get_value("beam", "ex");
     sprintf(c_dummy, "@ EX               %%le  %22.12g", tmp);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     tmp = get_value("beam", "ey");
     sprintf(c_dummy, "@ EY               %%le  %22.12g", tmp);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     tmp = get_value("beam", "et");
     sprintf(c_dummy, "@ ET               %%le  %22.12g", tmp);
     t->header->p[t->header->curr++] = tmpbuff(c_dummy);
     for (i = 0; i < 6; i++)
       {
        for (j = 0; j < 6; j++)
	  {
	   sprintf(c_dummy, "@ E%d%d              %%le  %22.12g",
                   i+1, j+1, eigen[6*j+i]);
           t->header->p[t->header->curr++] = tmpbuff(c_dummy);
	  }
       }
    }
}

void track_run(struct in_cmd* cmd)
{
  char rout_name[] = "track_run";
  int e_flag, flag = 1, izero = 0, npart = stored_track_start->curr;
  int *ibuf1, *ibuf2, *ibuf3;
  double orbit[6];
  double d_dummy, *buf1, *buf2, *buf_dxt, *buf_dyt, *buf3, *buf4, *buf5, 
    *buf6;
  struct table* t;
  int turns = command_par_value("turns", cmd->clone);
  if (track_is_on == 0)
    {
     warning("track_run: no TRACK command seen yet", "ignored");
     return;
    }
  if (npart == 0)
    {
     warning("track_run: no START command seen yet", "ignored");
     return;
    }
  adjust_beam();
  if (probe_beam) probe_beam = delete_command(probe_beam);
  probe_beam = clone_command(current_beam);
  adjust_probe(track_deltap); /* sets correct gamma, beta, etc. */
  adjust_rfc(); /* sets freq in rf-cavities from probe */
  zero_double(orbit0, 6);
  zero_double(oneturnmat, 36);
  if (get_option("onepass") == 0) 
    {
      tmrefo_(&izero,orbit0,orbit,oneturnmat); 
       /* closed orbit and one-turn linear transfer map */
    }
  track_tables_create(cmd);
  /* allocate buffers */
  ibuf1 = (int*) mymalloc(rout_name,npart*sizeof(int));
  ibuf2 = (int*) mymalloc(rout_name,npart*sizeof(int));
  ibuf3 = (int*) mymalloc(rout_name,current_sequ->n_nodes*sizeof(int));
  buf1 = (double*) mymalloc(rout_name,npart*sizeof(double));
  buf2 = (double*) mymalloc(rout_name,6*npart*sizeof(double));
  buf_dxt = (double*) mymalloc(rout_name,npart*sizeof(double));
  buf_dyt = (double*) mymalloc(rout_name,npart*sizeof(double));
  buf3 = (double*) mymalloc(rout_name,6*npart*sizeof(double));
  buf4 = (double*) mymalloc(rout_name,36*sizeof(double));
  buf5 = &d_dummy;
  buf6 = (double*) mymalloc(rout_name, current_sequ->n_nodes*sizeof(double));
  trrun_(&flag, &turns,orbit0, oneturnmat, ibuf1, ibuf2, buf1, buf2,
	 buf_dxt, buf_dyt, buf3, buf4, buf5, &e_flag, ibuf3, buf6);
  t = 
  table_register->tables[name_list_pos("tracksumm", table_register->names)];
  if (get_option("info"))  print_table(t);
  if (get_option("track_dump")) track_tables_dump();
  /* free buffers */
  free(ibuf1); free(ibuf2); free(ibuf3); 
  free(buf1); free(buf2); free(buf_dxt); free(buf_dyt); free(buf3); 
  free(buf4); free(buf6);
  fprintf(prt_file, "\n*****  end of trrun  *****\n");
}

void track_ripple(struct in_cmd* cmd)
{
  if (track_is_on == 0)
    {
     warning("track_ripple: no TRACK command seen yet", "ignored");
     return;
    }
  puts("entered track_ripple routine");
}

void track_start(struct command* comm)
{
  char name[FNAME_L];
  if (track_is_on == 0)
    {
     warning("track_start: no TRACK command seen yet", "ignored");
     return;
    }
  track_start_cnt++;
  strcpy(name, "start.");
  sprintf(c_dummy, "%d", track_start_cnt);
  strcat(name, c_dummy);
  add_to_command_list(name,comm,stored_track_start,1);
}

void track_tables_create(struct in_cmd* cmd)
{
  int i, j;
  char tab_name[NAME_L];
  struct table* t;
  t = make_table("tracksumm", "tracksumm", tracksumm_table_cols, 
  		 tracksumm_table_types, 2*stored_track_start->curr);
  add_to_table_list(t, table_register);
  if (get_option("onetable"))
    {
     t = make_table("trackone", "trackone", trackone_table_cols, 
  		 trackone_table_types, stored_track_start->curr*TRACK_ROWS);
     add_to_table_list(t, table_register);
    }
  else
    {
     for (i = 0; i < curr_obs_points; i++)
       {
        for (j = 0; j < stored_track_start->curr; j++) /* open tables */
          {
           sprintf(tab_name, "track.obs%04d.p%04d", i+1, j+1);
           t = make_table(tab_name, "trackobs", track_table_cols, 
  		     track_table_types, TRACK_ROWS);
           add_to_table_list(t, table_register);
	  }
       }
    }
}

void track_tables_dump()
{
  int j;
  for (j = 0; j < table_register->names->curr; j++)
    {
     if (strstr(table_register->names->names[j], "track.obs")
         || strcmp(table_register->names->names[j], "trackone") == 0)
       {
        strcpy(l_work, track_filename);
        strcat(l_work, &table_register->names->names[j][5]);
        strcat(l_work, track_fileext);
        out_table("track", table_register->tables[j], l_work);
       }
    }
}

void track_track(struct in_cmd* cmd)
{
  int k=0, pos, one = 1;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;

  if (current_sequ == NULL || current_sequ->ex_start == NULL)
    {
     warning("sequence not active,", "TRACK ignored");
     return;
    }
  if (attach_beam(current_sequ) == 0)
    fatal_error("TRACK - sequence without beam:", current_sequ->name);
  if (track_is_on)
    {
     warning("already inside TRACK command group,", "ignored");
     return;
    }
  track_is_on = 1;
  puts("enter TRACK module");
  if ((k = get_value(current_command->name,"onepass")) != 0)
     fprintf(prt_file, "one pass is on\n");
  set_option("onepass", &k);
  if ((k = get_value(current_command->name,"aperture")) != 0)
     fprintf(prt_file, "aperture tracking is on\n");
  set_option("aperture", &k);
  k = get_value(current_command->name,"dump");
  set_option("track_dump", &k);
  k = get_value(current_command->name,"onetable");
  set_option("onetable", &k);
  track_deltap=get_value(current_command->name,"deltap");
  set_variable("track_deltap", &track_deltap);
  if(track_deltap != 0) fprintf(prt_file, "track_deltap: %f\n",track_deltap);
  curr_obs_points = 1;  /* default: always observe at machine end */
  pos = name_list_pos("file", nl);
  if (nl->inform[pos]) set_option("track_dump", &one);
  if ((track_filename = pl->parameters[pos]->string) == NULL)
    {
     if (pl->parameters[pos]->call_def != NULL)
     track_filename = pl->parameters[pos]->call_def->string;
     else track_filename = permbuff("dummy");
    }
  track_filename = permbuff(track_filename);
  track_fileext = NULL;
  pos = name_list_pos("extension", nl);
  if ((track_fileext = pl->parameters[pos]->string) == NULL)
    {
     if (pl->parameters[pos]->call_def != NULL)
     track_fileext = pl->parameters[pos]->call_def->string;
     if (track_fileext == NULL)  track_fileext = permbuff("\0");
    }
  track_fileext = permbuff(track_fileext);
}

int twiss_input(struct command* tw)
     /* returns -1 if an invalid beta0 given,
        returns -1 if only betx or bety given,
        returns 1 if betx and bety are given,
        or a valid beta0 (which is then loaded), else 0 */
{
  struct name_list* nl = tw->par_names;
  struct command_parameter_list* pl = tw->par;
  struct command* beta;
  int i = -1, ret = 0, pos, sb = 0;
  char* name;
  double val;
  pos = name_list_pos("beta0", nl);
  if (nl->inform[pos] && (name = pl->parameters[pos]->string) != NULL)
    {
     if ((pos = name_list_pos(name, beta0_list->list)) > -1)
       {
	ret = 1;
	beta = beta0_list->commands[pos];
        do
	  {
	   i++;
           if (nl->inform[name_list_pos(nl->names[i], nl)] == 0) /* not read */
	     {
	      if (beta->par->parameters[i]->expr != NULL) 
	        val = expression_value(beta->par->parameters[i]->expr, 2);
              else val = beta->par->parameters[i]->double_value;
              pl->parameters[i]->double_value = val;
              nl->inform[name_list_pos(nl->names[i], nl)] = 1;
	     }
	  }
        while (strcmp(nl->names[i], "energy") != 0);
       }
     else ret = -1;
    }
  if (ret) return ret;
  /* if no beta0 given, betx and bety together set inval */
  if (nl->inform[name_list_pos("betx", nl)]) sb++;
  if (nl->inform[name_list_pos("bety", nl)]) sb++;
  if (sb)
    {
     if (sb < 2)  return -2;
     else         return 1;
    }
  else return 0;
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
   mass is the "nucleon" number times the average nucleon mass "nmass". 
*/
{
  struct name_list* nlc = comm->par_names;
  struct command_parameter_list* plc = comm->par;
  struct command_parameter_list* pl = current_beam->par;
  int pos, lp;
  char* name = blank;
  double energy = 0, beta = 0, gamma = 0, charge = 0, freq0 = 0, bcurrent = 0, 
         npart = 0, mass = 0, pc = 0, ex, exn, ey, eyn, alfa, circ = one, 
         arad = 0, nucleon = 1;
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
     else 
       {
        pos = name_list_pos("nucleon", nlc);
        if (nlc->inform[pos]) nucleon = command_par_value("nucleon", comm);
        else nucleon = command_par_value("nucleon", current_beam);
        mass = nucleon * get_variable("nmass");
       }
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
     if (strcmp(name, "ion") == 0) energy *= fabs(charge);
     if (energy <= mass) fatal_error("energy must be","> mass");
     pc = sqrt(energy*energy - mass*mass);
     gamma = energy / mass;
     beta = pc / energy;
    }
  else if((pos = name_list_pos("pc", nlc)) > -1 && nlc->inform[pos])
    {
     pc = command_par_value("pc", comm);
     if (strcmp(name, "ion") == 0) pc *= fabs(charge);
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
     if (strcmp(name, "ion") == 0) energy *= fabs(charge);
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
  store_comm_par_value("nucleon", nucleon, current_beam);
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
	e_par = e_pl->parameters[pos];
        par = pl->parameters[pos];
        switch (par->type)
          {
          case 0:
          case 1:
          case 2:
            e_par->double_value = par->double_value;
            e_par->expr = clone_expression(par->expr);
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

void update_node_constraints(struct node* c_node, struct constraint_list* cl)
{
  int i, j, k;
  k = 1; set_option("match_local", &k); /* flag */
  if (c_node->cl == NULL) c_node->cl = new_constraint_list(cl->curr);
  for (j = 0; j < cl->curr; j++)
    {
     k = -1;
     for (i = 0; i < c_node->cl->curr; i++)
       {
	if (strcmp(cl->constraints[j]->name, 
                   c_node->cl->constraints[i]->name) == 0) k = i;
       }
     if (k < 0)
       {
	if (c_node->cl->curr == c_node->cl->max) 
           grow_constraint_list(c_node->cl);
        c_node->cl->constraints[c_node->cl->curr++] = cl->constraints[j];
        total_const++;
       }
     else c_node->cl->constraints[k] = cl->constraints[j];
    }
}

void update_sequ_constraints(struct sequence* sequ, struct constraint_list* cl)
{
  int i, j, k;
  if (sequ->cl == NULL) sequ->cl = new_constraint_list(10);
  for (j = 0; j < cl->curr; j++)
    {
     k = -1;
     for (i = 0; i < sequ->cl->curr; i++)
       {
	if (strcmp(cl->constraints[j]->name, 
                   sequ->cl->constraints[i]->name) == 0) k = i;
       }
     if (k < 0)
       {
	if (sequ->cl->curr == sequ->cl->max) 
           grow_constraint_list(sequ->cl);
        sequ->cl->constraints[sequ->cl->curr++] = cl->constraints[j];
        total_const++;
       }
     else sequ->cl->constraints[k] = cl->constraints[j];
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

void use_sequ(struct in_cmd* cmd)
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int pos, lp;
  char* name;
  struct command* keep_beam = current_beam;
  if (sequ_is_on)
     fatal_error("no endsequence yet for sequence:", current_sequ->name);
  pos = name_list_pos("period", nl);
  if (nl->inform[pos] == 0) pos = name_list_pos("sequence", nl);
  if (nl->inform[pos])  /* parameter has been read */
    {
     if (current_range != NULL)
       {
	free(current_range); current_range = NULL;
       }
     name = pl->parameters[pos]->string;
     if ((pos = name_list_pos(name, line_list->list)) > -1) 
       make_sequ_from_line(name); 
     if ((lp = name_list_pos(name, sequences->list)) > -1)
       {
	current_sequ = sequences->sequs[lp];
        if (attach_beam(current_sequ) == 0)
           fatal_error("USE - sequence without beam:", current_sequ->name);
        current_sequ->beam = current_beam;
        pos = name_list_pos("range", nl);
        if (nl->inform[pos])  /* parameter has been read */
           current_range = tmpbuff(pl->parameters[pos]->string);
        expand_curr_sequ(0);        
       }
     else warning("unknown sequence skipped:", name);
    }
  current_beam = keep_beam;
}

double get_variable(char* name)
{
  char comm[NAME_L];
  char par[NAME_L];
  double val = zero;
  struct variable* var;
  struct element* el;
  struct command* cmd;
  char *p, *n = c_dummy, *q = comm;
  mycpy(c_dummy, name);
  if ((p = strstr(c_dummy, "->")) == NULL) /* variable */
    {
     if ((var = find_variable(c_dummy, variable_list)) != NULL)
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

void vector_to_table(char* table, char* col, int* nval, double* vals)
     /* puts nval values of array vals at the current line into columns
        starting with column whose name is in "col";
        The table count is increased separately with "augment_count" */
{
  int j, pos, c_pos, last;
  struct table* t;

  mycpy(c_dummy, table);
  if ((pos = name_list_pos(c_dummy, table_register->names)) > -1)
    t = table_register->tables[pos];
  else return;
  mycpy(c_dummy, col);
  if ((c_pos = name_list_pos(c_dummy, t->columns)) > -1)
  last = mymin(c_pos + *nval, t->num_cols);
  for (j = c_pos; j < last; j++)
     if (t->columns->inform[j] < 3) t->d_cols[j][t->curr] = vals[j-c_pos];
}

void warning(char* t1, char* t2) 
{
  if (get_option("warn")) printf("++++++ warning: %s %s\n",t1,t2);
}
