#include "madx.h"

static double
act_value(int pos, const struct name_list* chunks)
  /* returns the actual value of a variable, element, or command parameter */
{
  const char* name = chunks->names[pos];
  const char *n = name, *p;
  char comm[NAME_L];
  char par[NAME_L];
  char *q = comm;
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
    if (strncmp(comm, "beam", 4) == 0) {
      // LD 2017.04.11: remove current_beam side effect!
      struct command* cmd2 = cmd = find_command("default_beam", beam_list);
      if ((p = strchr(comm, '%')) != NULL) {
        if ((cmd = find_command(++p, beam_list)) == NULL)
          cmd = cmd2;
      }
      val = command_par_value(par, cmd);
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

static int
simple_logic_expr(int nit, char* toks[])
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
      val1 = polish_value(deco, join(&toks[l1_start], l1_end + 1 - l1_start));
    else return -1;
  }
  else  val1 = simple_double(toks, l1_start, l1_end);
  if ((t2 = loc_expr(toks, nit, l2_start, &l2_end)) == 0) return -1;
  if (t2 == 2)
  {
    if (polish_expr(l2_end + 1 - l2_start, &toks[l2_start]) == 0)
      val2 = polish_value(deco, join(&toks[l2_start], l2_end + 1 - l2_start));
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

static int
logic_expr(int nit, char* toks[])
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

// public interface

void
deco_init(void)
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

int
act_special(int type, char* statement)
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
  cnt_1 = rs+1; start_2 = rs + 1;
  mystrcpy(loc_buff, statement); loc_buff->c[re] =  '\0';
  while(aux_buff->max < cnt_1) grow_char_array(aux_buff);
  strncpy(aux_buff->c, statement, cnt_1); aux_buff->c[cnt_1] = '\0';
  switch (type)
  {
    case 1:  /* if */
      pro->buffers[level]->flag = 0;
      __attribute__ ((fallthrough));
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

void
pro_input(char* statement)
  /* processes one special (IF() etc.), or one normal statement after input */
{
  int type, code, nnb, ktmp;
  char* sem;
  int rs, re, start = 0, l = strlen(statement);

  clearerrorflag(); /*reset global error flag */

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
              break;
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
      if (sem > &statement[start]) /* skip empty ';' */
      {
        *sem = '\0';
        this_cmd = new_in_cmd(400);
        pre_split(&statement[start], work, 1);
        check_table(work->c); check_tabindex(work->c); check_tabstring(work->c);
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
                break;
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
      start = sem - statement;
      if (start < l)
      {
        if ((nnb = next_non_blank_pos(sem)) < 0)  start = l;
        else start += nnb;
      }
    }
  }
}

void
process(void)  /* steering routine: processes one command */
{
  int pos;
  char* name;
  struct element* el;

  if (this_cmd == NULL) return;

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
      this_cmd = buffer_in_cmd(this_cmd);
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
        if ((pos = name_list_pos(name, sequences->list)) < 0)
          enter_element(this_cmd);
        else
        {
          this_cmd->cmd_def = find_command("sequence", defined_commands);
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
              join_b(this_cmd->tok_list->p, this_cmd->tok_list->curr));
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
    this_cmd = buffer_in_cmd(this_cmd);
  }
}

int
polish_expr(int c_item, char** item)   /* split input */
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

double
polish_value(struct int_array* deco, char* expr_string)
  /* coded input (see below) */
  /* description see polish_expression */
{
  int i, k, kc, c_stack = -1;
  double stack[MAX_ITEM];
  char tmp[20];

  if (++polish_cnt > MAX_LOOP)
    fatal_error("circular call in expression", expr_string);
  stack[0] = 0;
  for (i = 0; i < deco->curr; i++)   /* decoding loop */
  {
    k = deco->i[i];
    if ( k < 5)     /* operator */
    {
      if (c_stack < 0)
      {
        fatal_error("stack underflow in expression:", expr_string);
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
            warning("division by zero, result set to zero, expr:", expr_string);
            stack[c_stack] = 0.0;
            break;
          }
          stack[c_stack] /= stack[c_stack+1];
          break;
        case 4:
          stack[c_stack] = pow(stack[c_stack],stack[c_stack+1]);
          break;
        default:
          fatal_error("illegal operator, Polish, expr.:", expr_string);
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
            case 18: /* function "exist" */
              continue; /* value in stack not changed */
              break;
            case 19:
              stack[c_stack] = floor(stack[c_stack]);
              break;
            case 20:
              stack[c_stack] = ceil(stack[c_stack]);
              break;
            case 21:
              stack[c_stack] = rint(stack[c_stack]);
              break;
            case 22: {
              double int_part;
              stack[c_stack] = modf(stack[c_stack], &int_part);
              } break;
            case 23:
              stack[c_stack] = erf(stack[c_stack]);
              break;
            case 24:
              stack[c_stack] = erfc(stack[c_stack]);
              break;
            case 25: {
                double x = stack[c_stack];
                stack[c_stack] = fabs(x) < 1e-12 ? 1.0 : sin(x)/x;
              }
              break;

            default:
              fatal_error("polish_value: illegal function in expr:",
                          expr_string);
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

void
print_value(struct in_cmd* cmd)
{
  char** toks = &cmd->tok_list->p[cmd->decl_start];
  int n = cmd->tok_list->curr - cmd->decl_start;
  int start=0, end, j;
  char* val_spec;
  char* val_format;
  char* val_exprstring;
  double val_value;

  while (start < n)
  {
    for (j = start; j < n; j++)
     {
      if (*toks[j] == ',')
       {
         break;
       }
     }

    if (loc_expr(toks, j, start, &end) > 0)
    {
      int nitem = end - start + 1;
      if (polish_expr(nitem, &toks[start]) == 0)
       {

        val_format = v_format("%S = %F ;\n");
        val_spec = spec_join(&toks[start], nitem);
        /*need to copy the string otherwise it is modified with the value extraction with Intel compiler*/
        int bufflen = strlen(val_spec)+1;
        char  lbuff_local_tmp__[bufflen < 8192 ? bufflen : 1];
        char *lbuff = bufflen < 8192 ? lbuff_local_tmp__ : mymalloc_atomic("print_value",bufflen);

        /*val_spec_cpy = mymalloc_atomic("print_value",(strlen(val_spec)+1) * sizeof(char));*/

        strcpy(lbuff, val_spec);

        val_exprstring = join(&toks[start], nitem); /*val_spec is modified by this line because it works on the same buffer as spec_join*/
        val_value = polish_value(deco, val_exprstring);

        fprintf(prt_file, val_format , lbuff,  val_value);

        if (lbuff != lbuff_local_tmp__) myfree("print_value",lbuff);

       }
      else
      {
        warning("invalid expression:", spec_join(&toks[start], nitem));
        return;
      }
      start = end+1;
      if (start < n-1 && *toks[start] == ',') start++;
    }
    else
    {
      warning("invalid expression:", spec_join(&toks[start], n-start));
      return;
    }
  }
}
