#include "madx.h"

static int
is_operand(char c) {
  return (isalnum(c) || c == '_' || c == '.');
}

static int
is_operator(char c) {
  return strchr(op_string, c) != 0;
}

static int
is_expr_start(char c) {
  return strchr("-+(",c) || is_operand(c);
}

#if 0 // not used...
/* combine two parameters using compound expression */
static struct expression*
comb_param(struct command_parameter* param1, char* op, struct command_parameter* param2)
{
  return compound_expr(param1->expr,param1->double_value,op,param2->expr,param2->double_value);
}
#endif

static double
combine_expr_expr(struct expression* exp1, const char* oper,
                  struct expression* exp2, struct expression** comb_exp)
{
  strcpy(c_dum->c, exp1->string);
  strcat(c_dum->c, oper);
  strcat(c_dum->c, exp2->string);
  mysplit(c_dum->c, tmp_p_array);
  *comb_exp = make_expression(tmp_p_array->curr, tmp_p_array->p);
  return expression_value(*comb_exp, 2);
}

static double
combine_expr_val(struct expression* exp1, const char* oper, double val2, struct expression** comb_exp)
{
  strcpy(c_dum->c, exp1->string);
  sprintf(aux_buff->c, "%.12g", val2);
  strcat(c_dum->c, oper);
  strcat(c_dum->c, aux_buff->c);
  mysplit(c_dum->c, tmp_p_array);
  *comb_exp = make_expression(tmp_p_array->curr, tmp_p_array->p);
  return expression_value(*comb_exp, 2);
}

static double
combine_val_expr(double val1, const char* oper, struct expression* exp2, struct expression** comb_exp)
{
  sprintf(c_dum->c, "%.12g", val1);
  strcat(c_dum->c, oper);
  strcat(c_dum->c, exp2->string);
  mysplit(c_dum->c, tmp_p_array);
  *comb_exp = make_expression(tmp_p_array->curr, tmp_p_array->p);
  return expression_value(*comb_exp, 2);
}

// public interface

struct expression*
make_expression(int n, char** toks)
  /* makes an expression from a list of tokens */
{
  struct expression* expr = NULL;

  if (polish_expr(n, toks) == 0)
    expr = new_expression(join_b(toks, n), deco);
  else warning("Invalid expression starting at:", join_b(toks, n));
  return expr;
}

struct expression*
new_expression(const char* in_string, struct int_array* polish)
{
  const char *rout_name = "new_expression";
  struct expression* ex = mycalloc(rout_name, 1, sizeof *ex);
  strcpy(ex->name, "expression");
  ex->stamp = 123456;
  int len = strlen(in_string)+1;
  ex->string = mymalloc_atomic(rout_name, len * sizeof *ex->string);
  strcpy(ex->string, in_string);
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", ex->name);
  if (polish != NULL) {
    ex->polish = new_int_array(polish->curr);
    ex->polish->curr = polish->curr;
    for (int j = 0; j < polish->curr; j++) ex->polish->i[j] = polish->i[j];
  }
  return ex;
}

struct expr_list*
new_expr_list(int length)
{
  const char *rout_name = "new_expr_list";
  struct expr_list* ell = mycalloc(rout_name, 1, sizeof *ell);
  strcpy(ell->name, "expr_list");
  ell->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", ell->name);
  ell->list  = mycalloc(rout_name, length, sizeof *ell->list);
  ell->max = length;
  return ell;
}

double
expression_value(struct expression* expr, int flag) /* recursive */
  /* returns the value of an expression if valid, else zero */
{
  double val = zero;
  if (expr->status == 0 || flag == 2)
  {
    if (expr->polish != NULL)
    {
      val = expr->value = polish_value(expr->polish, expr->string);
      expr->status = 1;
    }
  }
  else val = expr->value;
  return val;
}

struct expression*
clone_expression(struct expression* p)
{
  struct expression* clone;
  if (p == NULL) return NULL;
  clone = new_expression(p->string, p->polish);
  clone->status = p->status;
  clone->value = p->value;
  return clone;
}

struct expr_list*
clone_expr_list(struct expr_list* p)
{
  int i;
  struct expr_list* clone;
  if (p == NULL)  return NULL;
  clone = new_expr_list(p->curr);
  for (i = 0; i < p->curr; i++) clone->list[i] = clone_expression(p->list[i]);
  clone->curr = p->curr;
  return clone;
}

struct expression*
delete_expression(struct expression* expr)
{
  const char *rout_name = "delete_expression";
  if (expr == NULL) return NULL;
  if (stamp_flag && expr->stamp != 123456)
    fprintf(stamp_file, "d_ex double delete --> %s\n", expr->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", expr->name);
  if (expr->polish != NULL) expr->polish = delete_int_array(expr->polish);
  if (expr->string != NULL) myfree(rout_name, expr->string);
  myfree(rout_name, expr);
  return NULL;
}

struct expr_list*
delete_expr_list(struct expr_list* exprl)
{
  const char *rout_name = "delete_expr_list";
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

void
grow_expr_list(struct expr_list* p)
{
  const char *rout_name = "grow_expr_list";
  p->max *= 2;
  if (p->max == 0) p->max++;
  p->list = myrecalloc(rout_name, p->list, p->curr * sizeof *p->list, p->max * sizeof *p->list);
}

void
dump_expression(struct expression* ex)
{
  ex->value = expression_value(ex, 2);
  fprintf(prt_file, v_format("expression: %s :: value: %F\n"),
          ex->string, ex->value);
}

double
expr_combine(struct expression* exp1, double val1, const char* oper,
             struct expression*  exp2, double val2,
             struct expression** exp_comb)
{
  double val = 0;

  if (exp1 == NULL && exp2 == NULL)
  {
    *exp_comb = NULL;
    switch(oper[1])
    {
      case '+': val = val1 + val2; break;
      case '-': val = val1 - val2; break;
    }
  }
  else if(exp1 == NULL) val = combine_val_expr (val1, oper, exp2, exp_comb);
  else if(exp2 == NULL) val = combine_expr_val (exp1, oper, val2, exp_comb);
  else                  val = combine_expr_expr(exp1, oper, exp2, exp_comb);
  return val;
}

void
update_vector(struct expr_list* ell, struct double_array* da)
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

void
fill_expr_list(char** toks, int s_start, int s_end, struct expr_list* p)
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

void
fill_expr_var_list(struct el_list* ell, struct expression* expr, struct var_list* varl)
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

double
double_from_expr(char** toks, int s_start, int s_end)
  /* returns the value of an expression if valid, else INVALID */
{
  int end, nitem = s_end + 1;
  int type = loc_expr(toks, nitem, s_start, &end);
  if (type == 1) /* simple number */
    return simple_double(toks, s_start, end);
  else if (polish_expr(end + 1 - s_start, &toks[s_start]) == 0)
    return polish_value(deco, join( &toks[s_start],end + 1 - s_start) );
  else return INVALID;
}

int
loc_expr(char** items, int nit, int start, int* end)
  /* Returns the type and end of an expression, or 0 if illegal */
{
  char c;
  int i, e_type = 1, par_level = 0, ltog = -1;
  *end = start - 1;
  if (nit > start && is_expr_start(*items[start])) {
    for (i = start; i < nit; i++) {
      c = *items[i];
      if (c == '(')  { par_level++; e_type = 2; }
      else if (c == ')') {
        if (par_level == 0) return 0;
        par_level--; ltog = 0;
      }
      else if (par_level == 0) {
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

int
scan_expr(int c_item, char** item)   /* split input */

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
        if (strcmp("exist", functs[lp]) == 0  /* special function */
            && i+3 < c_item && *item[i+1] == '(' && *item[i+3] == ')')
        {
          if (find_variable(item[i+2], variable_list) == NULL)
            strcpy(item[i+2], "0");
          else strcpy(item[i+2], "1");
        }
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
      cat_doubles->a[l_cat++] = myatof(item[i]);
    }
    else if (is_operator(c))
    {
      if (l_cat == cat->max ) grow_int_array(cat);
      if (l_cat == oper->max) grow_int_array(oper);
      cat->i[l_cat] = 4;
      oper->i[l_cat++] = str_pos(op_string, c);
      /* oper->i[l_cat++] = (int)strchr(op_string, c) -(int)op_string;*/
    }
    else return 2;  /* illegal character */
  }
  if (level != 0)  return 1;  /* unclosed parentheses */
  cat->curr = l_cat;
  return 0;
}

struct expression*
compound_expr(struct expression* e1, double v1, const char* oper, struct expression* e2, double v2, int parentheses)
/* make one out of two expressions, using oper to connect them
 hbu 9/2005 moved from madxn.c to makethin.c as only used here
 and increased precision   sprintf(tmp, "%e"  ->   sprintf(tmp, "%.14g" */
{
   char lb[] = "(", rb[] = ")";
  char** toks = tmp_l_array->p;
  struct expression* expr = NULL;
  char tmp[30], op[30];
  int n;
  strcpy(op, oper);
  if(parentheses== 0) { //In this way if no parenthesis are chosen they
    lb[0] =  (char) 0;
    rb[0] =  (char) 0;
    if(e2 && e2->string[0]=='-') op[0] =' ';
  }
  if (e1 != NULL || e2 != NULL)
  {
    if (e1 != NULL)
    {
      if (e2 != NULL)
      {
        toks[0] = lb;
        toks[1] = e1->string; toks[2] = rb;
        toks[3] = op;
        toks[4] = lb; 
        toks[5] = e2->string; toks[6] = rb;
      }
      else
      {
        sprintf(tmp, "%.14g", v2); /* hbu */
        toks[0] = lb; 
        toks[1] = e1->string; toks[2] = rb;
        toks[3] = op;
        toks[4] = lb; 
        toks[5] = tmp; toks[6] = rb;
      }
    }
    else
    {
      sprintf(tmp, "%.14g", v1);  /* hbu */
      toks[0] = lb; 
      toks[1] = tmp; toks[2] = rb;
      toks[3] = op;
      toks[4] = lb; 
      toks[5] = e2->string; toks[6] = rb;
    }
    join(toks, 7);
    pre_split(c_join->c, l_wrk, 0);
    n = mysplit(l_wrk->c, tmp_l_array);
    expr = make_expression(n, toks);
  }
  return expr;
}

/* scale an expression by a number - or leave it NULL */
struct expression*
scale_expr(struct expression* expr,double scale)
{
  if (expr) return compound_expr(expr,0,"*",NULL,scale,1);
  return NULL;
}

