#include "madx.h"

static struct macro*
clone_macro(struct macro* org)
{
  int i;
  struct macro* clone = new_macro(org->n_formal, org->body->curr, org->tokens->curr);
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

static struct macro*
delete_macro(struct macro* macro)
{
  const char *rout_name = "delete_macro";
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

static void
grow_macro_list(struct macro_list* p)
{
  const char *rout_name = "grow_macro_list";
  p->max *= 2;
  p->macros = myrecalloc(rout_name, p->macros, p->curr * sizeof *p->macros, p->max * sizeof *p->macros);
}

#if 0 // not used...
static void
dump_macro(struct macro* m)
{
  fprintf(prt_file, "name: %s\n", m->name);
  if (m->formal != NULL) dump_char_p_array(m->formal);
  dump_char_array(m->body);
  if (m->tokens != NULL) dump_char_p_array(m->tokens);
}
#endif

#if 0 // not used...
static void
dump_macro_list(struct macro_list* ml)
{
  int i;
  puts("++++++ dump of macro list");
  for (i = 0; i < ml->curr; i++)
    dump_macro(ml->macros[i]);
}
#endif

void
disable_line( /* prevents line from further expansion by "use" */
  char* name, struct macro_list* nll)
{
  int pos;
  if ((pos = name_list_pos(name, nll->list)) > -1) 
   {
     nll->macros[pos]->dead = 1;
   }
}

#if 0 // not used
static void
remove_from_macro_list( /* removes macro from alphabetic macro list */
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
#endif

// public interface

struct macro*
new_macro(int n_formal, int length, int p_length)
{
  const char *rout_name = "new_macro";
  struct macro* m = mycalloc(rout_name, 1, sizeof *m);
  strcpy(m->name, "macro");
  m->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", m->name);
  if ((m->n_formal  = n_formal) > 0) m->formal = new_char_p_array(n_formal);
  if (p_length > 0) m->tokens = new_char_p_array(p_length);
  ++length;
  m->body = new_char_array(length + 1);
  return m;
}

struct macro_list*
new_macro_list(int length)
{
  const char *rout_name = "new_macro_list";
  struct macro_list* nll = mycalloc(rout_name, 1, sizeof *nll);
  strcpy(nll->name, "macro_list");
  nll->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", nll->name);
  nll->list = new_name_list(nll->name, length);
  nll->macros = mycalloc(rout_name, length, sizeof *nll->macros);
  nll->max = length;
  return nll;
}

int
make_macro(char* statement)
  /* makes a new macro from input command, stores name in name list */
{
  struct macro* m;
  char** toks = tmp_l_array->p;
  int i, n, rs, re, start_2;
  int len = strlen(statement);
  while(len >= aux_buff->max) grow_char_array(aux_buff);
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
    m->original = new_char_array(len);
  strcpy(m->original->c, statement);
  add_to_macro_list(m, macro_list);
  return 0;
}
void
save_macros2file(const char* fname){
  FILE *fptr;

    // opening file in writing mode
  fptr = fopen(fname, "w");
  for(int i=0; i<macro_list->curr; i++){
 
    fprintf(fptr, "%s \n", macro_list->macros[i]->original->c);
  }

  fclose(fptr);

}

void
exec_macro(struct in_cmd* cmd, int pos)
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
    for (i = 0; i < any; i++) {
      int len = strlen(toks[rs+i]);
      sum += len;
    }
    {
      int len = strlen(pro->buffers[level]->c_a->c);
      while (l_wrk->max < len+sum)
        grow_char_array(l_wrk);
    }
    for (i = 0; i < any; i++)
    {
      myrepl(macro_list->macros[pos]->formal->p[i], toks[rs+i],
             pro->buffers[level]->c_a->c, l_wrk->c);
      mystrcpy(pro->buffers[level]->c_a, l_wrk->c);
    }
  }

  if (get_option("echomacro")) {
    printf("=== echoing exec %s", macro_list->macros[pos]->name);
    if (macro_list->macros[pos]->n_formal > 0) {
      printf("(");
      for (i=3; i<cmd->tok_list->curr-2; i++) printf("%s,", cmd->tok_list->p[i]);
      printf("%s)", cmd->tok_list->p[i]);
    }
    printf("\n");
    puts(pro->buffers[level]->c_a->c);
    printf("=== end of echoing %s\n", macro_list->macros[pos]->name);
  }

  pro_input(pro->buffers[level]->c_a->c);
  pro->curr--;
}

void
add_to_macro_list( /* adds macro to alphabetic macro list */
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

void
replace_lines(struct macro* org, int replace, char** reps)
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
      line_buffer->p[line_buffer->curr++] = tmpbuff(p);
    }
  }
  delete_macro(line);
}


