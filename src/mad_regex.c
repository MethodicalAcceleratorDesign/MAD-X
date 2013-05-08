#include "madx.h"

// private types

struct reg_token
{
    int type;         /* 0: empty, 1: char, 2: simple, 3: list, 4: dot */
    int rep;          /* 0 no, 1 yes (including 0 repetitions) */
    int rep_cnt;      /* no. of repitions at current match */
    int rep_max;      /* max. no. of repitions at current match */
    int invert;       /* 0 no, 1 yes (for type 3 only) */
    int match_end;    /* 0 no, 1 yes (only last token for '$' */
    char c;           /* for type = 1 */
    struct r_char_array* simple;
    struct r_char_array* list;
    struct reg_token* next;
    struct reg_token* previous;
};

// note: very similar to char_array from mad_array.h

struct r_char_array             /* dynamic array of char for regex */
{
    int  max,                   /* max. array size */
      curr;                     /* current occupation */
    char* chars;
};

// forward declarations

// static struct reg_token* make_list(struct reg_token*, char*, int, int);

// private funtions (utils)

static int
char_count(char c, char* s)
{
  int n = 0, i = 0;
  while (*s != '\0')
  {
    if (*s++ == c)
    {
      i++;
      if (n < i) n = i;
    }
    else i = 0;
  }
  return n;
}

#if 0 // not used
static void
dump_tokens(struct reg_token* tk)
{
  while (tk)
  {
    printf("type: %d  rep = %d  invert = %d  max = %d  dollar = %d\n",
           tk->type, tk->rep, tk->invert, tk->rep_max, tk->match_end);
    tk = tk->next;
  }
}
#endif

static void
grow_r_char_array(struct r_char_array* a)
{
  const char *rout_name = "grow_r_char_array";
  char* loc = a->chars;
  a->max *= 2;
  a->chars = mymalloc_atomic(rout_name, a->max * sizeof *a->chars);
  strncpy(a->chars, loc, a->curr);
  myfree(rout_name, loc);
}

static struct reg_token*
make_dot(struct reg_token* rt)
{
  const char *rout_name = "make_dot";
  struct reg_token *rn = rt;
  if (rn->type != 0) rn = mycalloc(rout_name, 1, sizeof *rn);
  rn->type = 4;
  if (rn != rt) {
    rt->next = rn; rn->previous = rt;
  }
  return rn;
}

static void
fill_list(char s, char e, struct r_char_array* a)
{
  int i;
  while ((int)s <= (int)e)
  {
    if (a->curr == a->max) grow_r_char_array(a);
    a->chars[a->curr++] = s;
    i = (int)s + 1; s = (char) i;
  }
}

static struct reg_token*
make_list(struct reg_token* rt, char* pattern,
                            int is, int ie)
{
  const char *rout_name = "make_list";
  int i;
  struct reg_token *rn = rt;
  if (rn->type != 0) {
    rn = mycalloc(rout_name, 1, sizeof *rn);
    rt->next = rn; rn->previous = rt;
  }
  rn->type = 3;
  rn->list = mycalloc(rout_name, 1, sizeof *rn->list);
  rn->list->chars = mymalloc_atomic(rout_name, 100 * sizeof *rn->list->chars);
  rn->list->max = 100;
  if (pattern[is] == '^') { /* invert list */
    rn->invert = 1; is++;
  }
  for (i = is; i <= ie; i++) {
    if (i < ie && pattern[i+1] == '-') {
      fill_list(pattern[i], pattern[i+2], rn->list);
      i += 2;
    }
    else fill_list(pattern[i], pattern[i], rn->list);
  }
  rn->list->chars[rn->list->curr] = '\0';
  return rn;
}

static int
list_count(struct r_char_array* list, int invert, char* s)
{
  char* p;
  int n = 0, i = 0;
  while (*s != '\0')
  {
    p = strchr(list->chars, *s++);
    if ((p != NULL && invert == 0) || (p == NULL && invert != 0))
    {
      i++;
      if (n < i) n = i;
    }
    else i = 0;
  }
  return n;
}

static int
new_comb(struct reg_token* start)
{
  struct reg_token* rt = start;
  while (rt)
  {
    if (rt->rep_cnt++ < rt->rep_max) return 1;
    rt->rep_cnt = 0;
    rt = rt->next;
  }
  return 0;
}

// private funtions (regex)

static struct reg_token*
flag(struct reg_token* rt, int* err)
{
  const char *rout_name = "flag";
  struct reg_token *rn = rt;
  if (rt == NULL || rt->type == 0) {
    *err = 3; return rt; /* illegal '*' */
  }
  if (rt->rep != 0) {
    *err = 4; return rt;  /* double '*' */
  }
  *err = 0;
  if (rt->type == 2) {
    if (rt->simple->curr > 1) {
      rn = mycalloc(rout_name, 1, sizeof *rn);
      rn->c = rt->simple->chars[--rt->simple->curr];
      rt->simple->chars[rt->simple->curr] = '\0';
      rn->type = 1;
      rt->next = rn;
    }
    else {
      rt->type = 1;
      rt->c = rt->simple->chars[0];
      myfree(rout_name, rt->simple->chars);
      myfree(rout_name, rt->simple);
      rt->simple = NULL;
    }
  }
  rn->rep = 1;
  return rn;
}

static struct reg_token*
add_tok(char c, struct reg_token* rt)
{
  const char *rout_name = "add_tok";
  struct reg_token *rn = rt;
  if (rn->type == 1 || rn->type == 3 || rn->type == 4 || rn->rep == 1) {
    rn = mycalloc(rout_name, 1, sizeof *rn);
    rt->next = rn; rn->previous = rt;
  }
  if (rn->type == 0) {
    rn->type = 2;
    rn->simple = mycalloc(rout_name, 1, sizeof *rn->simple);
    rn->simple->chars = mymalloc_atomic(rout_name, 100*sizeof *rn->simple->chars);
    rn->simple->max = 100;
  }
  if (rn->simple->curr == rn->simple->max) grow_r_char_array(rn->simple);
  rn->simple->chars[rn->simple->curr++] = c;
  rn->simple->chars[rn->simple->curr] = '\0';
  return rn;
}

static struct reg_token*
convert_pattern(char* pattern, int dollar, int* error)
{
  const char *rout_name = "convert_pattern";
  struct reg_token *rt, *first_token;
  int i, j, k, toggle = 0, last = strlen(pattern);
  *error = 0;
  first_token = rt = mycalloc(rout_name, 1, sizeof *rt);
  if (pattern[0] == '*') {*error = 1; return NULL;}
  for (i = 0; i < last; i++) {
    if (pattern[i] == '\\') toggle = 1;
    else if (toggle == 0 && pattern[i] == '[') {
      k = 0;
      for (j = i+2; j < last; j++) {
        if (pattern[j] == ']') {
          k = j; break;
        }
      }
      if (k == 0)  {*error = 2; return NULL;} /* no closing right ']' */
      rt = make_list(rt, pattern, i+1, j-1);
      i = j;
    }
    else if (toggle == 0 && pattern[i] == '*') {
      rt = flag(rt, error); if (*error != 0)  return NULL;
    }
    else if (toggle == 0 && pattern[i] == '.') rt = make_dot(rt);
    else {
      rt = add_tok(pattern[i], rt);
      toggle = 0;
    }
  }
  rt->match_end = dollar;
  return first_token;
}

static void
myregend(char* mypat, struct reg_token* start)
{
  const char *rout_name = "myregend";
  struct reg_token *rp, *aux;
  if (mypat != NULL) myfree(rout_name, mypat); mypat = NULL;
  rp = start;
  while (rp != NULL)
  {
    if (rp->simple != NULL)
    {
      if (rp->simple->chars != NULL) myfree(rout_name, rp->simple->chars);
      myfree(rout_name, rp->simple);
    }
    if (rp->list != NULL)
    {
      if (rp->list->chars != NULL) myfree(rout_name, rp->list->chars);
      myfree(rout_name, rp->list);
    }
    aux = rp;
    rp = rp->next;
    myfree(rout_name, aux);
  }
}

static void
edit_tokens(struct reg_token* start, char* pattern, char* string, int dollar)
{
  struct reg_token *tk, *tp;
  (void)pattern;
  tk = tp = start;
  while (tk)
  {
    if (tk->rep)
    {
      if (tk->type == 1) tk->rep_max = char_count(tk->c, string);
      else if (tk->type == 3)
        tk->rep_max = list_count(tk->list, tk->invert, string);
      else if (tk->type == 4)  tk->rep_max = strlen(string);
    }
    tp = tk;
    tk = tk->next;
  }
  tp->match_end = dollar;
}

static char*
match_token(struct reg_token* rt, char* p)
{
  char *q;
  int j, n;
  if (*p == '\0')  return NULL;
  switch (rt->type)
  {
    case 1:      /* simple character match - exact repetition */
      for (j = 0; j < rt->rep_cnt; j++) if (*p++ != rt->c) return NULL;
      return p;
    case 2:      /* simple string */
      if (strncmp(rt->simple->chars, p, rt->simple->curr) == 0)
        return &p[rt->simple->curr];
      else return NULL;
    case 3:    /* match character from list */
      if (rt->rep == 0)  n = 1;
      else               n = rt->rep_cnt;
      for (j = 0; j < n; j++)
      {
        q = strchr(rt->list->chars, *p++);
        if ((q == NULL && rt->invert == 0)
            || (q != NULL && rt->invert != 0)) return NULL;
      }
      return p;
    case 4:      /* dot */
      if (rt->rep == 0)  n = 1;
      else               n = rt->rep_cnt;
      for (j = 0; j < n; j++) if (*p++ == '\0') return NULL;
      return p;
    default:
      return NULL;
  }
}

static int
(match_all)(struct reg_token* start, char* string)
{
  struct reg_token* rt;
  char *p;
  if (string[0] == '\0')  return 0;
  restart:
  rt = start; p = string;
  while (rt)
  {
    if (((p = match_token(rt, p)) == NULL)) break;
    else if (rt->match_end && *p != '\0')
    {
      if (rt->rep == 0)  break;
      else if (rt->rep_cnt != 0) break;
      else if (rt == start)  break;
    }
    rt = rt->next;
  }
  if (rt == NULL) return 0;
  if (new_comb(start))  goto restart;
  return 1;
}

// public interface

int
myregex(char* pattern, char* string)
{
  /* returns 0 if pattern = regex matches string, else 1 */
  const char *rout_name = "myregex";
  char *mypat;
  struct reg_token *first_token;
  int error, j, l, dollar = 0, res = 0;
  if (pattern == NULL || (l = strlen(pattern)) == 0)  return 0;
  int len = strlen(pattern)+5;
  mypat = mymalloc_atomic(rout_name, len * sizeof *mypat);
  strcpy(mypat, pattern);
/* $ at end ? */
  if (mypat[l-1] == '$')
  {
    if (l == 1)  return 0;
    else if (l > 2 && strncmp(&mypat[l-3], ".*", 2) == 0)
    {
      if (l == 3)  return 0;
      if (mypat[l-4] != '\\') mypat[l-3] = '\0';
    }
    else
    {
      mypat[l-1] = '\0';
      dollar = 1;
    }
  }
  else if (l > 1 && strncmp(&mypat[l-2], ".*", 2) == 0)
  {
    if (l == 2)  return 0;
    if (mypat[l-3] != '\\') mypat[l-2] = '\0';
  }
  l = strlen(mypat);
/* ^ at start ? */
  if (mypat[0] == '^')
  {
    if (l == 0)  return 0;
    else
    {
      for (j = 0; j < l; j++) mypat[j] = mypat[j+1];
    }
  }
  else if (strncmp(mypat, ".*", 2) != 0)
  {
    for (j = l; j >= 0; j--)  mypat[j+2] = mypat[j];
    mypat[0] = '.'; mypat[1] = '*';
  }
  l = strlen(mypat);
  /*   printf("mypat: %s\n", mypat); */
  first_token = convert_pattern(mypat, dollar, &error);
  if (error)
  {
    if (error == 1) puts("+++ illegal '*' in pattern");
    if (error == 2) puts("+++ missing ']' in pattern");
    if (error == 3) puts("+++ illegal '*' in pattern");
    if (error == 4) puts("+++ double '*' in pattern");
    myregend(mypat, NULL);
    return 1;
  }
  edit_tokens(first_token, mypat, string, dollar);
  /*  dump_tokens(first_token); */
  res = match_all(first_token, string);
  myregend(mypat, first_token);
  return res;
}

