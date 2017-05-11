#include "madx.h"

static struct constraint*
new_constraint(int type)
{
  const char *rout_name = "new_constraint";
  struct constraint* new = mycalloc(rout_name, 1, sizeof *new);
  strcpy(new->name, "constraint");
  new->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", new->name);
  new->type = type;
  return new;
}

static void
grow_constraint_list(struct constraint_list* p)
{
  const char *rout_name = "grow_constraint_list";
  p->max *= 2;
  p->constraints = myrecalloc(rout_name, p->constraints, p->curr * sizeof *p->constraints, p->max * sizeof *p->constraints);
}

static struct constraint*
make_constraint(int type, int find_npos, struct command_parameter* par)
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

  new->n_pos = find_npos ? next_constr_namepos(new->name) : 0;
  if (find_npos && new->n_pos == 0) {
    fatal_error(
        " +-+-+- fatal error\n"
        "match - collect: illegal name\n"
        "      - try with the \"slow\" option\n"
        "      - name = ", new->name
    );
  }

  return new;
}

static void
add_to_constraint_list(struct constraint* cs, struct constraint_list* cl)
  /* add constraint cs to the constraint list cl */
{
  if (cl->curr == cl->max) grow_constraint_list(cl);
  cl->constraints[cl->curr++] = cs;
}

// public interface

struct constraint_list*
new_constraint_list(int length)
{
  const char *rout_name = "new_constraint_list";
  struct constraint_list* il = mycalloc(rout_name, 1, sizeof *il);
  strcpy(il->name, "constraint_list");
  il->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", il->name);
  il->curr = 0;
  il->max = length;
  il->constraints = mycalloc(rout_name, length, sizeof *il->constraints);
  return il;
}

struct constraint*
delete_constraint(struct constraint* cst)
{
  const char *rout_name = "delete_constraint";
  if (cst == NULL)  return NULL;
  if (stamp_flag && cst->stamp != 123456)
    fprintf(stamp_file, "d_c double delete --> %s\n", cst->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", "constraint");
  myfree(rout_name, cst);
  return NULL;
}

struct constraint_list*
delete_constraint_list(struct constraint_list* cl)
{
  const char *rout_name = "delete_constraint_list";
  if (cl == NULL)  return NULL;
  if (stamp_flag && cl->stamp != 123456)
    fprintf(stamp_file, "d_c_l double delete --> %s\n", cl->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", "constraint_list");
  myfree(rout_name, cl);
  return NULL;
}

void
dump_constraint(struct constraint* c)
{
  fprintf(prt_file,
          v_format("name: %s type: %I value: %F min: %F max: %F weight: %F\n"),
          c->name, c->type, c->value, c->c_min, c->c_max, c->weight);
}

#if 0 // not used...
void
dump_constraint_list(struct constraint_list* cl)
{
  int i;
  for (i = 0; i < cl->curr; i++)
  {
    if (cl->constraints[i]) dump_constraint(cl->constraints[i]);
  }
}
#endif

#if 0 // not used...
int
constraint_name(char* name, int* name_l, int* index)
  /* returns the name of the constraint */
{
  int ncp, nbl, len;
  struct constraint* c_c;
  c_c = current_node->cl->constraints[*index];
  len = strlen(c_c->name);
  ncp = len < *name_l ? len : *name_l; // min(len, *name_l)
  nbl = *name_l - ncp;
  strncpy(name, c_c->name, ncp);
  return 1;
}
#endif

void
fill_constraint_list(int type /* 1 node, 2 global */,
                     struct command* cd, struct constraint_list* cl)
{
  struct command_parameter_list* pl = cd->par;
  struct name_list* nl = cd->par_names;
  struct constraint* l_cons;
  int j, find_npos = type == 1 && !get_option("slow");
  for (j = 0; j < pl->curr; j++)
  {
    if (nl->inform[j] && pl->parameters[j]->type == 4)
    {
      l_cons = make_constraint(type, find_npos, pl->parameters[j]);
      add_to_constraint_list(l_cons, cl);
    }
  }
}

// fetch info about one constraint and advance to the next
int
next_constraint(char* name, int* name_l, int* type, double* value,
                double* c_min, double* c_max, double* weight,
                int* slow_match, int* pos, double* val, char* node_name, int* nn_len)
  /* returns the parameters of the next constraint; 0 = none, else count */
  // NOTE: does NOT set `val` for match2: match2 seems to invoke this function
  // only via `jacob_print` where `val` is not needed.
{
  struct constraint* c_c;
  char s, takenextmacro, nomore; /* RDM fork */

  /* RDM fork */
  if (match_is_on==2) {
    int i=match2_cons_curr[0];
    int j=match2_cons_curr[1];
    int k=match2_cons_curr[2];

    if(match2_cons_name[i][j] == NULL) {
      j++;

      if(j >= MAX_MATCH_CONS)
        takenextmacro = 1;
      else if (match2_cons_name[i][j]==NULL) // if above is true this one can cause seg fault
        takenextmacro = 1;
      else
        takenextmacro = 0;

      if(takenextmacro) {
        i++; j=0;

        if(i>=MAX_MATCH_MACRO)
          nomore = 1;
        else if(match2_cons_name[i][j]==NULL)
          nomore = 1;
        else // i,j is the next constraint
          nomore = 0;

        if(nomore == 0) {
          name = match2_cons_name[i][j];
          *name_l = strlen(name)+1;
          *type = 2; /* to be changed according to s or <,=,>*/
          *value = match2_cons_value[i][j];
          s = match2_cons_sign[i][j];
               if (s == '>' && *value > 0) *value=0;
          else if (s == '<' && *value < 0) *value=0;
          c_min = value; /* unknown use */
          c_max = value; /* unknown use */
          *weight = 1; /*hardcode no weight with this interface */
          string_from_table_row("twiss ", "name ", pos, node_name);
          k++;
          match2_cons_curr[0]=i;
          match2_cons_curr[1]=j;
          match2_cons_curr[2]=k;
          return k;
        } else {
          match2_cons_curr[0]=0;
          match2_cons_curr[1]=0;
          match2_cons_curr[2]=0;
          return 0;
        }
      }
    }
  }
  else { /* RDM old match */
    if (current_node->cl == NULL) return 0;

    if (current_node->con_cnt == current_node->cl->curr) {
      current_node->con_cnt = 0;
      return 0;
    }

    c_c = current_node->cl->constraints[current_node->con_cnt];
    strfcpy(name, c_c->name, *name_l);
    *type = c_c->type;

    if (c_c->ex_value == NULL) *value = c_c->value;
    else                       *value = expression_value(c_c->ex_value,2);

    if (c_c->ex_c_min == NULL) *c_min = c_c->c_min;
    else                       *c_min = expression_value(c_c->ex_c_min,2);

    if (c_c->ex_c_max == NULL) *c_max = c_c->c_max;
    else                       *c_max = expression_value(c_c->ex_c_max,2);

    *weight = c_c->weight;

    if (*slow_match) {
      string_from_table_row("twiss ", "name ", pos, node_name);
      double_from_table_row("twiss ",  name,   pos, val);
    }
    else {
      *val = current_node->match_data[c_c->n_pos-1];
      current_node_name(node_name, nn_len);
    }

    return ++current_node->con_cnt;
  }
  /* RDM fork */
  return 0;
}

int
next_constr_namepos(char* name)
/* returns the Fortran (!) position of the named variable
   in the opt_fun array of twiss (LD: weak!) */
{
  int pos = 0;
  switch (*name) {
    case 'a':
           if (name[3] == 'x') pos = 4;       // a??x
      else if (name[3] == 'y') pos = 7;       // a??y
      break;

    case 'b':
           if (name[3] == 'x') pos = 3;       // b??x
      else if (name[3] == 'y') pos = 6;       // b??y
      break;

    case 'd':
      if      (name[1] == 'x') pos = 15;      // dx
      else if (name[1] == 'y') pos = 17;      // dy
      else if (name[1] == 'p') {
             if (name[2] == 'x') pos = 16;    // dpx
        else if (name[2] == 'y') pos = 18;    // dpy
      }
      else if (name[1] == 'm') {
             if (name[3] == 'x') pos = 21;    // dm?x
        else if (name[3] == 'y') pos = 24;    // dm?y
      }
      else if (name[1] == 'd') {
	           if (name[2] == 'x') pos = 25;    // ddx
	      else if (name[2] == 'y') pos = 27;    // ddy
	      else if (name[2] == 'p') {
	        if      (name[3] == 'x') pos = 26;  // ddpx
	        else if (name[3] == 'y') pos = 28;  // ddpy
	      }
	    }
      break;

    case 'e':
      pos = 33;                               // e
      break;

    case 'm':
           if (name[2] == 'x') pos = 5;       // m?x
      else if (name[2] == 'y') pos = 8;       // m?y
      break;

    case 'p':
           if (name[1] == 'x') pos = 10;      // px
      else if (name[1] == 'y') pos = 12;      // py
      else if (name[1] == 't') pos = 14;      // pt
      else if (name[1] == 'h') {              // was 14 (LD: BUG)
	      if      (name[3] == 'x') pos = 20;    // ph?x
	      else if (name[3] == 'y') pos = 23;    // ph?y
      	}
      break;

    case 'r':
      if (name[1] == '1') {
             if (name[2] == '1') pos = 29;    // r11
        else if (name[2] == '2') pos = 30;    // r12
	    }
      else if (name[1] == '2') {
          if      (name[2] == '1') pos = 31;  // r21
          else if (name[2] == '2') pos = 32;  // r22
	    }
      else if (name[1] == 'e')
        /* start mod HG 10.10.2010 - trying to be ascii-independent */
        pos = (name[2]-'1')*6 + name[3]-'1' + 34;  // re11-re66, no range check...
      break;

    case 't':
      pos = 13;
      break;

    case 'w':
           if (name[1] == 'x') pos = 19;      // wx
      else if (name[1] == 'y') pos = 22;      // wy
      break;

    case 'x':
      pos = 9;                                // x
      break;

    case 'y':                                 // y
      pos = 11;
      break;
    }

  return pos;
}

int
next_global(char* name, int* name_l, int* type, double* value, double* c_min, double* c_max, double* weight)
  /* returns the parameters of the next global constraint;
     0 = none, else count */
{
  struct constraint* c_c;
  if (current_sequ->cl == NULL) return 0;
  if (current_sequ->con_cnt == current_sequ->cl->curr)
  {
    current_sequ->con_cnt = 0; return 0;
  }
  c_c = current_sequ->cl->constraints[current_sequ->con_cnt];
  strfcpy(name, c_c->name, *name_l);
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

void
update_node_constraints(struct node* c_node, struct constraint_list* cl)
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

void
update_sequ_constraints(struct sequence* sequ, struct constraint_list* cl)
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
