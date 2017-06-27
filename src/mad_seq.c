#include "madx.h"

static void
make_sequ_node(struct sequence* sequ, int occ_cnt)
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

static struct sequence_list*
delete_sequence_list(struct sequence_list* sql)
{
  const char *rout_name = "delete_sequence_list";
  if (sql == NULL) return NULL;
  if (stamp_flag && sql->stamp != 123456)
    fprintf(stamp_file, "d_s_l double delete --> %s\n", sql->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", sql->name);
  if (sql->list != NULL) delete_name_list(sql->list);
  if (sql->sequs != NULL) myfree(rout_name, sql->sequs);
  myfree(rout_name, sql);
  return NULL;
}

static void
dump_exp_sequ(struct sequence* sequ, int level)
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

#if 0 // not used...
static void
dump_sequ(struct sequence* c_sequ, int level)
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
#endif

static void
grow_sequence_list(struct sequence_list* l)
{
  const char *rout_name = "grow_sequence_list";
  struct sequence** sloc = l->sequs;
  int new = 2*l->max;
  l->max = new;
  l->sequs = mycalloc(rout_name, new, sizeof *l->sequs);
  for (int j = 0; j < l->curr; j++) l->sequs[j] = sloc[j];
  myfree(rout_name, sloc);
}

static void
make_occ_list(struct sequence* sequ)
  /* makes the node occurrence list */
{
  struct node* c_node = sequ->start;
  int i;

  while (c_node != NULL) {
    if (c_node->p_elem != NULL) {
      if ((i = name_list_pos(c_node->p_elem->name, occ_list)) < 0)
        i = add_to_name_list(c_node->p_elem->name, 1, occ_list);
      else ++occ_list->inform[i];
    }
    if (c_node == sequ->end) break;
    c_node = c_node->next;
  }
}

static void
all_node_pos(struct sequence* sequ)
  /* calculates all node positions in an expanded sequence */
{
  struct node* node = sequ->start;

  while (node != NULL) {
    if (node->p_elem != NULL)
      node->length = node->p_elem->length
                   = element_value(node, "l");
    else if (node->p_sequ != NULL)
      node->length = sequence_length(node->p_sequ);
    else fatal_error("node is neither element nor sequence:",
                     node->name);

    // 2014-Mar-19  16:27:15  ghislain:
    // The following assumes that the sequence is circular!!!!
    // this also conflicts with the insertion of sequences where a refpos is given.
    //if ((node->position = get_node_pos(node, sequ)) < zero)
    //  node->position += sequence_length(sequ);
    node->position = get_node_pos(node, sequ);

    if (node == sequ->end) break;
    node = node->next;
  }
}

static struct sequence*
extract_sequence(char* name, struct sequence* sequ, struct node* from, struct node* to, char* refpos)
{
  struct element* el;
  int pos, marker_pos, from_cond;
  struct command* clone;
  struct sequence* keep_curr_sequ = current_sequ;
  struct sequence* new_sequ;
  struct node *q = from, *from_node;
  struct node_list* new_nodes;
  double start_value, end_value;
  char tmp_name[2*NAME_L];

  if (get_option("info")) printf("+++ extracting sequence %s from %s to %s\n",
                                 sequ->name, from->name, to->name);
  current_sequ = new_sequence(name, sequ->ref_flag);
  current_sequ->share = 1; /* make shared for recombination */
  current_sequ->refpos = refpos;
  current_sequ->cavities = new_el_list(100);
  new_nodes = new_node_list(1000);
/* fill new_node list for 'from' references */
  while (q)
  {
    add_to_node_list(q, 0, new_nodes);
    if (q == to)  break;
    q = q->next;
  }
/* sequence construction */
  q = from;
  start_value = get_node_pos(from, sequ);
  end_value = get_node_pos(to, sequ);
  current_sequ->l_expr = NULL;
  current_sequ->length = end_value - start_value;

  // 2014-Mar-21  20:42:19  ghislain:
  // this might be dangerous for non-circular sequences!!!
  // if (current_sequ->length < zero) current_sequ->length += sequ->length;

  marker_pos = name_list_pos("marker", defined_commands->list);
  clone = clone_command(defined_commands->commands[marker_pos]);
  sprintf(c_dum->c, "%s$start", name);
  el = make_element(c_dum->c, "marker", clone, 0);
  current_node = NULL;
  make_elem_node(el, 1);
  current_sequ->start = current_node;
  make_elem_node(from->p_elem, from->occ_cnt);
  current_node->at_value = 0;
  current_node->at_expr = NULL;
/* loop over sequ nodes from 'from' to 'to' */
  while (current_node != NULL)
  {
    while (strchr(q->next->name, '$')) q = q->next; /* suppress internal markers */
    if (q->next->p_elem)
      make_elem_node(q->next->p_elem, q->next->occ_cnt);
    else if(q->next->p_sequ)
      make_sequ_node(q->next->p_sequ, q->next->occ_cnt);
    else
      fatal_error("node has neither element nor sequence reference:", q->next->name);
    q = q->next;
    if (q->p_elem && strcmp(q->p_elem->base_type->name, "rfcavity") == 0 &&
        find_element(q->p_elem->name, current_sequ->cavities) == NULL)
      add_to_el_list(&q->p_elem, 0, current_sequ->cavities, 0);
    from_cond = 0;
    if (q->from_name)
    {
      strcpy(tmp_name, q->from_name);
      strcat(tmp_name, ":1");
      if ((pos = name_list_pos(tmp_name, sequ->nodes->list)) > -1)
      {
        from_node = sequ->nodes->nodes[pos];
        from_cond = 1;
        if (name_list_pos(tmp_name, new_nodes->list) > -1)
          from_cond = 2;
      }
    }
    if (from_cond == 2)
    {
      current_node->from_name = q->from_name;
      current_node->at_value = q->at_value;
      current_node->at_expr = q->at_expr;
    }
    else
    {
      if (from_cond == 1)
        current_node->at_value = expr_combine(q->at_expr, q->at_value, " + ",
                                              from_node->at_expr, from_node->at_value,
                                              &current_node->at_expr);
      else
      {
        current_node->at_expr = q->at_expr;
        current_node->at_value = q->at_value;
      }
      current_node->at_value = expr_combine(current_node->at_expr,
                                            current_node->at_value, " - ",
                                            NULL, start_value,
                                            &current_node->at_expr);
      if (current_node->at_value < zero)
        current_node->at_value = expr_combine(current_node->at_expr,
                                              current_node->at_value, " + ",
                                              sequ->l_expr, sequ->length,
                                              &current_node->at_expr);
    }
    if (q == to) break;
  }
  clone = clone_command(defined_commands->commands[marker_pos]);
  sprintf(c_dum->c, "%s$end", name);
  el = make_element(c_dum->c, "marker", clone, 0);
  make_elem_node(el, 1);
  current_node->at_expr = current_sequ->l_expr;
  current_node->at_value = current_sequ->length;
  current_sequ->end = current_node;
  current_node->next = current_sequ->start;
  current_sequ->start->previous = current_node;
  new_sequ = current_sequ;
  current_sequ = keep_curr_sequ;
  if (get_option("info"))
    printf("+++ new sequence: %s  with current length = %.12g\n\n",
           new_sequ->name, new_sequ->length);
  return new_sequ;
}

static void
fill_sequ_var_list(struct sequence_list* sql, struct el_list* ell, struct var_list* varl)
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

static void
seq_edit_ex(struct sequence* seq)
{
  edit_sequ = seq;
  edit_is_on = 1;
  seqedit_install = seqedit_move = seqedit_remove = seqedit_replace = 0;
  if (edit_sequ->ex_start != NULL) {
    edit_sequ->ex_nodes = delete_node_list(edit_sequ->ex_nodes);
    edit_sequ->ex_start = delete_node_ring(edit_sequ->ex_start);
  }
  if (occ_list == NULL)
    occ_list = new_name_list("occ_list", 10000);  /* for occurrence count */
  else occ_list->curr = 0;
  resequence_nodes(edit_sequ);
  all_node_pos(edit_sequ);
}

static void
seq_end_ex(void)
{
  occ_list->curr = 0;
  resequence_nodes(edit_sequ);
  selected_ranges->curr = 0;
  selected_ranges->list->curr = 0;
  edit_is_on = 0;
}

static void
expand_line(struct char_p_array* l_buff)
  /* expands a beam line, applies rep. count and inversion */
{

  int add=0, i=0, j=0, k=0, n=0, number=0, dummy=0, rep=-1, pos=0;
  int level = 0, l_max = 0, b_cnt = 0;
  char* p = blank;
  struct int_array* lbpos = new_int_array(l_buff->curr);
  struct int_array* rbpos = new_int_array(l_buff->curr);
  struct int_array* b_level = new_int_array(l_buff->curr);

  /* first get all bracket pairs with their level; keep max. level */
  for (i = 0; i < l_buff->curr; i++) {
    if (*l_buff->p[i] == '(') {
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
  for (level = l_max; level >=0; level--) {
    for (i = 0; i < b_cnt; i++) {
      if (b_level->i[i] == level && (pos = lbpos->i[i]) > 1) {
        if (*l_buff->p[pos-1] == '*') {
          sscanf(l_buff->p[pos-2], "%d", &rep);
	  add = rep - 1;
          number = rbpos->i[i] - pos - 1; /* inside bracket */
          n = number * add; /* extra tokens */
          while (l_buff->curr + n >= l_buff->max) grow_char_p_array(l_buff);

          for (j = l_buff->curr; j > pos + number; j--) /* shift upwards */
            l_buff->p[j+n] = l_buff->p[j];

          l_buff->curr += n;

          for (k = 1; k <= add; k++) {
            for (j = pos+1; j <= pos+number; j++)
              l_buff->p[j+k*number] = tmpbuff(l_buff->p[j]);
          }

          for (j = 0; j < b_cnt; j++) {  /* reset bracket pointers */
            if (lbpos->i[j] > pos + number) lbpos->i[j] += n;
            if (rbpos->i[j] > pos + number) rbpos->i[j] += n;
          }

          l_buff->p[pos-1] = l_buff->p[pos-2] = blank;
        }
      }
    }
  }

  /* loop over buffer, expand simple element repetition defined with '*' f.g. 10*myquad  */
  for (pos = 2; pos < l_buff->curr; pos++) {
    if (*l_buff->p[pos] == '*') {
      rep = -1;
      sscanf(l_buff->p[pos-1], "%d", &rep);
      if (rep < 0) fatal_error("expand_line","Problem with reading number of copies");
      n = add = rep - 1;
      while (l_buff->curr + n >= l_buff->max) grow_char_p_array(l_buff);

      for (j = l_buff->curr; j > pos + 1; j--) /* shift upwards */
        l_buff->p[j+n] = l_buff->p[j];

      l_buff->curr += n;

      j = pos+1;
      for (k = 1; k <= add; k++) l_buff->p[j+k] = tmpbuff(l_buff->p[j]);

      for (j = 0; j < b_cnt; j++) {  /* reset bracket pointers */
        if (lbpos->i[j] > pos + 1) lbpos->i[j] += n;
        if (rbpos->i[j] > pos + 1) rbpos->i[j] += n;
      }

      // 2014-Aug-18  19:56:34  ghislain: in case of single element with rep_count,
      // l_buff->p[pos-1] points to the rep_count
      // and l_buff->p[pos-2] the previous element which is then lost!!!
      // l_buff->p[pos-1] = l_buff->p[pos-2] = blank;
      l_buff->p[pos-1] = blank;
    }
  }

  /* get bracket pointers including new ones */
  while (b_level->max < l_buff->curr) grow_int_array(b_level);
  while (lbpos->max < l_buff->curr) grow_int_array(lbpos);
  while (rbpos->max < l_buff->curr) grow_int_array(rbpos);
  level = b_cnt = 0;
  for (i = 0; i < l_buff->curr; i++) {
    if (*l_buff->p[i] == '(') {
      lbpos->i[b_cnt] = i;
      b_level->i[b_cnt++] = level++;
    }
    else if (*l_buff->p[i] == ')')  level--;
  }

  for (i = 0; i < b_cnt; i++)
    get_bracket_t_range(l_buff->p, '(', ')', lbpos->i[i],
                        l_buff->curr-1, &dummy, &rbpos->i[i]);
  lbpos->curr = rbpos->curr = b_level->curr = b_cnt;

  /* now loop over level from highest down to zero, invert if '-' */
  for (level = l_max; level >= 0; level--) {
    for (i = 0; i < b_cnt; i++) {
      pos = lbpos->i[i];
      if (b_level->i[i] == level) {
        p = blank;
        for (j = pos - 1; j > 0; j--) {
          p = l_buff->p[j];
          if (*p != ' ')  break;
        }
        if (*p == '-') {
          number = rbpos->i[i] - pos - 1; // length of sequence
          n = number / 2;
          for (j = 0; j < n; j++) {
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
  for (i = 0; i < l_buff->curr; i++) {
    if (isalpha(*l_buff->p[i])) l_buff->p[n++] = l_buff->p[i];
  }
  l_buff->curr = n;
  lbpos = delete_int_array(lbpos);
  rbpos = delete_int_array(rbpos);
  b_level = delete_int_array(b_level);
}

static void
expand_sequence(struct sequence* sequ, int flag)
  /* expands a sequence into nodes, expands sequence nodes */
{
  /* Transfers errors from original nodes if flag != 0;
     this is needed for SXF input  */
  struct node *p, *q = sequ->start;

  int debug=get_option("debug");

  // 2015-Jun-09  14:14:24  ghislain: guard added for negative seq length
  if (sequ->length < 0)
    fatal_error("trying to expand sequence with negative length:", sequ->name);

  if (debug)
    printf("\n\nTOP Expand_sequence name %s with length %e, ref_flag %d\n",
	   sequ->name, sequ->length, sequ->ref_flag);

  p = sequ->ex_start = clone_node(sequ->start, 0);
  add_to_node_list(p, 0, sequ->ex_nodes);

  while (p != NULL) {
    if (q == sequ->end) break;
    p->next = clone_node(q->next, flag);
    p->next->previous = p;
    p = p->next;
    q = q->next;

    if (p->p_sequ == NULL) { // simple element, not a subsequence
      // 2015-Jun-09  14:14:24  ghislain: guard added for negative element length
      if (p->length < 0)
	fatal_error("trying to add node with negative length to current sequence:", p->name);

      add_to_node_list(p, 0, sequ->ex_nodes);
    }
    else { // subsequence
      if (p->p_sequ->refpos != NULL) { // REFPOS given for subsequence, ignore REFER of current sequence
	if (debug)
	  printf("\n\n Expand sub-sequence %s with initial position %e, final position %e, length %e, ref_flag %d, refpos '%s'\n",
		 p->name, p->position, p->position - sequ->ref_flag*p->p_sequ->length/2.,
		 p->length, sequ->ref_flag, p->p_sequ->refpos);
	p = expand_node(p, sequ, sequ, p->position - sequ->ref_flag*p->p_sequ->length/2. );
	if (debug) printf("\n\n");
      }
      else { // no REFPOS given
	if (debug)
	  printf("\n\n Expand sub-sequence %s with position %e, length %e, ref_flag %d\n",
		 p->name, p->position, p->length, sequ->ref_flag);
	p = expand_node(p, sequ, sequ, p->position);
	if (debug) printf("\n\n");
      }
    }
  }

  sequ->ex_end = p;
  sequ->ex_end->next = sequ->ex_start;
  sequ->ex_start->previous = sequ->ex_end;

  p = sequ->ex_start;
  while (p != sequ->ex_end) {
    if (strstr(p->base_name, "kicker") || strstr(p->base_name, "monitor"))
      p->enable = 1; /* flag for orbit correction module */
    p = p->next;
  }

  /*
    Attempt to discard attached twiss table not anymore valid...
    note: it cannot be done properly as table_register keep references
          on-pointer-to-table (address of address) not to-table itself...
          hence
          if (sequ->tw_table) sequ->tw_table = NULL;
          breaks the list of table_register used everywhere,
          leading to bus error in many places (pointers are almost never checked)
  */
}

static void
seq_flatten(struct sequence* sequ)
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

  sequ->start = sequ->ex_start;     sequ->ex_start = NULL;
  sequ->end = sequ->ex_end;         sequ->ex_end = NULL;

  c_node = sequ->start;
  while (c_node != NULL) {
    c_node->at_value = c_node->position;
    c_node->at_expr = NULL;
    c_node->from_name = NULL;
    if (c_node == sequ->end) break;
    c_node = c_node->next;
  }
  sequ->ref_flag = 0;
}

static void
insert_elem(struct sequence* sequ, struct node* node)
  /* inserts an element in a sequence as function of its position */
{
  struct node* c_node = sequ->start->next;

  while (c_node != NULL) {
    if (node->position <= c_node->position || c_node == sequ->end) break;
    c_node = c_node->next;
  }
  link_in_front(node, c_node);
}

static struct node*
install_one(struct element* el, char* from_name, double at_value, struct expression* at_expr, double position)
  /* adds an element to a sequence */
{
  struct node* node;
  int i, occ = 1;

  if (strcmp(el->base_type->name, "rfcavity") == 0 &&
      find_element(el->name, edit_sequ->cavities) == NULL)
    add_to_el_list(&el, 0, edit_sequ->cavities, 0);
  if ((i = name_list_pos(el->name, occ_list)) < 0)
    i = add_to_name_list(el->name, occ, occ_list);
  else occ = ++occ_list->inform[i];

  node = new_elem_node(el, occ);
  add_to_node_list(node, 0, edit_sequ->nodes);
  node->position = position;
  node->at_value = at_value;
  node->at_expr = at_expr;
  node->from_name = from_name;
  set_command_par_value("at", el->def, position);
  insert_elem(edit_sequ, node);
  return node;
}

static void
make_sequ_from_line(char* name)
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
  line->dead = 1;   /* prevent line from further conversion to sequence */
  line_buffer = new_char_p_array(1000);

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
    occ_list = new_name_list("occ_list", 10000);  /* for occurrence count */
  else occ_list->curr = 0;
  sprintf(c_dum->c, "%s$start", current_sequ->name);
  el = make_element(c_dum->c, "marker", clone, 0);
  current_node = NULL;
  make_elem_node(el, 1);
  current_sequ->start = current_node;


  current_sequ->length = line_nodes(line_buffer);

  sprintf(c_dum->c, "%s$end", current_sequ->name);
  el = make_element(c_dum->c, "marker", clone, 0);
  make_elem_node(el, 1);
  current_node->at_value = current_node->position = current_sequ->length;
  current_sequ->end = current_node;
  current_sequ->start->previous = current_sequ->end;
  current_sequ->end->next = current_sequ->start;
  current_sequ->line = 1; /* remember origin of sequence */

  if(line_buffer) delete_char_p_array(line_buffer,1);
}

static void
export_sequence(struct sequence* sequ, FILE* file, int noexpr)
  /* exports sequence in mad-X format */
{
  char num[2*NAME_L];
  struct element* el;
  struct sequence* sq;
  struct node* c_node = sequ->start;
  int exp_par_flag;
  int seqref = 0;
  char rpos[3][6] = {"exit", "centre", "entry"};
  seqref = sequ->ref_flag;  /* uncomment line to get entry or exit */
  *c_dum->c = '\0';
  if (sequ->share) strcat(c_dum->c, "shared ");
  strcat(c_dum->c, sequ->export_name);
  strcat(c_dum->c, ": sequence");
  if (seqref)
  {
    strcat(c_dum->c, ", refer = ");
    strcat(c_dum->c, rpos[seqref+1]);
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
      if (exp_par_flag) export_el_def(c_node->p_elem, c_dum->c, noexpr);
      write_nice(c_dum->c, file);
    }
    if (c_node == sequ->end)  break;
    c_node = c_node->next;
  }
  strcpy(c_dum->c, "endsequence");
  write_nice(c_dum->c, file);
}

static void
export_sequ_8(struct sequence* sequ, struct command_list* cl, FILE* file)
  /* exports sequence in mad-8 format */
  /* set refer = centre always (local var. seqref) HG 9.1.09 */
{
  char num[2*NAME_L];
  int exp_par_flag;
  int seqref = 0;
  struct element* el;
  struct sequence* sq;
  struct node* c_node = sequ->start;

  seqref = sequ->ref_flag;  /* uncomment line to get entry or exit */

  if (pass_select_list(sequ->name, cl) == 0)  return;

  *c_dum->c = '\0';
  strcat(c_dum->c, sequ->export_name);
  strcat(c_dum->c, ": sequence");
  if (seqref == 1)  strcat(c_dum->c, ", refer=entry");
  write_nice_8(c_dum->c, file);

  while(c_node != NULL) {
    exp_par_flag = 0;
    *c_dum->c = '\0';

    if (strchr(c_node->name, '$') == NULL
        && strcmp(c_node->base_name, "drift") != 0) {
      if ((el = c_node->p_elem) != NULL) {
        if (c_node->p_elem->def_type) {
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
      else {
        sprintf(num, v_format("%F"), c_node->at_value);
        strcat(c_dum->c, supp_tb(num));
      }

      if (c_node->from_name != NULL) {
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

static void
write_sequs(struct sequence_list* sql,struct command_list* cl, FILE* file, int noexpr)
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
          export_sequence(sql->sequs[i], file, noexpr);
      }
  }
}

static void
seq_cycle(struct in_cmd* cmd)
  /* cycles a sequence */
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  struct node *node, *clone;
  char* name = NULL;
  int pos = name_list_pos("start", nl);
  if (nl->inform[pos] && (name = pl->parameters[pos]->string) != NULL)
  {
    sprintf(c_dum->c, "%s:1", name);
    if ((pos = name_list_pos(c_dum->c, edit_sequ->nodes->list)) > -1)
    {
      node = edit_sequ->nodes->nodes[pos];
      sprintf(c_dum->c, "%s%s_p_", edit_sequ->name,strip(node->name));
      if (strstr(node->previous->name, "_p_") == NULL)
      {
        clone = clone_node(node, 0);
        clone->p_elem = clone_element(node->p_elem);
        strcpy(clone->p_elem->name, c_dum->c);

        /* IA 29.11.07 : fixes a bug with aperture module */
	/* Removed HG 11.10.2009 */
        /* sprintf(c_dum->c, " ");
           set_command_par_string("apertype", clone->p_elem->def,c_dum->c); */

        add_to_el_list(&clone->p_elem, node->p_elem->def->mad8_type,element_list,1);
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

static void
seq_edit(struct in_cmd* cmd)
  /* executes seqedit command */
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  char* name = NULL;
  int pos;

  pos = name_list_pos("sequence", nl);

  if (nl->inform[pos] && (name = pl->parameters[pos]->string) != NULL) {
    if ((pos = name_list_pos(name, sequences->list)) >= 0) {
      if (sequences->sequs[pos]->line)
        warning("sequence originates from line,","edit ignored");
      else  seq_edit_ex(sequences->sequs[pos]);
    }
    else warning("unknown sequence:", "ignored");
  }
  else warning("seqedit without sequence:", "ignored");
}

static void
seq_end(struct in_cmd* cmd)
  /* executes endedit command */
{
  char tmp[8];
  (void)cmd;
  sprintf(tmp, "%d", seqedit_install);
  put_info("seqedit - number of elements installed: ", tmp);
  sprintf(tmp, "%d", seqedit_move);
  put_info("seqedit - number of elements moved:     ", tmp);
  sprintf(tmp, "%d", seqedit_remove);
  put_info("seqedit - number of elements removed:   ", tmp);
  sprintf(tmp, "%d", seqedit_replace);
  put_info("seqedit - number of elements replaced:  ", tmp);
  seq_end_ex();
}

static void
seq_install(struct in_cmd* cmd)
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

  if (nl->inform[pos_e] && (name_e = pl->parameters[pos_e]->string) != NULL) {
    if (nl->inform[pos_c] && (name_c = pl->parameters[pos_c]->string) != NULL) {
      if ((cl = find_element(name_c, element_list)) == NULL) {
        warning("ignored, unknown class:", name_c);
        return;
      } else {
        el = clone_element(cl);
        strcpy(el->name, name_e);
        add_to_el_list(&el, cl->def->mad8_type, element_list, 2);
      }
    } else if ((el = find_element(name_e, element_list)) == NULL) {
      warning("ignored, unknown command or element:", name_c);
      return;
    }
  } else {
    warning("no element specified,","ignored");
    return;
  }

  if (nl->inform[name_list_pos("at", nl)] == 0) {
    warning("no 'at':", "ignored");
    return;
  }

  at = command_par_value("at", cmd->clone);
  expr = clone_expression(command_par_expr("at", cmd->clone));

  pos = name_list_pos("from", nl);
  if (nl->inform[pos]) {
    from_name = pl->parameters[pos]->string;

    if (strcmp(from_name, "selected") == 0) {
      if (seqedit_select->curr == 0) {
        warning("no active select commands:", "ignored");
	return;
      }
      if (get_select_ranges(edit_sequ, seqedit_select, selected_ranges) == 0) any = 1;
      c_node = edit_sequ->start;
      while (c_node != NULL) {
	if (any || name_list_pos(c_node->name, selected_ranges->list) > -1) {
	  for (k = 0; k < seqedit_select->curr; k++) {
	    myrepl(":", "[", c_node->name, name);
	    strcat(name, "]");
	    if (strchr(name, '$') == NULL && pass_select(c_node->name, seqedit_select->commands[k])) break;
	  }
	  if (k < seqedit_select->curr) {
	    from = get_node_pos(c_node, edit_sequ);
	    pname = permbuff(name);
	    install_one(el, pname, at, expr, at+from);
	    seqedit_install++;
	  }
	}
	if (c_node == edit_sequ->end) break;
	c_node = c_node->next;
      }

    } else {
      from_name = permbuff(pl->parameters[pos]->string);
      if ((from = hidden_node_pos(from_name, edit_sequ)) == INVALID) {
        warning("ignoring 'from' reference to unknown element:", from_name);
        return;
      }
      install_one(el, from_name, at, expr, at+from);
      seqedit_install++;
    }
  } else {
    install_one(el, from_name, at, expr, at);
    seqedit_install++;
  }
}

static void
seq_move(struct in_cmd* cmd)
  /* executes move command */
{
  char *name, *from_name;
  double at, by, to, from = zero;
  int any = 0, k;
  struct node *node, *next;
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
        while (node != edit_sequ->end)
        {
          node = node->next; node->moved = 0;
        }
        node = edit_sequ->start;
        while (node != NULL && node != edit_sequ->end)
        {
          next = node->next;
          if (node->moved == 0)
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
                  node = install_one(el, NULL, at, NULL, at);
                  node->moved = 1;
                  seqedit_move++;
                }
              }
            }
          }
          node = next;
        }
      }
    }
    else
    {
      strcpy(c_dum->c, name);
      square_to_colon(c_dum->c);
      if ((pos = name_list_pos(c_dum->c, edit_sequ->nodes->list)) > -1)
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

static void
seq_reflect(struct in_cmd* cmd)
  /* executes reflect command */
{
  struct node *tmp, *c_node;

  (void)cmd;
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
      = sequence_length(edit_sequ) - c_node->position;
    if (c_node == edit_sequ->end) break;
    c_node = c_node->next;
  }
}

static void
seq_remove(struct in_cmd* cmd)
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
      strcpy(c_dum->c, name);
      square_to_colon(c_dum->c);
      if ((pose = name_list_pos(c_dum->c, edit_sequ->nodes->list)) > -1)
      {
        seqedit_remove += remove_one(edit_sequ->nodes->nodes[pose]);
      }
      else warning("ignored because of unknown element:", name);
    }
  }
  else  warning("no element specified,","ignored");
}

static void
seq_replace(struct in_cmd* cmd)
  /* executes replace command */
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  struct node** rep_nodes = NULL;
  struct element** rep_els = NULL;
  struct node *node, *c_node;
  char* name;
  struct element* el;
  int count = count_nodes(edit_sequ);
  int any = 0, k, rep_cnt = 0, pos;

  pos = name_list_pos("element", nl);
  if ( !(nl->inform[pos]) || (name = pl->parameters[pos]->string) == NULL) {
    warning("no element specified, ","ignored");
    return;
  }

  if (strcmp(name, "selected") == 0) { // replace selected elements
    if (seqedit_select->curr == 0) {
      warning("no active select commands:", "ignored");
      return;
    }

    pos = name_list_pos("by", nl);
    if ( !(nl->inform[pos]) || (name = pl->parameters[pos]->string) == NULL) {
      warning("'by' missing, ","ignored");
      return;
    }

    if ((el = find_element(name, element_list)) == NULL) {
      warning("ignoring unknown 'by' element:",name);
      return;
    }

    rep_nodes = mymalloc("seq_replace", count*sizeof *rep_nodes);
    rep_els = mymalloc("seq_replace", count*sizeof *rep_els);

    if (get_select_ranges(edit_sequ, seqedit_select, selected_ranges) == 0) any = 1;

    c_node = edit_sequ->start;
    while (c_node != NULL) {
      if (any || name_list_pos(c_node->name, selected_ranges->list) > -1) {
	name = NULL;
	for (k = 0; k < seqedit_select->curr; k++) {
	  if (c_node->p_elem != NULL) name = c_node->p_elem->name;
	  if (name != NULL && strchr(name, '$') == NULL &&
	      pass_select(name, seqedit_select->commands[k]) ) break;
	}
	if (k < seqedit_select->curr) {
	  rep_els[rep_cnt] = el;
	  rep_nodes[rep_cnt++] = c_node;
	}
      }
      if (c_node == edit_sequ->end) break;
      c_node = c_node->next;
    }
  }

  else { // replace named elements
    rep_nodes = mymalloc("seq_replace", count*sizeof *rep_nodes);
    rep_els = mymalloc("seq_replace", count*sizeof *rep_els);

    strcpy(c_dum->c, name);
    square_to_colon(c_dum->c);

    if ((pos = name_list_pos(c_dum->c, edit_sequ->nodes->list)) > -1) {
      node = edit_sequ->nodes->nodes[pos];
      pos = name_list_pos("by", nl);
      if (nl->inform[pos] && (name = pl->parameters[pos]->string) != NULL) {
	if ((el = find_element(name, element_list)) != NULL) {
	  rep_els[rep_cnt] = el;
	  rep_nodes[rep_cnt++] = node;
	}
	else warning("ignoring unknown 'by' element: ",name);
      }
      else warning("'by' missing, ","ignored");
    }
    else warning("ignored because of unknown element: ", name);
  }

  for (k = 0; k < rep_cnt; k++)  replace_one(rep_nodes[k], rep_els[k]);

  seqedit_replace = rep_cnt;

  if (rep_nodes) myfree("seq_replace", rep_nodes);
  if (rep_els)   myfree("seq_replace", rep_els);

}

static struct element*
get_drift(double length)
  /* makes a drift space with the required length */
{
  const double tol = 1e-12; // length tolerance for sharing implicit drift
  struct element *p, *bt;
  struct command* clone;
  char key[NAME_L];

  for (int i = 0; i < drift_list->curr; i++) {
    p = drift_list->elem[i];
    if (fabs(p->length - length) < tol) return p;
  }

  sprintf(key, "drift_%d", drift_list->curr);
  bt = find_element("drift", base_type_list);
  clone = clone_command(bt->def);
  store_comm_par_value("l", length, clone);
  p = make_element(key, "drift", clone, 0);
  add_to_el_list(&p, 1, drift_list, 0);
  return p;
}

static int
add_drifts(struct node* c_node, struct node* end)
{
  const double tol = 1e-6;
  int cnt;

  char buf[256];

  int debug;
  debug = get_option("debug");

  if (!c_node) return 0;

  for (cnt=1; c_node != end && c_node->next; c_node = c_node->next, cnt++) {
    double drift_beg = c_node->position + c_node->length / 2;
    double drift_end = c_node->next->position - c_node->next->length / 2;
    double drift_len = drift_end-drift_beg;

    if (drift_len < -tol) {
      // implicit drift with negative length
      sprintf(buf, " %s and %s, length %e", c_node->name, c_node->next->name, drift_len);

      if (debug) {
	printf("\ncurrent node name %s position: %e length %e \n",
	       c_node->name, c_node->position, c_node->length);
	printf("next    node name %s position: %e length %e \n\n",
	       c_node->next->name, c_node->next->position, c_node->next->length);
      }

     fatal_error("negative drift between elements", buf);
    }
    else if (drift_len > tol) {
      // create or share 'long-enough' implicit drift
      struct element *drift = get_drift(drift_len);
      struct node *drift_node = new_elem_node(drift, 0);
      link_in_front(drift_node, c_node->next);
      drift_node->position = drift_beg + drift_len / 2;
      if (debug) printf("inserting a drift of length %e at position %e \n \n",
			drift_len,drift_beg + drift_len / 2);
      cnt++;
    }
    else
      // length in [-tol, tol], nothing to do (no drift inserted)
      // 2014-Feb-04  11:52:07  ghislain: tghought of adding a warning that a very short drift was ignored.
      // but the number of warnings explodes VERY quicky so leave it out for now.
      // sprintf(buf, " Length of drift: %e vs. tolerance: %e", drift_len, tol);
      // warning("Drift length below tolerance level was ignored.", buf);
      (void)0;
  }

  return cnt;
}


void
use_sequ(struct in_cmd* cmd)
{
  const char *rout_name = "use_sequ";
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  struct command* keep_beam = current_beam;
  int pos, lp;
  char* name;

  if (sequ_is_on)
    fatal_error("no endsequence yet for sequence:", current_sequ->name);

  // YIL: calling use screws up any twiss table calculated.. this is a "quick fix"
  pos = name_list_pos("period", nl);
  if (nl->inform[pos] == 0) pos = name_list_pos("sequence", nl);

  if (nl->inform[pos]) {  /* parameter has been read */
    if (current_range != NULL) {
      myfree(rout_name, current_range);
      current_range = NULL;
    }

    name = pl->parameters[pos]->string;
    if ((pos = name_list_pos(name, line_list->list)) > -1 && line_list->macros[pos]->dead == 0)
      make_sequ_from_line(name); /* only if not disabled */

    if ((lp = name_list_pos(name, sequences->list)) > -1) {
      current_sequ = sequences->sequs[lp];

      if (attach_beam(current_sequ) == 0)
        fatal_error("USE - sequence without beam:", current_sequ->name);

      current_sequ->beam = current_beam;
      pos = name_list_pos("range", nl);
      if (nl->inform[pos])  /* parameter has been read */
        current_range = tmpbuff(pl->parameters[pos]->string);

      current_sequ->tw_valid = 0;
      expand_curr_sequ(0);

      pos = name_list_pos("survey", nl);
      if (nl->inform[pos]) {  /* parameter has been read */
         pro_use_survey();
         pos = name_list_pos("survtest", nl);
         if (nl->inform[pos]) survtest_();
         exec_delete_table("survey");
     	}
    }
    else warning("unknown sequence skipped:", name);
  }

  current_beam = keep_beam;
}

int
sequ_check_valid_twiss(struct sequence * sequ)
{
  return sequ->tw_table != NULL && sequ->tw_valid;
}

int
set_cont_sequence(void)
{
  if (current_sequ->next_sequ != NULL)
    {
     set_sequence(current_sequ->next_sequ);
     return 1;
    }
  else return 0;
}

void
set_sequence(char* name)
{
  int lp;
  struct sequence* t_sequ;
  mycpy(c_dum->c, name);
  /* makes sequence "name" the current sequence */
  if ((lp = name_list_pos(c_dum->c, sequences->list)) > -1)
    t_sequ = sequences->sequs[lp];
  else
  {
    warning("unknown sequence ignored:", name);
    return;
  }
  if (attach_beam(t_sequ) == 0)
    fatal_error("USE - sequence without beam:", t_sequ->name);
  t_sequ->beam = current_beam;
  if (t_sequ == NULL || t_sequ->ex_start == NULL)
  {
    warning("sequence not active,", "SET ignored");
    return;
  }
  current_sequ = t_sequ;
}

// public interface

struct sequence*
new_sequence(const char* name, int ref)
{
  const char *rout_name = "new_sequence";
  struct sequence* s = mycalloc(rout_name, 1, sizeof *s);
  strcpy(s->name, name);
  s->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", s->name);
  s->ref_flag = ref;
  s->nodes = new_node_list(10000);
  return s;
}

struct sequence_list*
new_sequence_list(int length)
{
  const char *rout_name = "new_sequence_list";
  struct sequence_list* s = mycalloc(rout_name, length, sizeof *s);
  strcpy(s->name, "sequence_list");
  s->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", s->name);
  s->max = length;
  s->list = new_name_list(s->name, length);
  s->sequs = mycalloc(rout_name, length, sizeof *s->sequs);
  return s;
}

struct node*
new_sequ_node(struct sequence* sequ, int occ_cnt)
{
  struct node* p;
  p = new_node(compound(sequ->name, occ_cnt));
  p->p_sequ = sequ;
  p->length = sequence_length(sequ);
  p->base_name = permbuff("sequence");
  return p;
}

struct sequence*
delete_sequence(struct sequence* sequ)
{
  const char *rout_name = "delete_sequence";
  int lp;
  if (sequ->ex_start != NULL)
  {
    sequ->ex_nodes = delete_node_list(sequ->ex_nodes);
    sequ->ex_start = delete_node_ring(sequ->ex_start);
    sequ->orbits = delete_vector_list(sequ->orbits);
    myfree(rout_name, sequ->all_nodes);
  }
  if ((lp = name_list_pos(sequ->name, sequences->list)) > -1)
    remove_from_sequ_list(sequences->sequs[lp], sequences);
  if (sequ->l_expr) sequ->l_expr = delete_expression(sequ->l_expr);
  sequ->nodes = delete_node_list(sequ->nodes);
  sequ->start = delete_node_ring(sequ->start);
  if (sequ->cavities) sequ->cavities = delete_el_list(sequ->cavities);
  myfree(rout_name, sequ);
  return NULL;
}

void
remove_from_sequ_list(struct sequence* sequ, struct sequence_list* sql)
  /* removes sequence sequ from sequence list sql */
{
  int i;
  if ((i = remove_from_name_list(sequ->name, sql->list)) > -1)
    sql->sequs[i] = sql->sequs[--sql->curr];
  return;
}

double
sequence_length(struct sequence* sequ)
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

void
enter_sequence(struct in_cmd* cmd)
  /* handles sequence start and end on input */
{
  struct name_list* nl;
  struct command_parameter_list* pl;
  int i, k = 0, pos, aux_pos;
  char** toks = cmd->tok_list->p;
  struct element* el;
  struct command* clone;

  aux_pos = strcmp(toks[0], "shared") == 0 ? 1 : 0;

  if (strcmp(toks[0], "endsequence") == 0) {
    pos = name_list_pos("marker", defined_commands->list);
    clone = clone_command(defined_commands->commands[pos]);
    sprintf(c_dum->c, "%s$end", current_sequ->name);
    el = make_element(c_dum->c, "marker", clone, 0);
    make_elem_node(el, 1);
    current_node->at_expr = current_sequ->l_expr;
    current_node->at_value = current_sequ->length;
    current_sequ->end = current_node;
    current_sequ->start->previous = current_sequ->end;
    current_sequ->end->next = current_sequ->start;
  }

  else if (strcmp(toks[aux_pos+2], "sequence") == 0) {
    for (i = aux_pos+3; i < cmd->tok_list->curr; i++) {
      if (strcmp(toks[i], "refer") == 0) {
	if (i+2 < cmd->tok_list->curr) {
	  if      (strcmp(toks[i+2], "entry") == 0)  k = 1;
	  else if (strcmp(toks[i+2], "exit")  == 0)  k = -1;
	}
	break;
      }
    }

    if ((pos = name_list_pos(toks[aux_pos], sequences->list)) >= 0) {
      // sequence exists already; delete the old one
      /*printf("enter_sequence: removing %s\n", sequences->sequs[pos]->name);*/
      remove_from_sequ_list(sequences->sequs[pos], sequences);
      sequences->sequs[pos] = delete_sequence(sequences->sequs[pos]);
    }

    current_sequ = new_sequence(toks[aux_pos], k);
    add_to_sequ_list(current_sequ, sequences);
    /* prevent a line with this name from expansion */
    disable_line(current_sequ->name, line_list);

    cmd->clone = clone_command(cmd->cmd_def);
    scan_in_cmd(cmd);
    nl = cmd->clone->par_names;
    pl = cmd->clone->par;
    current_sequ->l_expr = command_par_expr("l", cmd->clone);
    current_sequ->length = command_par_value("l", cmd->clone);
    current_sequ->add_pass = command_par_value("add_pass", cmd->clone);

    if (current_sequ->l_expr == NULL && sequence_length(current_sequ) == zero)
      fatal_error("missing length for sequence:", toks[aux_pos]);

    pos = name_list_pos("refpos", nl);
    if (nl->inform[pos])
      current_sequ->refpos = permbuff(pl->parameters[pos]->string);

    pos = name_list_pos("next_sequ", nl);
    if (nl->inform[pos])
      current_sequ->next_sequ = permbuff(pl->parameters[pos]->string);

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

int
aperture_count(struct sequence* sequ)
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

void
enter_sequ_reference(struct in_cmd* cmd, struct sequence* sequ)
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

void
exec_dumpsequ(struct in_cmd* cmd)
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

void
exec_save(struct in_cmd* cmd)
  /* save a sequence with all necessary parameters and sub-sequences */
{
  int i, n = 0, pos, prev = 0,
    beam_save = log_val("beam", cmd->clone),
    mad8      = log_val("mad8", cmd->clone),
    bare      = log_val("bare", cmd->clone),
    noexpr    = log_val("noexpr", cmd->clone),
    all_sequ  = 0;
  char *name, *filename, *new_name = NULL;
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

  if (nl->inform[i] == 0) {
    warning("save without file:", "ignored");
    return;
  }

  filename = pl->parameters[i]->string;
  if ((out_file = fopen(filename, "w")) == NULL) {
    warning("cannot open output file:", filename);
    return;
  }

/* get export name for sequence (newname in SAVE) HG 15.10.07 */
  i = name_list_pos("newname", nl);
  if (nl->inform[i] != 0) new_name = pl->parameters[i]->string;
/* end -- export name for sequence (newname in SAVE) HG 15.10.07 */

  warning("SAVE makes all previous USE invalid !", " ");

  pos = name_list_pos("sequence", nl);
  clp = cmd->clone->par->parameters[pos];

  if (nl->inform[pos] == 0) {  /* no sequence given, all sequences saved */
    sqo = sequences; all_sequ = 1;
  }
  else {
    sqo = new_sequence_list(20);
    for (i = 0; i < clp->m_string->curr; i++) {
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

  for (pos = 0; pos < sqo->curr; pos++) {
    sequ = sqo->sequs[pos];
    add_to_sequ_list(sequ, sql);
    /* check for inserted sequences, flatten if necessary  - HG 23.3.04 */
    c_node = sequ->start;
    while(c_node != NULL) {
      if (c_node->p_sequ != NULL) {
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

    if (beam_save && bare == 0) {
      if (mad8 == 0) save_beam(sequ, out_file, noexpr); /* only mad-X */
      else warning("when mad-8 format requested,","beam not saved");
    }
  }

  for (i = sql->curr-1; i >= 0; i--) { /* loop over sequences, get elements */
    /* set export name for sequence (newname in SAVE) HG 15.10.07 */
    if (new_name == NULL)
      strcpy(sql->sequs[i]->export_name, sql->sequs[i]->name);
    else strcpy(sql->sequs[i]->export_name, new_name);
    /* end mod HG 15.10.07 */
    c_node = sql->sequs[i]->start;
    while (c_node != NULL) {
      if ((el = c_node->p_elem) != NULL && strchr(el->name, '$') == NULL
          && strcmp(el->base_type->name, "drift") != 0) {
        while (el->base_type != el) {
          add_to_el_list(&el, 0, ell, 0);
          el = el->parent;
        }
      }
      if (c_node == sql->sequs[i]->end) break;
      c_node = c_node->next;
    }
  }

  if (all_sequ == 0) {
    while (prev < ell->curr) { /* loop over elements, get variables -
                                  recursive, since elements may be added */
      prev = ell->curr;
      for (i = n; i < ell->curr; i++)
        fill_elem_var_list(ell->elem[i], ell, varl);
      n = prev;
    }
    fill_sequ_var_list(sql, ell, varl); /* variables for positions */
  }

  if (mad8) {
    if (bare == 0) {
      write_vars_8(varl, save_select, out_file);
      write_elems_8(ell, save_select, out_file);
    }
    for (pos = 0; pos < sql->curr; pos++) {
      sequ = sql->sequs[pos];
      all_node_pos(sequ);
      sequ->ex_nodes = new_node_list(2*sequ->nodes->curr);
      expand_sequence(sequ, 0);
      export_sequ_8(sequ, save_select, out_file);
      sequ->ex_nodes = delete_node_list(sequ->ex_nodes);
    }
  }
  else {
    if (bare == 0) {
      write_vars(varl, save_select, out_file, noexpr);
      write_elems(ell, save_select, out_file, noexpr);
    }
    write_sequs(sql, save_select, out_file, noexpr);
  }

  fclose(out_file);
  if (sqo != sequences) sqo = delete_sequence_list(sqo);
  sql = delete_sequence_list(sql);
  ell = delete_el_list(ell);
  varl = delete_var_list(varl);
  current_sequ = NULL;
}

void
exec_extract(struct in_cmd* cmd)
  /* "extract" command - extract part of  sequence from marker_1 to marker_2 */
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  struct sequence *split_sequ = NULL, *part;
  int i, j;
  char *name = NULL, *refpos = NULL;
  char newname[NAME_L];
  struct node *from, *to;
  /*
    start command decoding
  */
  i = name_list_pos("sequence", nl);
  if(nl->inform[i]) /* sequence specified */
  {
    name = pl->parameters[i]->string;
    if ((j = name_list_pos(name, sequences->list)) > -1)
      split_sequ = sequences->sequs[j];
    else
    {
      warning("unknown sequence ignored:", name);
      return;
    }
  }
  i = name_list_pos("newname", nl);
  if (nl->inform[i] == 0) sprintf(newname, "%s_1", name);
  else strcpy(newname, pl->parameters[i]->string);
  i = name_list_pos("from", nl);
  if (nl->inform[i] == 0)
  {
    warning("no 'from' marker given", " ");
    return;
  }
  sprintf(c_dum->c, "%s:1", pl->parameters[i]->string);
  if ((j = name_list_pos(c_dum->c, split_sequ->nodes->list)) > -1)
    from = split_sequ->nodes->nodes[j];
  else
  {
    warning("not in sequence:", pl->parameters[i]->string);
    return;
  }
  i = name_list_pos("to", nl);
  if (nl->inform[i] == 0)
  {
    warning("no 'to' marker given", " ");
    return;
  }
  if (strchr(pl->parameters[i]->string, '$'))
  {
    warning("extract: use of internal markers forbidden:",
            pl->parameters[i]->string);
    warning("sequence extraction aborted"," ");
    return;
  }
  sprintf(c_dum->c, "%s:1", pl->parameters[i]->string);
  if ((j = name_list_pos(c_dum->c, split_sequ->nodes->list)) > -1)
    to = split_sequ->nodes->nodes[j];
  else
  {
    warning("not in sequence:", pl->parameters[i]->string);
    return;
  }
  i = name_list_pos("refpos", nl);
  if (nl->inform[i]) refpos = pl->parameters[i]->string;
  /*
    end of command decoding - action!
  */
  part = extract_sequence(newname, split_sequ, from, to, refpos);
  add_to_sequ_list(part, sequences);
}

void
expand_curr_sequ(int flag)
  /* expands the current sequence, i.e. flattens it, inserts drifts etc. */
  /* The sequence length is updated - new feature HG 26.5.03 */
{
  const char *rout_name = "expand_curr_sequ";
  struct node* c_node;
  int j;

  current_sequ->end->at_value = current_sequ->end->position = sequence_length(current_sequ);

  if (current_sequ->ex_start != NULL) {
    current_sequ->ex_nodes = delete_node_list(current_sequ->ex_nodes);
    current_sequ->ex_start = delete_node_ring(current_sequ->ex_start);
    current_sequ->orbits = delete_vector_list(current_sequ->orbits);
  }

  if (current_sequ->ex_start == NULL) {
    use_count++;
    if (occ_list == NULL)
      occ_list = new_name_list("occ_list",10000);  /* for occurrence count */
    else
      occ_list->curr = 0;
    make_occ_list(current_sequ);
    all_node_pos(current_sequ);
    current_sequ->ex_nodes = new_node_list(2*current_sequ->nodes->curr);

    /* flatten the current sequence */
    expand_sequence(current_sequ, flag);
    /* add implict drifts in current sequence */
    current_sequ->n_nodes = add_drifts(current_sequ->ex_start, current_sequ->ex_end);

    if (current_sequ->all_nodes != NULL) myfree(rout_name, current_sequ->all_nodes);
    current_sequ->all_nodes = mymalloc(rout_name, current_sequ->n_nodes * sizeof *current_sequ->all_nodes);

    c_node = current_sequ->ex_start;
    for (j = 0; j < current_sequ->n_nodes; j++) {
      current_sequ->all_nodes[j] = c_node;
      c_node = c_node->next;
    }
  }

  set_node_bv(current_sequ); /* set bv factors for all nodes */

  if (current_range)
    set_range(current_range, current_sequ);
  else {
    current_sequ->range_start = current_sequ->ex_start;
    current_sequ->range_end = current_sequ->ex_end;
  }
}

void
reset_errors(struct sequence* sequ)
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

void
add_to_sequ_list(struct sequence* sequ, struct sequence_list* sql)
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

void
reset_sector(struct sequence* sequ, int val)
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

int
restart_sequ(void)
{
  if (current_sequ == 0x0)
   {
     warning("restart_sequ","Current sequence is not set");
     return -1;
   }
  current_node = current_sequ->range_start;
  return 1;
}

void
seq_edit_main(struct in_cmd* cmd)
  /* controls sequence editing */
{
  int k = cmd->decl_start - 1;
  char** toks = cmd->tok_list->p;
  if (strcmp(toks[k], "seqedit") == 0)  seq_edit(cmd);
  else if(edit_is_on) {
         if (strcmp(toks[k], "install") == 0)  seq_install(cmd);
    else if (strcmp(toks[k], "move") == 0)     seq_move(cmd);
    else if (strcmp(toks[k], "remove") == 0)   seq_remove(cmd);
    else if (strcmp(toks[k], "cycle") == 0)    seq_cycle(cmd);
    else if (strcmp(toks[k], "flatten") == 0)  seq_flatten(edit_sequ);
    else if (strcmp(toks[k], "reflect") == 0)  seq_reflect(cmd);
    else if (strcmp(toks[k], "replace") == 0)  seq_replace(cmd);
    else if (strcmp(toks[k], "endedit") == 0)  seq_end(cmd);
  }
  else warning("seqedit command outside edit", "ignored");
}

int
set_enable(const char* type, struct in_cmd* cmd)
{
  char* name;
  struct command_parameter* cp;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  struct sequence* sequ;
  struct node* nodes[2];
  struct node* c_node;
  int pos, n, status, count = 0; // k, // not used
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
    if (get_ex_range(name, sequ, nodes) == 0) // (k = // not used
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

struct sequence* find_sequence(const char* name, struct sequence_list* sequs)
{
  if (!sequs)
    sequs = sequences;
  int pos = name_list_pos(name, sequs->list);
  return pos >= 0 ? sequs->sequs[pos] : NULL;
}
