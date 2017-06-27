#include "madx.h"

// TG: start_iter_select/fetch_node_select can be used to replace/implement
// get_select_ex_ranges and the same pattern in a few other places.
struct select_iter
{
  int i_seq, full;
  struct command* cmd;
  struct sequence* sequ;
  struct sequence_list* sequs;
  const char* range;
  struct node *node, *range_end;
};

struct select_iter*
start_iter_select(struct command* cmd, struct sequence_list* sequs, struct sequence* sequ)
{
  const char* name;

  if (!sequs)
    sequs = sequences;

  if (sequs && !sequ && (name = command_par_string("sequence", cmd))) {
    sequ = find_sequence(name, sequs);
    sequs = NULL;
    if (!sequ) {
      warning("unknown sequence, skipped select: ", name);
      return NULL;
    }
  }

  struct select_iter* it = mycalloc("start_iter_select", 1, sizeof(struct select_iter));
  it->cmd = cmd;
  it->sequ = sequ ? sequ : sequs->sequs[0];
  it->full = log_val("full", cmd);    // `full` means "unfiltered"
  it->sequs = sequs;
  it->range = it->full ? NULL : command_par_string("range", cmd);

  return it;
}

int
fetch_node_select(struct select_iter* it, struct node** node, struct sequence** seq)
{
  struct node* nodes[2];

  if (!it)
    return 0;

  for (;;) {
    if (it->node) {
      it->node = it->node == it->range_end ? NULL : it->node->next;
    }
    else if (it->range) {
      if (get_ex_range(it->range, it->sequ, nodes)) {
        it->node = nodes[0];
        it->range_end = nodes[1];
      }
    }
    else {
      it->node = it->sequ->ex_start;
      it->range_end = it->sequ->ex_end;
    }

    if (!it->node && it->sequs && it->i_seq != it->sequs->curr) {
      it->sequ = it->sequs->sequs[++it->i_seq];
      continue;
    }

    if (!it->node || it->full || pass_select(it->node->p_elem->name, it->cmd)) {
      if (seq) *seq = it->sequ;
      if (node) *node = it->node;
      return it->node != NULL;
    }
  }
}

static int
get_select_ex_ranges(struct sequence* sequ, struct command_list* select, struct node_list* s_ranges)
  /* makes a list of nodes of an expanded sequence that pass the range
     selection */
{
  /*returns 0 if invalid sequence pointer
    1 if nodes in s_ranges (including 0) */
  struct node* c_node;
  if (sequ == NULL) return 0;
  s_ranges->curr = 0;
  s_ranges->list->curr = 0;
  for (int i = 0; i < select->curr; i++)
  {
    struct select_iter* it = start_iter_select(select->commands[i], NULL, sequ);
    while (fetch_node_select(it, &c_node, NULL)) {
      add_to_node_list(c_node, 0, s_ranges);
    }
    if (it->full) break;
  }
  return 1;
}

// public interface

int
pass_select(char* name, struct command* sc)
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

int
pass_select_list(char* name, struct command_list* cl)
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

int
get_select_ranges(struct sequence* sequ, struct command_list* select, struct node_list* s_ranges)
  /* makes a list of nodes of a sequence that pass the range selection */
{
  struct name_list* nl;
  struct command_parameter_list* pl;
  char* name;
  char full_range[] = "#s/#e";
  int i, pos; //, k // not used
  struct node* c_node;
  struct node* nodes[2];
  for (i = 0; i < select->curr; i++)
  {
    nl = select->commands[i]->par_names;
    pl = select->commands[i]->par;
    pos = name_list_pos("range", nl);
    if (pos > -1 && nl->inform[pos])  /* parameter has been read */
      name = pl->parameters[pos]->string;
    else name = full_range;
    if (get_range(name, sequ, nodes) > 0) // (k = not used
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
  return s_ranges->curr;
}

void
get_select_t_ranges(struct command_list* select, struct command_list* deselect, struct table* t)
  /* makes a list of table rows that pass the range selection and
     subsequent deselection */
{
  int rows[2];
  struct name_list* nl;
  struct command_parameter_list* pl;
  int i, pos;
  s_range->curr = 0; e_range->curr = 0;
  if (select != NULL)
  {
    for (i = 0; i < select->curr; i++)
    {
      nl = select->commands[i]->par_names;
      pl = select->commands[i]->par;
      pos = name_list_pos("range", nl);
      if (pos > -1 && nl->inform[pos]  /* parameter has been read */
          && get_table_range(pl->parameters[pos]->string, t, rows)
          && (rows[0] <= rows[1]))
      {
        if (s_range->max == s_range->curr) grow_int_array(s_range);
        if (e_range->max == e_range->curr) grow_int_array(e_range);
        s_range->i[s_range->curr++] = rows[0];
        e_range->i[e_range->curr++] = rows[1];
      }
      else
      {
        if (s_range->max == s_range->curr) grow_int_array(s_range);
        if (e_range->max == e_range->curr) grow_int_array(e_range);
        s_range->i[s_range->curr++] = 0;
        e_range->i[e_range->curr++] = t->curr - 1;
      }
    }
  }
  if (deselect != NULL)
  {
    for (i = 0; i < deselect->curr; i++)
    {
      nl = deselect->commands[i]->par_names;
      pl = deselect->commands[i]->par;
      pos = name_list_pos("range", nl);
      if (pos > -1 && nl->inform[pos]  /* parameter has been read */
          && get_table_range(pl->parameters[pos]->string, t, rows)
          && (rows[0] <= rows[1]))
      {
        if (sd_range->max == sd_range->curr) grow_int_array(sd_range);
        if (ed_range->max == ed_range->curr) grow_int_array(ed_range);
        sd_range->i[sd_range->curr++] = rows[0];
        ed_range->i[ed_range->curr++] = rows[1];
      }
    }
  }
}

void
set_selected_errors(void)
{
  int i;
  if (get_select_ex_ranges(current_sequ, error_select, selected_ranges) != 0)
    for (i = 0; i < selected_ranges->curr; i++)
      selected_ranges->nodes[i]->sel_err = 1;
}

void
set_range(char* range, struct sequence* sequ)
{
  struct node* nodes[2];
  current_sequ->range_start = current_sequ->ex_start;
  current_sequ->range_end = current_sequ->ex_end;
  if (get_ex_range(range, sequ, nodes) == 0) return;
  current_sequ->range_start = nodes[0];
  current_sequ->range_end = nodes[1];
}

void
set_sector(void)
{
  int i;
  if (sector_select->curr == 0) reset_sector(current_sequ, 1);
  else
  {
    sector_ranges->curr = 0; sector_ranges->list->curr = 0;
    if (get_select_ex_ranges(current_sequ, sector_select, sector_ranges) != 0)
      for (i = 0; i < sector_ranges->curr; i++)
        sector_ranges->nodes[i]->sel_sector = 1;
  }
}

static int
get_interval_sub_range(
             const char* range, struct node_list* sequ, struct node** nodes,
             struct node* start, struct node* stop)
  /* returns start and end node (nodes[0] and nodes[1])
     of a range in the full expanded sequence */
// LD: Same function as get_table_range
// TG: can be merged by making `table::node_nm` a `name_list` (needs to be
// non-unique unless shared drifts are reworked).
{
  if (sequ == NULL) return 0;

  char buf[5*NAME_L], *c[2];
  stolower(strcpy(buf, range));
  c[0] = strtok(buf, "/");
  int n = (c[1] = strtok(NULL,"/")) ? 2 : 1;

  for (int i = 0; i < n; i++) {
    if (!(nodes[i] = find_node_by_name(c[i], sequ, start, stop))) {
      warning("illegal expand range ignored:", range);
      return 0;
    }
  }

  if (n == 1) nodes[1] = nodes[0];
  return n;
}

int
get_ex_range(const char* range, struct sequence* sequ, struct node** nodes)
  /* returns start and end node (nodes[0] and nodes[1])
     of a range in the full expanded sequence */
{
  return sequ ? get_interval_sub_range(
      range, sequ->ex_nodes, nodes,
      sequ->ex_start, sequ->ex_end) : 0;
}

int
get_sub_range(const char* range, struct sequence* sequ, struct node** nodes)
{
  /* returns start and end node (nodes[0] and nodes[1])
     of a range between range_start and range_end of an expanded sequence */
  return sequ ? get_interval_sub_range(
      range, sequ->ex_nodes, nodes,
      sequ->range_start, sequ->range_end) : 0;
}

int
get_range(const char* range, struct sequence* sequ, struct node** nodes)
  /* returns start and end node (nodes[0] and nodes[1])
     of a range in the non-expanded sequence */
{
  return sequ ? get_interval_sub_range(
      range, sequ->nodes, nodes,
      sequ->start, sequ->end) : 0;
}

static int
par_defined(const char* str, const char* none)
{
  return str && strcmp(str, "") != 0 && strcmp(str, none) != 0;
}

static int
has_filter_condition(struct command* cmd)
{
  return par_defined(command_par_string("sequence", cmd), none)
      || par_defined(command_par_string("range",    cmd), "#s/#e")
      || par_defined(command_par_string("class",    cmd), none)
      || par_defined(command_par_string("pattern",  cmd), "any");
}

static int
store_select_command(struct in_cmd* cmd, struct command_list** plist, const char* name, int clear_all)
{
  // TG: `clear_all` forces unqualified clear (old behaviour). With `!clear_all`,
  // behave like a DESELECT.
  if (strcmp(cmd->clone->name, "select") == 0 &&        // !deselect
      log_val("clear", cmd->clone) &&
      (clear_all || !has_filter_condition(cmd->clone)))
  {
    delete_command_list(*plist);
    *plist = new_command_list(name, 10);
  }
  else {
    if ((*plist)->curr == (*plist)->max) {
      grow_command_list(*plist);
    }
    (*plist)->commands[(*plist)->curr++] = cmd->clone;
    cmd->clone_flag = 1;
  }
  return (*plist)->curr;
}

void
store_deselect(struct in_cmd* cmd)
{
  char* flag_name;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  struct command_list* dscl;
  int pos = name_list_pos("flag", nl);
  if (nl->inform[pos] == 0 ||
      (flag_name = pl->parameters[pos]->string) == NULL)
  {
    warning("no FLAG specified", "ignored");
    return;
  }
  if (strcmp(flag_name, "seqedit") == 0)
  {
  }
  else if (strcmp(flag_name, "error") == 0)
  {
  }
  else if (strcmp(flag_name, "makethin") == 0)
  {
  }
  else if (strcmp(flag_name, "save") == 0)
  {
  }
  else if (strcmp(flag_name, "sectormap") == 0)
  {
  }
  else if (strcmp(flag_name, "interpolate") == 0)
  {
    store_select_command(cmd, &interp_select, "interp_select", 0);
    select_interp(cmd->clone);
  }
  else /* store deselect for all tables */
  {
    if ((dscl = find_command_list(flag_name, table_deselect)) == NULL)
    {
      dscl = new_command_list("deselect", 10);
      add_to_command_list_list(flag_name, dscl, table_deselect);
    }
    if (log_val("clear", cmd->clone))
    {
      dscl = new_command_list("deselect", 10);
      add_to_command_list_list(flag_name, dscl, table_deselect);
    }
    else
    {
      if (dscl->curr == dscl->max) grow_command_list(dscl);
      dscl->commands[dscl->curr++] = cmd->clone;
      cmd->clone_flag = 1; /* do not drop */
    }
  }
}

void
store_select(struct in_cmd* cmd)
{
#ifdef _ONLINE
  char *sdds_pattern;
  char sdum[1000];
#endif
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
    if (store_select_command(cmd, &seqedit_select, "seqedit_select", 1) == 0) {
      selected_ranges->curr = 0;
      selected_ranges->list->curr = 0;
    }
  }
  else if (strcmp(flag_name, "error") == 0)
  {
    if (store_select_command(cmd, &error_select, "error_select", 1) == 0) {
      selected_ranges->curr = 0;
      selected_ranges->list->curr = 0;
      reset_errors(current_sequ);
    }
  }
  else if (strcmp(flag_name, "sdds") == 0)
#ifdef _ONLINE
  {
    if (log_val("clear", cmd->clone))
    {
      if (get_option("debug")) printf("here CLEAR all selected arrays\n");
      if(sdds_pat != NULL) sdds_pat = delete_char_p_array(sdds_pat, 1);
    }
    else
    {
      if(sdds_pat == NULL) sdds_pat = new_char_p_array(100);
      if (get_option("debug")) printf("current selections %d\n",sdds_pat->curr);
      sdds_pattern = command_par_string("pattern",cmd->clone);
      if(sdds_pattern != NULL) {
        if (get_option("debug")) printf("pattern in select: %s\n",sdds_pattern);
        strcpy(sdum,sdds_pattern);
        sdds_pat->p[sdds_pat->curr++] = tmpbuff(sdum);
      }
    }
  }
#else
  {
    warning("ignored, only available in ONLINE model:", "FLAG=SDDS");
  }
#endif

  else if (strcmp(flag_name, "makethin") == 0)
  {
    store_select_command(cmd, &slice_select, "slice_select", 1);
  }
  else if (strcmp(flag_name, "save") == 0)
  {
    store_select_command(cmd, &save_select, "save_select", 1);
  }
  else if (strcmp(flag_name, "sectormap") == 0)
  {
    if (sector_ranges == NULL)   sector_ranges = new_node_list(10000);
    if (store_select_command(cmd, &sector_select, "sector_select", 1) == 0) {
      sector_ranges->curr = 0;
      sector_ranges->list->curr = 0;
    }
  }
  else if (strcmp(flag_name, "interpolate") == 0)
  {
    store_select_command(cmd, &interp_select, "interp_select", 0);
    select_interp(cmd->clone);
  }
  else /* store select for all tables */
  {
    if ((scl = find_command_list(flag_name, table_select)) == NULL)
    {
      scl = new_command_list("select", 10);
      add_to_command_list_list(flag_name, scl, table_select);
    }
    if (log_val("clear", cmd->clone))
    {
      scl = new_command_list("select", 10);
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


