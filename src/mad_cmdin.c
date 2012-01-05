#include "madx.h"

#if 0 // not used...
static void
dump_in_cmd(struct in_cmd* p_inp)
{
  fprintf(prt_file, "%s: type =%d, sub_type = %d, decl_start = %d\n",
          p_inp->label, p_inp->type, p_inp->sub_type, p_inp->decl_start);
  if (p_inp->cmd_def != NULL)
  {
    fprintf(prt_file, "defining command: %s\n", p_inp->cmd_def->name);
    /* dump_command(p_inp->cmd_def); */
  }
}
#endif

// public interface

struct in_cmd*
buffered_cmd(struct in_cmd* cmd)
  /* returns a buffered command if found */
{
  int k;
  if ((k = name_list_pos(cmd->tok_list->p[cmd->decl_start],
                         buffered_cmds->labels)) > -1)
    return buffered_cmds->in_cmds[k];
  else return cmd;
}

void
buffer_in_cmd(struct in_cmd* cmd)
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

struct in_cmd*
new_in_cmd(int length)
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

struct in_cmd_list*
new_in_cmd_list(int length)
{
  char rout_name[] = "new_in_cmd_list";
  struct in_cmd_list* il =
    (struct in_cmd_list*) mycalloc(rout_name,1, sizeof(struct in_cmd_list));
  strcpy(il->name, "in_cmd_list");
  il->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", il->name);
  il->curr = 0;
  il->max = length;
  il->labels = new_name_list(il->name, length);
  il->in_cmds
    = (struct in_cmd**) mycalloc(rout_name,length, sizeof(struct in_cmd*));
  return il;
}

struct in_cmd*
delete_in_cmd(struct in_cmd* cmd)
{
  char rout_name[] = "delete_in_cmd";
  if (cmd == NULL) return NULL;
  if (stamp_flag && cmd->stamp != 123456)
    fprintf(stamp_file, "d_i_c double delete --> %s\n", cmd->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", cmd->name);
  if (cmd->tok_list != NULL)
    cmd->tok_list = delete_char_p_array(cmd->tok_list, 0);
  myfree(rout_name, cmd);
  return NULL;
}

void
grow_in_cmd_list(struct in_cmd_list* p)
{
  char rout_name[] = "grow_in_cmd_list";
  struct in_cmd** c_loc = p->in_cmds;
  int j, new = 2*p->max;

  p->max = new;
  p->in_cmds
    = (struct in_cmd**) mycalloc(rout_name,new, sizeof(struct in_cmd*));
  for (j = 0; j < p->curr; j++) p->in_cmds[j] = c_loc[j];
  myfree(rout_name, c_loc);
}

void
scan_in_cmd(struct in_cmd* cmd)
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
      if ((k = name_list_pos(cmd->tok_list->p[i], cmd->cmd_def->par_names)) < 0)  /* try alias */
      {
        k = name_list_pos(alias(cmd->tok_list->p[i]), cmd->cmd_def->par_names);
        if (k < 0)
          fatal_error("illegal keyword:", cmd->tok_list->p[i]);
        break;
      }
      else if ((i = decode_par(cmd, i, n, k, log)) < 0)
      {
        fatal_error("illegal format near:", cmd->tok_list->p[-i]);
        break;
      }
      cmd->clone->par_names->inform[k] = ++cnt; /* mark parameter as read */
      if (strcmp(cmd->tok_list->p[i], "true_") == 0
          || strcmp(cmd->tok_list->p[i], "false_") == 0)
         cmd->cmd_def->par->parameters[k]->double_value =
	   cmd->clone->par->parameters[k]->double_value;
    }
    i++;
  }
}

