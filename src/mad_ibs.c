#include "madx.h"

void
pro_ibs(struct in_cmd* cmd)
  /* control for IBS module */
{
  struct command* keep_beam = current_beam;
  struct name_list* nl = current_ibs->par_names;
  struct command_parameter_list* pl = current_ibs->par;
  char *filename = NULL, *table_name = NULL;
  int pos, w_file;

  (void)cmd;
  
  if (twiss_table == NULL)
    warning("no TWISS table present","IBS command ignored");

  else {
    if ((current_beam = find_command(twiss_table->org_sequ->name, beam_list)) == NULL)
      current_beam = find_command("default_beam", beam_list);

    pos = name_list_pos("file", nl);
    if (nl->inform[pos]) {
      if ((filename = pl->parameters[pos]->string) == NULL) {
        if (pl->parameters[pos]->call_def != NULL)
          filename = pl->parameters[pos]->call_def->string;
      }
      if (filename == NULL) filename = permbuff("dummy");
      w_file = 1;
    }
    else w_file = 0;

    set_option("ibs_table", &w_file); /* fill only if output */

    if (w_file) {
      table_name = permbuff("ibs");
      ibs_table = make_table(table_name, "ibs", ibs_table_cols,
                             ibs_table_types, current_sequ->n_nodes);
      add_to_table_list(ibs_table, table_register);
    }
    
    // LD 2016.02.18: START
    // adjust_beam();
    probe_beam = clone_command(current_beam);
    // adjust_rfc(); /* sets freq in rf-cavities from probe */

    // LD 2016.02.17: BUG, depends on the previous oneturnmap and disp0 -> alpha is wrong even with dp=0
    adjust_probe(0); /* sets correct gamma, beta, etc. */
    // adjust_rfc(); /* sets freq in rf-cavities from probe */
    // LD 2016.02.18: END

    ibs_();

    if (w_file) out_table(table_name, ibs_table, filename);
    probe_beam = delete_command(probe_beam); // LD: added...
    current_beam = keep_beam;
  }
}


