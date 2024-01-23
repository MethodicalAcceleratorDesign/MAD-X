#include "madx.h"

void
pro_ibs(struct in_cmd* cmd)
  /* control for IBS module */
{
  struct command* keep_beam = current_beam;
  char *filename = NULL, *table_name = NULL;
  int k, w_file;

  (void)cmd;
  
  if (twiss_table == NULL)
    warning("no TWISS table present","IBS command ignored");

  else {
    if ((current_beam = find_command(twiss_table->org_sequ->name, beam_list)) == NULL)
      current_beam = find_command("default_beam", beam_list);

    w_file = command_par_string_user2("file", current_ibs, &filename);
    if (w_file && !filename) {      // TG: should be impossible
      filename = permbuff("dummy");
    }

    /* declare and create the IBS table */
    table_name = permbuff("ibs");
    ibs_table = make_table(table_name, "ibs", ibs_table_cols,
                            ibs_table_types, current_sequ->n_nodes);
    add_to_table_list(ibs_table, table_register);

    // LD 2016.04.19
    adjust_beam();
    probe_beam = clone_command(current_beam);
    adjust_probe_fp(0); /* sets correct gamma, beta, etc. */

    ibs_();

    /* write the IBS table to file if asked by user */
    if (w_file) out_table(table_name, ibs_table, filename);
    probe_beam = delete_command(probe_beam); // LD: added...
    current_beam = keep_beam;
  }
}
