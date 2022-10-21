#include "madx.h"

static void
dynap_tables_create(struct in_cmd* cmd)
  /* creates the dynamic tables for DYNAP execution */
{
  int npart = stored_track_start->curr;
  struct table* t;

  (void)cmd;
  t = make_table("tracksumm", "tracksumm", tracksumm_table_cols,
                 tracksumm_table_types, 2*stored_track_start->curr);
  add_to_table_list(t, table_register);
  
  t = make_table("dynap", "dynap", dynap_table_cols, dynap_table_types, 10);
  add_to_table_list(t, table_register);
  
  t = make_table("dynaptune", "dynaptune", dynaptune_table_cols, dynaptune_table_types, npart);
  add_to_table_list(t, table_register);
  
  if (table_exists("mytracksumm")) {
  /*  hrr Sep 2021 table mytracksumm is not cleaned so pre-existence printf is not needed. 
      printf("Table mytracksumm does exist already\n"); hrr Sep 2021 */
  }
  else {
    t = make_table("mytracksumm", "mytracksumm", mytracksumm_table_cols,
                   mytracksumm_table_types, 2*stored_track_start->curr);
    add_to_table_list(t, table_register);
  }

}

// public interface

void
track_dynap(struct in_cmd* cmd)
{
  const char *rout_name = "track_dynap";
  int e_flag, flag = 2, izero = 0, turns = command_par_value("turns", cmd->clone);
  int npart = 2*stored_track_start->curr;
  int   *ibuf1, *ibuf2, *ibuf3;
  double *buf1, *buf2, *buf_dxt, *buf_dyt, *buf3, *buf4, *buf5, *buf6, *buf7, *buf8, *buf9, *buf10, *buf11;
  double orbit[6];
  struct table* t;
  int damp, quantum;

  damp = 0; quantum = 0;
  if (get_value("dynap","damp") == 0) {
    damp = 1;     fprintf(prt_file, "damp is on\n");
  }
  if (get_value("dynap","quantum") == 0) {
    quantum = 1;  fprintf(prt_file, "quantum is on\n");
  }
  set_option("damp", &damp);
  set_option("quantum", &quantum);

  if (track_is_on == 0) {
    warning("track_dynap: no TRACK command seen yet", "ignored");
    return;
  }
  if (npart == 0) {
    warning("track_dynap: no START command seen yet", "ignored");
    return;
  }
  if (turns < 64) {
    warning("track_dynap: turns cannot be < 64", "reset to 64");
    turns = 64;
  }

  // LD 2016.04.19
  zero_double(orbit0, 6);
  adjust_beam();
  probe_beam = clone_command(current_beam);
  adjust_probe_fp(track_deltap); /* sets correct gamma, beta, etc. */

  if (get_option("onepass") == 0)
    /* closed orbit and one-turn linear transfer map */
    tmrefo_(&izero,orbit0,orbit,oneturnmat);

  dynap_tables_create(cmd);

  /* allocate buffers */
  int nnode = current_sequ->n_nodes;
  ibuf1   = mymalloc_atomic(rout_name, npart   * sizeof *ibuf1);
  ibuf2   = mymalloc_atomic(rout_name, npart   * sizeof *ibuf2);
  ibuf3   = mymalloc_atomic(rout_name, nnode   * sizeof *ibuf3);
  buf_dxt = mymalloc_atomic(rout_name, npart   * sizeof *buf_dxt);
  buf_dyt = mymalloc_atomic(rout_name, npart   * sizeof *buf_dyt);
  buf1    = mymalloc_atomic(rout_name, npart   * sizeof *buf1);
  buf2    = mymalloc_atomic(rout_name, 6*npart * sizeof *buf2);
  buf3    = mymalloc_atomic(rout_name, 6*npart * sizeof *buf3);
  buf4    = mymalloc_atomic(rout_name, 6*6     * sizeof *buf4); /* eigenvectors */
  buf5    = mymalloc_atomic(rout_name, 6*npart * (turns+1) * sizeof *buf5);
  buf6    = mymalloc_atomic(rout_name, nnode   * sizeof *buf6);
  buf7    = mymalloc_atomic(rout_name, turns   * sizeof *buf7);
  buf8    = mymalloc_atomic(rout_name, 6*turns * sizeof *buf8);
  buf9    = mymalloc_atomic(rout_name, 2*turns * sizeof *buf9);
  buf10   = mymalloc_atomic(rout_name, turns   * sizeof *buf10);
  buf11   = mymalloc_atomic(rout_name, turns   * sizeof *buf11);

  trrun_(&flag, &turns,orbit0, oneturnmat, ibuf1, ibuf2, buf1, buf2,
         buf_dxt, buf_dyt, buf3, buf4, buf5, &e_flag, ibuf3, buf6);

  t = find_table("tracksumm");
  print_table(t);
  if (e_flag) {
    warning("track_dynap: particle lost before last turn,", "ignored");
    return;
  }
  
  trdynrun_(buf4, buf5, &turns, &npart, buf7, buf8, buf10, buf11, buf9);

  /* t = find_table("dynapsumm");
     print_table(t);
     if (get_option("dynap_dump")) dynap_tables_dump();  */

  probe_beam = delete_command(probe_beam);

  /* free buffers */
  myfree(rout_name, ibuf1);   myfree(rout_name, ibuf2);   myfree(rout_name, ibuf3);
  myfree(rout_name, buf_dxt); myfree(rout_name, buf_dyt);
  myfree(rout_name, buf1);    myfree(rout_name, buf2);
  myfree(rout_name, buf3);    myfree(rout_name, buf4);    myfree(rout_name, buf5);
  myfree(rout_name, buf6);    myfree(rout_name, buf7);    myfree(rout_name, buf8);
  myfree(rout_name, buf9);    myfree(rout_name, buf10);   myfree(rout_name, buf11);
}

