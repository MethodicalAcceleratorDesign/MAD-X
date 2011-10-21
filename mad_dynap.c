#include "madx.h"

void
dynap_tables_create(struct in_cmd* cmd)
  /* creates the dynamic tables for DYNAP execution */
{
  int npart = stored_track_start->curr;

  struct table* t;
  t = make_table("tracksumm", "tracksumm", tracksumm_table_cols,
                 tracksumm_table_types, 2*stored_track_start->curr);
  add_to_table_list(t, table_register);
  t = make_table("dynap", "dynap", dynap_table_cols, dynap_table_types, 10);
  add_to_table_list(t, table_register);
  t = make_table("dynaptune", "dynaptune", dynaptune_table_cols,
                 dynaptune_table_types, npart);
  add_to_table_list(t, table_register);
}

void
track_dynap(struct in_cmd* cmd)
{
  char rout_name[] = "track_dynap";
  int e_flag, flag = 2, izero = 0,
    turns = command_par_value("turns", cmd->clone),
    npart = 2*stored_track_start->curr;
  int *ibuf1, *ibuf2, *ibuf3;
  double orbit[6];
  double *buf1, *buf2, *buf_dxt, *buf_dyt, *buf3, *buf4, *buf5, *buf6,
    *buf7, *buf8, *buf9, *buf10, *buf11;
  struct table* t;
  int kopt01,kopt02;


  kopt02=0;
  kopt01 = get_value("dynap","damp");
  if (kopt01 == 0) {
    kopt02=1;
    fprintf(prt_file, "damp is on\n");}
  set_option("damp", &kopt02);



  kopt02=0;
  kopt01 = get_value("dynap","quantum");
  if (kopt01 == 0) {
    kopt02=1;
    fprintf(prt_file, "quantum is on\n");}
  set_option("quantum", &kopt02);




  if (track_is_on == 0)
  {
    warning("track_dynap: no TRACK command seen yet", "ignored");
    return;
  }
  if (npart == 0)
  {
    warning("track_dynap: no START command seen yet", "ignored");
    return;
  }
  if (turns < 64)
  {
    warning("track_dynap: turns cannot be < 64", "set to 64");
    turns = 64;
  }
  adjust_beam();
  if (probe_beam) probe_beam = delete_command(probe_beam);
  probe_beam = clone_command(current_beam);
  adjust_probe(track_deltap); /* sets correct gamma, beta, etc. */
  adjust_rfc(); /* sets freq in rf-cavities from probe */
  zero_double(orbit0, 6);
  zero_double(oneturnmat, 36);
  if (get_option("onepass") == 0)
  {
    tmrefo_(&izero,orbit0,orbit,oneturnmat);
    /* closed orbit and one-turn linear transfer map */
  }
  dynap_tables_create(cmd);
  /* allocate buffers */
  ibuf1 = (int*) mymalloc(rout_name,npart*sizeof(int));
  ibuf2 = (int*) mymalloc(rout_name,npart*sizeof(int));
  ibuf3 = (int*) mymalloc(rout_name, current_sequ->n_nodes*sizeof(int));
  buf1 = (double*) mymalloc(rout_name,npart*sizeof(double));
  buf2 = (double*) mymalloc(rout_name,6*npart*sizeof(double));
  buf_dxt = (double*) mymalloc(rout_name,npart*sizeof(double));
  buf_dyt = (double*) mymalloc(rout_name,npart*sizeof(double));
  buf3 = (double*) mymalloc(rout_name,6*npart*sizeof(double));
  buf4 = (double*) mymalloc(rout_name,36*sizeof(double)); /* eigenvectors */
  buf5 = (double*) mymalloc(rout_name,6*npart*(turns+1)*sizeof(double));
  buf6 = (double*) mymalloc(rout_name, current_sequ->n_nodes*sizeof(double));
  buf7 = (double*) mymalloc(rout_name, turns*sizeof(double));
  buf8 = (double*) mymalloc(rout_name, 6*turns*sizeof(double));
  buf9 = (double*) mymalloc(rout_name, 2*turns*sizeof(double));
  buf10 = (double*) mymalloc(rout_name, turns*sizeof(double));
  buf11 = (double*) mymalloc(rout_name, turns*sizeof(double));
  trrun_(&flag, &turns,orbit0, oneturnmat, ibuf1, ibuf2, buf1, buf2,
         buf_dxt, buf_dyt, buf3, buf4, buf5, &e_flag, ibuf3, buf6);
  t =
    table_register->tables[name_list_pos("tracksumm", table_register->names)];
  print_table(t);
  if (e_flag)
  {
    warning("track_dynap: particle lost before last turn,", "ignored");
    return;
  }
  dynap_(buf4, buf5, &turns, &npart, buf7, buf8, buf9, buf10, buf11);
  /*
    table_register->tables[name_list_pos("dynapsumm", table_register->names)];
    print_table(t);
    if (get_option("dynap_dump")) dynap_tables_dump();
  */
  /* free buffers */
  myfree(rout_name, ibuf1); myfree(rout_name, ibuf2);
  myfree(rout_name, ibuf3); myfree(rout_name, buf1); myfree(rout_name, buf2);
  myfree(rout_name, buf_dxt); myfree(rout_name, buf_dyt);
  myfree(rout_name, buf3); myfree(rout_name, buf4); myfree(rout_name, buf5);
  myfree(rout_name, buf6); myfree(rout_name, buf7); myfree(rout_name, buf8);
  myfree(rout_name, buf9); myfree(rout_name, buf10);
  myfree(rout_name, buf11);
}

