#include "madx.h"

// private functions

static void
track_observe(struct in_cmd* cmd)
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  struct node* nodes[2];
  int pos;
  if (track_is_on == 0)
  {
    warning("track_observe: no TRACK command seen yet,", "ignored");
    return;
  }
  pos = name_list_pos("place", nl);
  if (get_ex_range(pl->parameters[pos]->string, current_sequ, nodes))
  {
    nodes[0]->obs_point = ++curr_obs_points;
    nodes[0]->obs_orbit = new_double_array(6);
    nodes[0]->obs_orbit->curr = 6;

    // LD 2016.04.19
    zero_double(orbit0, 6);
    adjust_beam();
    probe_beam = clone_command(current_beam);
    adjust_probe_fp(track_deltap); /* sets correct gamma, beta, etc. */

    if (get_option("onepass") == 0)
    {
      tmrefo_(&curr_obs_points,orbit0,nodes[0]->obs_orbit->a,oneturnmat);
      /* closed orbit and one-turn linear transfer map */
    }

    probe_beam = delete_command(probe_beam); // LD: added...
  }
  else
  {
    warning("track_observe: unknown place,", "ignored");
    return;
  }
}

static void
track_run(struct in_cmd* cmd)
{
  const char *rout_name = "track_run";
  int e_flag, flag = 1, izero = 0, npart = stored_track_start->curr;
  int *ibuf1, *ibuf2, *ibuf3;
  double orbit[6];
  double *buf1, *buf2, *buf_dxt, *buf_dyt, *buf3, *buf4, buf5, *buf6;
  struct table* t;

  int turns = command_par_value("turns", cmd->clone);

  if (track_is_on == 0) {
    warning("track_run: no TRACK command seen yet", "ignored");
    return;
  }

  if (npart == 0) {
    warning("track_run: no START command seen yet", "ignored");
    return;
  }

  // LD 2016.04.19
  zero_double(orbit0, 6);
  adjust_beam();
  probe_beam = clone_command(current_beam);
  adjust_probe_fp(track_deltap); /* sets correct gamma, beta, etc. */

  if (get_option("onepass") == 0) {
    tmrefo_(&izero,orbit0,orbit,oneturnmat);
    /* closed orbit and one-turn linear transfer map */
  }

  if (!command_par_value("keeptrack", cmd->clone)) { // LD: added 2016.10.14
    track_tables_delete(); /* deleting all track related tables,
                              emptying does not work because different number of particles*/
    track_tables_create(cmd);
  }

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
  buf4    = mymalloc_atomic(rout_name, 36      * sizeof *buf4);
  buf6    = mymalloc_atomic(rout_name, nnode   * sizeof *buf6);

  // run track rountine
  trrun_(&flag, &turns,orbit0, oneturnmat, ibuf1, ibuf2, buf1, buf2,
         buf_dxt, buf_dyt, buf3, buf4, &buf5, &e_flag, ibuf3, buf6);

  // summary
  t = table_register->tables[name_list_pos("tracksumm", table_register->names)];
  if (get_option("info")) print_table(t);
  if (get_option("track_dump")) track_tables_dump();

  probe_beam = delete_command(probe_beam); // LD: added 2016.02.17

  /* free buffers */
  myfree(rout_name, ibuf1);   myfree(rout_name, ibuf2);  myfree(rout_name, ibuf3);
  myfree(rout_name, buf_dxt); myfree(rout_name, buf_dyt);
  myfree(rout_name, buf1);    myfree(rout_name, buf2);  myfree(rout_name, buf3);
  myfree(rout_name, buf4);    myfree(rout_name, buf6);
}

static void
track_end(struct in_cmd* cmd)
{
  int i;
  struct node* c_node;

  (void)cmd;
  if (track_is_on == 0)
  {
    warning("track_end: no TRACK command seen yet", "ignored");
    return;
  }
  for (i = 0; i < stored_track_start->curr; i++)
    stored_track_start->commands[i] =
      delete_command(stored_track_start->commands[i]);
  stored_track_start->curr = 0;
  c_node = current_sequ->ex_start;
  while(c_node != NULL) /* clean observation points */
  {
    c_node->obs_point = 0;
    c_node->obs_orbit = delete_double_array(c_node->obs_orbit);
    if (c_node == current_sequ->ex_end)  break;
    c_node = c_node->next;
  }
  track_is_on = 0;
  fprintf(prt_file, "exit TRACK module\n\n");
}

static void
track_ripple(struct in_cmd* cmd)
{
  (void)cmd;

  warning("track_ripple routine is not implemented", "ignored");

  if (track_is_on == 0)
  {
    warning("track_ripple: no TRACK command seen yet", "ignored");
    return;
  }
}

static void
track_track(struct in_cmd* cmd)
{
  int k=0, pos, one = 1;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;

  if (current_sequ == NULL || current_sequ->ex_start == NULL)
  {
    warning("sequence not active,", "TRACK ignored");
    return;
  }
  if (attach_beam(current_sequ) == 0)
    fatal_error("TRACK - sequence without beam:", current_sequ->name);
  if (track_is_on)
  {
    warning("already inside TRACK command group,", "ignored");
    return;
  }
  track_is_on = 1;
  puts("enter TRACK module");
  if ((k = get_value(current_command->name,"onepass")) != 0)
    fprintf(prt_file, "one pass is on\n");
  set_option("onepass", &k);

  if ((k = get_value(current_command->name,"update")) != 0)
    fprintf(prt_file, "update is on\n");
  set_option("update", &k);


  if ((k = get_value(current_command->name,"damp")) != 0)
    fprintf(prt_file, "damp is on\n");
  set_option("damp", &k);

  if ((k = get_value(current_command->name,"quantum")) != 0) {
    if ((pos = name_list_pos("seed", nl)) > -1) {
      if (nl->inform[pos]) {
        int seed = command_par_value("seed", cmd->clone);
        init55(seed);
        fprintf(prt_file, "quantum is on with seed %d\n", seed);
      } else
      fprintf(prt_file, "quantum is on\n");
    } else
      fprintf(prt_file, "quantum is on\n");
  }
  set_option("quantum", &k);

  if ((k = get_value(current_command->name,"aperture")) != 0)
    fprintf(prt_file, "aperture tracking is on\n");
  set_option("aperture", &k);
  if ((k = get_value(current_command->name,"recloss")) != 0)
    fprintf(prt_file, "losses recorded\n");
  set_option("recloss", &k);
  k = get_value(current_command->name,"dump");
  set_option("track_dump", &k);
  k = get_value(current_command->name,"onetable");
  set_option("onetable", &k);
  track_deltap=get_value(current_command->name,"deltap");
  set_variable("track_deltap", &track_deltap);
  if(track_deltap != 0) fprintf(prt_file, v_format("track_deltap: %F\n"),
                                track_deltap);
  curr_obs_points = 1;  /* default: always observe at machine end */
  pos = name_list_pos("file", nl);
  if (nl->inform[pos]) set_option("track_dump", &one);
  if ((track_filename = pl->parameters[pos]->string) == NULL)
  {
    if (pl->parameters[pos]->call_def != NULL)
      track_filename = pl->parameters[pos]->call_def->string;
    else track_filename = permbuff("dummy");
  }
  track_filename = permbuff(track_filename);
  track_fileext = NULL;
  pos = name_list_pos("extension", nl);
  if ((track_fileext = pl->parameters[pos]->string) == NULL)
  {
    if (pl->parameters[pos]->call_def != NULL)
      track_fileext = pl->parameters[pos]->call_def->string;
    if (track_fileext == NULL)  track_fileext = permbuff("\0");
  }
  track_fileext = permbuff(track_fileext);
}

/**
 * copies track positions from commands to array
 * returns number of copied tracks, value <= 0 in case of error
 *
 * Used by gettrack
 */
static int
copytrackstoarray(void)
{
  int ntracks = 0;/*number of tracks : returned value */
  int n = 0; /*interator over tracks*/
  struct command* comm;
  if (trackstrarpositions)
    deletetrackstrarpositions();

  ntracks = getnumberoftracks();
  if (ntracks <= 0) {
    printf("ERROR: copytrackstoarray: number of tracks is 0! Nothing to copy!");
    return 0;
  }
  trackstrarpositions = mymalloc("copytrackstoarray", ntracks * sizeof *trackstrarpositions);

  for (n = 0; n < ntracks; n++) {
    trackstrarpositions[n] = mymalloc_atomic("copytrackstoarray", 6 * sizeof *trackstrarpositions[0]);

    comm = stored_track_start->commands[n];
    trackstrarpositions[n][0] = command_par_value("x",  comm);
    trackstrarpositions[n][1] = command_par_value("px", comm);
    trackstrarpositions[n][2] = command_par_value("y",  comm);
    trackstrarpositions[n][3] = command_par_value("py", comm);
    trackstrarpositions[n][4] = command_par_value("t",  comm);
    trackstrarpositions[n][5] = command_par_value("pt", comm);
  }

  return ntracks;
}

// public interface

int
next_start(double* x,double* px,double* y,double* py,double* t,
           double* deltae,double* fx,double* phix,double* fy,double* phiy,
           double* ft,double* phit)
  /* returns the parameters of the next particle to track;
     0 = none, else count */
{
  struct command* comm;

  if (start_cnt == stored_track_start->curr)
  {
    start_cnt = 0; return 0;
  }
  comm = stored_track_start->commands[start_cnt];
  *x  = command_par_value("x", comm);
  *px = command_par_value("px", comm);
  *y  = command_par_value("y", comm);
  *py = command_par_value("py", comm);
  *t  = command_par_value("t", comm);
  *deltae = command_par_value("pt", comm);
  *fx   = command_par_value("fx", comm);
  *phix = command_par_value("phix", comm);
  *fy   = command_par_value("fy", comm);
  *phiy = command_par_value("phiy", comm);
  *ft   = command_par_value("ft", comm);
  *phit = command_par_value("phit", comm);
  return ++start_cnt;
}

void
pro_track(struct in_cmd* cmd)
  /* controls track module */
{
  if (current_sequ == NULL || current_sequ->ex_start == NULL) {
    warning("TRACK, but no active sequence:", "ignored");
    return;
  }

       if (strcmp(cmd->tok_list->p[0], "track"   ) == 0) track_track  (cmd);
  else if (strcmp(cmd->tok_list->p[0], "dynap"   ) == 0) track_dynap  (cmd);
  else if (strcmp(cmd->tok_list->p[0], "endtrack") == 0) track_end    (cmd);
  else if (strcmp(cmd->tok_list->p[0], "observe" ) == 0) track_observe(cmd);
  else if (strcmp(cmd->tok_list->p[0], "run"     ) == 0) track_run    (cmd);
  else if (strcmp(cmd->tok_list->p[0], "ripple"  ) == 0) track_ripple (cmd);
  else if (strcmp(cmd->tok_list->p[0], "start"   ) == 0) track_start  (cmd->clone), cmd->clone_flag = 1;
}

void
track_pteigen(double* eigen)
{
  int i, j, pos;
  struct table* t;

  if ((pos = name_list_pos("trackone", table_register->names)) > -1) {
    t = table_register->tables[pos];

    if (t->header == NULL)
      t->header = new_char_p_array(45);
    else {
      // if (t->header->max - t->header->curr < 45)
      warning("Table trackone: multiple runs of track detected", "header values not updated");
      return;
    }

    table_add_header(t, "@ XC               %%le  %F", orbit0[0]);
    table_add_header(t, "@ PXC              %%le  %F", orbit0[1]);
    table_add_header(t, "@ YC               %%le  %F", orbit0[2]);
    table_add_header(t, "@ PYC              %%le  %F", orbit0[3]);
    table_add_header(t, "@ TC               %%le  %F", orbit0[4]);
    table_add_header(t, "@ PTC              %%le  %F", orbit0[5]);
    table_add_header(t, "@ EX               %%le  %F", get_value("beam", "ex"));
    table_add_header(t, "@ EY               %%le  %F", get_value("beam", "ey"));
    table_add_header(t, "@ ET               %%le  %F", get_value("beam", "et"));
    for (i = 0; i < 6; i++) {
      for (j = 0; j < 6; j++) {
        table_add_header(t, "@ E%d%d              %%le  %F", i+1, j+1, eigen[6*j+i]);
      }
    }
  }
}

void
track_start(struct command* comm)
{
  char name[FNAME_L];
  if (track_is_on == 0)
  {
    warning("track_start: no TRACK command seen yet", "ignored");
    return;
  }
  track_start_cnt++;
  strcpy(name, "start.");
  sprintf(c_dum->c, "%d", track_start_cnt);
  strcat(name, c_dum->c);
  add_to_command_list(name,comm,stored_track_start,1);
}

void
track_tables_create(struct in_cmd* cmd)
{
  int i, j, pos;
  char tab_name[NAME_L];
  struct table* t;
  int t_size;
  int turns = command_par_value("turns", cmd->clone);
  int ffile = command_par_value("ffile", cmd->clone);
  if (ffile <= 0) ffile = 1;
  t_size = turns / ffile + 10;

  if ((pos = name_list_pos("tracksumm", table_register->names)) > -1) {
    printf("Table tracksumm does exist already\n");

  }
  else {
    t = make_table("tracksumm", "tracksumm", tracksumm_table_cols,
                   tracksumm_table_types, 2*stored_track_start->curr);
    add_to_table_list(t, table_register);
  }
  if (get_option("recloss"))
  {
    if ((pos = name_list_pos("trackloss", table_register->names)) > -1) {
      printf("Table trackloss does exist already\n");
    }
    else {
      t = make_table("trackloss", "trackloss", trackloss_table_cols,
                     trackloss_table_types, stored_track_start->curr*t_size);
      add_to_table_list(t, table_register);
    }
  }
  if (get_option("onetable"))
  {
    if ((pos = name_list_pos("trackone", table_register->names)) > -1) {
      printf("Table trackone does exist already\n");
    }
    else {
      t = make_table("trackone", "trackone", trackone_table_cols,
                     trackone_table_types, stored_track_start->curr*t_size);
      add_to_table_list(t, table_register);
    }
  }
  else
  {
    for (i = 0; i < curr_obs_points; i++)
    {
      for (j = 0; j < stored_track_start->curr; j++) /* open tables */
      {
        sprintf(tab_name, "track.obs%04d.p%04d", i+1, j+1);
        t = make_table(tab_name, "trackobs", track_table_cols,
                       track_table_types, t_size);
        add_to_table_list(t, table_register);
      }
    }
  }
}

void
track_tables_delete(void)
{
  int j;
  /*
  tracksumm
  */
  exec_delete_table("tracksumm");
  for (j = table_register->names->curr - 1; j >= 0; j--)
  {

    if (   strstr(table_register->names->names[j], "track.obs")
        || (strcmp(table_register->names->names[j], "trackone") == 0)
        || (strcmp(table_register->names->names[j], "trackloss") == 0))
    {
      exec_delete_table(table_register->names->names[j]);
    }
  }
}


void
track_tables_dump(void)
{
  int j;
  for (j = 0; j < table_register->names->curr; j++)
  {
    if (strstr(table_register->names->names[j], "track.obs")
        || strcmp(table_register->names->names[j], "trackone") == 0)
    {
      strcpy(l_wrk->c, track_filename);
      strcat(l_wrk->c, &table_register->names->names[j][5]);
      strcat(l_wrk->c, track_fileext);
      out_table("track", table_register->tables[j], l_wrk->c);
    }
  }
}

int
getnumberoftracks(void)
{
/*returns number of input tracks */
  if (stored_track_start == 0x0)
  {
    return 0;
  }

  return stored_track_start->curr;

}

int
getcurrentcmdname(char* string)
{
  if (current_command == 0x0)
  {
    return 0;
  }

  strcpy(string, current_command->name);
  return strlen(current_command->name);

}

const char* getcurrentelementname()
{
/*returns name of the current element 
  Used in rviewer plugin
*/

  if (current_node == 0x0)
  {
    return 0x0;
  }

  return current_node->name;

}

int
gettrack(int* nt, double* x,double* px,double* y,double* py,double* t,double* pt)
{
  /* returns the parameters of track n;
     0 = none, else count */
  int n = *nt - 1;

  if ( trackstrarpositions == 0x0 )
  {
    copytrackstoarray();
  }
  if ( (n<0) || (n >= stored_track_start->curr) )
  {
    printf("gettrack: track number %d out of range",n);
    return 1;
  }


  *x      = trackstrarpositions[n][0];
  *px     = trackstrarpositions[n][1];
  *y      = trackstrarpositions[n][2];
  *py     = trackstrarpositions[n][3];
  *t      = trackstrarpositions[n][4];
  *pt     = trackstrarpositions[n][5];
  return 0;
}

void
deletetrackstrarpositions(void)
{
  /* deletes the array with track positions */
  int i;
  for ( i = 0; i < stored_track_start->curr; i++)
  {
    myfree("deletetrackstrarpositions",trackstrarpositions[i]);
  }
  myfree("deletetrackstrarpositions",trackstrarpositions);

  trackstrarpositions = 0x0;
}

