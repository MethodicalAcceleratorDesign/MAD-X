#include "madx.h"

static void
fill_twiss_header_ptc(struct table* t, double ptc_deltap)
  /* puts beam parameters etc. at start of twiss table */
{
  int i, h_length = 39+3+1+1+6+4; /* change when adding header lines ! - last 6 for the closed orbit */
  double dtmp;
  /*  struct table* s; */
  char tmp[16];

  // int returnStatus; // not used
  int row;

  if (t == NULL) return;
  /* ATTENTION: if you add header lines, augment h_length accordingly */
  if (t->header == NULL)  t->header = new_char_p_array(h_length);
  strcpy(tmp, t->org_sequ->name);
  sprintf(c_dum->c, v_format("@ SEQUENCE         %%%02ds \"%s\""),
          strlen(tmp),stoupper(tmp));
  t->header->p[t->header->curr++] = tmpbuff(c_dum->c);
  i = get_string("beam", "particle", tmp);
  sprintf(c_dum->c, v_format("@ PARTICLE         %%%02ds \"%s\""),
          i, stoupper(tmp));
  t->header->p[t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = get_value("beam", "mass");
  sprintf(c_dum->c, v_format("@ MASS             %%le  %F"), dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = get_value("beam", "charge");
  sprintf(c_dum->c, v_format("@ CHARGE           %%le  %F"), dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = get_value("beam", "energy");
  sprintf(c_dum->c, v_format("@ ENERGY           %%le  %F"), dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = get_value("beam", "pc");
  sprintf(c_dum->c, v_format("@ PC               %%le  %F"), dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = get_value("beam", "gamma");
  sprintf(c_dum->c, v_format("@ GAMMA            %%le  %F"), dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = get_value("beam", "kbunch");
  sprintf(c_dum->c, v_format("@ KBUNCH           %%le  %F"), dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = get_value("beam", "bcurrent");
  sprintf(c_dum->c, v_format("@ BCURRENT         %%le  %F"), dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = get_value("beam", "sige");
  sprintf(c_dum->c, v_format("@ SIGE             %%le  %F"), dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = get_value("beam", "sigt");
  sprintf(c_dum->c, v_format("@ SIGT             %%le  %F"), dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = get_value("beam", "npart");
  sprintf(c_dum->c, v_format("@ NPART            %%le  %F"), dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = get_value("beam", "ex");
  sprintf(c_dum->c, v_format("@ EX               %%le  %F"), dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = get_value("beam", "ey");
  sprintf(c_dum->c, v_format("@ EY               %%le  %F"), dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = get_value("beam", "et");
  sprintf(c_dum->c, v_format("@ ET               %%le  %F"), dtmp);
  t->header->p[t->header->curr++] = tmpbuff(c_dum->c);
  sprintf(c_dum->c, v_format("@ DELTAP           %%le  %F"), ptc_deltap);
  t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

  /* one-turn information gets computed iff ptc_twiss_summary set to 1 in madx_ptc_twiss.f90 */
  if (get_option("ptc_twiss_summary") != zero){

    /* now retreive all pieces of information from the ptc_twiss*/

    row = 1; /* this particular table has only one row filled-in */

    /* length of the machine */
    double_from_table_row("ptc_twiss_summary","length",&row,&dtmp); // returnStatus = not used
    /* returnStatus should always be equal to zero */
    sprintf(c_dum->c, v_format("@ LENGTH           %%le  %F"), dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);


    /* momentum compaction factor, phase-slip factor and energy transition */
    double_from_table_row("ptc_twiss_summary","alpha_c", &row, &dtmp); // returnStatus = not used
    /* returnStatus should always be equal to zero */
    sprintf(c_dum->c, v_format("@ ALPHA_C          %%le  %F"), dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);
    
    /* momentum compaction factor first order derivative w.r.t delta-p/p */
    double_from_table_row("ptc_twiss_summary","alpha_c_p", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ ALPHA_C_P        %%le  %F"),dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);
    
    
     /* WARNING when restoring the following two lines don't forget to replace 39+1  by 39+2 for h_length */
     /* momentum compaction factor second order derivative w.r.t delta-p/p */
     /* uncomment the following once computation of alpha_c_p2 is reliable */
    double_from_table_row("ptc_twiss_summary","alpha_c_p2", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ ALPHA_C_P2       %%le  %F"),dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c); 
    
     /* WARNING when restoring the following two lines don't forget to replace 39+2  by 39+3 for h_length */
    /* momentum compaction factor third order derivative w.r.t delta-p/p */
    /* uncomment the following once computation of alpha_c_p3 is reliable */
    double_from_table_row("ptc_twiss_summary","alpha_c_p3", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ ALPHA_C_P3       %%le  %F"),dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);
       
    double_from_table_row("ptc_twiss_summary","eta_c", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ ETA_C            %%le  %F"), dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

    double_from_table_row("ptc_twiss_summary","gamma_tr", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ GAMMA_TR         %%le  %F"), dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);
    
    /* tunes and chromaticities */
    double_from_table_row("ptc_twiss_summary","q1", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ Q1               %%le  %F"), dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

    double_from_table_row("ptc_twiss_summary","q2", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ Q2               %%le  %F"), dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

    double_from_table_row("ptc_twiss_summary","dq1", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ DQ1              %%le  %F"), dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

    double_from_table_row("ptc_twiss_summary","dq2", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ DQ2              %%le  %F"), dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

    /* 26 november 2009 */
    double_from_table_row("ptc_twiss_summary","qs", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ QS               %%le  %F"), dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);    


    /* extremas of the beta-function */
    double_from_table_row("ptc_twiss_summary","beta_x_min", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ BETA_X_MIN       %%le  %F"), dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

    double_from_table_row("ptc_twiss_summary","beta_x_max", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ BETA_X_MAX       %%le  %F"), dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

    double_from_table_row("ptc_twiss_summary","beta_y_min", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ BETA_Y_MIN       %%le  %F"), dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

    double_from_table_row("ptc_twiss_summary","beta_y_max", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ BETA_Y_MAX       %%le  %F"), dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

    /* now for the 6 closed orbits */
    double_from_table_row("ptc_twiss_summary","orbit_x", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ ORBIT_X          %%le  %F"),dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

    double_from_table_row("ptc_twiss_summary","orbit_px", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ ORBIT_PX         %%le  %F"),dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

    double_from_table_row("ptc_twiss_summary","orbit_y", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ ORBIT_Y          %%le  %F"),dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

    double_from_table_row("ptc_twiss_summary","orbit_py", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ ORBIT_PY         %%le  %F"),dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

    double_from_table_row("ptc_twiss_summary","orbit_pt", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ ORBIT_PT         %%le  %F"),dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

    double_from_table_row("ptc_twiss_summary","orbit_-cT", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ ORBIT_-CT        %%le  %F"),dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

/* orbits RMS */
    double_from_table_row("ptc_twiss_summary","xcorms", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ XCORMS           %%le  %F"),dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

    double_from_table_row("ptc_twiss_summary","pxcorms", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ PXCORMS          %%le  %F"),dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

    double_from_table_row("ptc_twiss_summary","ycorms", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ YCORMS           %%le  %F"),dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

    double_from_table_row("ptc_twiss_summary","pycorms", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ PYCORMS          %%le  %F"),dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

/* orbits MAX */
    double_from_table_row("ptc_twiss_summary","xcomax", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ XCOMAX           %%le  %F"),dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

    double_from_table_row("ptc_twiss_summary","pxcomax", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ PXCOMAX          %%le  %F"),dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

    double_from_table_row("ptc_twiss_summary","ycomax", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ YCOMAX           %%le  %F"),dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

    double_from_table_row("ptc_twiss_summary","pycomax", &row, &dtmp); // returnStatus = not used
    sprintf(c_dum->c, v_format("@ PYCOMAX          %%le  %F"),dtmp);
    t->header->p[t->header->curr++] = tmpbuff(c_dum->c);

  }
}

static int
pro_ptc_select_checkpushtable(struct in_cmd* cmd, struct int_array** tabnameIA, struct int_array** colnameIA)
{
  struct command_parameter_list* c_parameters= cmd->clone->par;
  struct name_list*              c_parnames  = cmd->clone->par_names;
  struct table*                  aTable      = 0x0;
  int                            pos         = 0;
  char*                          tablename   = 0x0;
  char*                          columnname  = 0x0;


  /*extracts column specified by the user*/
  pos        = name_list_pos("column", c_parnames);
  if (pos < 0)
  {
    printf("mad_ptc.c: pro_ptc_select: column parameter does not exist.\n");
    return 5;
  }

  columnname  = c_parameters->parameters[pos]->string;
  if ( columnname == 0x0 )
  {
/*    warning("mad_ptc.c: pro_ptc_select: Column name is empty: ", "ignored");*/
    return 6;
  }

  *colnameIA = new_int_array(1+strlen(columnname));
  conv_char(columnname,*colnameIA);


  /*extracts table specified by the user*/
  pos   = name_list_pos("table", c_parnames);
  if (pos < 0)
  {
    printf("mad_ptc.c: pro_ptc_select: table parameter does not exist.\n");
    return 1;
  }

  tablename  = c_parameters->parameters[pos]->string;
  if ( tablename == 0x0 )
  {
    return -1;/*This means that table name was not specified at all*/
  }
  if ( tablename[0] == 0 )
  {
    return -1; /*This means that table name was not specified at all*/
  }
  pos = name_list_pos(tablename, table_register->names);
  if (pos < 0)
  {
    printf("mad_ptc.c: pro_ptc_select: table <<%s>> does not exist: Create table first\n",tablename);
    return 3;
  }

  aTable = table_register->tables[pos];
  if (aTable == 0x0)
  {
    printf("mad_ptc.c: pro_ptc_select: table <<%s>> is NULL: \n",tablename);
    return 4;
  }


  /*checks if the specified column exists*/
  pos = name_list_pos(columnname,aTable->columns);
  if (pos < 0)
  {
    error("mad_ptc.c: pro_ptc_select","Can not find column named <<%s>> in table <<%s>>.",
          columnname,aTable->name);
    return 7;
  }

  pos = name_list_pos("name",aTable->columns);
  if (pos < 0)
  {
    warning("mad_ptc.c: pro_ptc_selectaTable->name: There is no column named <<name>> in table <<%s>>.",aTable->name);
    return 8;
  }

  /*so none of the columns is filled */
  aTable->org_cols = aTable->num_cols;


  *tabnameIA = new_int_array(1+strlen(tablename));
  conv_char(tablename,*tabnameIA);

  return 0;
}

// public interface

int
minimum_acceptable_order(void)
{
  return min_order;
}

int
select_ptc_idx(void)
{
  struct table* t;
  int pos;

  if ((pos = name_list_pos("normal_results", table_register->names)) > -1)
  {
    t = table_register->tables[pos];
    return t->curr;
  }
  else
    return pos;
}

void
ptc_track_end(void)
{
  int i;
  struct node* c_node;
  if (track_is_on == 0)
  {
    warning("ptc_track_end: no PTC_TRACK command seen yet", "ignored");
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
  curr_obs_points = 1;
  track_is_on = 0;
}

void
ptc_track_observe(struct in_cmd* cmd)
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  struct node* nodes[2];
  int pos;

  pos = name_list_pos("place", nl);
  if (get_ex_range(pl->parameters[pos]->string, current_sequ, nodes))
  {
    nodes[0]->obs_point = ++curr_obs_points;
    printf("obs_points: %d \n",curr_obs_points);
  }
  else
  {
    warning("ptc_track_observe: unknown place,", "ignored");
    return;
  }
}
/*________________________________________________________________*/
void
ptc_putbeambeam(struct in_cmd* cmd)
{
/*Installs beam beam interaction on a given integration step*/
/*might be defined by global s or element name and local s*/

  (void)cmd;
  w_ptc_putbeambeam_();
}


/*________________________________________________________________*/
void
ptc_dumpmaps(struct in_cmd* cmd)
/*Dumps PTC map for each element in the current sequence*/
{
  (void)cmd;
  w_ptc_dumpmaps_();
}

void
ptc_oneturnmap(struct in_cmd* cmd)
{
  (void)cmd;
}

void
pro_ptc_twiss(void)
  /* controls ptc_twiss module */
{
  struct command* keep_beam = current_beam;
  struct name_list* nl = current_twiss->par_names;
  struct command_parameter_list* pl = current_twiss->par;
  struct int_array* tarr;
  int ptc_twiss_summary = 0 ; /* set to 1 when a summary-table is filled-in */
  struct int_array* summary_tarr; /* to pass summary-table name to Fortran */
  struct node *nodes[2], *use_range[2];
  double ptc_deltap;
  char *filename = NULL, *table_name;
  char *summary_filename = NULL, *summary_table_name; /* for summary table */
  int j,l ,pos, w_file,beta_def;
  int w_file_summary; /* toggle to write the summary table into a file */
  /*
    start command decoding
  */
  use_range[0] = current_sequ->range_start;
  use_range[1] = current_sequ->range_end;

  if ((pos = name_list_pos("range", nl)) > -1 && nl->inform[pos])
  {
    if (get_sub_range(pl->parameters[pos]->string, current_sequ, nodes))
    {
      current_sequ->range_start = nodes[0];
      current_sequ->range_end = nodes[1];
    }
    else warning("illegal range ignored:", pl->parameters[pos]->string);
  }
  for (j = 0; j < current_sequ->n_nodes; j++)
  {
    if (current_sequ->all_nodes[j] == current_sequ->range_start) break;
  }

  if (attach_beam(current_sequ) == 0)
    fatal_error("PTC_TWISS - sequence without beam:", current_sequ->name);

  pos = name_list_pos("table", nl);
  if(nl->inform[pos]) /* table name specified - overrides save */
  {
    table_name = pl->parameters[pos]->string;
    if (table_name == NULL)
    {
      table_name = pl->parameters[pos]->call_def->string;
    }
  }
  else
  {
    /*strcpy(table_name,"ptc_twiss");*/
    table_name = "ptc_twiss";
  }

  /* --- */
  /* do the same as above for the table holding summary data after one-turn */
  pos = name_list_pos("summary_table",nl);
  if (nl->inform[pos]){ /* summary-table's name specified */
    summary_table_name = pl->parameters[pos]->string;
    if (summary_table_name == NULL){
      summary_table_name = pl->parameters[pos]->call_def->string;
    }
  }
  else {
    summary_table_name = "ptc_twiss_summary";
  }
  /* --- */

  pos = name_list_pos("file", nl);
  if (nl->inform[pos])
  {
    if ((filename = pl->parameters[pos]->string) == NULL)
    {
      if (pl->parameters[pos]->call_def != NULL)
        filename = pl->parameters[pos]->call_def->string;
    }
    if (filename == NULL) filename = permbuff("dummy");
    w_file = 1;
  }
  else w_file = 0;

  /* --- */
  /* do the same as above for the file to hold the summary-table */
  pos = name_list_pos("summary_file",nl);
  if (nl->inform[pos]){
    if ((summary_filename = pl->parameters[pos]->string) == NULL){
      if (pl->parameters[pos]->call_def != NULL)
        summary_filename = pl->parameters[pos]->call_def->string;
    }
    if (summary_filename == NULL) summary_filename = permbuff("dummy");
    w_file_summary = 1;
  }
  else w_file_summary = 0;

  /*
    end of command decoding
  */
  if ((beta_def = twiss_input(current_twiss)) < 0)
  {
    if (beta_def == -1) warning("unknown beta0,", "Twiss ignored");
    else if (beta_def == -2)
      warning("betx or bety missing,", "Twiss ignored");
    /*      set_variable("twiss_tol", &tol_keep); */
    return;
  }

  set_option("twiss_inval", &beta_def);
  adjust_beam();
  probe_beam = clone_command(current_beam);
  ptc_deltap = get_value(current_command->name,"deltap");
  adjust_probe(ptc_deltap); /* sets correct gamma, beta, etc. */
  adjust_rfc(); /* sets freq in rf-cavities from probe */
  l = strlen(table_name);
  tarr = new_int_array(l+1);
  conv_char(table_name, tarr);

  twiss_table = make_table(table_name, "twiss", twiss_table_cols,
                           twiss_table_types, current_sequ->n_nodes);

  twiss_table->dynamic = 1;
  add_to_table_list(twiss_table, table_register);
  current_sequ->tw_table = twiss_table;
  twiss_table->org_sequ = current_sequ;
  twiss_table->curr= 0;
  current_node = current_sequ->ex_start;
  /* w_ptc_twiss_(tarr->i); */

  /* --- */
  /* create additional table to hold summary data after one-turn */
  /* such as momentum compaction factor, tune and chromaticities */
  l = strlen(summary_table_name); /* reuse of l */
  summary_tarr = new_int_array(l+1);
  conv_char(summary_table_name, summary_tarr);

  ptc_twiss_summary_table = make_table(summary_table_name, "twiss summary",
                                       ptc_twiss_summary_table_cols, ptc_twiss_summary_table_types,
                                       5); /* only one  row would be enough  */
  /*  18 january 2010 static table - ptc_twiss_summary_table->dynamic = 1; */
 /* actually static for time-being */
  add_to_table_list( ptc_twiss_summary_table, table_register );
  /* --- */

  w_ptc_twiss_(tarr->i,summary_tarr->i);

  /* upon completion of the Fortran ptc_twiss call ... */

  fill_twiss_header_ptc(twiss_table,ptc_deltap);
  if (w_file) out_table(table_name, twiss_table, filename);

  /* --- */
  if (w_file_summary) out_table(summary_table_name, ptc_twiss_summary_table,
                                summary_filename);
  /* --- */

  /* cleanup */
  current_beam = keep_beam;
  probe_beam = delete_command(probe_beam);
  current_sequ->range_start = use_range[0];
  current_sequ->range_end = use_range[1];
  delete_int_array(tarr);

  /* --- */
  delete_int_array(summary_tarr);
  /* --- */

  /* For the time-being, summary data are only available in case of a closed machine */
  ptc_twiss_summary = get_option("ptc_twiss_summary");
  if (ptc_twiss_summary) {
    print_table(ptc_twiss_summary_table);
  }
}

void
pro_ptc_create_layout(void)
  /* controls ptc_create_layout module */
{
//  int pos; not used
  struct command* keep_beam = current_beam;
  if (attach_beam(current_sequ) == 0)
    fatal_error("ptc_create_layout - sequence without beam:", current_sequ->name);
  adjust_beam();
  probe_beam = clone_command(current_beam);

  if (name_list_pos("errors_dipole", table_register->names) <= -1) // (pos = not used
  {
    errors_dipole = make_table("errors_dipole", "efield", efield_table_cols,
                               efield_table_types, 10000);
    add_to_table_list(errors_dipole, table_register);
  }
  else
  {
    reset_count("errors_dipole");
  }

  if (name_list_pos("errors_field", table_register->names) <= -1) // (pos = not used
  {
    errors_field = make_table("errors_field", "efield", efield_table_cols,
                              efield_table_types, 10000);
    add_to_table_list(errors_field, table_register);
  }
  else
  {
    reset_count("errors_field");
  }

  if (name_list_pos("errors_total", table_register->names) <= -1) // (pos = not used
  {
    errors_total = make_table("errors_total", "efield", efield_table_cols,
                              efield_table_types, 10000);
    add_to_table_list(errors_total, table_register);
  }
  else
  {
    reset_count("errors_total");
  }

  w_ptc_create_layout_();
  /* cleanup */
  current_beam = keep_beam;
  probe_beam = delete_command(probe_beam);
}

void
pro_ptc_read_errors(void)
  /* controls ptc_read_errors module */
{
  struct command* keep_beam = current_beam;
  if (attach_beam(current_sequ) == 0)
    fatal_error("ptc_read_errors - sequence without beam:", current_sequ->name);
  adjust_beam();
  probe_beam = clone_command(current_beam);

  w_ptc_read_errors_();
  /* cleanup */
  current_beam = keep_beam;
  probe_beam = delete_command(probe_beam);
}

void
pro_ptc_refresh_k(void)
  /* controls ptc_refresh_k module */
{
  struct command* keep_beam = current_beam;
  if (attach_beam(current_sequ) == 0)
    fatal_error("ptc_refresh_k - sequence without beam:", current_sequ->name);
  adjust_beam();
  probe_beam = clone_command(current_beam);

  w_ptc_refresh_k_();
  /* cleanup */
  current_beam = keep_beam;
  probe_beam = delete_command(probe_beam);
}

void
select_ptc_normal(struct in_cmd* cmd)
  /* sets up all columns of the table normal_results except the last one (value) */
{
  struct name_list* nl;
  struct command_parameter_list* pl;
  struct table* t;
  int pos;
  int i, j, jj, curr;
  int skew, mynorder,myn1,myn2,mynres,indexa[4][1000];
  char* order_list;
  int min_req_order;
  double order[4],n1,n2,n3,n4;

  nl = this_cmd->clone->par_names;
  pl = this_cmd->clone->par;
  if (log_val("clear", cmd->clone))
  {
    min_order = 1;
    min_req_order = 1;
    mynres = 0;
    skew = 0;
    reset_count("normal_results");
/*    if ((pos = name_list_pos("normal_results", table_register->names)) > -1) delete_table(table_register->tables[pos]);*/
    return;
  }
  if ((pos = name_list_pos("normal_results", table_register->names)) <= -1)
  {
    /* initialise table */
    normal_results = make_table("normal_results", "normal_res", normal_res_cols,
                                normal_res_types, MAX_ROWS);
    normal_results->dynamic = 1;
    add_to_table_list(normal_results, table_register);
    reset_count("normal_results");
    pos = name_list_pos("normal_results", table_register->names);
    min_order = 1;
    min_req_order = 1;
  }
  t = table_register->tables[pos];

  /* initialise order array */
  order[0] = zero;
  order[1] = zero;
  order[2] = zero;
  order[3] = zero;
  if (t->curr == t->max) grow_table(t);

  for (j = 0; j < PTC_NAMES_L; j++)
  {
    /* Treat each ptc variable */

    pos = name_list_pos(names[j], nl);
    if (pos > -1 && nl->inform[pos])
    {
      curr = pl->parameters[pos]->m_string->curr;
      if (curr > 4)
        printf("Too many values for the attribute %s. Only the first four are retained.\n",names[j]);
      for (i = 0; i < curr; i++)
      {
        order_list = pl->parameters[pos]->m_string->p[i];
        order[i] = atoi(order_list);
      }

      if (j == 10 || j == 11)
      {
        min_req_order = order[0]+order[1]+order[2];
        mynres = 0;
        skew = 0;
        mynorder = (int)order[0];
        if (mynorder < 0) skew = 1;
        mynorder = abs(mynorder);
        myn1 = (int)order[1];
        myn2 = (int)order[2];
        min_req_order = mynorder;
        res_index_(&skew, &mynorder, &myn1, &myn2, &indexa[0][0], &mynres);
	if (mynres > 0)
        {
          if (j == 10)
          {
            for (jj = 0; jj < mynres; jj++)
            {
              n1 = (double)indexa[0][jj];
              n2 = (double)indexa[1][jj];
              n3 = (double)indexa[2][jj];
              n4 = (double)indexa[3][jj];
              string_to_table_curr("normal_results", "name", "hamc");
              double_to_table_curr("normal_results", "order1", &n1);
              double_to_table_curr("normal_results", "order2", &n2);
              double_to_table_curr("normal_results", "order3", &n3);
              double_to_table_curr("normal_results", "order4", &n4);
              augment_count("normal_results");
              string_to_table_curr("normal_results", "name", "hams");
              double_to_table_curr("normal_results", "order1", &n1);
              double_to_table_curr("normal_results", "order2", &n2);
              double_to_table_curr("normal_results", "order3", &n3);
              double_to_table_curr("normal_results", "order4", &n4);
              augment_count("normal_results");
              string_to_table_curr("normal_results", "name", "hama");
              double_to_table_curr("normal_results", "order1", &n1);
              double_to_table_curr("normal_results", "order2", &n2);
              double_to_table_curr("normal_results", "order3", &n3);
              double_to_table_curr("normal_results", "order4", &n4);
              augment_count("normal_results");
            }
            string_to_table_curr("normal_results", "name", "haml");
            double_to_table_curr("normal_results", "order1", &order[0]);
            double_to_table_curr("normal_results", "order2", &order[1]);
            double_to_table_curr("normal_results", "order3", &order[2]);
            double_to_table_curr("normal_results", "order4", &order[3]);
            n1 = (double)mynres;
            double_to_table_curr("normal_results", "value", &n1);
            augment_count("normal_results");
          }
          if (j == 11)
          {
            for (jj = 0; jj < mynres; jj++)
            {
              n1 = (double)indexa[0][jj];
              n2 = (double)indexa[1][jj];
              n3 = (double)indexa[2][jj];
              n4 = (double)indexa[3][jj];
              string_to_table_curr("normal_results", "name", "gnfc");
              double_to_table_curr("normal_results", "order1", &n1);
              double_to_table_curr("normal_results", "order2", &n2);
              double_to_table_curr("normal_results", "order3", &n3);
              double_to_table_curr("normal_results", "order4", &n4);
              augment_count("normal_results");
              string_to_table_curr("normal_results", "name", "gnfs");
              double_to_table_curr("normal_results", "order1", &n1);
              double_to_table_curr("normal_results", "order2", &n2);
              double_to_table_curr("normal_results", "order3", &n3);
              double_to_table_curr("normal_results", "order4", &n4);
              augment_count("normal_results");
              string_to_table_curr("normal_results", "name", "gnfa");
              double_to_table_curr("normal_results", "order1", &n1);
              double_to_table_curr("normal_results", "order2", &n2);
              double_to_table_curr("normal_results", "order3", &n3);
              double_to_table_curr("normal_results", "order4", &n4);
              augment_count("normal_results");
            }
            string_to_table_curr("normal_results", "name", "gnfu");
            double_to_table_curr("normal_results", "order1", &order[0]);
            double_to_table_curr("normal_results", "order2", &order[1]);
            double_to_table_curr("normal_results", "order3", &order[2]);
            double_to_table_curr("normal_results", "order4", &order[3]);
            n1 = (double)mynres;
            double_to_table_curr("normal_results", "value", &n1);
            augment_count("normal_results");
          }
        }
      }
      else
      {
        if(((strcmp(names[j], "dq1") == zero)|| (strcmp(names[j], "dq2") == zero)) && order[0] == zero) order[0]=one;
        string_to_table_curr("normal_results", "name", names[j]);
        double_to_table_curr("normal_results", "order1", &order[0]);
        double_to_table_curr("normal_results", "order2", &order[1]);
        double_to_table_curr("normal_results", "order3", &order[2]);
        double_to_table_curr("normal_results", "order4", &order[3]);
        augment_count("normal_results");
        if(j == 12)
        {
          min_req_order = 1;
        }
        else
        {
          min_req_order = order[0]+order[1]+order[2];
          if (j >= 8) min_req_order += order[0]+order[1];
          if (j >= 6) min_req_order += 1;
        }
      }
      if (min_order < min_req_order) min_order = min_req_order;
    }
  }
 if (debuglevel > 2) 
  {
    printf("The minimum required order is %d \n--------------------------------\n",min_order);
  }  
}

void
pro_ptc_trackline(struct in_cmd* cmd)
{
  /*Does PTC tracking taking to the account acceleration */
  /*it is basically wrapper to subroutine ptc_trackline() in madx_ptc_trackline.f90*/

  int pos, one;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int parexist = -1;
  double value = 0;
  int ivalue = 0;

  struct command* keep_beam = current_beam;

  if (attach_beam(current_sequ) == 0)
    fatal_error("PTC_TRACKLINE - sequence without beam:", current_sequ->name);
  

  pos = name_list_pos("file", nl);

  if (nl->inform[pos])
  {
    set_option("track_dump", &one);
  }

  if ((track_filename = pl->parameters[pos]->string) == NULL)
  {
    if (pl->parameters[pos]->call_def != NULL)
    {
      track_filename = pl->parameters[pos]->call_def->string;
    }
    else
    {
      track_filename = permbuff("dummy");
    }
  }
  track_filename = permbuff(track_filename);
  track_fileext = NULL;
  pos = name_list_pos("extension", nl);

  if ((track_fileext = pl->parameters[pos]->string) == NULL)
  {
    if (pl->parameters[pos]->call_def != NULL)
    {
      track_fileext = pl->parameters[pos]->call_def->string;
    }
    if (track_fileext == NULL)
    {
      track_fileext = permbuff("\0");
    }
  }

  track_fileext = permbuff(track_fileext);

  if (command_par_value("everystep",cmd->clone) != 0)
  {
    printf("Enforcing onetable=true, current is %f\n", command_par_value("onetable",cmd->clone));
    set_command_par_value("onetable", cmd->clone, 1.0);
    printf("Now is %f\n", command_par_value("onetable",cmd->clone));
  }

  
  parexist = command_par_value2("onetable",cmd->clone,&value);
  
  if (parexist)
   { 
     ivalue = (int)value;
     set_option("onetable", &ivalue);
   }
   

  adjust_beam();
  probe_beam = clone_command(current_beam);
  adjust_rfc(); /* sets freq in rf-cavities from probe */


  track_tables_create(cmd);
  
  
  if (command_par_value("everystep",cmd->clone) != 0)
  {
    /*printf("Calling PTC track line every step\n");*/
    w_ptc_track_everystep_(&curr_obs_points);
  }
  else
  {
    /*printf("Calling STD PTC track line\n");*/
    w_ptc_trackline_(&curr_obs_points);
  }

  track_tables_dump();

  /* cleanup */
  current_beam = keep_beam;
  probe_beam = delete_command(probe_beam);
  
}

void
pro_ptc_enforce6d(struct in_cmd* cmd)
{
  /*Does PTC tracking taking to the account acceleration */
  /*it is basically wrapper to subroutine ptc_trackline() in madx_ptc_trackline.f90*/
  double switchvalue;
  struct name_list* nl;
  int flag;

  if (cmd == 0x0)
  {
    warning("pro_ptc_enforce6d:","Command is null!!!");
    return;
  }

  if (cmd->clone == 0x0)
  {
    error("pro_ptc_enforce6d","Command Definintion is null!!!");
    return;
  }

  nl = cmd->clone->par_names;

  /*DEBUG LEVEL SWITCH*/
  if ( name_list_pos("flag", nl) >=0 )
  {
    command_par_value2("flag", cmd->clone, &switchvalue);
    flag = (int)switchvalue;
    w_ptc_enforce6d_(&flag);
  }
  else
  {
    printf("flag is not present\n");
  }
}

void
pro_ptc_setswitch(struct in_cmd* cmd)
{
  /*Does PTC tracking taking to the account acceleration */
  /*it is basically wrapper to subroutine ptc_trackline() in madx_ptc_trackline.f90*/
  int i;
  double switchvalue;
  struct name_list* nl;

  if (cmd == 0x0)
  {
    warning("pro_ptc_setswitch:","Command is null!!!");
    return;
  }

  if (cmd->clone == 0x0)
  {
    printf("pro_ptc_setswitch: Command Definintion is null!!!\n");
    return;
  }

  if (match_is_on == kMatch_PTCknobs)
  {
    madx_mpk_setsetswitch(cmd);
    return;
  }

  nl = cmd->clone->par_names;

  /*DEBUG LEVEL SWITCH*/
  if ( name_list_pos("debuglevel", nl) >=0 )
  {
    command_par_value2("debuglevel", cmd->clone, &switchvalue);
    debuglevel = (int)switchvalue;
    w_ptc_setdebuglevel_(&debuglevel);
  }
  else
  {
    printf("debuglevel is not present\n");
  }


  /*ACCELERATION SWITCH*/
  if ( name_list_pos("maxacceleration", nl) >=0 )
  {
    command_par_value2("maxacceleration", cmd->clone, &switchvalue);
    if (debuglevel > 0) printf("maxaccel is found and its value is %f\n", switchvalue);
    i = (int)switchvalue;
    w_ptc_setaccel_method_(&i);
  }
  else
  {
    if (debuglevel > 0) printf("maxaccel is not present\n");
  }


  /*EXACT SWITCH*/
  if ( name_list_pos("exact_mis", nl) >=0 )
  {
    command_par_value2("exact_mis", cmd->clone, &switchvalue);
    if (debuglevel > 0) printf("exact_mis is found and its value is %f\n", switchvalue);
    i = (int)switchvalue;
    w_ptc_setexactmis_(&i);
  }
  else
  {
    if (debuglevel > 0)  printf("exact_mis is not present\n");
  }


  /*radiation SWITCH*/
  if ( name_list_pos("radiation", nl) >=0 )
  {
    command_par_value2("radiation", cmd->clone, &switchvalue);
    if (debuglevel > 0) printf("radiation is found and its value is %f\n", switchvalue);
    i = (int)switchvalue;
    w_ptc_setradiation_(&i);
  }
  else
  {
    if (debuglevel > 0) printf("radiation is not present\n");
  }

  /*fringe SWITCH*/
  if ( name_list_pos("fringe", nl) >=0 )
  {
    command_par_value2("fringe", cmd->clone, &switchvalue);
    if (debuglevel > 0) printf("fringe is found and its value is %f\n", switchvalue);
    i = (int)switchvalue;
    w_ptc_setfringe_(&i);
  }
  else
  {
    if (debuglevel > 0) printf("fringe is not present\n");
  }



  /*totalpath SWITCH*/
  if ( name_list_pos("totalpath", nl) >=0 )
  {
    command_par_value2("totalpath", cmd->clone, &switchvalue);
    if (debuglevel > 0) printf("totalpath is found and its value is %f\n", switchvalue);
    i = (int)switchvalue;
    w_ptc_settotalpath_(&i);
  }
  else
  {
    if (debuglevel > 0) printf("totalpath is not present\n");
  }


  /*TIME SWITCH*/
  if ( name_list_pos("time", nl) >=0 )
  {
    command_par_value2("time", cmd->clone, &switchvalue);
    if (debuglevel > 0) printf("time is found and its value is %f\n", switchvalue);
    i = (int)switchvalue;
    w_ptc_settime_(&i);
  }
  else
  {
    if (debuglevel > 0) printf("time is not present\n");
  }

  /*NOCAVITY SWITCH*/
  if ( name_list_pos("nocavity", nl) >=0 )
  {
    command_par_value2("nocavity", cmd->clone, &switchvalue);
    if (debuglevel > 0) printf("nocavity is found and its value is %f\n", switchvalue);
    i = (int)switchvalue;
    w_ptc_setnocavity_(&i);
  }
  else
  {
    if (debuglevel > 0) printf("nocavity is not present\n");
  }

  if (debuglevel > 0) printf("obs_points pro_ptc_setswitch Done\n");
}

void
pro_ptc_printparametric(struct in_cmd* cmd)
{
  struct command_parameter_list* c_parameters= cmd->clone->par;
  struct name_list*              c_parnames  = cmd->clone->par_names;
  int                            pos         = 0;

  char*                          filename    = 0x0;
  struct int_array*              filenameIA      = 0x0;
  static int                     zeroint = 0;
  int*                           filep = 0x0;

  pos   = name_list_pos("filename", c_parnames);
  if (pos < 0)
  {

    filep = &zeroint;
  }
  else
  {
    filename  = c_parameters->parameters[pos]->string;
    if ( filename == 0x0 )
    {
      filep = &zeroint;
    }
    else
    {
      filenameIA = new_int_array(1+strlen(filename));
      conv_char(filename,filenameIA);
      filep = filenameIA->i;
    }
  }

  pos   = name_list_pos("format", c_parnames);
  if (pos < 0)
  {
    printf("mad_ptc.c: pro_ptc_printparametric: format parameter does not exist.\n");
    return;
  }

  w_ptc_writeparresults_(filep);

  delete_int_array(filenameIA);
}

void
pro_ptc_printframes(struct in_cmd* cmd)
{
  struct command_parameter_list* c_parameters= cmd->clone->par;
  struct name_list*              c_parnames  = cmd->clone->par_names;
  int                            pos         = 0;

  char*                          filename    = 0x0;
  struct int_array*              filenameIA  = 0x0;
  char*                          format      = 0x0;
/* 
 * Piotr.Skowronski@cern.ch
 *  Routine that writes coordinates of magnets from PTC
 *  either as text or ROOT macro which executed gives 3d image of the layout
 *  Requires one parameter: filename 
 */
  pos   = name_list_pos("file", c_parnames);
  if (pos < 0)
  {
    printf("mad_ptc.c: pro_ptc_printframes: file parameter does not exist.\n");
    return;
  }

  filename  = c_parameters->parameters[pos]->string;
  if ( filename == 0x0 )
  {
    warning("mad_ptc.c: pro_ptc_printframes: no file name: ", "ignored");
    return;
  }


  pos   = name_list_pos("format", c_parnames);
  if (pos < 0)
  {
    printf("mad_ptc.c: pro_ptc_printframes: format parameter does not exist.\n");
    return;
  }

  format  = c_parameters->parameters[pos]->string;
  printf("mad_ptc.c: pro_ptc_printframes: format is %s.\n", format);

  filenameIA = new_int_array(1+strlen(filename));

  conv_char(filename,filenameIA);


  if (strcmp(format,"rootmacro") == 0)
  {
    w_ptc_printlayout_rootm_(filenameIA->i);
  }
  else
  {
    w_ptc_printframes_(filenameIA->i);
  }

  delete_int_array(filenameIA);
}

void
pro_ptc_export_xml(struct in_cmd* cmd)
{
  struct command_parameter_list* c_parameters= cmd->clone->par;
  struct name_list*              c_parnames  = cmd->clone->par_names;
  int                            pos         = 0;
  char*                          filename    = 0x0;
  struct int_array*              filenameIA      = 0x0;
/*  char*                          format    = 0x0; */

  pos   = name_list_pos("file", c_parnames);
  if (pos < 0)
  {
    /* should never enter here */
    printf("mad_ptc.c: pro_ptc_export_xml: file parameter does not exist.\n");
    return;
  }

  filename  = c_parameters->parameters[pos]->string;
  printf("will write to file %s\n",filename);

  if ( filename == 0x0 )
  {
    warning("mad_ptc.c: pro_ptc_export_xml: no file name: ", "ignored");
    return;
  }

  filenameIA = new_int_array(1+strlen(filename));

  conv_char(filename,filenameIA);

  w_ptc_export_xml_(filenameIA->i);

  delete_int_array(filenameIA);
}

void
pro_ptc_eplacement(struct in_cmd* cmd)
{/*
   Sets a parameter
 */
  struct command_parameter_list* c_parameters= cmd->clone->par;
  struct name_list*              c_parnames  = cmd->clone->par_names;
  int                            pos         = 0;
  int                            k         = 0;
  struct node*                   nodes[2]={0x0,0x0};
  struct node*                   anode=0x0;
  char*                          element;
  int                            refframe=0;/*0 global, 1 current position, 2 end face if the previous element*/


  pos   = name_list_pos("refframe", c_parnames);
  if (pos < 0)
  {
    printf("mad_ptc.c: pro_ptc_eplacement: refframe parameter does not exist.\n");
    return;
  }

  if ( c_parnames->inform[pos] != 0 )
  {
    /*if it is zero it is not specified*/

    if ( c_parameters->parameters[pos]->string == 0x0 )
    {
      warning("mad_ptc.c: pro_ptc_eplacement: string describing refframe is null: ", "using default");
      refframe = 0;
    }
    else
    {
      /*printf("refframe is %s.\n", c_parameters->parameters[pos]->string );*/

      if ( strcmp(c_parameters->parameters[pos]->string,"current")  == 0 )
      {
        refframe = 1;
      }

      if ( strcmp(c_parameters->parameters[pos]->string,"previouselement") == 0 )
      {
        refframe = 2;
      }
    }
  }


  pos   = name_list_pos("range", c_parnames);
  if (pos < 0)
  {
    printf("mad_ptc.c: pro_ptc_eplacement: range parameter does not exist.\n");
    return;
  }

  if ( c_parnames->inform[pos] == 0 )
  {
    printf("mad_ptc.c: pro_ptc_eplacement: inform for range is 0.\n");
    return;
  }

  element  = c_parameters->parameters[pos]->string;
  if ( element == 0x0 )
  {
    warning("mad_ptc.c: pro_ptc_eplacement: no element name: ", "ignored");
    return;
  }


  k = get_ex_range(element, current_sequ, nodes);
  if ( k != 1)
  {
    if (k > 1)
    {
      warningnew("pro_ptc_eplacement","More then one element correstponds to the range <<%s>>.",element);
      seterrorflag(1,"pro_ptc_eplacement","More then one element correstponds to the range");
      return;
    }
    else
    {
      warningnew("pro_ptc_eplacement","Element <<%s>> not found",element);
      seterrorflag(1,"pro_ptc_eplacement","Element not found");
      return;
    }
  }



  pos = 0;
  anode=current_sequ->range_start;
  while(anode)
  {
    /*printf("%d: Comparing %#x %s with  %#x %s \n",pos, nodes[0], nodes[0]->name, anode, anode->name);*/
    if ( nodes[0]  == anode  )
    {
      /*printf("Element is at pos %d !\n",pos);*/
      break;
    }

    if (anode == current_sequ->range_end)
    {
      warningnew("pro_ptc_eplacement","Reached the end of sequence - Element <<%s>> not found",element);
      return;
    }

    anode = anode->next;
    pos++; /*before if because fortran numerates from 1*/

  }
  w_ptc_eplacement_(&pos,&refframe);
}

void
pro_ptc_varyknob(struct in_cmd* cmd)
{/*
   Sets a variable based on parameter
 */

  if (match_is_on != kMatch_PTCknobs)
  {
    warningnew("pro_ptc_varyknob","Match with ptcknobs is not active, command ignored");
    return;
  }

  madx_mpk_addvariable(cmd);

}

void
pro_ptc_knob(struct in_cmd* cmd)
{/*
   Sets a parameter
 */
  struct command_parameter_list* c_parameters= cmd->clone->par;
  struct name_list*              c_parnames  = cmd->clone->par_names;
  int                            pos         = 0;

  char*                          element     = 0x0;
  struct int_array*              elementIA   = 0x0;
  char*                          initialp    = 0x0;
  struct int_array*              initialpIA  = 0x0;
  char*                                   p  = 0x0;

  pos   = name_list_pos("element", c_parnames);
  if (pos < 0)
  {
    printf("mad_ptc.c: pro_ptc_knob: element parameter does not exist.\n");
    return;
  }
  element  = c_parameters->parameters[pos]->string;

  pos   = name_list_pos("initial", c_parnames);
  if (pos < 0)
  {
    printf("mad_ptc.c: pro_ptc_knob: initial parameter does not exist.\n");
    return;
  }

  initialp  = c_parameters->parameters[pos]->string;

  if ( (element == 0x0) && (initialp == 0x0) )
  {
    warning("mad_ptc.c: pro_ptc_knob: no element name neither initial pareter specified: ",
            "command ignored");
    return;
  }

  if (initialp)
  {
    mycpy(c_dum->c, initialp);

    stolower(c_dum->c);

    initialpIA = new_int_array(1+strlen(c_dum->c));

    conv_char(c_dum->c,initialpIA);

    w_ptc_addknob_i_(initialpIA->i);

    delete_int_array(initialpIA);

  }
  else
  {
    mycpy(c_dum->c, element);

    if (command_par_value("exactmatch",cmd->clone) != 0)
    {
      p = strstr(c_dum->c,"[");
      if (p)
      {
        *p = ':';
        p = strstr(c_dum->c,"]");
        if (p == 0x0)
        {
          warningnew("mad_ptc.c: pro_ptc_knob:","element %s is bady defned. Commnad ignored.",element);
          return;
        }
        *p=0;
      }
      else
      { /*assume it is the first element*/
        p = &(c_dum->c[strlen(c_dum->c)]);
        p[0]=':';
        p[1]='1';
        p[2]= 0;
      }

    }
    stoupper(c_dum->c);

    elementIA = new_int_array(1+strlen(c_dum->c));

    conv_char(c_dum->c,elementIA);

    w_ptc_addknob_(elementIA->i);

    delete_int_array(elementIA);
  }
}

void
pro_ptc_setknobvalue(struct in_cmd* cmd)
{/*
   Sets a parameter value
 */
  struct command_parameter_list* c_parameters= cmd->clone->par;
  struct name_list*              c_parnames  = cmd->clone->par_names;
  int                            pos         = 0;

  char*                          element    = 0x0;
  struct int_array*              elementIA      = 0x0;

  pos   = name_list_pos("element", c_parnames);
  if (pos < 0)
  {
    printf("mad_ptc.c: pro_ptc_knob: element parameter does not exist.\n");
    return;
  }

  element  = c_parameters->parameters[pos]->string;
  if ( element == 0x0 )
  {
    warning("mad_ptc.c: pro_ptc_knob: no element name: ", "ignored");
    return;
  }
  mycpy(c_dum->c, element);

  stoupper(c_dum->c);

  elementIA = new_int_array(1+strlen(c_dum->c));

  conv_char(c_dum->c,elementIA);

  w_ptc_setknobvalue_(elementIA->i);

  delete_int_array(elementIA);


}

void
pro_ptc_setfieldcomp(struct in_cmd* cmd)
{/*
   Sets a parameter value
 */
  struct command_parameter_list* c_parameters= cmd->clone->par;
  struct name_list*              c_parnames  = cmd->clone->par_names;
  int                            pos         = 0;
  int                            k         = 0;
  struct node*                   nodes[2]={0x0,0x0};
  struct node*                   anode=0x0;
  char*                          element;



  pos   = name_list_pos("element", c_parnames);
  if (pos < 0)
  {
    printf("mad_ptc.c: pro_ptc_setfieldcomp: range parameter does not exist.\n");
    return;
  }

  if ( c_parnames->inform[pos] == 0 )
  {
    printf("mad_ptc.c: pro_ptc_setfieldcomp: inform for range is 0.\n");
    return;
  }

  element  = c_parameters->parameters[pos]->string;
  if ( element == 0x0 )
  {
    warning("mad_ptc.c: pro_ptc_setfieldcomp: no element name: ", "ignored");
    return;
  }


  k = get_range(element, current_sequ, nodes);
  if ( k != 1)
  {
    if (k > 1)
    {
      warningnew("pro_ptc_setfieldcomp","More then one element correstponds to the range <<%s>>.",element);
      seterrorflag(1,"pro_ptc_setfieldcomp","More then one element correstponds to the range");
      return;
    }
    else
    {
      warningnew("pro_ptc_setfieldcomp","Element <<%s>> not found",element);
      seterrorflag(1,"pro_ptc_setfieldcomp","Element not found");
      return;
    }
  }



  pos = 0;
  anode=current_sequ->range_start;
  while( anode != 0x0  )
  {
    if ( nodes[0]  == current_sequ->nodes->nodes[pos]  )
    {
      /* printf("Element is at pos %d !\n",pos);*/
      break;
    }

    if (current_sequ->nodes->nodes[pos] == current_sequ->range_end)
    {
      warningnew("pro_ptc_setfieldcomp","Reached the end of sequence - Element <<%s>> not found",element);
      return;
    }
    pos++; /*before if because fortran numerates from 1*/
  }

  w_ptc_setfieldcomp_(&pos);

}

void
pro_ptc_select(struct in_cmd* cmd)
{/*
   processes ptc_select command
   it directs ptc_twiss to store given QUANTITY in a named TABLE's COLUMN
   Then, it these values are accessible for other MAD-X modules for calculations.
   The most important one is the matching module.
 */

  char*                          monomial    = 0x0;

  int                            element     = 0;
  int*                           tablep      = 0;
  int*                           columnp     = 0;
  static int                     zeroint     = 0;/*if there is no column name or table name these are passed as null strings */
  struct int_array*              tabnameIA   = 0x0;/*string passing to fortran is tricky*/
  struct int_array*              colnameIA   = 0x0;/*and is done via integer arrays*/
  struct int_array*              monoIA      = 0x0;

  monomial = command_par_string("monomial",cmd->clone);

  if (monomial == 0x0)
  {
    warning("mad_ptc.c: pro_ptc_select: monomial is NULL ", "ignored");
    return;
  }

  monoIA = new_int_array(1+strlen(monomial));
  conv_char(monomial,monoIA);

  element = command_par_value("polynomial",cmd->clone);

  pro_ptc_select_checkpushtable(cmd,&tabnameIA,&colnameIA);

  if ( tabnameIA )
  {
    tablep = tabnameIA->i;
  }
  else
  {
    tablep = &zeroint;
  }

  if ( colnameIA )
  {
    columnp = colnameIA->i;
  }
  else
  {
    columnp = &zeroint;
  }

  w_ptc_addpush_(tablep,columnp,&element,monoIA->i);

  delete_int_array(tabnameIA);
  delete_int_array(colnameIA);
  delete_int_array(monoIA);

}

int
pro_ptc_moments(struct in_cmd* cmd)
{
  int no = command_par_value("no", cmd->clone);

  w_ptc_moments_(&no);

  return 1;
}

int
pro_ptc_select_moment(struct in_cmd* cmd)
{
  /*adds a moment or more moments to the list
    if parametric switch is present they will be stored also as Taylor Series*/

  int pos, tablepos;
  int i, j;
  int mdefi[6];
  char* mdefin, *pchar;
  char  tablename[48];
  char  colname[9];
  int   clen = 0;
  struct int_array*              tabIA      = 0x0;
  struct int_array*              mdefIA      = 0x0;
  struct command_parameter_list* c_parameters= cmd->clone->par;
  struct name_list*              c_parnames  = cmd->clone->par_names;
  int                            parametric = 0;
  int                            int_arr[100];


  tablepos = name_list_pos("table", c_parnames);
  if ( tablepos < 0)
  {
    printf("Weired: table parameter is not defined\n");
    return 1;
  }

  pchar  = c_parameters->parameters[tablepos]->string;
  if ( pchar == 0x0 )
  {
    strcpy(tablename,"moments");
  }
  else if ( pchar[0] == 0 )
  {
    strcpy(tablename,"moments");
  }
  else
  {
    strcpy(tablename,pchar);
  }

  tabIA = new_int_array(1+strlen(tablename));
  conv_char(tablename,tabIA);

  pos = name_list_pos("moment_s", c_parnames);
  if ( pos < 0)
  {
    printf("Weired: moments parameter is not defined\n");
    return 1;
  }

  if (c_parnames->inform[pos])
  {
    for (j = 0; j < c_parameters->parameters[pos]->m_string->curr; j++)
    {
      strcpy(colname,"mu000000");

      mdefin = c_parameters->parameters[pos]->m_string->p[j];

      /*printf("String no %d is %s\n", j, mdefin);*/

      clen = strlen(mdefin);

      /*the loop below decodes string monomial to integer monomial */
      for (i = 0; i< 6; i++)
      {
        if (clen > i)
        {
          mdefi[i] = mdefin[i] - '0';
          colname[2+i] = mdefin[i];
        }
        else
        {
          mdefi[i] = 0;
        }
      }


      mdefIA = new_int_array(9);/*size + "mu" + 6charnumbers*/
      conv_char(colname,mdefIA);

      w_ptc_addmoment_(&(mdefi[0]),&(mdefi[1]),&(mdefi[2]),&(mdefi[3]),&(mdefi[4]),&(mdefi[5]),
                       tabIA->i, mdefIA->i, &parametric);

      delete_int_array(mdefIA);

    }
  }


  pos = name_list_pos("moment", c_parnames);
  /*i is dummy... we know there is no strings or doubles */
  comm_para_("moment", &pos, &i, &i, int_arr, 0x0, 0x0, 0x0);

  if (pos > 6) pos = 6;
/*if there is something and it is not only one zero */
  if ( (pos >= 0)       && !((pos == 1) && (int_arr[0] == 0))  )
  {

    for (i=0;i<pos;i++)
    {
      if (int_arr[i] < 0) break;

      mdefi[i] = int_arr[i];
    }

    for (i=pos; i<6;i++)
    {
      mdefi[i] = 0;
    }

    tablepos = name_list_pos("column", c_parnames);
    if ( tablepos < 0)
    {
      printf("Weired: column parameter is not defined\n");
      return 1;
    }

    tablename[0] = 0;
    pchar  = c_parameters->parameters[tablepos]->string;
    if ( pchar != 0x0 )
    {
      if ( pchar[0] != 0 )
      {
        strcpy(tablename,pchar);
      }

    }

    if (tablename[0] == 0)
    {
      sprintf(tablename,"mu_%d_%d_%d_%d_%d_%d",
              mdefi[0],mdefi[1],mdefi[2],mdefi[3],mdefi[4],mdefi[5]);
      printf("pro_ptc_select_moment: Column name not provied, generated one is %s",tablename);
    }

    mdefIA = new_int_array(1+strlen(tablename));

    conv_char(tablename,mdefIA);
    w_ptc_addmoment_(&(mdefi[0]),&(mdefi[1]),&(mdefi[2]),&(mdefi[3]),&(mdefi[4]),&(mdefi[5]),
                     tabIA->i, mdefIA->i, &parametric);

    delete_int_array(mdefIA);
  }



  delete_int_array(tabIA);

  return 0;
}

int
makemomentstables(void)
{
  static const int maxtables = 100;
  char*             tables[100];           /*tables[maxtables];*/
  struct name_list* cols[100];  /*cols[maxtables][maxcols];*/
  struct table*     t;
  char              tabname[20];
  char              colname[17];
  int               nmom;
  int i,j; // ,k; not used


  memset(tables,0x0,maxtables*sizeof(char*));


  nmom = w_ptc_getnmoments_();
  for (i = 1; i <= nmom; i++)
  {
    w_ptc_getmomentstabcol_(&i, tabname, colname);
    /*printf(" mom %d: %s %s\n",i, tabname, colname);*/

    for(j=0; tables[j] != 0x0 ;j++)
    {
      if ((strcmp(tables[j],tabname) == 0)) break;
    }
    /*printf(" index of this table is %d \n",j);*/

    if (tables[j] == 0x0)
    {
      tables[j] = (char*)mycalloc("makemomentstables",strlen(tabname) + 1, sizeof(char));
      strcpy(tables[j],tabname);
      cols[j] = new_name_list("columns", 15);
      add_to_name_list(permbuff("name"),3,cols[j]);
      add_to_name_list(permbuff("s"),2,cols[j]);
    }

    add_to_name_list(permbuff(colname),2,cols[j]); // k = not used

  }

  if (moments_tables)
  {
    myfree("",moments_tables->tables);
    delete_name_list(moments_tables->names);
    myfree("",moments_tables);
  }

  moments_tables = new_table_list(10);

  for(j=0; tables[j] != 0x0 ;j++)
  {
    /*printf("Making table %s\n",tables[j]);*/

    t = new_table(tables[j], "usermoments", current_sequ->n_nodes, cols[j]);
    t->org_cols = cols[j]->curr;
    /*print_table(t);*/
    add_to_table_list(t, table_register);
    add_to_table_list(t, moments_tables);
  }


/*  make_table("moments", "twiss", twiss_table_cols,
    twiss_table_types, current_sequ->n_nodes); */
  return 0;
}

void
augmentcountmomtabs(double* s)
{
  int i;
  struct table* t;

  if (moments_tables == 0x0)
  {
    warning("augmentcountmomtabs","moments_tables is NULL\n");
    return;
  }

  for ( i = 0; i <  moments_tables->curr; i++)
  {
    t = moments_tables->tables[i];
    t->s_cols[0][t->curr] = tmpbuff(current_node->name);
    t->d_cols[1][t->curr] = *s;
    if (t->num_cols > t->org_cols)  add_vars_to_table(t);
    if (++t->curr == t->max) grow_table(t);
  }
}

void
pro_ptc_script(struct in_cmd* cmd)
{/*
   processes ptc_script command
   it directs ptc_twiss to store given QUANTITY in a named TABLE's COLUMN
   Then, it these values are accessible for other MAD-X modules for calculations.
   The most important one is the matching module.
 */

  struct command_parameter_list* c_parameters= cmd->clone->par;
  struct name_list*              c_parnames  = cmd->clone->par_names;
  int                            pos         = 0;
  char*                          scriptname   = 0x0;
  struct int_array*              scriptnameIA = 0x0;/*string passing to fortran is tricky*/

  /*extracts table specified by the user*/
  pos   = name_list_pos("file", c_parnames);
  if (pos < 0)
  {
    printf("mad_ptc.c: pro_ptc_script: file parameter does not exist.\n");
    return;
  }

  scriptname  = c_parameters->parameters[pos]->string;
  if ( scriptname == 0x0 )
  {
    warning("mad_ptc.c: pro_ptc_script: no script name: ", "ignored");
    return;
  }

  scriptnameIA = new_int_array(1+strlen(scriptname));
  conv_char(scriptname,scriptnameIA);

  w_ptc_script_(scriptnameIA->i);/*calls the fortran*/

  delete_int_array(scriptnameIA);

}

void
pro_ptc_open_gino(struct in_cmd* cmd)
{
/*
  processes ptc_open_gino command
*/

  struct command_parameter_list* c_parameters= cmd->clone->par;
  struct name_list*              c_parnames  = cmd->clone->par_names;
  int                            pos         = 0;
  char*                          scriptname   = 0x0;
  struct int_array*              scriptnameIA = 0x0;/*string passing to fortran is tricky*/

  pos   = name_list_pos("command", c_parnames);
  scriptname  = c_parameters->parameters[pos]->string;
  if ( scriptname == 0x0 )
  {
    warning("mad_ptc.c: pro_ptc_open_gino: no script name: ", "ignored");
    return;
  }

  scriptnameIA = new_int_array(1+strlen(scriptname));
  conv_char(scriptname,scriptnameIA);

  w_ptc_open_gino_(scriptnameIA->i);/*calls the fortran*/

  delete_int_array(scriptnameIA);

}

void
pro_ptc_track(struct in_cmd* cmd)
{
  int k=0, pos, one = 1;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
/*  char rout_name[] = "ptc_track"; */
  int npart = stored_track_start->curr;
  struct table* t;
/*  int turns = command_par_value("turns", cmd->clone); */

  track_is_on = 1;
  puts("enter PTC_TRACK module");
  if (current_sequ == NULL || current_sequ->ex_start == NULL)
  {
    warning("sequence not active,", "TRACK ignored");
    return;
  }
  if (attach_beam(current_sequ) == 0)
    fatal_error("TRACK - sequence without beam:", current_sequ->name);
  if ((k = get_value(current_command->name,"onepass")) != 0)
    fprintf(prt_file, "one pass is on\n");
  /*
    if ((k = get_value(current_command->name,"damp")) != 0)
    fprintf(prt_file, "damp is on\n");
    set_option("damp", &k);
    if ((k = get_value(current_command->name,"quantum")) != 0)
    fprintf(prt_file, "quantum is on\n");
    set_option("quantum", &k);
  */
  set_option("onepass", &k);
  if ((k = get_value(current_command->name,"aperture")) != 0)
    fprintf(prt_file, "aperture tracking is on\n");
  set_option("aperture", &k);
  k = get_value(current_command->name,"dump");
  set_option("track_dump", &k);
  k = get_value(current_command->name,"onetable");
  set_option("onetable", &k);
  track_deltap=get_value(current_command->name,"deltap");
  set_variable("track_deltap", &track_deltap);
  if(track_deltap != 0) fprintf(prt_file, v_format("track_deltap: %F\n"),
                                track_deltap);
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

  if (npart == 0)
  {
    warning("track_run: no START command seen yet", "ignored");
    return;
  }
  track_tables_create(cmd);
  printf("obs_points ptc_track: %d \n",curr_obs_points);
  w_ptc_track_(&curr_obs_points);
  t = table_register->tables[name_list_pos("tracksumm", table_register->names)];
  if (get_option("info"))  print_table(t);
  if (get_option("track_dump")) track_tables_dump();
  fprintf(prt_file, "\n*****  end of ptc_run  *****\n");
}

