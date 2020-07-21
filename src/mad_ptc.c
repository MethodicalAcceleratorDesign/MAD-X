#include "madx.h"

void  makerdtstwisstable(void);
void printpoly(int*, int );

static void
fill_twiss_header_ptc(struct table* t, double ptc_deltap)
  /* puts beam parameters etc. at start of twiss table */
{
  int i, h_length = 100; /*39+3+1+1+6+4;  change when adding header lines ! - last 6 for the closed orbit */
  double dtmp;
  /*  struct table* s; */
  static const int tmplen=16;
  char tmp[tmplen];
  int row;

  static const char * const beampars[] = {"mass", "charge", "energy", "pc",
                                          "gamma", "kbunch", "bcurrent",
		  "sige","sigt",
		  "npart",
		  "ex","ey","et"};

  int Nbeampars = sizeof(beampars)/sizeof(beampars[0]);

  static const char * const ptcpars[] = {"length", "alpha_c", "alpha_c_p", "alpha_c_p2",
                                         "alpha_c_p3", "eta_c", "gamma_tr",
		 "q1","q2","dq1","dq2","qs",
		 "beta_x_min","beta_x_max","beta_y_min","beta_y_max",
		 "beta11min","beta11max", "beta12min","beta12max", "beta13min","beta13max",
		 "beta21min","beta21max", "beta22min","beta22max", "beta23min","beta23max",
		 "beta31min","beta31max", "beta32min","beta32max", "beta33min","beta33max",
		 "disp1min","disp1max",
		 "disp2min","disp2max",
		 "disp3min","disp3max",
		 "disp4min","disp4max",
		 "orbit_x", "orbit_px","orbit_y","orbit_py","orbit_pt","orbit_t",
		 "xcorms",  "pxcorms","ycorms","pycorms","tcorms","ptcorms",
		 "xcomin",  "xcomax",
		 "pxcomin", "pxcomax",
		 "ycomin",  "ycomax",
		 "pycomin", "pycomax",
		 "tcomin",  "tcomax",
		 "ptcomin", "ptcomax"};

  int Nptcpars = sizeof(ptcpars)/sizeof(ptcpars[0]);

  /*printf("There is %d beam parameters and %d PTC pars to add\n",Nbeampars,Nptcpars);*/

  if (t == NULL) return;
  /* ATTENTION: if you add header lines, augment h_length accordingly */
  if (t->header == NULL)  t->header = new_char_p_array(h_length);
  strncpy(tmp, t->org_sequ->name,tmplen);
  sprintf(c_dum->c, v_format("@ SEQUENCE         %%%02ds \"%s\""),
          strlen(tmp),stoupper(tmp));
  addto_char_p_array(t->header,c_dum);

  i = get_string("beam", "particle", tmp);
  sprintf(c_dum->c, v_format("@ PARTICLE         %%%02ds \"%s\""),
          i, stoupper(tmp));
  addto_char_p_array(t->header,c_dum);


  for (i=0; i<Nbeampars;i++)
   {
      strncpy(tmp,beampars[i],tmplen);  /*we have to copy for stoupper that can not change a constant string*/
      dtmp = get_value("beam", tmp);
      sprintf(c_dum->c, v_format("@ %-16.16s %%le  %F"),stoupper(tmp) ,dtmp);
      addto_char_p_array(t->header,c_dum);
   }

  sprintf(c_dum->c, v_format("@ DELTAP           %%le  %F"), ptc_deltap);
  addto_char_p_array(t->header,c_dum);

  if (get_option("ptc_twiss_summary") != zero){
    /* one-turn information gets computed iff ptc_twiss_summary set to 1 in madx_ptc_twiss.f90 */
    /* retreive all pieces of information from the ptc_twiss*/
    row = 1; /* this particular table has only one row filled-in */

    for (i=0; i<Nptcpars;i++)
     {
        strncpy(tmp,ptcpars[i],tmplen);
        double_from_table_row("ptc_twiss_summary",tmp,&row,&dtmp);
        sprintf(c_dum->c, v_format("@ %-16.16s %%le  %F"),stoupper(tmp) ,dtmp);
        addto_char_p_array(t->header,c_dum);
     }

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
  aTable = find_table(tablename);
  if (!aTable)
  {
    printf("mad_ptc.c: pro_ptc_select: table <<%s>> does not exist: Create table first\n",tablename);
    return 3;
  }


  /*checks if the specified column exists*/
  pos = name_list_pos(columnname,aTable->columns);
  if (pos < 0)
  {
    mad_error("mad_ptc.c: pro_ptc_select","Can not find column named <<%s>> in table <<%s>>.",
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
  struct table* t = find_table("normal_results");;
  return t ? t->curr : -1;
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

  if (current_sequ)
   {
     c_node = current_sequ->ex_start;
   }
  else
   {
     c_node = 0x0;
   }

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
  struct node* nodes[2];

  const char* place = command_par_string("place", cmd->clone);
  if (get_ex_range(place, current_sequ, nodes))
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
pro_ptc_normal(void)
  /* controls ptc_normal module */
{
  // need to set probe_beam to read beam parameters for beam-beam, for example
  struct command* keep_beam = current_beam;
  double ptc_deltap;  
  ptc_deltap = get_value(current_command->name,"deltap");
  
  adjust_beam();
  
  probe_beam = clone_command(current_beam);
  
  adjust_probe_fp(ptc_deltap); /* sets correct gamma, beta, etc. */
  
  w_ptc_normal_();

  /* cleanup */
  current_beam = keep_beam;
  probe_beam = delete_command(probe_beam);
  
  
}
void
pro_ptc_twiss(void)
  /* controls ptc_twiss module */
{
  struct command* keep_beam = current_beam;
  struct int_array* tarr;
  int ptc_twiss_summary = 0 ; /* set to 1 when a summary-table is filled-in */
  struct int_array* summary_tarr; /* to pass summary-table name to Fortran */
  struct node *nodes[2], *use_range[2];
  double ptc_deltap;
  const char *table_name = NULL, *summary_table_name = NULL;
  char *filename = NULL, *summary_filename = NULL; /* for summary table */
  int j, w_file,beta_def;
  int w_file_summary; /* toggle to write the summary table into a file */
  struct table* nonlin_table = 0;


  /*
    start command decoding
  */
  use_range[0] = current_sequ->range_start;
  use_range[1] = current_sequ->range_end;

  char* range = command_par_string_user("range", current_twiss);
  if (range)
  {
    if (get_sub_range(range, current_sequ, nodes))
    {
      current_sequ->range_start = nodes[0];
      current_sequ->range_end = nodes[1];
    }
    else warning("illegal range ignored:", range);
  }
  for (j = 0; j < current_sequ->n_nodes; j++)
  {
    if (current_sequ->all_nodes[j] == current_sequ->range_start) break;
  }

  if (attach_beam(current_sequ) == 0)
    fatal_error("PTC_TWISS - sequence without beam:", current_sequ->name);

  table_name = command_par_string_user("table", current_twiss);
  if(!table_name)
    table_name = "ptc_twiss";

  /* --- */
  /* do the same as above for the table holding summary data after one-turn */
  summary_table_name = command_par_string_user("summary_table", current_twiss);
  if (!summary_table_name)
    summary_table_name = "ptc_twiss_summary";

  /* --- */

  w_file = command_par_string_user2("file", current_twiss, &filename);
  if (w_file && !filename)                      // TG: should be impossible?
    filename = permbuff("dummy");

  /* --- */
  /* do the same as above for the file to hold the summary-table */
  w_file_summary = command_par_string_user2("summary_file", current_twiss, &summary_filename);
  if (w_file_summary && !summary_filename)      // TG: should be impossible?
    summary_filename = permbuff("dummy");

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
  ptc_deltap = get_value(current_command->name,"deltap");

  // LD 2016.04.19
  adjust_beam();
  probe_beam = clone_command(current_beam);

  adjust_probe_fp(ptc_deltap); /* sets correct gamma, beta, etc. */


  nonlin_table = make_table("nonlin", "nonlin", nonlin_table_cols,
                           nonlin_table_types, current_sequ->n_nodes);
  /*nonlin_table->dynamic = 1;*/
  add_to_table_list(nonlin_table, table_register);

  tarr = new_int_array(strlen(table_name)+1);
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


  if (command_par_value("trackrdts",current_twiss) != 0)
   {
     makerdtstwisstable();
   }
  /* --- */
  /* create additional table to hold summary data after one-turn */
  /* such as momentum compaction factor, tune and chromaticities */
  summary_tarr = new_int_array(strlen(summary_table_name)+1);
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

  // LD 2016.04.19
  adjust_beam();
  probe_beam = clone_command(current_beam);
  adjust_probe_fp(0);

  if (!table_exists("errors_dipole"))
  {
    errors_dipole = make_table("errors_dipole", "efield", efield_table_cols,
                               efield_table_types, 10000);
    add_to_table_list(errors_dipole, table_register);
  }
  else
  {
    reset_count("errors_dipole");
  }

  if (!table_exists("errors_field"))
  {
    errors_field = make_table("errors_field", "efield", efield_table_cols,
                              efield_table_types, 10000);
    add_to_table_list(errors_field, table_register);
  }
  else
  {
    reset_count("errors_field");
  }

  if (!table_exists("errors_total"))
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

  // LD 2016.04.19
  adjust_beam();
  probe_beam = clone_command(current_beam);
  adjust_probe_fp(0);

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

  // LD 2016.04.19
  adjust_beam();
  probe_beam = clone_command(current_beam);
  adjust_probe_fp(0);

  w_ptc_refresh_k_();
  /* cleanup */
  current_beam = keep_beam;
  probe_beam = delete_command(probe_beam);
}

void
select_ptc_normal(struct in_cmd* cmd)
  /* sets up all columns of the table normal_results except the last one (value) */
{
  struct table* t;
  struct command_parameter* cp;
  int i, j, jj, curr;
  int skew, mynorder,myn1,myn2,mynres,indexa[4][1000];
  char* order_list;
  int min_req_order;
  double order[4],n1,n2,n3,n4;

  if (log_val("clear", cmd->clone))
  {
    min_order = 1;
    min_req_order = 1;
    mynres = 0;
    skew = 0;
    reset_count("normal_results");
/*    if (t = find_table("normal_results")) delete_table(t);*/
    return;
  }
  if (!(t = find_table("normal_results")))
  {
    /* initialise table */
    normal_results = make_table("normal_results", "normal_res", normal_res_cols,
                                normal_res_types, MAX_ROWS);
    normal_results->dynamic = 1;
    add_to_table_list(normal_results, table_register);
    reset_count("normal_results");
    t = find_table("normal_results");
    min_order = 1;
    min_req_order = 1;
  }

  /* initialise order array */
  order[0] = zero;
  order[1] = zero;
  order[2] = zero;
  order[3] = zero;
  if (t->curr == t->max) grow_table(t);

  for (j = 0; j < PTC_NAMES_L; j++)
  {
    /* Treat each ptc variable */

    if (command_par(names[j], this_cmd->clone, &cp))
    {
      curr = cp->m_string->curr;
      if (curr > 4)
        printf("Too many values for the attribute %s. Only the first four are retained.\n",names[j]);
      for (i = 0; i < curr; i++)
      {
        order_list = cp->m_string->p[i];
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

  int one = 1;
  int parexist = -1;
  double value = 0;
  int ivalue = 0;

  struct command* keep_beam = current_beam;

  if (attach_beam(current_sequ) == 0)
    fatal_error("PTC_TRACKLINE - sequence without beam:", current_sequ->name);

  if (command_par_string_or_calldef("file", cmd->clone, &track_filename))
    set_option("track_dump", &one);
  if (!track_filename)
    track_filename = permbuff("dummy");
  track_filename = permbuff(track_filename);

  command_par_string_or_calldef("extension", cmd->clone, &track_fileext);
  if (!track_fileext)
    track_fileext = permbuff("\0");
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

  // LD 2016.04.19
  adjust_beam();
  probe_beam = clone_command(current_beam);
  adjust_probe_fp(0);

  track_tables_delete(); /* deleting all track related tables*/

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
    mad_error("pro_ptc_enforce6d","Command Definintion is null!!!");
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
  int found;

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

  /*DEBUG DEBUG LEVEL*/
  if ( name_list_pos("debuglevel", nl) >=0 )
  {
    found = command_par_value2("debuglevel", cmd->clone, &switchvalue);
    debuglevel = (int)switchvalue;
    w_ptc_setdebuglevel_(&debuglevel);
  }
  else
  {
    printf("debuglevel is not present (keeping current value)\n");
  }

  /*MAPDUMP LEVEL*/
  if ( name_list_pos("mapdump", nl) >=0 )
  {
    found = command_par_value2("mapdump", cmd->clone, &switchvalue);
    int mapdump = (int)switchvalue;
    w_ptc_setmapdumplevel_(&mapdump);
  }

  /*ACCELERATION SWITCH*/
  found = command_par_value_user2("maxacceleration", cmd->clone, &switchvalue);
  if (found)
   {
     if (debuglevel > 0) printf("maxaccel is found and its value is %f\n", switchvalue);
     i = (int)switchvalue;
     w_ptc_setaccel_method_(&i);
   }
  else
   {
     if (debuglevel > 0) printf("maxaccel is not present (keeping current value)\n");
   }

  /*EXACT_MIS SWITCH*/
  found = command_par_value_user2("exact_mis", cmd->clone, &switchvalue);
  if (found)
   {
     if (debuglevel > 0) printf("exact_mis is found and its value is %f\n", switchvalue);
     i = (int)switchvalue;
     w_ptc_setexactmis_(&i);
   }
  else
   {
     if (debuglevel > 0)  printf("exact_mis is not present (keeping current value)\n");
   }


  /*radiation SWITCH*/
  found = command_par_value_user2("radiation", cmd->clone, &switchvalue);
  if (found)
   {
    if (debuglevel > 0) printf("radiation is found and its value is %f\n", switchvalue);
    i = (int)switchvalue;
    w_ptc_setradiation_(&i);
   }
  else
   {
    if (debuglevel > 0) printf("radiation is not present (keeping current value)\n");
   }


  /*modulation SWITCH*/
  found = command_par_value_user2("modulation", cmd->clone, &switchvalue);
  if (found)
   {
    if (debuglevel > 0) printf("modulation is found and its value is %f\n", switchvalue);
    i = (int)switchvalue;
    w_ptc_setmodulation_(&i);
   }
  else
   {
    if (debuglevel > 0) printf("modulation is not present (keeping current value)\n");
   }

  /*stochastic SWITCH*/
  found = command_par_value_user2("stochastic", cmd->clone, &switchvalue);
  if (found)
   {
    if (debuglevel > 0) printf("stochastic is found and its value is %f\n", switchvalue);
    i = (int)switchvalue;
    w_ptc_setstochastic_(&i);
   }
  else
   {
    if (debuglevel > 0) printf("stochastic is not present (keeping current value)\n");
   }

  /*envelope SWITCH*/
  found = command_par_value_user2("envelope", cmd->clone, &switchvalue);
  if (found)
   {
    if (debuglevel > 0) printf("envelope is found and its value is %f\n", switchvalue);
    i = (int)switchvalue;
    w_ptc_setenvelope_(&i);
   }
  else
   {
    if (debuglevel > 0) printf("envelope is not present (keeping current value)\n");
   }

  /*fringe SWITCH*/
  found = command_par_value_user2("fringe", cmd->clone, &switchvalue);
  if (found)
   {
    if (debuglevel > 0) printf("fringe is found and its value is %f\n", switchvalue);
    i = (int)switchvalue;
    w_ptc_setfringe_(&i);
   }
  else
   {
    if (debuglevel > 0) printf("fringe is not present (keeping current value)\n");
   }


  /*totalpath SWITCH*/
  found = command_par_value_user2("totalpath", cmd->clone, &switchvalue);
  if (found)
   {
    if (debuglevel > 0) printf("totalpath is found and its value is %f\n", switchvalue);
    i = (int)switchvalue;
    w_ptc_settotalpath_(&i);
   }
  else
   {
    if (debuglevel > 0) printf("totalpath is not present (keeping current value)\n");
   }


  /*TIME SWITCH*/
  found = command_par_value_user2("time", cmd->clone, &switchvalue);
  if (found)
   {
    if (debuglevel > 0) printf("time is found and its value is %f\n", switchvalue);
    i = (int)switchvalue;
    w_ptc_settime_(&i);
   }
  else
   {
    if (debuglevel > 0) printf("time is not present (keeping current value)\n");
   }

  /*NOCAVITY SWITCH*/
  found = command_par_value_user2("nocavity", cmd->clone, &switchvalue);
  if (found)
   {
    if (debuglevel > 0) printf("nocavity is found and its value is %f\n", switchvalue);
    i = (int)switchvalue;
    w_ptc_setnocavity_(&i);
   }
  else
   {
    if (debuglevel > 0) printf("nocavity is not present (keeping current value)\n");
   }

  /*SEED SETTING*/
  found = command_par_value_user2("seed", cmd->clone, &switchvalue);
  if (found)
   {
    if (debuglevel > 0) printf("seed is found and its value is %f\n", switchvalue);
    i = (int)switchvalue;
    w_ptc_setseed_(&i);
   }
  else
   {
    if (debuglevel > 0) printf("seed is not present (keeping current value)\n");
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
  int                            pos         = 0;
  int                            k         = 0;
  struct node*                   nodes[2]={0x0,0x0};
  struct node*                   anode=0x0;
  char*                          element;
  int                            refframe=-1;/*0 global, 1 current position, 2 end face if the previous element*/

  char* s_refframe = command_par_string_user("refframe", cmd->clone);
  if (s_refframe)
  {
    /*if it is zero it is not specified*/
    /*printf("refframe is %s.\n", s_refframe );*/

    if ( strcmp(s_refframe,"current")  == 0 )
     {
       refframe = 1;
     }

    if ( strcmp(s_refframe,"previouselement") == 0 )
     {
       refframe = 2;
     }
     
     if ( strcmp(s_refframe,"gcs") == 0 )
     {
       refframe = 0;
     }

    if (refframe < 0)
     {
       warning("mad_ptc.c: pro_ptc_eplacement: did not recognize string describing refframe, using default  ", s_refframe);
       refframe = 0;
     }

  }
    else
    {
      warning("mad_ptc.c: pro_ptc_eplacement: string describing refframe is null: ", "using default");
      refframe = 0;
    }


  element = command_par_string_user("range", cmd->clone);
  if ( !element )
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
  int                            pos         = 0;
  int                            k         = 0;
  struct node*                   nodes[2]={0x0,0x0};
  struct node*                   anode=0x0;
  char*                          element;

  element  = command_par_string_user("element", cmd->clone);
  if ( !element )
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
  static const int tabnamelen = 48;
  char  tablename[tabnamelen];
  char  colname[9];
  int   clen = 0;
  struct int_array*              tabIA      = 0x0;
  struct int_array*              mdefIA      = 0x0;
  struct command_parameter_list* c_parameters= cmd->clone->par;
  struct name_list*              c_parnames  = cmd->clone->par_names;
  struct command_parameter*      cp;
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
    strncpy(tablename,"moments",tabnamelen);
  }
  else if ( pchar[0] == 0 )
  {
    strncpy(tablename,"moments",tabnamelen);
  }
  else
  {
    strncpy(tablename,pchar,tabnamelen);
  }

  tabIA = new_int_array(1+strlen(tablename));
  conv_char(tablename,tabIA);

  pos = name_list_pos("moment_s", c_parnames);
  if ( pos < 0)
  {
    printf("Weired: moments parameter is not defined\n");
    return 1;
  }

  if (command_par("moment_s", cmd->clone, &cp))
  {
    for (j = 0; j < cp->m_string->curr; j++)
    {
      strcpy(colname,"mu000000");

      mdefin = cp->m_string->p[j];

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
        strncpy(tablename,pchar,tabnamelen);
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
  enum { MAXTABLES = 100 };
  char*             tables[MAXTABLES];
  struct name_list* cols[MAXTABLES];
  struct table*     t;
  static const int  tabnamelen = 48;
  char              tabname[tabnamelen];
  char              colname[17];
  int               nmom;
  int i,j; // ,k; not used

  memset(tables, 0, MAXTABLES*sizeof(char*));

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

    if (tables[j] == 0x0) {
      int len = strlen(tabname)+1;
      tables[j] = mymalloc_atomic("makemomentstables", len * sizeof *tables[0]);
      strcpy(tables[j],tabname);
      cols[j] = new_name_list("columns", 15);
      add_to_name_list(permbuff("name"),3,cols[j]);
      add_to_name_list(permbuff("s"),2,cols[j]);
    }

    add_to_name_list(permbuff(colname),2,cols[j]); // k = not used
  }

  if (moments_tables) {
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
    if (t->num_cols > t->org_cols)  add_vars_to_table(t,1);
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
  int k=0, one = 1;
/*  const char *rout_name = "ptc_track"; */
  int npart = stored_track_start->curr;
  struct table* t;
/*  int turns = command_par_value("turns", cmd->clone); */

  track_is_on = 1;
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
  k = get_value(current_command->name,"recloss");
  set_option("recloss", &k);
  k = get_value(current_command->name,"onetable");
  set_option("onetable", &k);
  track_deltap=get_value(current_command->name,"deltap");
  set_variable("track_deltap", &track_deltap);
  if(track_deltap != 0) fprintf(prt_file, v_format("track_deltap: %F\n"),
                                track_deltap);

  if (command_par_string_or_calldef("file", cmd->clone, &track_filename))
      set_option("track_dump", &one);
  if (!track_filename) track_filename = permbuff("dummy");
  track_filename = permbuff(track_filename);

  command_par_string_or_calldef("extension", cmd->clone, &track_fileext);
  if (!track_fileext)
    track_fileext = permbuff("\0");
  track_fileext = permbuff(track_fileext);

  if (npart == 0)
  {
    warning("track_run: no START command seen yet", "ignored");
    return;
  }


  track_tables_delete(); /* deleting all track related tables*/

  track_tables_create(cmd);
  if (debuglevel > 2)
   {
     printf("obs_points ptc_track: %d \n",curr_obs_points);
   }

  w_ptc_track_(&curr_obs_points);
  t = find_table("tracksumm");
  if (get_option("info"))  print_table(t);
  if (get_option("track_dump")) track_tables_dump();

 if (debuglevel > 1)
  {
    fprintf(prt_file, "\n*****  end of ptc_run  *****\n");
  }
}
/*_______________________________________________________*/

void printpoly(int p[6], int dim )
{
 int i;

 printf("f"); /*icase*/

 for (i=0; i<dim; i++)
  {
    printf("%1d",p[i]); /*icase*/
  }

 printf("\n");

}

/*_______________________________________________________*/
void makerdtstwisstable()
{
  int i;

  struct table* rdts_table;
  char** table_cols;
  int*  table_type;


  table_cols = mymalloc_atomic("",9*sizeof(char*));
  table_type = mymalloc_atomic("",9*sizeof(int));

  for (i=0; i<9; i++)
  {
    table_cols[i] = mymalloc_atomic("",10*sizeof(char));
    table_type[i] = 2;
  }
  table_type[0] = 3;

  strcpy(table_cols[0], "name");
  strcpy(table_cols[1], "s"); /*can not be s becuase it will not plot than*/
  strcpy(table_cols[2], "k1l");
  strcpy(table_cols[3], "k1sl");
  strcpy(table_cols[4], "k2l");
  strcpy(table_cols[5], "k2sl");
  strcpy(table_cols[6], "k3l");
  strcpy(table_cols[7], "k3sl");
  strcpy(table_cols[8], " ");


  char name[] = "twissrdt";

  rdts_table = make_table2(name, name, table_cols,
                           table_type, current_sequ->n_nodes);

  rdts_table->dynamic = 1;
  add_to_table_list(rdts_table, table_register);
  rdts_table->org_sequ = current_sequ;
  rdts_table->curr= 0;


}
