#include "madx.h"

// private functions

static double
mult_par(const char* par, struct element* el)
  /* returns multipole parameter for par = "k0l" or "k0sl" etc. */
{
  double val = zero;
  char tmp[12] ,*p;
  strcpy(tmp, par);
  if (*tmp == 'k' && (p = strchr(tmp, 'l')) != NULL)
  {
    *p = '\0';  /* suppress trailing l */
    int skew = 0;
    if ((p = strchr(tmp, 's')) != NULL)
    {
      skew = 1; *p = '\0';
    }
    int k = 0;
    sscanf(&tmp[1], "%d", &k);
    double vect[FIELD_MAX];
    int l;
    if (skew) l = element_vector(el, "ksl", vect);
    else      l = element_vector(el, "knl", vect);
    if (k < l) val = vect[k];
  }
  return val;
}


static void
fill_beta0(struct command* beta0, struct node* node)
{
  /* makes uses of fact that beta0 and twiss_table have the same variables
     at the same place (+2)  up to energy inclusive */
  struct command_parameter_list* pl = beta0->par;
  struct name_list* nl = beta0->par_names;
  int i = -1, pos;
  if (twiss_table == NULL) return;
  for (pos = 0; pos < twiss_table->curr; pos++)
  {
    if (twiss_table->p_nodes[pos] == node)  break;
  }
  if (pos < twiss_table->curr)
  {
    do
    {
      i++;
      pl->parameters[i]->double_value = twiss_table->d_cols[i+3][pos];
      /* if (strstr(nl->names[i], "mu")) pl->parameters[i]->double_value *= twopi; frs 19.10.2006*/
    }
    while (strcmp(nl->names[i], "energy") != 0);
  }
}

static void
exec_savebeta(void)
  /* stores twiss values in a beta0 structure */
{
  struct node* nodes[2];
  struct command* beta0;
  char* label;
  int i, pos;
  for (i = 0; i < savebeta_list->curr; i++)
  {
    struct command* cmd = savebeta_list->commands[i];
    label = command_par_string("label", cmd);

    if (find_command(label, beta0_list) == NULL) { /* fill only once */
      const char* sequence = command_par_string_user("sequence", cmd);
      if (!sequence || strcmp(sequence, current_sequ->name) == 0) {
        if (get_ex_range(command_par_string("place", cmd), current_sequ, nodes)) {
          pos = name_list_pos("beta0", defined_commands->list);
          beta0 = clone_command(defined_commands->commands[pos]);
          fill_beta0(beta0, nodes[0]);
          add_to_command_list(label, beta0, beta0_list, 0);
        }
      }
    }
  }
}

static double*
col_data(struct table* t, const char* name)
{
    int pos = name_list_pos(name, t->columns);
    return t->d_cols[pos];
}

static void
fill_twiss_header(struct table* t)
  /* puts beam parameters etc. at start of twiss table */
{
  int i, h_length = 40; /* change adding header lines ! */
  struct table* s;
  char tmp[NAME_L];

  if (t == NULL) return;
  /* ATTENTION: if you add header lines, augment h_length accordingly */
  if (t->header == NULL)  t->header = new_char_p_array(h_length);

  strncpy(tmp, t->org_sequ->name, NAME_L);
  table_add_header(t, "@ SEQUENCE         %%%02ds \"%s\"", strlen(tmp),stoupper(tmp));
  i = get_string("beam", "particle", tmp);
  table_add_header(t, "@ PARTICLE         %%%02ds \"%s\"", i, stoupper(tmp));

  table_add_header(t, "@ MASS             %%le  %F", get_value("beam", "mass"));
  table_add_header(t, "@ CHARGE           %%le  %F", get_value("beam", "charge"));
  table_add_header(t, "@ ENERGY           %%le  %F", get_value("beam", "energy"));
  table_add_header(t, "@ PC               %%le  %F", get_value("beam", "pc"));
  table_add_header(t, "@ GAMMA            %%le  %F", get_value("beam", "gamma"));
  table_add_header(t, "@ KBUNCH           %%le  %F", get_value("beam", "kbunch"));
  table_add_header(t, "@ BCURRENT         %%le  %F", get_value("beam", "bcurrent"));
  table_add_header(t, "@ SIGE             %%le  %F", get_value("beam", "sige"));
  table_add_header(t, "@ SIGT             %%le  %F", get_value("beam", "sigt"));
  table_add_header(t, "@ NPART            %%le  %F", get_value("beam", "npart"));
  table_add_header(t, "@ EX               %%le  %F", get_value("beam", "ex"));
  table_add_header(t, "@ EY               %%le  %F", get_value("beam", "ey"));
  table_add_header(t, "@ ET               %%le  %F", get_value("beam", "et"));
  table_add_header(t, "@ BV_FLAG          %%le  %F", get_value("beam", "bv"));

  if ((s = find_table("summ")))
  {
    table_add_header(t, "@ LENGTH           %%le  %F", col_data(s, "length")[0]);
    table_add_header(t, "@ ALFA             %%le  %F", col_data(s, "alfa")[0]);
    table_add_header(t, "@ ORBIT5           %%le  %F", col_data(s, "orbit5")[0]);
    table_add_header(t, "@ GAMMATR          %%le  %F", col_data(s, "gammatr")[0]);
    table_add_header(t, "@ Q1               %%le  %F", col_data(s, "q1")[0]);
    table_add_header(t, "@ Q2               %%le  %F", col_data(s, "q2")[0]);
    table_add_header(t, "@ DQ1              %%le  %F", col_data(s, "dq1")[0]);
    table_add_header(t, "@ DQ2              %%le  %F", col_data(s, "dq2")[0]);
    table_add_header(t, "@ DXMAX            %%le  %F", col_data(s, "dxmax")[0]);
    table_add_header(t, "@ DYMAX            %%le  %F", col_data(s, "dymax")[0]);
    table_add_header(t, "@ XCOMAX           %%le  %F", col_data(s, "xcomax")[0]);
    table_add_header(t, "@ YCOMAX           %%le  %F", col_data(s, "ycomax")[0]);
    table_add_header(t, "@ BETXMAX          %%le  %F", col_data(s, "betxmax")[0]);
    table_add_header(t, "@ BETYMAX          %%le  %F", col_data(s, "betymax")[0]);
    table_add_header(t, "@ XCORMS           %%le  %F", col_data(s, "xcorms")[0]);
    table_add_header(t, "@ YCORMS           %%le  %F", col_data(s, "ycorms")[0]);
    table_add_header(t, "@ DXRMS            %%le  %F", col_data(s, "dxrms")[0]);
    table_add_header(t, "@ DYRMS            %%le  %F", col_data(s, "dyrms")[0]);
    table_add_header(t, "@ DELTAP           %%le  %F", col_data(s, "deltap")[0]);
    table_add_header(t, "@ SYNCH_1          %%le  %F", col_data(s, "synch_1")[0]);
    table_add_header(t, "@ SYNCH_2          %%le  %F", col_data(s, "synch_2")[0]);
    table_add_header(t, "@ SYNCH_3          %%le  %F", col_data(s, "synch_3")[0]);
    table_add_header(t, "@ SYNCH_4          %%le  %F", col_data(s, "synch_4")[0]);
    table_add_header(t, "@ SYNCH_5          %%le  %F", col_data(s, "synch_5")[0]);
    
  }
}


static void
pro_embedded_twiss(struct command* current_global_twiss)
  /* controls twiss embedded module */
{
  struct command* keep_beam = current_beam;
  struct command* keep_twiss;
  struct name_list* nl = current_twiss->par_names;
  struct table* twiss_tb;
  struct table* keep_table = NULL;
  char *filename = NULL, *name, *sector_name;
  double tol,tol_keep;
  double betx,bety,alfx,mux,alfy,muy,x,px,y,py,t,pt,dx,dpx,dy,dpy,wx,
    phix,dmux,wy,phiy,dmuy,ddx,ddpx,ddy,ddpy,
    r11,r12,r21,r22,s;
  int i, jt=0, lp, k_orb = 0, u_orb = 0, k = 1;
  int w_file, beta_def, inval = 1, chrom_flg; // ks, err not used
  int keep_info;
  int orbit_input = 0; // counter of number of elements of initial orbit given on command line

  /* Set embedded_flag */

  embedded_flag = 1;

  keep_info = get_option("info");
  i = keep_info * get_option("twiss_print");
  set_option("info", &i);

  if(get_option("twiss_print"))
    fprintf(prt_file, "enter Twiss module\n");

  /*
    start command decoding
  */
  name = command_par_string_user("sequence", current_twiss);
  if(name) { /* sequence specified */
    if ((lp = name_list_pos(name, sequences->list)) > -1)
      current_sequ = sequences->sequs[lp];
    else {
      warning("unknown sequence ignored:", name);
      return;
    }
  }

  if (current_sequ == NULL || current_sequ->ex_start == NULL) {
    warning("sequence not active,", "Twiss ignored");
    return;
  }

  if (attach_beam(current_sequ) == 0)
    fatal_error("TWISS - sequence without beam:", current_sequ->name);

  if (!current_sequ->tw_table) {
    warning("no TWISS table present", "TWISS command ignored");
    return;
  }

  char *table_name = current_sequ->tw_table->name;

  if (get_value(current_command->name,"sectormap") != 0) { // (ks = not used
    set_option("twiss_sector", &k);

    command_par_string_or_calldef("sectorfile", current_twiss, &sector_name);

    if ((sec_file = fopen(sector_name, "w")) == NULL)
      fatal_error("cannot open output file:", sector_name);
  }

  /* Find index to the twiss table */

  if((twiss_tb = find_table(table_name))) {
    if (twiss_tb->origin ==1) return; /* table is read, has no node pointers */
    for (jt = 0; jt < twiss_tb->curr; jt++) {
      // need `->previous` since the start node itself may be replaced in
      // `interpolate_node`:
      if (twiss_tb->p_nodes[jt]->previous == current_sequ->range_start->previous) break;
    }
  }

  name = command_par_string_user("useorbit", current_twiss);
  if(name) { /* orbit specified */
    if (current_sequ->orbits == NULL)
      warning("orbit not found, ignored: ", name);
    else {
      if ((u_orb = name_list_pos(name, current_sequ->orbits->names)) < 0)
        warning("orbit not found, ignored: ", name);
      else set_option("useorbit", &k);
    }
  }

  name = command_par_string_user("keeporbit", current_twiss);
  if(name) { /* orbit specified tw_table*/
    if (current_sequ->orbits == NULL)
      current_sequ->orbits = new_vector_list(10);
    else if (current_sequ->orbits->curr == current_sequ->orbits->max)
      grow_vector_list(current_sequ->orbits);
    if ((k_orb = name_list_pos(name, current_sequ->orbits->names)) < 0) {
      k_orb = add_to_name_list(permbuff(name), 0, current_sequ->orbits->names);
      current_sequ->orbits->vectors[k_orb] = new_double_array(6);
    }
    set_option("keeporbit", &k);
  }

  w_file = command_par_string_user2("file", current_twiss, &filename);
  if (w_file && !filename)
    filename = permbuff("dummy");

  tol_keep = get_variable("twiss_tol");
  if (par_present("tolerance", current_twiss)) {
    tol = command_par_value("tolerance", current_twiss);
    set_variable("twiss_tol", &tol);
  }

  chrom_flg = command_par_value("chrom", current_twiss);

  /*
    end of command decoding
  */

  /* Initialise Twiss parameters */

  keep_twiss = current_twiss;

  if ((beta_def = twiss_input(current_twiss)) < 0) {
    if (beta_def == -1) warning("unknown beta0,", "Twiss ignored");
    else if (beta_def == -2) warning("betx or bety missing,", "Twiss ignored");
    set_variable("twiss_tol", &tol_keep);
    return;
  }

  set_option("twiss_inval", &beta_def);
  set_option("twiss_summ", &k);
  set_option("twiss_chrom", &chrom_flg);
  set_option("twiss_save", &k);

  /* Read Twiss parameters */

  current_twiss = current_global_twiss;

  /*jt is the row number of previous element*/
  if (jt <= 0) warning("Embedded Twiss failed: ", "MAD-X continues"); //jump to cleanup
  else {
    double_from_table_row(table_name, "betx", &jt, &betx);
    double_from_table_row(table_name, "bety", &jt, &bety);
    double_from_table_row(table_name, "alfx", &jt, &alfx);
    double_from_table_row(table_name, "mux", &jt, &mux);
    /* mux = mux*twopi; frs 19.10.2006 */
    double_from_table_row(table_name, "alfy", &jt, &alfy);
    double_from_table_row(table_name, "muy", &jt, &muy);
    /* muy = muy*twopi; frs 19.10.2006 */
    double_from_table_row(table_name, "x", &jt, &x);
    double_from_table_row(table_name, "px", &jt, &px);
    double_from_table_row(table_name, "y", &jt, &y);
    double_from_table_row(table_name, "py", &jt, &py);
    double_from_table_row(table_name, "t", &jt, &t);
    double_from_table_row(table_name, "pt", &jt, &pt);
    double_from_table_row(table_name, "dx", &jt, &dx);
    double_from_table_row(table_name, "dpx", &jt, &dpx);
    double_from_table_row(table_name, "dy", &jt, &dy);
    double_from_table_row(table_name, "dpy", &jt, &dpy);
    double_from_table_row(table_name, "wx", &jt, &wx);
    double_from_table_row(table_name, "phix", &jt, &phix);
    double_from_table_row(table_name, "dmux", &jt, &dmux);
    double_from_table_row(table_name, "wy", &jt, &wy);
    double_from_table_row(table_name, "phiy", &jt, &phiy);
    double_from_table_row(table_name, "dmuy", &jt, &dmuy);
    double_from_table_row(table_name, "ddx", &jt, &ddx);
    double_from_table_row(table_name, "ddpx", &jt, &ddpx);
    double_from_table_row(table_name, "ddy", &jt, &ddy);
    double_from_table_row(table_name, "ddpy", &jt, &ddpy);
    double_from_table_row(table_name, "r11",&jt, &r11);
    double_from_table_row(table_name, "r12",&jt, &r12);
    double_from_table_row(table_name, "r21",&jt, &r21);
    double_from_table_row(table_name, "r22",&jt, &r22);
    double_from_table_row(table_name, "s",&jt, &s);

    /* Store these Twiss parameters as initial values */

    current_twiss = keep_twiss;
    set_value("twiss", "betx" , &betx);    nl->inform[name_list_pos("betx",nl)] = 1;
    set_value("twiss", "bety" , &bety);    nl->inform[name_list_pos("bety",nl)] = 1;
    set_value("twiss", "alfx" , &alfx);    nl->inform[name_list_pos("alfx",nl)] = 1;
    set_value("twiss", "mux", &mux);       nl->inform[name_list_pos("mux",nl)] = 1;
    set_value("twiss", "alfy", &alfy);     nl->inform[name_list_pos("alfy",nl)] = 1;
    set_value("twiss", "muy", &muy);       nl->inform[name_list_pos("muy",nl)] = 1;
    set_value("twiss", "x", &x);           nl->inform[name_list_pos("x",nl)] = 1;
    set_value("twiss", "px", &px);         nl->inform[name_list_pos("px",nl)] = 1;
    set_value("twiss", "y", &y);           nl->inform[name_list_pos("y",nl)] = 1;
    set_value("twiss", "py", &py);         nl->inform[name_list_pos("py",nl)] = 1;
    set_value("twiss", "t", &t);           nl->inform[name_list_pos("t",nl)] = 1;
    set_value("twiss", "pt", &pt);         nl->inform[name_list_pos("pt",nl)] = 1;
    set_value("twiss", "dx", &dx);         nl->inform[name_list_pos("dx",nl)] = 1;
    set_value("twiss", "dpx", &dpx);       nl->inform[name_list_pos("dpx",nl)] = 1;
    set_value("twiss", "dy", &dy);         nl->inform[name_list_pos("dy",nl)] = 1;
    set_value("twiss", "dpy", &dpy);       nl->inform[name_list_pos("dpy",nl)] = 1;
    set_value("twiss", "wx", &wx);         nl->inform[name_list_pos("wx",nl)] = 1;
    set_value("twiss", "phix", &phix);     nl->inform[name_list_pos("phix",nl)] = 1;
    set_value("twiss", "dmux", &dmux);     nl->inform[name_list_pos("dmux",nl)] = 1;
    set_value("twiss", "wy", &wy);         nl->inform[name_list_pos("wy",nl)] = 1;
    set_value("twiss", "phiy", &phiy);     nl->inform[name_list_pos("phiy",nl)] = 1;
    set_value("twiss", "dmuy", &dmuy);     nl->inform[name_list_pos("dmuy",nl)] = 1;
    set_value("twiss", "ddx", &ddx);       nl->inform[name_list_pos("ddx",nl)] = 1;
    set_value("twiss", "ddpx", &ddpx);     nl->inform[name_list_pos("ddpx",nl)] = 1;
    set_value("twiss", "ddy", &ddy);       nl->inform[name_list_pos("ddy",nl)] = 1;
    set_value("twiss", "ddpy", &ddpy);     nl->inform[name_list_pos("ddpy",nl)] = 1;
    set_value("twiss", "r11", &r11);       nl->inform[name_list_pos("r11",nl)] = 1;
    set_value("twiss", "r12", &r12);       nl->inform[name_list_pos("r12",nl)] = 1;
    set_value("twiss", "r21", &r21);       nl->inform[name_list_pos("r21",nl)] = 1;
    set_value("twiss", "r22", &r22);       nl->inform[name_list_pos("r22",nl)] = 1;

    summ_table = make_table("summ", "summ", summ_table_cols, summ_table_types, twiss_deltas->curr+1);
    add_to_table_list(summ_table, table_register);

    if (get_option("twiss_sector")) {
      reset_sector(current_sequ, 0);
      set_sector();
    }

    // 2014-May-30  12:33:48  ghislain: suppressed
    /* if (get_option("useorbit")) */
    /*   copy_double(current_sequ->orbits->vectors[u_orb]->a, orbit0, 6); */
    /* else if (guess_flag) { */
    /*   for (i = 0; i < 6; i++) { */
    /*     if (guess_orbit[i] != zero) orbit0[i] = guess_orbit[i]; */
    /*   } */
    /* } */

    // 2014-May-30  12:33:48  ghislain: modified order of priority
    //              and added input for values given on command line
    zero_double(orbit0, 6);

    if (guess_flag) {
      if (get_option("info"))
	      printf(" Found initial orbit vector from coguess values. \n");
      copy_double(guess_orbit,orbit0,6);
    }
    // if given, useorbit overrides coguess
    if (get_option("useorbit")) {
      if (get_option("info"))
	      printf(" Found initial orbit vector from twiss useorbit values. \n");
      copy_double(current_sequ->orbits->vectors[u_orb]->a, orbit0, 6);
    }
    // if given, orbit0 values from twiss command line modify individual values
    if (par_present("x",  current_twiss)) { orbit0[0] = command_par_value("x",  current_twiss); orbit_input++;}
    if (par_present("px", current_twiss)) { orbit0[1] = command_par_value("px", current_twiss); orbit_input++;}
    if (par_present("y",  current_twiss)) { orbit0[2] = command_par_value("y",  current_twiss); orbit_input++;}
    if (par_present("py", current_twiss)) { orbit0[3] = command_par_value("py", current_twiss); orbit_input++;}
    if (par_present("t",  current_twiss)) { orbit0[4] = command_par_value("t",  current_twiss); orbit_input++;}
    if (par_present("pt", current_twiss)) { orbit0[5] = command_par_value("pt", current_twiss); orbit_input++;}

    if (orbit_input > 0 && get_option("info"))
      printf(" Found %d initial orbit vector values from twiss command. \n", orbit_input);

    if (get_option("debug"))
      printf(" Initial orbit: %e %e %e %e %e %e\n", orbit0[0], orbit0[1], orbit0[2], orbit0[3], orbit0[4], orbit0[5]);
    // 2014-May-30  12:33:48  ghislain: end of modifications


    if(twiss_deltas->curr <= 0)
      fatal_error("PRO_TWISS_EMBEDDED "," - No twiss deltas");

    // LD 2016.04.19
    adjust_beam();
    probe_beam = clone_command(current_beam);

    for (i = 0; i < twiss_deltas->curr; i++) {

      const char *table_embedded_name = "embedded_twiss_table";
      struct int_array* tarr;
      struct int_array* dummy_arr; /* for the new signature of the twiss() Fortran function*/

      tarr = new_int_array(strlen(table_embedded_name)+1);
      conv_char(table_embedded_name, tarr);
      dummy_arr = new_int_array(5+1);
      conv_char("dummy",dummy_arr);

      twiss_table = make_table(table_embedded_name, "twiss", twiss_table_cols,
                               twiss_table_types, current_sequ->n_nodes + current_sequ->num_interp);
      twiss_table->dynamic = 1; /* flag for table row access to current row */
      add_to_table_list(twiss_table, table_register);

      keep_table = current_sequ->tw_table;
      current_sequ->tw_table = twiss_table;

      twiss_table->org_sequ = current_sequ;

      current_node = current_sequ->range_start;
      set_option("twiss_inval", &inval);

      // LD 2016.04.19
      adjust_probe_fp(twiss_deltas->a[i]); /* sets correct gamma, beta, etc. */

      // CALL TWISS
      twiss_(oneturnmat, disp0, tarr->i, dummy_arr->i); /* different call */

      if ((twiss_success = get_option("twiss_success"))) {
        if (get_option("keeporbit"))
	        copy_double(orbit0, current_sequ->orbits->vectors[k_orb]->a, 6);
        fill_twiss_header(twiss_table);
        if (i == 0) exec_savebeta(); /* fill beta0 at first delta_p only */
        if (w_file) out_table(table_embedded_name, twiss_table, filename);
      }
      else warning("Twiss failed: ", "MAD-X continues");

      tarr = delete_int_array(tarr);
      dummy_arr = delete_int_array(dummy_arr);
    }

    if (sec_file) {
      fclose(sec_file);
      sec_file = NULL;
    }

    if (twiss_success && get_option("twiss_print")) print_table(summ_table);
  }


  /* cleanup */
  current_beam = keep_beam;
  probe_beam = delete_command(probe_beam);
  set_option("twiss_print", &k);
  k = 0;
  set_option("couple", &k);
  set_option("chrom", &k);
  set_option("rmatrix", &k);
  /* set_option("centre", &k); */
  set_option("twiss_sector", &k);
  set_option("keeporbit", &k);
  set_option("useorbit", &k);
  set_option("info", &keep_info);
  set_variable("twiss_tol", &tol_keep);
  current_sequ->tw_table = keep_table;

  /* Reset embedded_flag */
  embedded_flag = 0;
}

static void
set_twiss_deltas(struct command* comm)
{
  char* string;
  int i, k = 0, n = 0;
  double s, sign = one, ar[3];
  twiss_deltas->curr = 1;
  twiss_deltas->a[0] = 0;
  string = command_par_string_user("deltap", comm);
  if (string)
  {
    pre_split(string, c_dum, 0);
    mysplit(c_dum->c, tmp_p_array);
    while (k < tmp_p_array->curr)
    {
      for (i = k; i < tmp_p_array->curr; i++)
        if (*tmp_p_array->p[i] == ':') break;
      ar[n++] = double_from_expr(tmp_p_array->p, k, i-1);
      k = i + 1;
    }
    if (n == 1) twiss_deltas->a[0] = ar[0];
    else  /* there is a range given - fill array */
    {
      if (n == 2) ar[n++] = ar[1] - ar[0];
      if (ar[2] == zero) twiss_deltas->a[0] = ar[0];
      else if (ar[2] * (ar[1] - ar[0]) < zero)
        warning("illegal deltap range ignored:", string);
      else
      {
        twiss_deltas->a[0] = ar[0];
        if (ar[2] < zero) sign = -sign;
        for (s = sign * (ar[0] + ar[2]);
             s <= sign * ar[1]; s+= sign * ar[2])
        {
          if (twiss_deltas->curr == twiss_deltas->max)
          {
            sprintf(c_dum->c, "%d values", twiss_deltas->max);
            warning("deltap loop cut at", c_dum->c);
            break;
          }
          twiss_deltas->a[twiss_deltas->curr]
            = twiss_deltas->a[twiss_deltas->curr-1] + ar[2];
          twiss_deltas->curr++;
        }
      }
    }
  }
}

// public interface

void
copy_twiss_data(double* twiss_data, int* offset, int* nval, int* interp_index)
{
  if (!current_node->cl) return;
  for (int i = 0; i < current_node->cl->curr; ++i) {
    struct constraint* c = current_node->cl->constraints[i];
    if (c->index != *interp_index-1)
      continue;
    int n_pos = c->n_pos - 1;
    if (n_pos >= *offset && n_pos < *offset + *nval)
      c->evaluated = twiss_data[n_pos];
  }
}

void
get_disp0(double* disp)
{
  copy_double(disp0, disp, 6);
}

void
complete_twiss_table(struct table* t)
  /* fills all items missing after "twiss" into twiss table */
{
  int i, j, mult, n, myrbend;
  double el, val;
  struct node* c_node;
  char tmp[16];

  if (t == NULL) return;
  i = t->curr;
  c_node = current_node;
  mult = strcmp(c_node->base_name, "multipole") == 0 ? 1 : 0;
  t->s_cols[0][i] = tmpbuff(c_node->name);
  t->s_cols[1][i] = tmpbuff(c_node->base_name);
  t->s_cols[twiss_fill_end+1][i] = tmpbuff(c_node->p_elem->parent->name);
  t->s_cols[twiss_fill_end+2][i] = tmpbuff(command_par_string("comments",c_node->p_elem->def));
  for (j = twiss_opt_end+1; j<= twiss_fill_end; j++)
  {
    el = c_node->length;
    strcpy(tmp, twiss_table_cols[j]);
    myrbend = (strcmp(c_node->p_elem->base_type->name, "rbend") == 0);
    if (strcmp(twiss_table_cols[j], "l") == 0) val = el;
    else if (strcmp(tmp, "slot_id") == 0) val =  el_par_value(tmp, c_node->p_elem);
    else if (strcmp(tmp, "e1") == 0 || strcmp(tmp, "e2") == 0)
    {
      if(myrbend)
      {
        val =  el_par_value(tmp, c_node->p_elem) +
          c_node->other_bv * el_par_value("angle", c_node->p_elem) / two;
      }
      else
      {
        val =  el_par_value(tmp, c_node->p_elem);
      }
    }
    else if (strcmp(tmp, "assembly_id") == 0) val =  el_par_value(tmp, c_node->p_elem);
    else if (strcmp(tmp, "mech_sep") == 0) val =  el_par_value(tmp, c_node->p_elem);
    /*== jln 11.11.2010 dealt with the new property v_pos as for mech_sep */
    else if (strcmp(tmp, "v_pos") == 0 ) val = el_par_value(tmp, c_node->p_elem);
    /*==*/
    else if (strcmp(tmp, "lrad") == 0) val =  el_par_value(tmp, c_node->p_elem);
    else if (strcmp(tmp, "h1") == 0 &&
             strcmp(c_node->base_name, "dipedge") == 0) val = el_par_value("h", c_node->p_elem);
    else if (strcmp(tmp, "bbcharge") == 0) val =  el_par_value("charge", c_node->p_elem);
    else if (strcmp(tmp, "xma") == 0) val =  el_par_value(tmp, c_node->p_elem);
    else if (strcmp(tmp, "yma") == 0) val =  el_par_value(tmp, c_node->p_elem);
    else if (strcmp(tmp, "sigx") == 0) val =  el_par_value(tmp, c_node->p_elem);
    else if (strcmp(tmp, "sigy") == 0) val =  el_par_value(tmp, c_node->p_elem); 
    else if(mult)
    {
      if(j<=twiss_mult_end)
      {
        val = mult_par(twiss_table_cols[j], c_node->p_elem);
        val *= c_node->other_bv; /* dipole_bv kill initiative SF TR FS */
      }
      else
      {
        val = el_par_value(tmp, c_node->p_elem);
      }
    }
    else
    {
      strcpy(tmp, twiss_table_cols[j]);
      n = strlen(tmp) - 1;
      if (n > 1 && tmp[0] == 'k' && isdigit(tmp[1]) && tmp[n] == 'l')
        tmp[n] = '\0'; /* suppress trailing l in k0l etc. */
      if (el != zero && n > 1 && tmp[0] == 'k' && tmp[1] == 's' && tmp[n] == 'i')
        tmp[n] = '\0'; /* suppress trailing i in ksi */
      val = el_par_value(tmp, c_node->p_elem);
      if (n > 1 && tmp[0] == 'k' && isdigit(tmp[1]))
        val *= c_node->other_bv; /* dipole_bv kill initiative SF TR FS */
      else if (strstr(tmp, "kick") || strcmp(tmp, "angle") == 0 ||
               strcmp(tmp, "ks") == 0 || strcmp(tmp, "ksi") == 0 ||
               strcmp(tmp, "volt") == 0 )
        val *= c_node->other_bv; /* dipole_bv kill initiative SF TR FS */
      if (el != zero)
      {
        if (strstr(tmp,"kick") == NULL && strcmp(tmp, "angle")
            && strcmp(tmp, "tilt") && strcmp(tmp, "e1") && strcmp(tmp, "e2")
            && strcmp(tmp, "h1") && strcmp(tmp, "h2") && strcmp(tmp, "hgap")
            && strcmp(tmp, "fint") && strcmp(tmp, "fintx")
            && strcmp(tmp, "volt") && strcmp(tmp, "lag")
            && strcmp(tmp, "freq") && strcmp(tmp, "harmon")) val *= el;
      }
    }


    t->d_cols[j][i] = val;
  }
}

#if 0 // for debugging
void disp_beam_(void);
void disp_beam_(void)
{
  static int i = 0;
  i++;
  fprintf(stdout, "%d:current_sequ = %p '%s'\n"    , i, (void*)current_sequ, current_sequ->name);
  fprintf(stdout, "%d:current_beam = %p '%s%%%s'\n", i, (void*)current_beam, current_beam->name, current_beam->par->parameters[name_list_pos("sequence", current_beam->par_names)]->string);
}
#endif

void
pro_twiss(void)
  /* controls twiss module */
{
  struct command* keep_beam = current_beam;
//  struct int_array* tarr;
//  struct int_array* tarr_sector;
  struct node *nodes[2], *use_range[2];
  const char *table_name = NULL;
  char *filename = NULL, *name, *sector_name = NULL;
  char dummy[NAME_L] = "dummy", *sector_table_name = dummy; /* second string required by) */
  /* will be set to a proper string in case twiss_sector option selected */
  double tol,tol_keep, q1_val_p = 0, q2_val_p = 0, q1_val, q2_val, dq1, dq2;
  int i, j, lp, k_orb = 0, u_orb = 0, pos, k_save = 1, k = 1, k_sect, w_file, beta_def;
  int chrom_flg;
  int orbit_input = 0; // counter of number of elements of initial orbit given on command line
  int debug = get_option("debug");
  int keep_info = get_option("info");
  i = keep_info * get_option("twiss_print");
  set_option("info", &i);

  if (current_twiss == NULL) {
    seterrorflag(2,"pro_twiss","No twiss command seen yet");
    warning("pro_twiss","No twiss command seen yet!");
    return;
  }

  if (current_twiss->par_names == NULL) {
    seterrorflag(3,"pro_twiss","Last twiss has NULL par_names pointer. Cannot proceed further.");
    warning("pro_twiss","Last twiss has NULL par_names pointer. Cannot proceed further.");
    return;
  }

  if (current_twiss->par == NULL) {
    seterrorflag(4,"pro_twiss","Last twiss has NULL par pointer. Cannot proceed further.");
    warning("pro_twiss","Last twiss has NULL par pointer. Cannot proceed further.");
    return;
  }

  if (match_is_on) k_save = 0;  /* match gets its own variable transfer -
				    this can be overridden with option "slow"
                                    on match command */

  /*    start command decoding  */
  name = command_par_string_user("sequence", current_twiss);
  if (name) { /* sequence specified */
    if ((lp = name_list_pos(name, sequences->list)) > -1)
      current_sequ = sequences->sequs[lp];
    else {
      warning("unknown sequence ignored:", name);
      return;
    }
  }

  if (current_sequ == NULL || current_sequ->ex_start == NULL) {
    warning("sequence not active,", "Twiss ignored");
    return;
  }

  if(get_option("twiss_print"))
    fprintf(prt_file, "enter Twiss module\n");

  if (attach_beam(current_sequ) == 0)
    fatal_error("TWISS - sequence without beam:", current_sequ->name);

  if ((table_name = command_par_string_user("table", current_twiss))) { }
  else if ((table_name = command_par_string_user("save", current_twiss))) {
    k_save = 1;
  }
  else table_name = "twiss";

  if ((k_sect = get_value(current_command->name,"sectormap")) != 0) {
    set_option("twiss_sector", &k_sect);
    command_par_string_or_calldef("sectortable", current_twiss, &sector_table_name);
    command_par_string_or_calldef("sectorfile", current_twiss, &sector_name);
  }

  use_range[0] = current_sequ->range_start;
  use_range[1] = current_sequ->range_end;
  if ((name = command_par_string_user("range", current_twiss))) {
    if (get_sub_range(name, current_sequ, nodes)) {
      current_sequ->range_start = nodes[0];
      current_sequ->range_end = nodes[1];
    }
    else warning("illegal range ignored:", name);
  }
  for (j = 0; j < current_sequ->n_nodes; j++) {
    if (current_sequ->all_nodes[j] == current_sequ->range_start) break;
  }

  if((name = command_par_string_user("useorbit", current_twiss))) {
    /* orbit specified */
    if (current_sequ->orbits == NULL)
      warning("orbit not found, ignored: ", name);
    else {
      if ((u_orb = name_list_pos(name, current_sequ->orbits->names)) < 0)
        warning("orbit not found, ignored: ", name);
      else set_option("useorbit", &k);
    }
  }

  if(par_present("centre", current_twiss)) {
    set_option("centre", &k);
    // Special check that it really is set to TRUE (not only included in the argument)
    if(command_par_value("centre", current_twiss)>0){
      current_sequ->tw_centre=1;
    }
  }
  else {
    k = 0;
    current_sequ->tw_centre=0;
    set_option("centre", &k);
    k = 1;

  }

  name = command_par_string_user("keeporbit", current_twiss);
  if (name) { /* orbit specified */
    if (current_sequ->orbits == NULL)
      current_sequ->orbits = new_vector_list(10);
    else if (current_sequ->orbits->curr == current_sequ->orbits->max)
      grow_vector_list(current_sequ->orbits);
    if ((k_orb = name_list_pos(name, current_sequ->orbits->names)) < 0) {
      k_orb = add_to_name_list(permbuff(name), 0, current_sequ->orbits->names);
      current_sequ->orbits->vectors[k_orb] = new_double_array(6);
    }
    set_option("keeporbit", &k);
  }

  w_file = command_par_string_user2("file", current_twiss, &filename);
  if (w_file && !filename)
    filename = permbuff("dummy");

  tol_keep = get_variable("twiss_tol");
  if (par_present("tolerance", current_twiss)) {
    tol = command_par_value("tolerance", current_twiss);
    set_variable("twiss_tol", &tol);
  }

  chrom_flg = command_par_value("chrom", current_twiss);

  if(par_present("notable", current_twiss)) k_save = 0;

  /*    end of command decoding  */

  if ((beta_def = twiss_input(current_twiss)) < 0) {
    if (beta_def == -1) warning("unknown beta0,", "Twiss ignored");
    else if (beta_def == -2) warning("betx or bety missing,", "Twiss ignored");
    set_variable("twiss_tol", &tol_keep);
    return;
  }

  set_option("twiss_inval", &beta_def);
  set_option("twiss_summ", &k);
  set_option("twiss_chrom", &chrom_flg);
  set_option("twiss_save", &k_save);

  set_twiss_deltas(current_twiss);

  summ_table = make_table("summ", "summ", summ_table_cols, summ_table_types, twiss_deltas->curr+1);
  add_to_table_list(summ_table, table_register);

  /* now create the sector table */
  struct int_array* tarr_sector = NULL;
  tarr_sector = new_int_array(strlen(sector_table_name)+1);
  conv_char(sector_table_name, tarr_sector);

  if (get_option("twiss_sector")) {
    reset_sector(current_sequ, 0);
    set_sector();
    twiss_sector_table = make_table(sector_table_name, "sectormap",
                                    twiss_sector_table_cols,
                                    twiss_sector_table_types,
                                    1); /* start with one row at first */
    twiss_sector_table->dynamic = 1; /* flag for access to current row */
    add_to_table_list(twiss_sector_table, table_register);
  }

  // 2014-May-30  12:33:48  ghislain: modified order of priority
  //              and added input for values given on command line
  zero_double(orbit0, 6);

  if (guess_flag) {
    if (get_option("info"))
      printf(" Found initial orbit vector from coguess values. \n");
    copy_double(guess_orbit,orbit0,6);
  }
  // if given, useorbit overrides coguess
  if (get_option("useorbit")) {
    if (get_option("info"))
      printf(" Found initial orbit vector from twiss useorbit values. \n");
    copy_double(current_sequ->orbits->vectors[u_orb]->a, orbit0, 6);
  }
  // if given, orbit0 values from twiss command line modify individual values
  if (par_present("x",  current_twiss)) { orbit0[0] = command_par_value("x",  current_twiss); orbit_input++;}
  if (par_present("px", current_twiss)) { orbit0[1] = command_par_value("px", current_twiss); orbit_input++;}
  if (par_present("y",  current_twiss)) { orbit0[2] = command_par_value("y",  current_twiss); orbit_input++;}
  if (par_present("py", current_twiss)) { orbit0[3] = command_par_value("py", current_twiss); orbit_input++;}
  if (par_present("t",  current_twiss)) { orbit0[4] = command_par_value("t",  current_twiss); orbit_input++;}
  if (par_present("pt", current_twiss)) { orbit0[5] = command_par_value("pt", current_twiss); orbit_input++;}

  if (orbit_input > 0 && get_option("info"))
    printf(" Found %d initial orbit vector values from twiss command. \n", orbit_input);

  if (debug)
    printf(" Initial orbit: %e %e %e %e %e %e\n", orbit0[0], orbit0[1], orbit0[2], orbit0[3], orbit0[4], orbit0[5]);
  // 2014-May-30  12:33:48  ghislain: end of modifications

  // LD 2016.04.19
  adjust_beam();
  probe_beam = clone_command(current_beam);

  twiss_success = 1;
  set_option("twiss_success", &twiss_success);
  struct int_array* tarr = NULL;
  for (i = 0; i < twiss_deltas->curr; i++) {

    tarr = new_int_array(strlen(table_name)+1);
    conv_char(table_name, tarr);

    if (chrom_flg) { /* calculate chromaticity from tune difference - HG 6.2.09*/
      twiss_table = make_table(table_name, "twiss", twiss_table_cols, twiss_table_types, current_sequ->n_nodes + current_sequ->num_interp);
      twiss_table->dynamic = 1; /* flag for table row access to current row */
      add_to_table_list(twiss_table, table_register);
      current_sequ->tw_table = twiss_table;
      current_sequ->tw_valid = 1;
      twiss_table->org_sequ = current_sequ;

      // LD 2016.04.19
      adjust_probe_fp(twiss_deltas->a[i]+DQ_DELTAP);
      current_node = current_sequ->ex_start;
      twiss_(oneturnmat, disp0, tarr->i, tarr_sector->i); // CALL TWISS
      if ((twiss_success = get_option("twiss_success")) == 0) break;

      pos = name_list_pos("q1", summ_table->columns);
      q1_val_p = summ_table->d_cols[pos][i];
      pos = name_list_pos("q2", summ_table->columns);
      q2_val_p = summ_table->d_cols[pos][i];
    }

    if (k_save) {
      twiss_table = make_table(table_name, "twiss", twiss_table_cols, twiss_table_types, current_sequ->n_nodes + current_sequ->num_interp);
      twiss_table->dynamic = 1; /* flag for table row access to current row */
      add_to_table_list(twiss_table, table_register);
      current_sequ->tw_table = twiss_table;
      current_sequ->tw_valid = 1;
      twiss_table->org_sequ = current_sequ;
    }

    // LD 2016.04.19
    adjust_probe_fp(twiss_deltas->a[i]); /* sets correct gamma, beta, etc. */
    current_node = current_sequ->ex_start;
    twiss_(oneturnmat, disp0, tarr->i, tarr_sector->i);
    if ((twiss_success = get_option("twiss_success")) == 0) break;
    augment_count("summ ");

    if (chrom_flg) { /* calculate chromaticity from tune difference - HG 6.2.09*/
      pos = name_list_pos("q1", summ_table->columns);
      q1_val = summ_table->d_cols[pos][i];
      pos = name_list_pos("q2", summ_table->columns);
      q2_val = summ_table->d_cols[pos][i];

      double beta0   = get_value("probe","beta");
      dq1 = (q1_val_p - q1_val) / (beta0*DQ_DELTAP);
      dq2 = (q2_val_p - q2_val) / (beta0*DQ_DELTAP);

      pos = name_list_pos("dq1", summ_table->columns);
      summ_table->d_cols[pos][i] = dq1;
      pos = name_list_pos("dq2", summ_table->columns);
      summ_table->d_cols[pos][i] = dq2;
    }

    if (get_option("keeporbit"))
      copy_double(orbit0, current_sequ->orbits->vectors[k_orb]->a, 6);

    if (k_save) fill_twiss_header(twiss_table);

    if (i == 0) exec_savebeta(); /* fill beta0 at first delta_p only */

    if (k_save && w_file) out_table(table_name, twiss_table, filename);

    if ((twiss_deltas->curr > 1) && (i < twiss_deltas->curr-1)) {
      struct table *t = detach_table_from_table_list(table_name, table_register);
      if (t) {
        char new_table_name[NAME_L];
        sprintf(new_table_name, "%s_%d", table_name, i+1);
        rename_table(t, new_table_name);
        add_to_table_list(t, table_register);
      }
    }
    tarr = delete_int_array(tarr);
  } // i = 0 .. twiss_deltas->curr-1

  if (!twiss_success) {
    seterrorflag(1,"pro_twiss","TWISS failed");
    warning("Twiss failed: ", "MAD-X continues");
  } else {
    if (get_option("twiss_sector") && 0 != strcmp("", sector_name))
      out_table( sector_table_name, twiss_sector_table, sector_name );
    if (get_option("twiss_print"))
      print_table(summ_table);
  }

  /* cleanup */
  tarr = delete_int_array(tarr);
  tarr_sector = delete_int_array(tarr_sector);
  current_beam = keep_beam;
  probe_beam = delete_command(probe_beam);
  k = 0;
  set_option("couple", &k);
  set_option("chrom", &k);
  set_option("rmatrix", &k);
  set_option("twiss_sector", &k);
  set_option("keeporbit", &k);
  set_option("useorbit", &k);
  set_option("info", &keep_info);
  set_variable("twiss_tol", &tol_keep);
  current_sequ->range_start = use_range[0];
  current_sequ->range_end = use_range[1];
}

void print_eigenvectors_(double *eigenvectors){
  int i,p;
  char *filename;
  FILE *fptr;
  
  filename = command_par_string_user("eigenfile", current_twiss);
  if ((fptr = fopen(filename, "w")) == NULL)
      fatal_error("cannot open output file:", filename);
  
  for(i=0; i < 6; i ++){
    p=i*6;
    fprintf(fptr, "%.9e %.9e %.9e %.9e %.9e %.9e \n", eigenvectors[p],eigenvectors[p+1],
    eigenvectors[p+2],eigenvectors[p+3],eigenvectors[p+4],eigenvectors[p+5]);
  }
  fclose(fptr);
}
int
embedded_twiss(void)
  /* controls twiss module to create a twiss table for interpolated nodes
     between two elements */
{
  struct name_list* tnl; /* OB 31.1.2002: local name list for TWISS input definition */
  struct in_cmd* cmd;
  struct command* current_global_twiss;
  struct command_parameter* cp;
  char* embedded_twiss_beta[2];
  int j, tpos;
  int izero = 0;

  cmd = embedded_twiss_cmd;
  keep_tw_print = get_option("twiss_print");
  set_option("twiss_print", &izero);

  /* START defining a TWISS input command for default sequence */

  local_twiss[0] = new_in_cmd(10);
  local_twiss[0]->type = 0;
  local_twiss[0]->clone = local_twiss[0]->cmd_def
    = clone_command(find_command("twiss", defined_commands));
  tnl = local_twiss[0]->cmd_def->par_names;
  tpos = name_list_pos("sequence", tnl);
  local_twiss[0]->cmd_def->par->parameters[tpos]->string = current_sequ->name;
  local_twiss[0]->cmd_def->par_names->inform[tpos] = 1;

  /* END defining a TWISS input command for default sequence */

  if (current_sequ == NULL || current_sequ->ex_start == NULL)
  {
    warning("Command called without active sequence,", "ignored");
    return 1;
  }
  /* END CHK-SEQ; OB 1.2.2002 */

  for (j = 0; j < local_twiss[0]->cmd_def->par->curr; j++)
  {
    tnl = local_twiss[0]->cmd_def->par_names;
    tpos = name_list_pos("sequence", tnl);
    if (j != tpos) local_twiss[0]->cmd_def->par_names->inform[j] = 0;
  }

  /* START CHK-BETA-INPUT; OB 1.2.2002 */
  /* START CHK-BETA0; OB 23.1.2002 */
  if (command_par("beta0", cmd->clone, &cp))  /* parameter has been read */
  {
    /* beta0 specified */
    embedded_twiss_beta[0] = buffer(cp->m_string->p[0]);


    /* START defining a TWISS input command for the sequence */
    tnl = local_twiss[0]->cmd_def->par_names;
    tpos = name_list_pos("beta0", tnl);
    local_twiss[0]->cmd_def->par_names->inform[tpos] = 1;
    local_twiss[0]->cmd_def->par->parameters[tpos]->string = embedded_twiss_beta[0];
    /* END defining a TWISS input command for the sequence */
  }

  /* END CHK-BETA0; OB 23.1.2002 */

  /* END CHK-RANGE; OB 12.11.2002 */

  /* START CHK-USEORBIT; HG 28.1.2003 */
  if (command_par("useorbit", cmd->clone, &cp))  /* parameter has been read */
  {
    /* useorbit specified */
    /* START adding useorbit to TWISS input command for each sequence */
    tnl = local_twiss[0]->cmd_def->par_names;
    tpos = name_list_pos("useorbit", tnl);
    local_twiss[0]->cmd_def->par_names->inform[tpos] = 1;
    local_twiss[0]->cmd_def->par->parameters[tpos]->string
      = buffer(cp->m_string->p[0]);
    /* END adding range to TWISS input command for each sequence */
  }
  /* END CHK-USEORBIT; HG 28.1.2003 */

  /* START CHK-KEEPORBIT; HG 28.1.2003 */
  if (command_par("keeporbit", cmd->clone, &cp))  /* parameter has been read */
  {
    /* keeporbit specified */
    /* START adding keeporbit to TWISS input command for each sequence */
    tnl = local_twiss[0]->cmd_def->par_names;
    tpos = name_list_pos("keeporbit", tnl);
    local_twiss[0]->cmd_def->par_names->inform[tpos] = 1;
    local_twiss[0]->cmd_def->par->parameters[tpos]->string
      = buffer(cp->m_string->p[0]);
    /* END adding range to TWISS input command for each sequence */
  }
  /* END CHK-KEEPORBIT; HG 28.1.2003 */

  /* END CHK-BETA-INPUT; OB 1.2.2002 */

  /* START generating a TWISS table via 'pro_twiss'; OB 1.2.2002 */

  current_global_twiss = current_twiss;
  current_twiss = local_twiss[0]->clone;
  pro_embedded_twiss(current_global_twiss);

  /* END generating a TWISS table via 'pro_twiss' */
  current_twiss = current_global_twiss;

  return 0;
}

void
store_beta0(struct in_cmd* cmd)
{
  int k = cmd->decl_start - 1;
  if (k == 0) warning("beta0 without label:", "ignored");
  else {
    cmd->clone_flag = 1; /* do not delete */
    add_to_command_list(cmd->tok_list->p[0], cmd->clone, beta0_list, 0);
  }
}

void
store_savebeta(struct in_cmd* cmd)
{
  char* name = NULL;
//  struct command* comm; // not used
  if (log_val("clear", cmd->clone)) {
    delete_command_list(savebeta_list);
    savebeta_list = new_command_list("savebeta_list", 10);
    delete_command_list(beta0_list);
    beta0_list = new_command_list("beta0_list", 10);
  }
  else {
    if (!par_present("place", cmd->clone)) {
      warning("savebeta without place:", "ignored");
      return;
    }
    if (!(name = command_par_string_user("label", cmd->clone))) {
      warning("savebeta without label:", "ignored");
      return;
    }

    cmd->clone_flag = 1; /* do not delete */
    if ( find_command(name, beta0_list) ) remove_from_command_list(name, beta0_list);
    add_to_command_list(permbuff(name), cmd->clone, savebeta_list, 0);
  }
}

int
twiss_input(struct command* tw)
  /* returns -1 if an invalid beta0 given,
     returns -1 if only betx or bety given,
     returns 1 if betx and bety are given,
     or a valid beta0 (which is then loaded), else 0 */
{
  struct name_list* nl = tw->par_names;
  struct command_parameter_list* pl = tw->par;
  struct command* beta;
  int i = -1, ret = 0, pos, sb = 0;
  char* name;
  double val;
  if ((name = command_par_string_user("beta0", tw)))
  {
    if ((pos = name_list_pos(name, beta0_list->list)) > -1)
    {
      ret = 1;
      beta = beta0_list->commands[pos];
      do
      {
        i++;
        if (!par_present(nl->names[i], tw)) /* not read */
        {
          if (beta->par->parameters[i]->expr != NULL)
            val = expression_value(beta->par->parameters[i]->expr, 2);
          else val = beta->par->parameters[i]->double_value;
          pl->parameters[i]->double_value = val;
          nl->inform[name_list_pos(nl->names[i], nl)] = 1;
        }
      }
      while (strcmp(nl->names[i], "energy") != 0);
    }
    else ret = -1;
  }
  if (ret) return ret;
  /* if no beta0 given, betx and bety together set inval */
  if (par_present("betx", tw)) sb++;
  if (par_present("bety", tw)) sb++;
  if (sb)
  {
    if (sb < 2)  return -2;
    else         return 1;
  }
  else return 0;
}


