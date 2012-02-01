#include <signal.h>

#include "madx.h"

// private functions

static void
mad_init_c(void)
  /* initializes program */
{
  struct variable* var;
  int ione = 1;

  // make stdout unbuffered for proper I/O synchronization with fortran
  setvbuf(stdout, 0, _IONBF, 0);

  interactive = intrac();
  init55(123456789);          /* random generator */
  if (watch_flag == 1)  debug_file = fopen("madx.debug", "w");
  else if (watch_flag == 2)  debug_file = stdout;
  if (stamp_flag == 1)  stamp_file = fopen("madx.stamp", "w");
  else if (stamp_flag == 2)  stamp_file = stdout;
  in = new_in_buff_list(100); /* list of input buffers, dynamic */
  in->input_files[0] = stdin;
  prt_file = stdout;
  pro = new_in_buff_list(100); /* list of process buffers, dynamic */
  pro->buffers[0] = new_in_buffer(IN_BUFF_SIZE);
  pro->curr = 1;
  c_dum = new_char_array(AUX_LG);
  c_join = new_char_array(AUX_LG);
  work = new_char_array(AUX_LG);
  l_wrk = new_char_array(AUX_LG);
  char_buff = new_char_array_list(100); /* list of character arrays, dynamic */
  char_buff->ca[char_buff->curr++] = new_char_array(CHAR_BUFF_SIZE);
  aux_buff = new_char_array(AUX_LG);  /* dynamic temporary buffer */
  drift_list = new_el_list(1000); /* dynamic list for internal drifts */
  variable_list = new_var_list(2000); /* dynamic list of variables */
  comm_constraints = new_constraint_list(10); /* dynamic constraint list */
  beam_list = new_command_list("beam_list", 10); /* dynamic beam list */
  stored_track_start = new_command_list("track_start", 100); /* dynamic */
  table_deselect = new_command_list_list(10); /* dynamic */
  table_select = new_command_list_list(10); /* dynamic */
  defined_commands = new_command_list("defined_commands", 100); /* dynamic */
  stored_commands = new_command_list("stored_commands", 500); /* dynamic */
  line_list = new_macro_list(100); /* dynamic */
  macro_list = new_macro_list(100); /* dynamic */
  base_type_list = new_el_list(60); /* dynamic */
  element_list = new_el_list(20000); /* dynamic */
  buffered_cmds = new_in_cmd_list(10000); /* dynamic */
  sequences = new_sequence_list(20); /* dynamic */
  match_sequs = new_sequence_list(2);
  selected_ranges = new_node_list(10000); /* dynamic */
  selected_elements = new_el_list(10000); /* dynamic */
  tmp_p_array = new_char_p_array(1000); /* dynamic */
  tmp_l_array = new_char_p_array(1000); /* dynamic */
  sxf_list = new_name_list("sxf_list", 50); /* dynamic */
  all_table_lists = new_table_list_list(10); /* dynamic */
  deco_init();
  get_defined_constants();
  get_defined_commands();
  get_sxf_names();
  pi = get_variable("pi");
  twopi = two * pi;
  degrad = 180 / pi;
  raddeg = pi / 180;
  e = get_variable("e");
  clight = get_variable("clight");
  hbar = get_variable("hbar");
  var = new_variable("twiss_tol", 1.e-6, 1, 1, NULL, NULL);
  add_to_var_list(var, variable_list, 1);
  title = permbuff("no-title");
  set_defaults("option");
  set_defaults("beam");
  add_to_command_list("default_beam", current_beam, beam_list, 0);
  set_defaults("set");
  set_defaults("setplot");
  set_defaults("threader");
  table_register = new_table_list(10); /* dynamic */
  beta0_list = new_command_list("beta0_list", 10); /* dynamic */
  savebeta_list = new_command_list("savebeta_list", 10); /* dynamic */
  seqedit_select = /* dynamic - for "select seqedit" commands */
    new_command_list("seqedit_select", 10);
  error_select = /* dynamic - for "select error" commands */
    new_command_list("error-select", 10);
  save_select = /* dynamic - for "select save" commands */
    new_command_list("save_select", 10);
  slice_select = /* dynamic - for "select makethin" commands */
    new_command_list("slice_select", 10);
  sector_select = /* dynamic - for "select sectormap" commands */
    new_command_list("sector_select", 10);
  s_range = new_int_array(10); /* dynamic */
  e_range = new_int_array(10); /* dynamic */
  sd_range = new_int_array(10); /* dynamic */
  ed_range = new_int_array(10); /* dynamic */
  zero_double(orbit0, 6);
  zero_double(disp0, 6);
  zero_double(guess_orbit,6);
  zero_double(oneturnmat, 36);
  set_option("twiss_print", &ione);
}

static void
set_sigterm(void)
{
  /* provide a termination routine for access to memory outside scope */
  if (signal(SIGSEGV, mad_mem_handler) == SIG_IGN)
    signal(SIGSEGV, SIG_IGN);
}

// public functions

void
madx_start(void)
  /* prints start message after having read madxdict.h */
{
  struct tm* tm;

  // set signal handlers (temporary)
  set_sigterm();

  // init C and Fortran
  mad_init_c();
  mad_init_f_();

  /*  setbuf(stdout,(char *)0); */ /* no buffering - for debugging */
  time(&start_time); /* initialize timing */
  tm = localtime(&start_time); /* split system time */
  last_time = start_time;

  // compute padding of OSTYPE
  const char *pad[] = { "", " ", "  ", "    " };
  const int pad_sz = sizeof pad/sizeof *pad;
  int pad_idx = strlen("Windows")-strlen(version_ostype);
  if (pad_idx >= pad_sz) pad_idx = pad_sz-1; 

  printf("\n  +++++++++++++++++++++++++++++++++++++++++++\n");
  printf("  +    %s  (%s bit, %s) %s    +\n", version_name, version_arch, version_ostype, pad[pad_idx]);
  printf("  +    %s     +\n", version_name[strlen(version_name)-1]=='0' ? version_type_pro : version_type_dev);
  printf("  + %s      +\n", version_date_mod);
  printf("  + Execution Time Stamp: %02d.%02d.%02d %02d.%02d.%02d +\n",
         tm->tm_mday, tm->tm_mon+1, tm->tm_year%100,
         tm->tm_hour, tm->tm_min, tm->tm_sec);
  printf("  +++++++++++++++++++++++++++++++++++++++++++\n");
}

void
madx_input(int top)
  /* loops over input until end of execution */
{
  while (in_stop == 0)
  {
    if (interactive && in->curr == 0) puts("X: ==>");
    if (return_flag || get_stmt(in->input_files[in->curr], 0) == 0)
    {
      if (in->curr == 0) return;
      fclose(in->input_files[in->curr--]);
      return_flag = 0;
      if (in->curr == top) return;
    }
    else
    {
      stolower_nq(in->buffers[in->curr]->c_a->c);
      pro_input(in->buffers[in->curr]->c_a->c);
      if (stop_flag)  return;
    }
  }
}

void
madx_finish(void)
  /* write the termination message, pass madx.ps through ps2ps */
{
  int warn_numb, warn_numbf, nwarnings;
  
  /* should work with Lahey on windows 24.03.2004 */

  match2_delete_expressions();
  match2_delete_arrays();

  if (final_message == 0)
  {
    final_message = 1;
    if (plots_made)
    {
      gxterm_();
      if(system("which ps2ps > tmp_plot.ps") == 0)
	    {
         system("cp madx.ps tmp_plot.ps");
         system("ps2ps tmp_plot.ps madx.ps");
	    }
      system("rm tmp_plot.ps");
    }
    mad_err_getwarn(&warn_numb, &warn_numbf);
    nwarnings = warn_numb + warn_numbf;
    printf("\n  Number of warnings: %d\n",nwarnings);
    if (nwarnings > 0)
    {
      printf("%d in C and %d in Fortran\n",warn_numb,warn_numbf);
    }
    if (get_option("trace")) time_stamp("end");
    printf("\n  ++++++++++++++++++++++++++++++++++++++++++++\n");
    printf("  + %s %s finished normally +\n", version_name, version_arch);
    printf("  ++++++++++++++++++++++++++++++++++++++++++++\n");
  }
}


