#include "madx.h"

void
exec_sodd(struct in_cmd* cmd)
{
  int ierr;

  /* use correct beam for sequence to be plotted - HG 031127 */
  struct command* keep_beam = current_beam;


  if (attach_beam(current_sequ) == 0)
    fatal_error("TWISS - sequence without beam:", current_sequ->name);
  /* end part1 of HG 031127 */

  /* <JMJ 7/11/2002> The following ifndef exclusion is a quick fix so that
     the WIN32 version
     does not try to do X11 graphics. However this has the consequence that
     the program will not make Postscript files.  HG needs to separate these things.
     </JMJ 7/11/2002> */
  /*FS 27.03.2004 works now on Windows using gxx11ps.F and gxx11psc.c courtesy HG */

  /* get nosixtrack */

  if (!this_cmd || !this_cmd->clone)
    fatal_error("SODD "," - No existing command");

  if(!par_present("nosixtrack", this_cmd->clone))
  {
    printf("Build-up of input file fc.34 by call to program sixtrack. \n");
    conv_sixtrack(cmd);
    // fclose(f34);
    printf("input file fc.34 is ready. \n");
  }
  sodd_table_70 = make_table("detune_1_end", "sodd_detune_5", sodd_detune_5_cols,
                             sodd_detune_5_types, 2);
  sodd_table_70->dynamic = 1;
  add_to_table_list(sodd_table_70, table_register);
  sodd_table_71 = make_table("detune_1_all", "sodd_detune_5", sodd_detune_5_cols,
                             sodd_detune_5_types, 2);
  sodd_table_71->dynamic = 1;
  add_to_table_list(sodd_table_71, table_register);
  sodd_table_72 = make_table("detune_2_hor", "sodd_detune_5", sodd_detune_5_cols,
                             sodd_detune_5_types, 2);
  sodd_table_72->dynamic = 1;
  add_to_table_list(sodd_table_72, table_register);
  sodd_table_73 = make_table("detune_2_ver", "sodd_detune_5", sodd_detune_5_cols,
                             sodd_detune_5_types, 2);
  sodd_table_73->dynamic = 1;
  add_to_table_list(sodd_table_73, table_register);
  sodd_table_74 = make_table("distort_1_f_end", "sodd_distort1_8", sodd_distort1_8_cols,
                             sodd_distort1_8_types, 2);
  sodd_table_74->dynamic = 1;
  add_to_table_list(sodd_table_74, table_register);
  sodd_table_75 = make_table("distort_1_h_end", "sodd_distort1_8", sodd_distort1_8_cols,
                             sodd_distort1_8_types, 2);
  sodd_table_75->dynamic = 1;
  add_to_table_list(sodd_table_75, table_register);
  sodd_table_76 = make_table("distort_1_f_all", "sodd_distort1_11", sodd_distort1_11_cols,
                             sodd_distort1_11_types, 2);
  sodd_table_76->dynamic = 1;
  add_to_table_list(sodd_table_76, table_register);
  sodd_table_77 = make_table("distort_1_h_all", "sodd_distort1_11", sodd_distort1_11_cols,
                             sodd_distort1_11_types, 2);
  sodd_table_77->dynamic = 1;
  add_to_table_list(sodd_table_77, table_register);
  sodd_table_78 = make_table("distort_2_f_end", "sodd_distort2_9", sodd_distort2_9_cols,
                             sodd_distort2_9_types, 2);
  sodd_table_78->dynamic = 1;
  add_to_table_list(sodd_table_78, table_register);
  sodd_table_79 = make_table("distort_2_h_end", "sodd_distort2_9", sodd_distort2_9_cols,
                             sodd_distort2_9_types, 2);
  sodd_table_79->dynamic = 1;
  add_to_table_list(sodd_table_79, table_register);
  soddin_(&ierr);

  /* part 2 of HG 031127 */
  current_beam = keep_beam;
  /* end of part 2 of HG 031127 */
}

