#include "madx.h"

// public interface

void
pro_taper(struct in_cmd* cmd)
  /* calls the taper module */
{
  struct command* taper = cmd->clone;
  double stepsize = 0.0;
  int iterate = 0, error = 0, t_file;
  char* filename = "no_taper_file";
  int debug = get_option("debug");
  
  if (current_sequ == NULL || current_sequ->ex_start == NULL) {
    warning("sequence not active,", "TAPER ignored");
    return;
  }

  if (attach_beam(current_sequ) == 0)
    fatal_error("TAPER - sequence without beam:", current_sequ->name);
  
  if (command_par_value("reset", taper) != 0) {
    taperreset_(&error); /* reset */
    fprintf(prt_file, "All taper values have been reset to zero \n\n"); 
    return;
  }

  iterate   = command_par_value("iterate", taper);
  stepsize  = command_par_value("stepsize", taper);

  if (par_present("file",taper)) {
    t_file = command_par_string_user2("file", taper, &filename);
  }
    
  if (iterate < 0) {
    warning("negative value for ITERATE, ", "reset to absolute value");
    iterate = abs(iterate);
  }

  if (iterate > 10) {
    warning("ITERATE value larger than 10 is probably useless, ", "reset to 10.");
    iterate = 10;
  }

  if (stepsize < 0) {
    warning("negative value for STEPSIZE, ", "reset to absolute value");    
    stepsize = fabs(stepsize);
  }  
  
  if (debug) {
    fprintf(prt_file, "\n Taper parameters: iterate = %d   stepsize = %e   file = %s\n\n" , iterate, stepsize, filename);
  }

  zero_double(orbit0, 6);
  
  if (guess_flag) { // guess_flag is set by COGUESS command
    copy_double(guess_orbit,orbit0,6);
    if (get_option("info")) { 
      fprintf(prt_file, " Found initial orbit vector from COGUESS values. \n");
    }
  }  
  
  taperreset_(&error); /* start from bare sequence for adjust_probe */
  adjust_beam();
  probe_beam = clone_command(current_beam);
  adjust_probe_fp(0);

  taper_(orbit0, &iterate, &stepsize, filename, &error); /* call taper module */
  
  probe_beam = delete_command(probe_beam);
}
