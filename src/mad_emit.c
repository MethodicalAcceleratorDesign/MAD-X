#include "madx.h"

// public interface

void
pro_emit(struct in_cmd* cmd)
  /* calls the emit module */
{
  struct command* emit = cmd->clone;
  double e_deltap, e_tol, u0;
  int j, error, keep;
  double emit_v[3], nemit_v[3], bmax[9], gmax[9], dismax[4], tunes[3],
    sig_v[4], pdamp[3], r0mat[4];
  char tmp[100];
  int updatebeam;

  fprintf(prt_file, "enter EMIT module\n");

  if (current_sequ == NULL || current_sequ->ex_start == NULL) {
    warning("sequence not active,", "EMIT ignored");
    return;
  }

  if (attach_beam(current_sequ) == 0)
    fatal_error("EMIT - sequence without beam:", current_sequ->name);

  e_deltap = command_par_value("deltap", emit);
  e_tol = command_par_value("tol", emit);

  keep = get_option("twiss_print");
  j = 0;
  set_option("twiss_print", &j);

  zero_double(orbit0, 6);
  zero_double(disp0, 6);
  zero_double(oneturnmat, 6*6);

  adjust_beam();
  probe_beam = clone_command(current_beam);

  tmrefe_(oneturnmat); /* one-turn linear transfer map */
  twcpin_(oneturnmat,disp0,r0mat,&error); /* added for disp0 computation */

  adjust_probe(e_deltap); /* sets correct gamma, beta, etc. */
  print_probe();
  adjust_rfc(); /* sets freq in rf-cavities from probe */

  // guess_flag is set by COGUESS command
  if (guess_flag) copy_double(guess_orbit, orbit0, 6);

  double tt[6*6*6] = {0};
  getclor_(orbit0, oneturnmat, tt, &error); /* closed orbit */

  if (error == 0) {
    current_node = current_sequ->ex_start;
    emit_(&e_deltap, &e_tol, orbit0, disp0, oneturnmat, &u0, emit_v, nemit_v,
          bmax, gmax, dismax, tunes, sig_v, pdamp, &updatebeam);

    if (e_deltap != zero) { 
      sprintf(tmp, v_format("%F"), e_deltap);
      warning("EMIT: beam not updated, non-zero deltap: ", tmp);
    }
    else if (updatebeam) {      
      store_comm_par_value("ex", emit_v[0], current_beam);
      store_comm_par_value("exn", nemit_v[0], current_beam);
      store_comm_par_value("ey", emit_v[1], current_beam);
      store_comm_par_value("eyn", nemit_v[1], current_beam);
      store_comm_par_value("et", emit_v[2], current_beam);
      store_comm_par_value("sigt", sig_v[2], current_beam);
      store_comm_par_value("sige", sig_v[3], current_beam);
      store_comm_par_value("u0", u0, current_beam);
      store_comm_par_value("qs", tunes[2], current_beam);
      store_comm_par_vector("pdamp", pdamp, current_beam);
      printf("\n EMIT: beam parameters have been updated.\n"); 
    }
    else warning("EMIT: beam not updated, RADIATE is false or longitudinal stability not ensured.", ""); 

    print_rfc();
  }

  set_option("twiss_print", &keep);
}

