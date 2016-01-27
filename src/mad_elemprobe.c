#include "madx.h"

void
adjust_probe(double delta_p)
  /* adjusts beam parameters to the current deltap */
{
  int j;
  double etas, slope, qs, fact, tmp, ds = oneturnmat[34];
  double alfa, beta, gamma, dtbyds, circ, freq0; // , deltat // not used
  double betas, gammas, et, sigt, sige;

  et = command_par_value("et", current_beam);
  sigt = command_par_value("sigt", current_beam);
  sige = command_par_value("sige", current_beam);
  beta = command_par_value("beta", current_beam);
  gamma = command_par_value("gamma", current_beam);
  circ = command_par_value("circ", current_beam);

  /* assume oneturnmap and disp0 already computed (see pro_emit) */ 
  for (j = 0; j < 4; j++) ds += oneturnmat[4 + 6*j] * disp0[j];
  tmp = - beta * beta * ds / circ;
  freq0 = (clight * ten_m_6 * beta) / (circ * (one + tmp * delta_p));
  etas = beta * gamma * (one + delta_p);
  gammas = sqrt(one + etas * etas);
  betas = etas / gammas;
  tmp = - betas * betas * ds / circ;
  alfa = one / (gammas * gammas) + tmp;
  dtbyds = delta_p * tmp / betas;
  // deltat = circ * dtbyds;
  store_comm_par_value("freq0", freq0, probe_beam);
  store_comm_par_value("alfa", alfa, probe_beam);
  store_comm_par_value("beta", betas, probe_beam);
  store_comm_par_value("gamma", gammas, probe_beam);
  store_comm_par_value("dtbyds", dtbyds, probe_beam);
  store_comm_par_value("deltap", delta_p, probe_beam);
  slope = -rfc_slope();
  qs = sqrt(fabs((tmp * slope) / (twopi * betas)));

  if (qs != zero)
  {
    fact = (tmp * circ) / (twopi * qs);
    if (et > zero)
    {
      sigt = sqrt(fabs(et * fact));
      sige = sqrt(fabs(et / fact));
    }
    else if (sigt > zero)
    {
      sige = sigt / fact;
      et = sige * sigt;
    }
    else if (sige > zero)
    {
      sigt = sige * fact;
      et = sige * sigt;
    }
  }
  if (sigt < ten_m_15)
  {
    put_info("Zero value of SIGT", "replaced by 1.");
    sigt = one;
  }
  if (sige < ten_m_15)
  {
    put_info("Zero value of SIGE", "replaced by 1/1000.");
    sigt = ten_m_3;
  }
  store_comm_par_value("qs", qs, probe_beam);
  store_comm_par_value("et", et, probe_beam);
  store_comm_par_value("sigt", sigt, probe_beam);
  store_comm_par_value("sige", sige, probe_beam);
}

