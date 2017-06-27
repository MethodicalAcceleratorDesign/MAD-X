#ifndef MAD_EXTRN_F_H
#define MAD_EXTRN_F_H

/*
 * Fortran types in C
 */

#ifndef MAD_TYPES_F_H
#include "mad_types_f.h"
#endif

/*
 * Called by Fortran
 */

// from mad_cmdpar.c
#define comm_para comm_para_
#define get_string get_string_
#define get_value get_value_ // **
#define get_vector get_vector_

// from mad_const.c
#define next_constraint next_constraint_
#define next_constr_namepos next_constr_namepos_
#define next_global next_global_

// from mad_bbeam.c
#define make_bb6d_ixy make_bb6d_ixy_

// from mad_elem.c
#define element_name element_name_
#define el_par_vector el_par_vector_
#define get_node_vector get_node_vector_

// from mad_elemerr.c
#define node_al_errors node_al_errors_
#define node_fd_errors node_fd_errors_

// from mad_err.c
#define augmentfwarn augmentfwarn_
#define geterrorflag geterrorflag_
#define seterrorflagfort seterrorflagfort_

// from mad_eval.c
#define pro_input pro_input_

// from mad_gxx11c.c
#define wopen    wopen_
#define wclose   wclose_
#define wclrwk   wclrwk_
#define wpl      wpl_
#define wfa      wfa_
#define wswn     wswn_
#define wtx      wtx_
#define wwait    wwait_
#define wsetci   wsetci_
#define wsetls   wsetls_
#define wstring  wstring_
#define cbyt     cbyt_
#define mydtime  mydtime_

// from mad_inter.c
#define interpolate_node interpolate_node_
#define reset_interpolation reset_interpolation_
#define start_interp_node start_interp_node_
#define fetch_interp_node fetch_interp_node_

// from mad_match.c
#define mtcond mtcond_
#define mtputconsname mtputconsname_

// from mad_node.c
#define advance_node advance_node_ // *
#define advance_to_pos advance_to_pos_
#define current_node_name current_node_name_
#define node_name node_name_
#define node_string node_string_
#define node_value node_value_ // **
#define retreat_node retreat_node_
#define store_node_value store_node_value_
#define store_node_vector store_node_vector_
#define store_no_fd_err store_no_fd_err_

// from mad_option.c
#define get_option get_option_ // *
#define set_option set_option_
#define set_cont_sequence set_cont_sequence_
#define set_sequence set_sequence_

// from mad_orbit.c
#define f_ctof f_ctof_

// from mad_plot.c
#define get_title get_title_
#define get_version get_version_
#define plot_option plot_option_

// from mad_ptc.c
#define augmentcountmomtabs augmentcountmomtabs_
#define makemomentstables makemomentstables_
#define minimum_acceptable_order minimum_acceptable_order_
#define select_ptc_idx select_ptc_idx_

// from mad_rand.c
#define frndm frndm_
#define grndm grndm_

// from mad_rplot.cpp
#define newrplot newrplot_
#define plottrack plottrack_
#define plottwiss plottwiss_
#define rplotfinish rplotfinish_
#define rviewer rviewer_
#define madxv_setfctnname madxv_setfctnname_
#define madxv_setknobname madxv_setknobname_
#define madxv_setfunctionat madxv_setfunctionat_

// from mad_seq.c
#define restart_sequ restart_sequ_

// from mad_table.c
// warning:augment_counts is provided by madx_ptc_knobs.f90
#define augment_count augment_count_
#define augmentcountonly augmentcountonly_
#define table_length table_length_
#define table_exists table_exists_
#define table_cell_exists table_cell_exists_
#define table_column_exists table_column_exists_
#define table_header_exists table_header_exists_
#define double_from_table_header double_from_table_header_
#define double_from_table_row double_from_table_row_ // *
#define string_from_table_row string_from_table_row_
#define double_to_table_row double_to_table_row_
#define string_to_table_row string_to_table_row_
#define double_to_table_curr double_to_table_curr_ // **
#define vector_to_table_curr vector_to_table_curr_
#define string_to_table_curr string_to_table_curr_
#define comment_to_table_curr comment_to_table_curr_
#define double_table double_table_
#define make_map_table make_map_table_
#define reset_count reset_count_
#define sector_out sector_out_
#define table_length table_length_
#define table_range table_range_

// from mad_track.c
#define deletetrackstrarpositions deletetrackstrarpositions_
#define getcurrentcmdname getcurrentcmdname_
#define getnumberoftracks getnumberoftracks_
#define gettrack gettrack_
#define next_start next_start_
#define track_pteigen track_pteigen_

// from mad_twiss.c
#define copy_twiss_data copy_twiss_data_
#define embedded_twiss embedded_twiss_
#define get_disp0 get_disp0_

// from mad_util.c
#define intrac intrac_

// from mad_var.c
#define set_variable set_variable_
#define get_variable get_variable_
#define next_vary next_vary_

// from mad_vec.c
#define vdot vdot_ // *
#define vmod vmod_

/*
 * Provided by Fortran
 */
 
// from mad_init_f.F90
void mad_init_f_(void);

// from dynap.f90
//void dynap_(F_DOUBLE eigen, F_DOUBLE coords, F_INTEGER turns, F_INTEGER npart, F_DOUBLE distvect,
//	    F_DOUBLE zn, F_DOUBLE dq, F_DOUBLE onelog, F_DOUBLE turnnumber);
void trdynrun_(F_DOUBLE eigen, F_DOUBLE coords, F_INTEGER turns, F_INTEGER npart, F_DOUBLE distvect,
	       F_DOUBLE zn, F_DOUBLE onelog, F_DOUBLE turnnumber, F_DOUBLE dq);

// from emit.f90
void emit_(F_DOUBLE deltap, F_DOUBLE tol, F_DOUBLE orbit0, F_DOUBLE disp0, F_DOUBLE rt,
	   F_DOUBLE u0, F_DOUBLE emit_v, F_DOUBLE nemit_v, F_DOUBLE bmax, F_DOUBLE gmax,
	   F_DOUBLE dismax, F_DOUBLE tunes, F_DOUBLE sig_v, F_DOUBLE pdamp, F_LOGICAL updatebeam);

// from gxx11.f90(Unix)
// or gxx11ps.f90(Windows)
void gxterm_(void);

// from ibsdb.f90
void ibs_(void);

// from match.f90
void collect_(F_INTEGER ncon, F_DOUBLE fsum, F_DOUBLE fvect);
void mtlmdf_(F_INTEGER ncon, F_INTEGER nvar, F_DOUBLE tol, F_INTEGER calls,
	     F_INTEGER call_lim, F_DOUBLE vect, F_DOUBLE dvect, F_DOUBLE fun_vec,
	     F_DOUBLE diag, F_DOUBLE w_ifjac, F_DOUBLE w_ipvt, F_DOUBLE w_qtf, 
	     F_DOUBLE w_iwa1, F_DOUBLE w_iwa2, F_DOUBLE w_iwa3, F_DOUBLE w_iwa4, 
	     F_DOUBLE xold);
void mtmigr_(F_INTEGER ncon, F_INTEGER nvar, F_INTEGER strategy, F_DOUBLE tol,
	     F_INTEGER calls, F_INTEGER call_lim, F_DOUBLE vect, F_DOUBLE dvect,
	     F_DOUBLE fun_vect, F_DOUBLE w_iwa1, F_DOUBLE w_iwa2, F_DOUBLE w_iwa3, 
	     F_DOUBLE w_iwa4, F_DOUBLE w_iwa5, F_DOUBLE w_iwa6, F_DOUBLE w_iwa7, F_DOUBLE w_iwa8);
void mtsimp_(F_INTEGER ncon, F_INTEGER nvar, F_DOUBLE tol, F_INTEGER calls, F_INTEGER call_lim,
	     F_DOUBLE vect, F_DOUBLE dvect, F_DOUBLE fun_vect, F_DOUBLE w_iwa1, F_DOUBLE w_iwa2,
	     F_DOUBLE w_iwa3);

// from matchjc.f90
void mtsvd_(F_INTEGER M, F_INTEGER N, F_DOUBLE fjac, F_DOUBLE SV, F_DOUBLE U, F_DOUBLE VT);
void mtjac_(F_INTEGER ncon, F_INTEGER nvar, F_INTEGER strategy, F_DOUBLE cool, F_DOUBLE balance, 
	    F_DOUBLE random, F_INTEGER nrep, F_INTEGER bisec, F_DOUBLE cond, F_INTEGER match_mode,
	    F_DOUBLE tol, F_INTEGER calls, F_INTEGER call_lim, F_DOUBLE vect, F_DOUBLE dvect,
	    F_DOUBLE fun_vec, F_DOUBLE w_ifjac, F_DOUBLE w_iwa4, F_DOUBLE fval, F_DOUBLE xstart,
	    F_DOUBLE xold);

// from matchsa.f90
void mtsa_(F_INTEGER ncon, F_INTEGER nvar, F_DOUBLE tol, F_INTEGER calls, F_INTEGER call_lim,
	   F_DOUBLE vect, F_DOUBLE fun_vect, F_INTEGER iseed, F_INTEGER iprint, F_DOUBLE lb, 
	   F_INTEGER nacp, F_DOUBLE ub, F_DOUBLE xopt, F_DOUBLE c, F_DOUBLE vm, F_DOUBLE xp);

// from orbf.f90
void setup_(F_DOUBLE resp, F_DOUBLE a, F_INTEGER im, F_INTEGER ic, F_INTEGER nm, F_INTEGER nc);
void setupi_(F_INTEGER resp, F_INTEGER a, F_INTEGER im, F_INTEGER ic, F_INTEGER nm, F_INTEGER nc);
void prdmat_(F_DOUBLE a, F_INTEGER nc, F_INTEGER nm);
void primat_(F_INTEGER a, F_INTEGER nc, F_INTEGER nm);
void micit_(F_DOUBLE a, F_CHARACTER conm, F_DOUBLE xin, F_DOUBLE cin, F_DOUBLE res, F_INTEGER nx,
	    F_DOUBLE rms, F_INTEGER im, F_INTEGER ic, F_INTEGER iter, F_INTEGER ny, F_DOUBLE ax,
	    F_DOUBLE cinx, F_DOUBLE xinx, F_DOUBLE resx, F_DOUBLE rho, F_DOUBLE ptop, F_DOUBLE rmss, F_DOUBLE xrms,
	    F_DOUBLE xptp, F_DOUBLE xiter, F_INTEGER ifail);
void haveit_(F_DOUBLE a, F_DOUBLE xin, F_DOUBLE cin, F_DOUBLE res, F_INTEGER nx, F_INTEGER im,
	     F_INTEGER ic, F_DOUBLE cb, F_DOUBLE xmeas, F_DOUBLE xres, F_DOUBLE y, F_DOUBLE z,
	     F_DOUBLE xd);
void svddec_(F_DOUBLE a, F_DOUBLE svdmat, F_DOUBLE umat, F_DOUBLE vmat, F_DOUBLE ws, F_DOUBLE wvec, 
	     F_INTEGER sortw, F_DOUBLE sngcut, F_DOUBLE sngval, F_INTEGER im, F_INTEGER ic, 
	     F_INTEGER iflag, F_INTEGER sing, F_INTEGER dbg);
void svdcorr_(F_DOUBLE a, F_DOUBLE svdmat, F_DOUBLE umat, F_DOUBLE vmat, F_DOUBLE wmat, F_DOUBLE utmat,
	      F_DOUBLE vtmat, F_DOUBLE wtmat, F_DOUBLE xin, F_DOUBLE xc, F_DOUBLE xout, 
	      F_DOUBLE xpred, F_DOUBLE ws, F_DOUBLE wvec, F_INTEGER sortw, F_INTEGER nx,
	      F_INTEGER im, F_INTEGER ic, F_INTEGER iflag, F_INTEGER dbg);

// from plot.f90
void pesopt_(F_INTEGER ierr);
void pefill_(F_INTEGER ierr);
void pemima_(void);
void plotit_(F_INTEGER initfl);

// from resindex.f90
void res_index_(F_LOGICAL skew, F_INTEGER mynorder, F_INTEGER myn1, F_INTEGER myn2, F_INTEGER indexa,
		F_INTEGER mynres);

// from sodd.f90
void soddin_(F_INTEGER ierr);

// from survey.f90
void survey_(void);
void survtest_(void);

// from toucheck.f90
void touschek_(void);

// from trrun.f90
void trrun_(F_INTEGER switch_, F_INTEGER turns, F_DOUBLE orbit0, F_DOUBLE rt, F_INTEGER part_id,
	    F_INTEGER last_turn, F_DOUBLE last_pos, F_DOUBLE z, F_DOUBLE dxt, F_DOUBLE dyt, 
	    F_DOUBLE last_orbit, F_DOUBLE eigen, F_DOUBLE coords, F_INTEGER e_flag,
	    F_INTEGER code_buf, F_DOUBLE l_buf);

// from twiss.f90
void tmrefe_(F_DOUBLE rf);
void tmrefo_(F_INTEGER kobs, F_DOUBLE orbit0, F_DOUBLE orbit, F_DOUBLE rt);
void twiss_(F_DOUBLE rt, F_DOUBLE disp0, F_INTEGER tab_name, F_INTEGER sector_tab_name);
void twcpin_(F_DOUBLE rt, F_DOUBLE disp0, F_DOUBLE r0mat, F_INTEGER error);
void twdisp_ini_(F_DOUBLE rt, F_DOUBLE disp0);

// from util.f90
void getclor_(F_DOUBLE orbit0, F_DOUBLE rt, F_DOUBLE tt, F_INTEGER error);

// from wrap.f90
void w_ptc_addknob_(F_INTEGER fibre);
void w_ptc_addknob_i_(F_INTEGER paramn);
void w_ptc_addmoment_(F_INTEGER x, F_INTEGER px, F_INTEGER y, F_INTEGER py, F_INTEGER t, F_INTEGER dp,
		      F_INTEGER tableIA, F_INTEGER columnIA, F_INTEGER parametric);
void w_ptc_addpush_(F_INTEGER tabname, F_INTEGER colname, F_INTEGER polinomial, F_INTEGER monomial);
void w_ptc_align_(void);
void w_ptc_putbeambeam_(void);
void w_ptc_create_layout_(void);
void w_ptc_create_universe_(void);
void w_ptc_dumpmaps_(void);
void w_ptc_end_(void);
void w_ptc_enforce6d_(F_INTEGER level);
void w_ptc_eplacement_(F_INTEGER elementidx, F_INTEGER rf);
void w_ptc_export_xml_(F_INTEGER filename);
void w_ptc_getnfieldcomp_(F_INTEGER fibreidx, F_INTEGER ncomp, F_DOUBLE nval);
void w_ptc_getsfieldcomp_(F_INTEGER fibreidx, F_INTEGER ncomp, F_DOUBLE nval);
void w_ptc_moments_(F_INTEGER no);
void w_ptc_move_to_layout_(void);
void w_ptc_normal_(void);
void w_ptc_open_gino_(F_INTEGER scriptname);
void w_ptc_printframes_(F_INTEGER scriptname);
void w_ptc_printlayout_rootm_(F_INTEGER filename);
void w_ptc_read_errors_(void);
void w_ptc_refresh_k_(void);
void w_ptc_refreshtables_(void);
void w_ptc_script_(F_INTEGER scriptname);
void w_ptc_setaccel_method_(F_INTEGER method);
void w_ptc_setdebuglevel_(F_INTEGER level);
void w_ptc_setexactmis_(F_INTEGER method);
void w_ptc_setfieldcomp_(F_INTEGER fibreidx);
void w_ptc_setfringe_(F_INTEGER method);
void w_ptc_setknobvalue_(F_INTEGER fible);
void w_ptc_setnocavity_(F_INTEGER method);
void w_ptc_setradiation_(F_INTEGER method);
void w_ptc_settime_(F_INTEGER method);
void w_ptc_settotalpath_(F_INTEGER method);
void w_ptc_track_(F_INTEGER max_obs);
void w_ptc_track_everystep_(F_INTEGER max_obs);
void w_ptc_trackline_(F_INTEGER max_obs);
void w_ptc_twiss_(F_INTEGER tab_name, F_INTEGER summary_name);
void w_ptc_writeparresults_(F_INTEGER filename);

// from madx_ptc_distrib.f90
void w_ptc_getmomentstabcol_(F_INTEGER n, F_CHARACTER tabn, F_CHARACTER coln);
int  w_ptc_getnmoments_(void);

// from madx_ptc_knobs.f90
void w_ptc_rviewer_(void);

#endif
