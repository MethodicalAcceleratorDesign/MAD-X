#ifndef MAD_EXTRN_F_H
#define MAD_EXTRN_F_H

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

// from mad_option.c
#define get_option get_option_ // *
#define set_option set_option_

// from mad_orbit.c
#define f_ctof f_ctof_

// from mad_plot.c
#define get_title get_title_
#define get_version get_version_
#define interp_node interp_node_
#define plot_option plot_option_
#define reset_interpolation reset_interpolation_

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
// warning: augment_counts is provided by madx_ptc_knobs.f90
#define augment_count augment_count_
#define augmentcountonly augmentcountonly_
#define char_from_table char_from_table_
#define comment_to_table comment_to_table_
#define double_from_table double_from_table_ // *
#define double_table double_table_
#define double_to_table double_to_table_ // **
#define double_to_table_row double_to_table_row_
#define headvalue headvalue_
#define make_map_table make_map_table_
#define reset_count reset_count_
#define sector_out sector_out_
#define string_from_table string_from_table_
#define string_to_table string_to_table_
#define table_length table_length_
#define table_range table_range_
#define vector_to_table vector_to_table_

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
#define get_twiss_data get_twiss_data_

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
void dynap_();

// from emit.f90
void emit_();

// from gxx11.f90   (Unix)
// or   gxx11ps.f90 (Windows)
void gxterm_(void);

// from ibsdb.f90
void ibs_(void);

// from match.f90
void collect_();
void mtlmdf_();
void mtmigr_();
void mtsimp_();

// from matchjc.f90
void mtsvd_();
void mtjac_();

// from matchsa.f90
void mtsa_();

// from orbif.f90
void setup_();
void micit_();
void haveit_();
void svddec_m_();
void svddec_c_();
void svdcorr_m_();
void svdcorr_c_();

// from plot.f90
void pesopt_();
void pefill_();
void pemima_(void);
void plotit_();

// from resindex.f90
void res_index_();

// from sodd.f90
void soddin_();

// from survey.f90
void survey_(void);

// from toucheck.f90
void touschek_(void);

// from trrun.f90
void trrun_();

// from twiss.f90
void tmrefe_();
void tmrefo_();
void twiss_();

// from util.f90
void getclor_();

// from wrap.f90
void w_ptc_addknob_();
void w_ptc_addknob_i_();
void w_ptc_addmoment_();
void w_ptc_addpush_();
void w_ptc_align_(void);
void w_ptc_create_layout_(void);
void w_ptc_create_universe_(void);
void w_ptc_dumpmaps_(void);
void w_ptc_end_(void);
void w_ptc_enforce6d_();
void w_ptc_eplacement_();
void w_ptc_export_xml_();
void w_ptc_getnfieldcomp_();
int  w_ptc_getnmoments_();
void w_ptc_getsfieldcomp_();
void w_ptc_moments_();
void w_ptc_move_to_layout_(void);
void w_ptc_normal_(void);
void w_ptc_open_gino_();
void w_ptc_printframes_();
void w_ptc_printlayout_rootm_();
void w_ptc_read_errors_(void);
void w_ptc_refresh_k_(void);
void w_ptc_refreshtables_(void);
void w_ptc_script_();
void w_ptc_setaccel_method_();
void w_ptc_setdebuglevel_();
void w_ptc_setexactmis_();
void w_ptc_setfieldcomp_();
void w_ptc_setfringe_();
void w_ptc_setknobvalue_();
void w_ptc_setnocavity_();
void w_ptc_setradiation_();
void w_ptc_settime_();
void w_ptc_settotalpath_();
void w_ptc_track_();
void w_ptc_track_everystep_();
void w_ptc_trackline_();
void w_ptc_twiss_();
void w_ptc_writeparresults_();

// from madx_ptc_distrib.f90
void w_ptc_getmomentstabcol_();

// from madx_ptc_knobs.f90
void w_ptc_rviewer_(void);

/*
 * Windows specific (should be removed)
 */

// cancel type_ofCall
#define type_ofCall

#endif
