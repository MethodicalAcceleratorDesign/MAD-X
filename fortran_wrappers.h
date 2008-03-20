#ifndef _FORTRAN_WRAPPERS_H
#define _FORTRAN_WRAPPERS_H
/* redirect FORTRAN calls to wrappers that synchronize FORTRAN and C stdout buffering */
/* when crossing the border upon calling FORTRAN from C. */
#include "fortran_wrappers_prototypes.h"
#define dynap_ dynap_wrapper
#define mtgetc_ mtgetc_wrapper
#define collect_ collect_wrapper
#define mtsvd_ mtsvd_wrapper
#define setup_ setup_wrapper
#define setupi_ setupi_wrapper
#define micit_ micit_wrapper
#define haveit_ haveit_wrapper
#define svddec_m_ svddec_m_wrapper
#define svddec_c_ svddec_c_wrapper
#define svdcorr_m_ svdcorr_m_wrapper
#define svdcorr_c_ svdcorr_c_wrapper
#define primat_ primat_wrapper
#define prdmat_ prdmat_wrapper
#define pefill_ pefill_wrapper
#define pesopt_ pesopt_wrapper
#define plotit_ plotit_wrapper
#define w_ptc_create_universe_ w_ptc_create_universe_wrapper
#define w_ptc_create_layout_ w_ptc_create_layout_wrapper
#define w_ptc_move_to_layout_ w_ptc_move_to_layout_wrapper
#define w_ptc_input_ w_ptc_input_wrapper
#define w_ptc_align_ w_ptc_align_wrapper
#define w_ptc_track_ w_ptc_track_wrapper
#define w_ptc_twiss_ w_ptc_twiss_wrapper
#define w_ptc_normal_ w_ptc_normal_wrapper
#define w_ptc_moments_ w_ptc_moments_wrapper
#define w_ptc_initmoments_ w_ptc_initmoments_wrapper
#define w_ptc_trackline_ w_ptc_trackline_wrapper
#define w_ptc_track_everystep_ w_ptc_track_everystep_wrapper
#define w_ptc_addmoment_ w_ptc_addmoment_wrapper
#define w_ptc_getnfieldcomp_ w_ptc_getnfieldcomp_wrapper
#define w_ptc_getsfieldcomp_ w_ptc_getsfieldcomp_wrapper
#define w_ptc_setfieldcomp_ w_ptc_setfieldcomp_wrapper
#define w_ptc_eplacement_ w_ptc_eplacement_wrapper
#define w_ptc_end_ w_ptc_end_wrapper
#define w_ptc_rviewer_ w_ptc_rviewer_wrapper
#define res_index_ res_index_wrapper
#define soddin_ soddin_wrapper
#define twiss_ twiss_wrapper
#define tmrefe_ tmrefe_wrapper
#define tmrefo_ tmrefo_wrapper
#define getclor_ getclor_wrapper
#define timest_ timest_wrapper
#define timex_ timex_wrapper
#define w_ptc_dumpmaps_ w_ptc_dumpmaps_wrapper
#define w_ptc_setdebuglevel_ w_ptc_setdebuglevel_wrapper
#define w_ptc_setaccel_method_ w_ptc_setaccel_method_wrapper
#define w_ptc_setexactmis_ w_ptc_setexactmis_wrapper
#define w_ptc_setradiation_ w_ptc_setradiation_wrapper
#define w_ptc_settotalpath_ w_ptc_settotalpath_wrapper
#define w_ptc_settime_ w_ptc_settime_wrapper
#define w_ptc_setnocavity_ w_ptc_setnocavity_wrapper
#define w_ptc_setfringe_ w_ptc_setfringe_wrapper
#define w_ptc_setknobvalue_ w_ptc_setknobvalue_wrapper
#define w_ptc_refreshtables_ w_ptc_refreshtables_wrapper
#define w_ptc_addknob_ w_ptc_addknob_wrapper
#define w_ptc_addknob_i_ w_ptc_addknob_i_wrapper
#define w_ptc_writeparresults_ w_ptc_writeparresults_wrapper
#define w_ptc_printframes_ w_ptc_printframes_wrapper
#define w_ptc_printlayout_rootm_ w_ptc_printlayout_rootm_wrapper
#define w_ptc_addpush_ w_ptc_addpush_wrapper
#define w_ptc_script_ w_ptc_script_wrapper
#define w_ptc_open_gino_ w_ptc_open_gino_wrapper

#endif
