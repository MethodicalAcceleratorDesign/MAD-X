/* to avoid warnings of implicit declaration from code calling the functions below */
#ifndef _FORTRAN_WRAPPERS_PROTOTYPES_H
#define _FORTRAN_WRAPPERS_PROTOTYPES_H
/* Wrap 'dynap' defined in 'dynap.F'*/
void dynap_wrapper(double* eigen,double* coords,int* turns,int* npart,double* distvect,double* zn,double* dq,double* onelog,double* turnnumber);
/* Wrap 'mtgetc' defined in 'match.F'*/
void mtgetc_wrapper(double* vect,double* dvect);
/* Wrap 'collect' defined in 'match.F'*/
void collect_wrapper(int* ncon,double* fsum,double* fvect);
/* Wrap 'mtsvd' defined in 'matchjc.F'*/
void mtsvd_wrapper(int* M,int* N,double* fjac,double* SV,double* U,double* VT);
/* Wrap 'setup' defined in 'orbf.F'*/
void setup_wrapper(double* resp,double* a,int* im,int* ic,int* nm,int* nc);
/* Wrap 'setupi' defined in 'orbf.F'*/
void setupi_wrapper(int* resp,int* a,int* im,int* ic,int* nm,int* nc);
/* Wrap 'micit' defined in 'orbf.F'*/
void micit_wrapper(double* a,char* conm,double* xin,double* cin,double* res,int* nx,float* rms,int* im,int* ic,int* iter,int* ny,float* ax,float* cinx,float* xinx,float* resx,float* rho,float* ptop,float* rmss,float* xrms,float* xptp,float* xiter,int* ifail);
/* Wrap 'haveit' defined in 'orbf.F'*/
void haveit_wrapper(double* a,double* xin,double* cin,double* res,int* nx,int* im,int* ic,double* cb,double* xmeas,double* xres,double* y,double* z,double* xd);
/* Wrap 'svddec_m' defined in 'orbf.F'*/
void svddec_m_wrapper(double* a,double* svdmat,double* umat,double* vmat,double* wmat,double* utmat,double* vtmat,double* wtmat,double* ws,double* wvec,int* sortw,int* im,int* ic,int* iflag,int* sing,int* dbg);
/* Wrap 'svddec_c' defined in 'orbf.F'*/
void svddec_c_wrapper(double* a,double* svdmat,double* umat,double* vmat,double* wmat,double* utmat,double* vtmat,double* wtmat,double* ws,double* wvec,int* sortw,int* im,int* ic,int* iflag,int* sing,int* dbg);
/* Wrap 'svdcorr_m' defined in 'orbf.F'*/
void svdcorr_m_wrapper(double* a,double* svdmat,double* umat,double* vmat,double* wmat,double* utmat,double* vtmat,double* wtmat,double* xin,double* xc,double* xout,double* xa,double* xb,double* xpred,double* ws,double* wvec,int* sortw,int* nx,int* im,int* ic,int* iflag,int* dbg);
/* Wrap 'svdcorr_c' defined in 'orbf.F'*/
void svdcorr_c_wrapper(double* a,double* svdmat,double* umat,double* vmat,double* wmat,double* utmat,double* vtmat,double* wtmat,double* xin,double* xc,double* xout,double* xa,double* xb,double* xpred,double* ws,double* wvec,int* sortw,int* nx,int* im,int* ic,int* iflag,int* dbg);
/* Wrap 'primat' defined in 'orbf.F'*/
void primat_wrapper(int* a,int* nc,int* nm);
/* Wrap 'prdmat' defined in 'orbf.F'*/
void prdmat_wrapper(double* a,int* nc,int* nm);
/* Wrap 'pefill' defined in 'plot.F'*/
void pefill_wrapper(int* ierr);
/* Wrap 'pesopt' defined in 'plot.F'*/
void pesopt_wrapper(int* ierr);
/* Wrap 'plotit' defined in 'plot.F'*/
void plotit_wrapper(int* initfl);
/* Wrap 'w_ptc_create_universe' defined in 'wrap.f90'*/
void w_ptc_create_universe_wrapper();
/* Wrap 'w_ptc_create_layout' defined in 'wrap.f90'*/
void w_ptc_create_layout_wrapper();
/* Wrap 'w_ptc_move_to_layout' defined in 'wrap.f90'*/
void w_ptc_move_to_layout_wrapper();
/* Wrap 'w_ptc_input' defined in 'wrap.f90'*/
void w_ptc_input_wrapper();
/* Wrap 'w_ptc_align' defined in 'wrap.f90'*/
void w_ptc_align_wrapper();
/* Wrap 'w_ptc_track' defined in 'wrap.f90'*/
void w_ptc_track_wrapper(int* max_obs);
/* Wrap 'w_ptc_twiss' defined in 'wrap.f90'*/
void w_ptc_twiss_wrapper(int* tab_name);
/* Wrap 'w_ptc_normal' defined in 'wrap.f90'*/
void w_ptc_normal_wrapper();
/* Wrap 'w_ptc_moments' defined in 'wrap.f90'*/
void w_ptc_moments_wrapper(int* no);
/* Wrap 'w_ptc_initmoments' defined in 'wrap.f90'*/
void w_ptc_initmoments_wrapper();
/* Wrap 'w_ptc_trackline' defined in 'wrap.f90'*/
void w_ptc_trackline_wrapper(int* max_obs);
/* Wrap 'w_ptc_track_everystep' defined in 'wrap.f90'*/
void w_ptc_track_everystep_wrapper(int* max_obs);
/* Wrap 'w_ptc_addmoment' defined in 'wrap.f90'*/
void w_ptc_addmoment_wrapper(int* x,int* px,int* y,int* py,int* t,int* dp,int* tableIA,int* columnIA,int* parametric);
/* Wrap 'w_ptc_getnfieldcomp' defined in 'wrap.f90'*/
void w_ptc_getnfieldcomp_wrapper(int* fibreidx,int* ncomp,double* nval);
/* Wrap 'w_ptc_getsfieldcomp' defined in 'wrap.f90'*/
void w_ptc_getsfieldcomp_wrapper(int* fibreidx,int* ncomp,double* nval);
/* Wrap 'w_ptc_setfieldcomp' defined in 'wrap.f90'*/
void w_ptc_setfieldcomp_wrapper(int* fibreidx);
/* Wrap 'w_ptc_eplacement' defined in 'wrap.f90'*/
void w_ptc_eplacement_wrapper(int* elementidx,int* rf);
/* Wrap 'w_ptc_end' defined in 'wrap.f90'*/
void w_ptc_end_wrapper();
/* Wrap 'w_ptc_rviewer' defined in 'madx_ptc_knobs.f90'*/
void w_ptc_rviewer_wrapper();
/* Wrap 'res_index' defined in 'resindex.F'*/
void res_index_wrapper(int* skew,int* mynorder,int* myn1,int* myn2,int* indexa,int* mynres);
/* Wrap 'soddin' defined in 'sodd.F'*/
void soddin_wrapper(int* ierr);
/* Wrap 'twiss' defined in 'twiss.F'*/
void twiss_wrapper(double* rt,double* disp0,int* tab_name);
/* Wrap 'tmrefe' defined in 'twiss.F'*/
void tmrefe_wrapper(double* rt);
/* Wrap 'tmrefo' defined in 'twiss.F'*/
void tmrefo_wrapper(int* kobs,double* orbit0,double* orbit,double* rt);
/* Wrap 'getclor' defined in 'util.F'*/
void getclor_wrapper(double* orbit0,double* rt,double* tt,int* error);
/* Wrap 'timest' defined in 'timest.f90'*/
void timest_wrapper(float* r1);
/* Wrap 'timex' defined in 'timex.f90'*/
void timex_wrapper(float* r1);
/* Wrap 'w_ptc_dumpmaps' defined in 'wrap.f90'*/
void w_ptc_dumpmaps_wrapper();
/* Wrap 'w_ptc_setdebuglevel' defined in 'wrap.f90'*/
void w_ptc_setdebuglevel_wrapper(int* level);
/* Wrap 'w_ptc_setaccel_method' defined in 'wrap.f90'*/
void w_ptc_setaccel_method_wrapper(int* method);
/* Wrap 'w_ptc_setexactmis' defined in 'wrap.f90'*/
void w_ptc_setexactmis_wrapper(int* method);
/* Wrap 'w_ptc_setradiation' defined in 'wrap.f90'*/
void w_ptc_setradiation_wrapper(int* method);
/* Wrap 'w_ptc_settotalpath' defined in 'wrap.f90'*/
void w_ptc_settotalpath_wrapper(int* method);
/* Wrap 'w_ptc_settime' defined in 'wrap.f90'*/
void w_ptc_settime_wrapper(int* method);
/* Wrap 'w_ptc_setnocavity' defined in 'wrap.f90'*/
void w_ptc_setnocavity_wrapper(int* method);
/* Wrap 'w_ptc_setfringe' defined in 'wrap.f90'*/
void w_ptc_setfringe_wrapper(int* method);
/* Wrap 'w_ptc_setknobvalue' defined in 'wrap.f90'*/
void w_ptc_setknobvalue_wrapper(int* fibre);
/* Wrap 'w_ptc_refreshtables' defined in 'wrap.f90'*/
void w_ptc_refreshtables_wrapper();
/* Wrap 'w_ptc_addknob' defined in 'wrap.f90'*/
void w_ptc_addknob_wrapper(int* fibre);
/* Wrap 'w_ptc_addknob_i' defined in 'wrap.f90'*/
void w_ptc_addknob_i_wrapper(int* paramn);
/* Wrap 'w_ptc_writeparresults' defined in 'wrap.f90'*/
void w_ptc_writeparresults_wrapper(int* filename);
/* Wrap 'w_ptc_printframes' defined in 'wrap.f90'*/
void w_ptc_printframes_wrapper(int* filename);
/* Wrap 'w_ptc_printlayout_rootm' defined in 'wrap.f90'*/
void w_ptc_printlayout_rootm_wrapper(int* filename);
/* Wrap 'w_ptc_addpush' defined in 'wrap.f90'*/
void w_ptc_addpush_wrapper(int* tabname,int* colname,int* polinomial,int* monomial);
/* Wrap 'w_ptc_script' defined in 'wrap.f90'*/
void w_ptc_script_wrapper(int* scriptname);
/* Wrap 'w_ptc_open_gino' defined in 'wrap.f90'*/
void w_ptc_open_gino_wrapper(int* scriptname);
#endif
