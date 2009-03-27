/* to avoid warnings of implicit declarations from fortran_wrappers.c */
#ifndef _FORTRAN_PROTOTYPES_H
#define _FORTRAN_PROTOTYPES_H
/* Wrap 'dynap' defined in 'dynap.f90'*/
void dynap_(double* eigen,double* coords,int* turns,int* npart,double* distvect,double* zn,double* dq,double* onelog,double* turnnumber);
/* Wrap 'w_ptc_rviewer' defined in 'madx_ptc_knobs.f90'*/
void w_ptc_rviewer_();
/* Wrap 'mtgetc' defined in 'match.f90'*/
void mtgetc_(double* vect,double* dvect);
/* Wrap 'collect' defined in 'match.f90'*/
void collect_(int* ncon,double* fsum,double* fvect);
/* Wrap 'mtsvd' defined in 'matchjc.f90'*/
void mtsvd_(int* M,int* N,double* fjac,double* SV,double* U,double* VT);
/* Wrap 'setup' defined in 'orbf.f90'*/
void setup_(double* resp,double* a,int* im,int* ic,int* nm,int* nc);
/* Wrap 'setupi' defined in 'orbf.f90'*/
void setupi_(int* resp,int* a,int* im,int* ic,int* nm,int* nc);
/* Wrap 'micit' defined in 'orbf.f90'*/
void micit_(double* a,char* conm,double* xin,double* cin,double* res,int* nx,float* rms,int* im,int* ic,int* iter,int* ny,float* ax,float* cinx,float* xinx,float* resx,float* rho,float* ptop,float* rmss,float* xrms,float* xptp,float* xiter,int* ifail);
/* Wrap 'haveit' defined in 'orbf.f90'*/
void haveit_(double* a,double* xin,double* cin,double* res,int* nx,int* im,int* ic,double* cb,double* xmeas,double* xres,double* y,double* z,double* xd);
/* Wrap 'svddec_m' defined in 'orbf.f90'*/
void svddec_m_(double* a,double* svdmat,double* umat,double* vmat,double* wmat,double* utmat,double* vtmat,double* wtmat,double* ws,double* wvec,int* sortw,double* sngcut,double* sngval,int* im,int* ic,int* iflag,int* sing,int* dbg);
/* Wrap 'svddec_c' defined in 'orbf.f90'*/
void svddec_c_(double* a,double* svdmat,double* umat,double* vmat,double* wmat,double* utmat,double* vtmat,double* wtmat,double* ws,double* wvec,int* sortw,double* sngcut,double* sngval,int* im,int* ic,int* iflag,int* sing,int* dbg);
/* Wrap 'svdcorr_m' defined in 'orbf.f90'*/
void svdcorr_m_(double* a,double* svdmat,double* umat,double* vmat,double* wmat,double* utmat,double* vtmat,double* wtmat,double* xin,double* xc,double* xout,double* xa,double* xb,double* xpred,double* ws,double* wvec,int* sortw,int* nx,int* im,int* ic,int* iflag,int* dbg);
/* Wrap 'svdcorr_c' defined in 'orbf.f90'*/
void svdcorr_c_(double* a,double* svdmat,double* umat,double* vmat,double* wmat,double* utmat,double* vtmat,double* wtmat,double* xin,double* xc,double* xout,double* xa,double* xb,double* xpred,double* ws,double* wvec,int* sortw,int* nx,int* im,int* ic,int* iflag,int* dbg);
/* Wrap 'primat' defined in 'orbf.f90'*/
void primat_(int* a,int* nc,int* nm);
/* Wrap 'prdmat' defined in 'orbf.f90'*/
void prdmat_(double* a,int* nc,int* nm);
/* Wrap 'pefill' defined in 'plot.f90'*/
void pefill_(int* ierr);
/* Wrap 'pesopt' defined in 'plot.f90'*/
void pesopt_(int* ierr);
/* Wrap 'plotit' defined in 'plot.f90'*/
void plotit_(int* initfl);
/* Wrap 'res_index' defined in 'resindex.F'*/
void res_index_(int* skew,int* mynorder,int* myn1,int* myn2,int indexa[4][1000],int* mynres);
/* Wrap 'soddin' defined in 'sodd.f90'*/
void soddin_(int* ierr);
/* Wrap 'timest' defined in 'timest.f90'*/
void timest_(float* r1);
/* Wrap 'timex' defined in 'timex.f90'*/
void timex_(float* r1);
/* Wrap 'twiss' defined in 'twiss.f90'*/
void twiss_(double* rt,double* disp0,int* tab_name,int* sector_tab_name);
/* Wrap 'tmrefe' defined in 'twiss.f90'*/
void tmrefe_(double* rt);
/* Wrap 'tmrefo' defined in 'twiss.f90'*/
void tmrefo_(int* kobs,double* orbit0,double* orbit,double* rt);
/* Wrap 'getclor' defined in 'util.f90'*/
void getclor_(double* orbit0,double* rt,double* tt,int* error);
/* Wrap 'w_ptc_create_universe' defined in 'wrap.f90'*/
void w_ptc_create_universe_();
/* Wrap 'w_ptc_create_layout' defined in 'wrap.f90'*/
void w_ptc_create_layout_();
/* Wrap 'w_ptc_export_xml' defined in 'wrap.f90'*/
void w_ptc_export_xml_(int* filename);
/* Wrap 'w_ptc_move_to_layout' defined in 'wrap.f90'*/
void w_ptc_move_to_layout_();
/* Wrap 'w_ptc_input' defined in 'wrap.f90'*/
void w_ptc_input_();
/* Wrap 'w_ptc_align' defined in 'wrap.f90'*/
void w_ptc_align_();
/* Wrap 'w_ptc_twiss' defined in 'wrap.f90'*/
void w_ptc_twiss_(int* tab_name,int* summary_name);
/* Wrap 'w_ptc_normal' defined in 'wrap.f90'*/
void w_ptc_normal_();
/* Wrap 'w_ptc_moments' defined in 'wrap.f90'*/
void w_ptc_moments_(int* no);
/* Wrap 'w_ptc_initmoments' defined in 'wrap.f90'*/
void w_ptc_initmoments_();
/* Wrap 'w_ptc_dumpmaps' defined in 'wrap.f90'*/
void w_ptc_dumpmaps_();
/* Wrap 'w_ptc_track' defined in 'wrap.f90'*/
void w_ptc_track_(int* max_obs);
/* Wrap 'w_ptc_trackline' defined in 'wrap.f90'*/
void w_ptc_trackline_(int* max_obs);
/* Wrap 'w_ptc_track_everystep' defined in 'wrap.f90'*/
void w_ptc_track_everystep_(int* max_obs);
/* Wrap 'w_ptc_setdebuglevel' defined in 'wrap.f90'*/
void w_ptc_setdebuglevel_(int* level);
/* Wrap 'w_ptc_setaccel_method' defined in 'wrap.f90'*/
void w_ptc_setaccel_method_(int* method);
/* Wrap 'w_ptc_setexactmis' defined in 'wrap.f90'*/
void w_ptc_setexactmis_(int* method);
/* Wrap 'w_ptc_setradiation' defined in 'wrap.f90'*/
void w_ptc_setradiation_(int* method);
/* Wrap 'w_ptc_settotalpath' defined in 'wrap.f90'*/
void w_ptc_settotalpath_(int* method);
/* Wrap 'w_ptc_settime' defined in 'wrap.f90'*/
void w_ptc_settime_(int* method);
/* Wrap 'w_ptc_setnocavity' defined in 'wrap.f90'*/
void w_ptc_setnocavity_(int* method);
/* Wrap 'w_ptc_setfringe' defined in 'wrap.f90'*/
void w_ptc_setfringe_(int* method);
/* Wrap 'w_ptc_end' defined in 'wrap.f90'*/
void w_ptc_end_();
/* Wrap 'w_ptc_getnfieldcomp' defined in 'wrap.f90'*/
void w_ptc_getnfieldcomp_(int* fibreidx,int* ncomp,double* nval);
/* Wrap 'w_ptc_getsfieldcomp' defined in 'wrap.f90'*/
void w_ptc_getsfieldcomp_(int* fibreidx,int* ncomp,double* nval);
/* Wrap 'w_ptc_setfieldcomp' defined in 'wrap.f90'*/
void w_ptc_setfieldcomp_(int* fibreidx);
/* Wrap 'w_ptc_setknobvalue' defined in 'wrap.f90'*/
void w_ptc_setknobvalue_(int* fibre);
/* Wrap 'w_ptc_refreshtables' defined in 'wrap.f90'*/
void w_ptc_refreshtables_();
/* Wrap 'w_ptc_addknob' defined in 'wrap.f90'*/
void w_ptc_addknob_(int* fibre);
/* Wrap 'w_ptc_addknob_i' defined in 'wrap.f90'*/
void w_ptc_addknob_i_(int* paramn);
/* Wrap 'w_ptc_addmoment' defined in 'wrap.f90'*/
void w_ptc_addmoment_(int* x,int* px,int* y,int* py,int* t,int* dp,int* tableIA,int* columnIA,int* parametric);
/* Wrap 'w_ptc_writeparresults' defined in 'wrap.f90'*/
void w_ptc_writeparresults_(int* filename);
/* Wrap 'w_ptc_printframes' defined in 'wrap.f90'*/
void w_ptc_printframes_(int* filename);
/* Wrap 'w_ptc_printlayout_rootm' defined in 'wrap.f90'*/
void w_ptc_printlayout_rootm_(int* filename);
/* Wrap 'w_ptc_eplacement' defined in 'wrap.f90'*/
void w_ptc_eplacement_(int* elementidx,int* rf);
/* Wrap 'w_ptc_addpush' defined in 'wrap.f90'*/
void w_ptc_addpush_(int* tabname,int* colname,int* polinomial,int* monomial);
/* Wrap 'w_ptc_script' defined in 'wrap.f90'*/
void w_ptc_script_(int* scriptname);
/* Wrap 'w_ptc_open_gino' defined in 'wrap.f90'*/
void w_ptc_open_gino_(int* scriptname);
#endif
