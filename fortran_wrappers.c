/* set of 64 wrappers to synchronize FORTRAN and C stdout buffers */
/* when crossing the border upon calling FORTRAN from C. */

#include <stdio.h>
#include "fortran_prototypes.h"

extern void call_fortran_flush_();

/* Wrap 'dynap' defined in 'dynap.f90' */
void dynap_wrapper(double* eigen,double* coords,int* turns,int* npart,double* distvect,double* zn,double* dq,double* onelog,double* turnnumber){
	fflush(stdout);
	dynap_(eigen,coords,turns,npart,distvect,zn,dq,onelog,turnnumber);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_rviewer' defined in 'madx_ptc_knobs.f90' */
void w_ptc_rviewer_wrapper(){
	fflush(stdout);
	w_ptc_rviewer_();
	call_fortran_flush_();
}
/* Wrap 'mtgetc' defined in 'match.f90' */
void mtgetc_wrapper(double* vect,double* dvect){
	fflush(stdout);
	mtgetc_(vect,dvect);
	call_fortran_flush_();
}
/* Wrap 'collect' defined in 'match.f90' */
void collect_wrapper(int* ncon,double* fsum,double* fvect){
	fflush(stdout);
	collect_(ncon,fsum,fvect);
	call_fortran_flush_();
}
/* Wrap 'mtsvd' defined in 'matchjc.f90' */
void mtsvd_wrapper(int* M,int* N,double* fjac,double* SV,double* U,double* VT){
	fflush(stdout);
	mtsvd_(M,N,fjac,SV,U,VT);
	call_fortran_flush_();
}
/* Wrap 'setup' defined in 'orbf.f90' */
void setup_wrapper(double* resp,double* a,int* im,int* ic,int* nm,int* nc){
	fflush(stdout);
	setup_(resp,a,im,ic,nm,nc);
	call_fortran_flush_();
}
/* Wrap 'setupi' defined in 'orbf.f90' */
void setupi_wrapper(int* resp,int* a,int* im,int* ic,int* nm,int* nc){
	fflush(stdout);
	setupi_(resp,a,im,ic,nm,nc);
	call_fortran_flush_();
}
/* Wrap 'micit' defined in 'orbf.f90' */
void micit_wrapper(double* a,char* conm,double* xin,double* cin,double* res,int* nx,float* rms,int* im,int* ic,int* iter,int* ny,float* ax,float* cinx,float* xinx,float* resx,float* rho,float* ptop,float* rmss,float* xrms,float* xptp,float* xiter,int* ifail){
	fflush(stdout);
	micit_(a,conm,xin,cin,res,nx,rms,im,ic,iter,ny,ax,cinx,xinx,resx,rho,ptop,rmss,xrms,xptp,xiter,ifail);
	call_fortran_flush_();
}
/* Wrap 'haveit' defined in 'orbf.f90' */
void haveit_wrapper(double* a,double* xin,double* cin,double* res,int* nx,int* im,int* ic,double* cb,double* xmeas,double* xres,double* y,double* z,double* xd){
	fflush(stdout);
	haveit_(a,xin,cin,res,nx,im,ic,cb,xmeas,xres,y,z,xd);
	call_fortran_flush_();
}
/* Wrap 'svddec_m' defined in 'orbf.f90' */
void svddec_m_wrapper(double* a,double* svdmat,double* umat,double* vmat,double* wmat,double* utmat,double* vtmat,double* wtmat,double* ws,double* wvec,int* sortw,double* sngcut,double* sngval,int* im,int* ic,int* iflag,int* sing,int* dbg){
	fflush(stdout);
	svddec_m_(a,svdmat,umat,vmat,wmat,utmat,vtmat,wtmat,ws,wvec,sortw,sngcut,sngval,im,ic,iflag,sing,dbg);
	call_fortran_flush_();
}
/* Wrap 'svddec_c' defined in 'orbf.f90' */
void svddec_c_wrapper(double* a,double* svdmat,double* umat,double* vmat,double* wmat,double* utmat,double* vtmat,double* wtmat,double* ws,double* wvec,int* sortw,double* sngcut,double* sngval,int* im,int* ic,int* iflag,int* sing,int* dbg){
	fflush(stdout);
	svddec_c_(a,svdmat,umat,vmat,wmat,utmat,vtmat,wtmat,ws,wvec,sortw,sngcut,sngval,im,ic,iflag,sing,dbg);
	call_fortran_flush_();
}
/* Wrap 'svdcorr_m' defined in 'orbf.f90' */
void svdcorr_m_wrapper(double* a,double* svdmat,double* umat,double* vmat,double* wmat,double* utmat,double* vtmat,double* wtmat,double* xin,double* xc,double* xout,double* xa,double* xb,double* xpred,double* ws,double* wvec,int* sortw,int* nx,int* im,int* ic,int* iflag,int* dbg){
	fflush(stdout);
	svdcorr_m_(a,svdmat,umat,vmat,wmat,utmat,vtmat,wtmat,xin,xc,xout,xa,xb,xpred,ws,wvec,sortw,nx,im,ic,iflag,dbg);
	call_fortran_flush_();
}
/* Wrap 'svdcorr_c' defined in 'orbf.f90' */
void svdcorr_c_wrapper(double* a,double* svdmat,double* umat,double* vmat,double* wmat,double* utmat,double* vtmat,double* wtmat,double* xin,double* xc,double* xout,double* xa,double* xb,double* xpred,double* ws,double* wvec,int* sortw,int* nx,int* im,int* ic,int* iflag,int* dbg){
	fflush(stdout);
	svdcorr_c_(a,svdmat,umat,vmat,wmat,utmat,vtmat,wtmat,xin,xc,xout,xa,xb,xpred,ws,wvec,sortw,nx,im,ic,iflag,dbg);
	call_fortran_flush_();
}
/* Wrap 'primat' defined in 'orbf.f90' */
void primat_wrapper(int* a,int* nc,int* nm){
	fflush(stdout);
	primat_(a,nc,nm);
	call_fortran_flush_();
}
/* Wrap 'prdmat' defined in 'orbf.f90' */
void prdmat_wrapper(double* a,int* nc,int* nm){
	fflush(stdout);
	prdmat_(a,nc,nm);
	call_fortran_flush_();
}
/* Wrap 'pefill' defined in 'plot.f90' */
void pefill_wrapper(int* ierr){
	fflush(stdout);
	pefill_(ierr);
	call_fortran_flush_();
}
/* Wrap 'pesopt' defined in 'plot.f90' */
void pesopt_wrapper(int* ierr){
	fflush(stdout);
	pesopt_(ierr);
	call_fortran_flush_();
}
/* Wrap 'plotit' defined in 'plot.f90' */
void plotit_wrapper(int* initfl){
	fflush(stdout);
	plotit_(initfl);
	call_fortran_flush_();
}
/* Wrap 'res_index' defined in 'resindex.F' */
void res_index_wrapper(int* skew,int* mynorder,int* myn1,int* myn2,int indexa[4][1000],int* mynres){
	fflush(stdout);
	res_index_(skew,mynorder,myn1,myn2,indexa,mynres);
	call_fortran_flush_();
}
/* Wrap 'soddin' defined in 'sodd.f90' */
void soddin_wrapper(int* ierr){
	fflush(stdout);
	soddin_(ierr);
	call_fortran_flush_();
}
/* Wrap 'timest' defined in 'timest.f90' */
void timest_wrapper(float* r1){
	fflush(stdout);
	timest_(r1);
	call_fortran_flush_();
}
/* Wrap 'timex' defined in 'timex.f90' */
void timex_wrapper(float* r1){
	fflush(stdout);
	timex_(r1);
	call_fortran_flush_();
}
/* Wrap 'twiss' defined in 'twiss.f90' */
void twiss_wrapper(double* rt,double* disp0,int* tab_name,int* sector_tab_name){
	fflush(stdout);
	twiss_(rt,disp0,tab_name,sector_tab_name);
	call_fortran_flush_();
}
/* Wrap 'tmrefe' defined in 'twiss.f90' */
void tmrefe_wrapper(double* rt){
	fflush(stdout);
	tmrefe_(rt);
	call_fortran_flush_();
}
/* Wrap 'tmrefo' defined in 'twiss.f90' */
void tmrefo_wrapper(int* kobs,double* orbit0,double* orbit,double* rt){
	fflush(stdout);
	tmrefo_(kobs,orbit0,orbit,rt);
	call_fortran_flush_();
}
/* Wrap 'getclor' defined in 'util.f90' */
void getclor_wrapper(double* orbit0,double* rt,double* tt,int* error){
	fflush(stdout);
	getclor_(orbit0,rt,tt,error);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_create_universe' defined in 'wrap.f90' */
void w_ptc_create_universe_wrapper(){
	fflush(stdout);
	w_ptc_create_universe_();
	call_fortran_flush_();
}
/* Wrap 'w_ptc_create_layout' defined in 'wrap.f90' */
void w_ptc_create_layout_wrapper(){
	fflush(stdout);
	w_ptc_create_layout_();
	call_fortran_flush_();
}
/* Wrap 'w_ptc_export_xml' defined in 'wrap.f90' */
void w_ptc_export_xml_wrapper(int* filename){
	fflush(stdout);
	w_ptc_export_xml_(filename);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_move_to_layout' defined in 'wrap.f90' */
void w_ptc_move_to_layout_wrapper(){
	fflush(stdout);
	w_ptc_move_to_layout_();
	call_fortran_flush_();
}
/* Wrap 'w_ptc_input' defined in 'wrap.f90' */
void w_ptc_input_wrapper(){
	fflush(stdout);
	w_ptc_input_();
	call_fortran_flush_();
}
/* Wrap 'w_ptc_align' defined in 'wrap.f90' */
void w_ptc_align_wrapper(){
	fflush(stdout);
	w_ptc_align_();
	call_fortran_flush_();
}
/* Wrap 'w_ptc_twiss' defined in 'wrap.f90' */
void w_ptc_twiss_wrapper(int* tab_name,int* summary_name){
	fflush(stdout);
	w_ptc_twiss_(tab_name,summary_name);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_normal' defined in 'wrap.f90' */
void w_ptc_normal_wrapper(){
	fflush(stdout);
	w_ptc_normal_();
	call_fortran_flush_();
}
/* Wrap 'w_ptc_moments' defined in 'wrap.f90' */
void w_ptc_moments_wrapper(int* no){
	fflush(stdout);
	w_ptc_moments_(no);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_initmoments' defined in 'wrap.f90' */
void w_ptc_initmoments_wrapper(){
	fflush(stdout);
	w_ptc_initmoments_();
	call_fortran_flush_();
}
/* Wrap 'w_ptc_dumpmaps' defined in 'wrap.f90' */
void w_ptc_dumpmaps_wrapper(){
	fflush(stdout);
	w_ptc_dumpmaps_();
	call_fortran_flush_();
}
/* Wrap 'w_ptc_track' defined in 'wrap.f90' */
void w_ptc_track_wrapper(int* max_obs){
	fflush(stdout);
	w_ptc_track_(max_obs);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_trackline' defined in 'wrap.f90' */
void w_ptc_trackline_wrapper(int* max_obs){
	fflush(stdout);
	w_ptc_trackline_(max_obs);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_track_everystep' defined in 'wrap.f90' */
void w_ptc_track_everystep_wrapper(int* max_obs){
	fflush(stdout);
	w_ptc_track_everystep_(max_obs);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_setdebuglevel' defined in 'wrap.f90' */
void w_ptc_setdebuglevel_wrapper(int* level){
	fflush(stdout);
	w_ptc_setdebuglevel_(level);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_setaccel_method' defined in 'wrap.f90' */
void w_ptc_setaccel_method_wrapper(int* method){
	fflush(stdout);
	w_ptc_setaccel_method_(method);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_setexactmis' defined in 'wrap.f90' */
void w_ptc_setexactmis_wrapper(int* method){
	fflush(stdout);
	w_ptc_setexactmis_(method);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_setradiation' defined in 'wrap.f90' */
void w_ptc_setradiation_wrapper(int* method){
	fflush(stdout);
	w_ptc_setradiation_(method);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_settotalpath' defined in 'wrap.f90' */
void w_ptc_settotalpath_wrapper(int* method){
	fflush(stdout);
	w_ptc_settotalpath_(method);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_settime' defined in 'wrap.f90' */
void w_ptc_settime_wrapper(int* method){
	fflush(stdout);
	w_ptc_settime_(method);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_setnocavity' defined in 'wrap.f90' */
void w_ptc_setnocavity_wrapper(int* method){
	fflush(stdout);
	w_ptc_setnocavity_(method);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_setfringe' defined in 'wrap.f90' */
void w_ptc_setfringe_wrapper(int* method){
	fflush(stdout);
	w_ptc_setfringe_(method);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_end' defined in 'wrap.f90' */
void w_ptc_end_wrapper(){
	fflush(stdout);
	w_ptc_end_();
	call_fortran_flush_();
}
/* Wrap 'w_ptc_getnfieldcomp' defined in 'wrap.f90' */
void w_ptc_getnfieldcomp_wrapper(int* fibreidx,int* ncomp,double* nval){
	fflush(stdout);
	w_ptc_getnfieldcomp_(fibreidx,ncomp,nval);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_getsfieldcomp' defined in 'wrap.f90' */
void w_ptc_getsfieldcomp_wrapper(int* fibreidx,int* ncomp,double* nval){
	fflush(stdout);
	w_ptc_getsfieldcomp_(fibreidx,ncomp,nval);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_setfieldcomp' defined in 'wrap.f90' */
void w_ptc_setfieldcomp_wrapper(int* fibreidx){
	fflush(stdout);
	w_ptc_setfieldcomp_(fibreidx);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_setknobvalue' defined in 'wrap.f90' */
void w_ptc_setknobvalue_wrapper(int* fibre){
	fflush(stdout);
	w_ptc_setknobvalue_(fibre);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_refreshtables' defined in 'wrap.f90' */
void w_ptc_refreshtables_wrapper(){
	fflush(stdout);
	w_ptc_refreshtables_();
	call_fortran_flush_();
}
/* Wrap 'w_ptc_addknob' defined in 'wrap.f90' */
void w_ptc_addknob_wrapper(int* fibre){
	fflush(stdout);
	w_ptc_addknob_(fibre);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_addknob_i' defined in 'wrap.f90' */
void w_ptc_addknob_i_wrapper(int* paramn){
	fflush(stdout);
	w_ptc_addknob_i_(paramn);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_addmoment' defined in 'wrap.f90' */
void w_ptc_addmoment_wrapper(int* x,int* px,int* y,int* py,int* t,int* dp,int* tableIA,int* columnIA,int* parametric){
	fflush(stdout);
	w_ptc_addmoment_(x,px,y,py,t,dp,tableIA,columnIA,parametric);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_writeparresults' defined in 'wrap.f90' */
void w_ptc_writeparresults_wrapper(int* filename){
	fflush(stdout);
	w_ptc_writeparresults_(filename);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_printframes' defined in 'wrap.f90' */
void w_ptc_printframes_wrapper(int* filename){
	fflush(stdout);
	w_ptc_printframes_(filename);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_printlayout_rootm' defined in 'wrap.f90' */
void w_ptc_printlayout_rootm_wrapper(int* filename){
	fflush(stdout);
	w_ptc_printlayout_rootm_(filename);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_eplacement' defined in 'wrap.f90' */
void w_ptc_eplacement_wrapper(int* elementidx,int* rf){
	fflush(stdout);
	w_ptc_eplacement_(elementidx,rf);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_addpush' defined in 'wrap.f90' */
void w_ptc_addpush_wrapper(int* tabname,int* colname,int* polinomial,int* monomial){
	fflush(stdout);
	w_ptc_addpush_(tabname,colname,polinomial,monomial);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_script' defined in 'wrap.f90' */
void w_ptc_script_wrapper(int* scriptname){
	fflush(stdout);
	w_ptc_script_(scriptname);
	call_fortran_flush_();
}
/* Wrap 'w_ptc_open_gino' defined in 'wrap.f90' */
void w_ptc_open_gino_wrapper(int* scriptname){
	fflush(stdout);
	w_ptc_open_gino_(scriptname);
	call_fortran_flush_();
}
