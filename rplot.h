#ifndef RPLOT_H
#define RPLOT_H
/*Piotr Skowronski, CERN*/
#ifdef ROOT_PLOT

#include <MadxPlotter.h>
#include <MadxViewer.h>
#define type_OfExtern "C"

#else

#define type_OfExtern

#endif

/*define function names depending if we are in MS WIN or elsewhere*/
#ifndef WIN32


# define newrplot newrplot_
# define plottrack plottrack_
# define plottwiss plottwiss_
# define rplotfinish rplotfinish_
# define rviewer rviewer_
# define madxv_setfctnname madxv_setfctnname_
# define madxv_setknobname madxv_setknobname_
# define madxv_setfunctionat madxv_setfunctionat_

# define type_ofCall


#else

# define newrplot NEWRPLOT
# define plottrack PLOTTRACK
# define plottwiss PLOTTWISS
# define rplotfinish RPLOTFINISH
# define rviewer RVIEWER
# define madxv_setfctnname MADXV_SETFCTNNAME
# define madxv_setknobname MADXV_SETKNOBNAME
# define madxv_setfunctionat MADXV_SETFUNCTIONAT

# define type_ofCall  _stdcall
#endif


void *rplot_handle = 0x0;

/*pointer to rplotter function, C intrtface to MadxPlotter::Fill, see MadxPlotter for details */
typedef void (*rplot_plottrack_fctn)(int,int,int,double, double,double,double,double,double,double);/*function type definition*/
rplot_plottrack_fctn rplot_plottrack = 0x0; /*pointer to function*/



void loadrplotlib();
 
extern type_OfExtern
void type_ofCall plottrack(int* particleno, int* obspoint,  int* turn, double* x, double* xp,
                           double* y, double* yp, double* dpOverP,  double* p,
                           double* length); 

  
extern type_OfExtern
void type_ofCall  plottwiss(int* obspoint, 
                            double* betax, double* alfax,
                            double* betay, double* alfay,
                            double* betaz, double* alfaz,
                            double* length);
  
extern type_OfExtern void type_ofCall rplotfinish();
extern type_OfExtern void type_ofCall rviewer();
extern type_OfExtern void type_ofCall madxv_setfctnname(int* n, const char* name);
extern type_OfExtern void type_ofCall madxv_setfunctionat(int* n, int* el, const char* name);





#endif 
