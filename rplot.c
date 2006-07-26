#include "rplot.h"

#include <stdio.h>




extern type_OfExtern
void type_ofCall plottrack(int* particleno, int* obspoint,
                           double* x, double* xp,
                           double* y, double* yp,
                           double* dpOverP,  double* p,
                           double* length)
{

#ifdef ROOT_PLOT

  if (MadxPlotter::GetDebug() >9)
  {
    printf("rlot.c: plottrack: %d      %d   %10.7f %10.7f %10.7f %10.7f %10.7f %10.7f  %10.7f \n",
           *particleno, *obspoint, *x,     *xp,    *y,     *yp,    *dpOverP, *p, *length);
  }

  MadxPlotter::Instance()->Fill(*particleno, *obspoint, *x,  *xp,  *y,  *yp,  *dpOverP, *p ,*length);
#endif
}
/*_________________________________________________________________________________________*/

extern type_OfExtern
void type_ofCall  plottwiss(int* obspoint,
                            double* betax, double* alfax,
                            double* betay, double* alfay,
                            double* betaz, double* alfaz,
                            double* length)
{

#ifdef ROOT_PLOT

  if (MadxPlotter::GetDebug() >9)
  {
    printf("rlot.c: plottwiss:   %d   %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n",
           *obspoint,  *betax, *alfax, *betay, *alfay, *betaz, *alfaz, *length);
  }

  MadxPlotter::Instance()->Fill(*obspoint, *betax, *alfax, *betay, *alfay, *betaz, *alfaz, *length);
#endif
}
/*_________________________________________________________________________________________*/



extern type_OfExtern void type_ofCall rplotfinish()
{
/*terminates plotter*/
#ifdef ROOT_PLOT
  MadxPlotter::Instance()->Finish();/*writes and deletes the current plot*/
#endif

}
/*_________________________________________________________________________________________*/

extern type_OfExtern void type_ofCall newrplot()
{
/*adds new plotter*/
#ifdef ROOT_PLOT
  MadxPlotter::Instance()->NewPlot();
#endif
}
/*_________________________________________________________________________________________*/

extern type_OfExtern void type_ofCall plotter()
{
/*adds new plotter*/
#ifdef ROOT_PLOT
  MadxPlotter::Instance()->Plotter();
#endif
}
/*_________________________________________________________________________________________*/


extern type_OfExtern void type_ofCall print()
{
/*adds new plotter*/
#ifdef ROOT_PLOT
  MadxPlotter::Instance()->Plotter();
#endif
}
/*_________________________________________________________________________________________*/


extern type_OfExtern void type_ofCall rviewer()
{
/*adds new plotter*/
#ifdef ROOT_PLOT
  
  MadxViewer::Instance();
#endif
}
/*_________________________________________________________________________________________*/

extern type_OfExtern 
void type_ofCall madxv_setknobname(int* n, const char* name)
{
#ifdef ROOT_PLOT
  MadxViewer::Instance()->SetKnobName(*n,name);
#endif
}
/*_________________________________________________________________________________________*/


extern type_OfExtern 
void type_ofCall madxv_setfctnname(int* n, const char* name)
{
#ifdef ROOT_PLOT
  MadxViewer::Instance()->SetFunctionName(*n,name);
#endif
}
/*_________________________________________________________________________________________*/


extern type_OfExtern 
void type_ofCall madxv_setfunctionat(int* el, int* n,  const char* name)
{
#ifdef ROOT_PLOT
  MadxViewer::Instance()->SetFunctionAt(*el,*n,name);
#endif
}

