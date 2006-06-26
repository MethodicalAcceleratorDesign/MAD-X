#include "rplot.h"

#include <stdio.h>

/*
  #define ROOT_PLOT
*/

#ifdef ROOT_PLOT
#include <MadxPlotter.h>
#define type_OfExtern "C"
#else
#define type_OfExtern
#endif


#ifndef WIN32
# define newrplot newrplot_
# define plottrack plottrack_
# define plottwiss plottwiss_
# define rplotfinish rplotfinish_
# define type_ofCall
#else
# define newrplot NEWRPLOT
# define plottrack PLOTTRACK
# define plottwiss PLOTTWISS
# define rplotfinish RPLOTFINISH
# define type_ofCall  _stdcall
#endif


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



extern type_OfExtern void type_ofCall rplotfinish()
{
/*terminates plotter*/
#ifdef ROOT_PLOT
  MadxPlotter::Instance()->Finish();/*writes and deletes the current plot*/
#endif

}

extern type_OfExtern void type_ofCall newrplot()
{
/*adds new plotter*/
#ifdef ROOT_PLOT
  MadxPlotter::Instance()->NewPlot();
#endif
}



extern type_OfExtern void type_ofCall plotter()
{
/*adds new plotter*/
#ifdef ROOT_PLOT
  MadxPlotter::Instance()->Plotter();
#endif
}
