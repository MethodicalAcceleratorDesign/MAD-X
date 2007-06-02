#include "rplot.h"
/*Piotr Skowronski, CERN*/

#include <stdio.h>

#define WIN32 _WIN32

#ifndef _WIN32
#include <dlfcn.h>
#endif


extern type_OfExtern void type_ofCall warning(const char*, const char*);

void loadrplotlib()
{
#if PLUGIN_SUPPORT
  void *handle;
  /*simple check if the library was linked dynamically*/
  handle = dlopen( 0,  RTLD_GLOBAL | RTLD_LAZY);
  if (handle == 0x0)
  {
    warning("rplot.c: loadrplotlib()",
            "Most probobly yur binary is statically linked, what disables plugin support.\n");
    return;
  }
/*   printf("0 loading Done\n");*/
  if (!handle) {
    fprintf (stderr, "%s\n", dlerror());
    return;
  }


  handle = dlopen( "libCint.so",   RTLD_GLOBAL | RTLD_LAZY);
  if (!handle) {
    fprintf (stderr, "%s\n", dlerror());
    return;
  }

  handle = dlopen( "libCore.so",   RTLD_GLOBAL | RTLD_LAZY);
  if (!handle) {
    fprintf (stderr, "%s\n", dlerror());
    return;
  }


  handle = dlopen( "libTree.so",   RTLD_GLOBAL | RTLD_LAZY);
  if (!handle) {
    fprintf (stderr, "%s\n", dlerror());
    return;
  }

  handle = dlopen( "libGui.so",   RTLD_GLOBAL | RTLD_LAZY);
  if (!handle) {
    fprintf (stderr, "%s\n", dlerror());
    return;
  }

  handle = dlopen( "libMatrix.so",   RTLD_GLOBAL | RTLD_LAZY);
  if (!handle) {
    fprintf (stderr, "%s\n", dlerror());
    return;
  }

  handle = dlopen( "libHist.so",   RTLD_GLOBAL | RTLD_LAZY);
  if (!handle) {
    fprintf (stderr, "%s\n", dlerror());
    return;
  }

  handle = dlopen( "libGraf.so",   RTLD_GLOBAL | RTLD_LAZY);
  if (!handle) {
    fprintf (stderr, "%s\n", dlerror());
    return;
  }

  handle = dlopen( "libGpad.so",   RTLD_GLOBAL | RTLD_LAZY);
  if (!handle) {
    fprintf (stderr, "%s\n", dlerror());
    return;
  }

  handle = dlopen( "libGed.so",   RTLD_GLOBAL | RTLD_NOW);
  if (!handle) {
    fprintf (stderr, "%s\n", dlerror());
    return;
  }

  printf("Loading librplot.so\n");
  rplot_handle = dlopen ("librplot.so", RTLD_GLOBAL | RTLD_LAZY );
  if (!rplot_handle) {
    fprintf (stderr, "%s\n", dlerror());
    return;
  }

#else
  warning("rplot.c: loadrplotlib()",
          "Plugin support is not enabled in this version. Unable to use rplot.");
#endif
}
/*_________________________________________________________________________________________*/

void unloadrplotlib()
{
#if PLUGIN_SUPPORT
  dlclose(rplot_handle);
  rplot_handle = 0x0;
  rplot_plottrack = 0x0;
#endif
}
/*_________________________________________________________________________________________*/

extern type_OfExtern
void type_ofCall plottrack(int* particleno, int* obspoint, int* turn,
                           double* x, double* xp,
                           double* y, double* yp,
                           double* dpOverP,  double* p,
                           double* length)
{

#ifdef ROOT_PLOT

  if (MadxPlotter::GetDebug() >9)
  {
    printf("rlot.c: plottrack: %d    %d %d %10.7f %10.7f %10.7f %10.7f %10.7f %10.7f  %10.7f \n",
           *particleno, *obspoint,*turn, *x,     *xp,    *y,     *yp,    *dpOverP, *p, *length);
  }

  MadxPlotter::Instance()->Fill(*particleno, *obspoint, *turn, *x,  *xp,  *y,  *yp,  *dpOverP, *p ,*length);
#elif defined PLUGIN_SUPPORT

  if (rplot_plottrack == 0x0) return;

  (*rplot_plottrack)(*particleno, *obspoint, *turn, *x,  *xp,  *y,  *yp,  *dpOverP, *p ,*length);

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
#elif defined PLUGIN_SUPPORT
  typedef void (*rvfun)();
  rvfun rfinish;
  char *error;

  if (rplot_handle == 0x0) return;

  rfinish = (rvfun)dlsym(rplot_handle, "madxplotter_rplotfinish");
  if ((error = dlerror()) != NULL)
  {
    fprintf (stderr, "%s\n", error);
    return;
  }

  (*rfinish)();

#endif

}
/*_________________________________________________________________________________________*/

extern type_OfExtern void type_ofCall newrplot()
{
/*adds new plotter*/
#ifdef ROOT_PLOT
  MadxPlotter::Instance()->NewPlot();
#elif defined PLUGIN_SUPPORT

  typedef void (*rvfun)();
  rvfun newrplot;
  char *error;


  if (rplot_handle == 0x0)
  {
    loadrplotlib();
    if (rplot_handle == 0x0)
    {
      /*It means that library was not loaded*/
      return;
    }
  }

  newrplot = (rvfun)dlsym(rplot_handle, "madxplotter_newrplot");
  if ((error = dlerror()) != NULL)
  {
    fprintf (stderr, "%s\n", error);
    return;
  }

  rplot_plottrack = (rplot_plottrack_fctn)dlsym(rplot_handle, "madxplotter_plottrack");
  if ((error = dlerror()) != NULL)
  {
    fprintf (stderr, "%s\n", error);
    return;
  }

  (*newrplot)();
#else
  warning("rplot.c: newrplot()",
          "Plugin support is not enabled in this version. Unable to use rplot.");
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

  printf("Compiled with librplot.so support\n");

  MadxViewer::Instance();

#elif defined PLUGIN_SUPPORT

  typedef void (*rvfun)();
  rvfun viewer;
  void* fctn = 0x0;
  char *error;

  loadrplotlib();

  if (rplot_handle == 0x0) return;

  printf("Looking for runrviewer() \n");
  fctn = dlsym(rplot_handle, "runrviewer");
  if ((error = dlerror()) != NULL)  {
    fprintf (stderr, "%s\n", error);
    return;
  }

  viewer = (rvfun)fctn;

  (*viewer)();

  unloadrplotlib();
#else
  warning("rplot.c: rviewer()",
          "Plugin support is not enabled in this version. Unable to use rviewer.");
#endif
}
/*_________________________________________________________________________________________*/

extern type_OfExtern
void type_ofCall madxv_setknobname(int* n, const char* name)
{
#ifdef ROOT_PLOT
  MadxViewer::Instance()->SetKnobName(*n,name);
#elif defined PLUGIN_SUPPORT
  typedef void (*sknfctn)(int*,const char*);
  sknfctn fctn;
  char *error;

  if(rplot_handle == 0x0) return;

  fctn = (sknfctn)dlsym(rplot_handle, "madxviewer_setknobname");
  if ((error = dlerror()) != NULL)
  {
    fprintf (stderr, "%s\n", error);
    return;
  }

  (*fctn)(n,name);

#endif
}
/*_________________________________________________________________________________________*/


extern type_OfExtern
void type_ofCall madxv_setfctnname(int* n, const char* name)
{
#ifdef ROOT_PLOT
  MadxViewer::Instance()->SetFunctionName(*n,name);
#elif defined PLUGIN_SUPPORT

  typedef void (*sknfctn)(int*,const char*);
  sknfctn fctn;
  char *error;

  if(rplot_handle == 0x0) return;

  fctn = (sknfctn)dlsym(rplot_handle, "madxviewer_setfctnname");
  if ((error = dlerror()) != NULL)
  {
    fprintf (stderr, "%s\n", error);
    return;
  }

  (*fctn)(n,name);

#endif
}
/*_________________________________________________________________________________________*/


extern type_OfExtern
void type_ofCall madxv_setfunctionat(int* el, int* n,  const char* name)
{
#ifdef ROOT_PLOT
  MadxViewer::Instance()->SetFunctionAt(*el,*n,name);
#elif defined PLUGIN_SUPPORT

  typedef void (*sknfctn)(int*,int*,const char*);
  sknfctn fctn;
  char *error;

  if(rplot_handle == 0x0) return;

  fctn = (sknfctn)dlsym(rplot_handle, "madxviewer_setfunctionat");
  if ((error = dlerror()) != NULL)
  {
    fprintf (stderr, "%s\n", error);
    return;
  }

  (*fctn)(el,n,name);

#endif
}

/*_________________________________________________________________________________________*/
