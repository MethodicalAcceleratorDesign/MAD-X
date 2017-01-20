/* Piotr Skowronski, CERN */

#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifdef __cplusplus
extern "C" {
#endif

#include "mad_extrn_f.h"
#include "mad_err.h"
#include "mad_rplot.h"

#ifdef __cplusplus
}
#endif

#ifndef _WIN32
#include <dlfcn.h>
#endif

#ifdef ROOT_PLOT

#include <MadxPlotter.h>
#include <MadxViewer.h>

#endif

#ifdef _PLUGIN

/* pointer to rplotter function, C intrtface to MadxPlotter::Fill, see MadxPlotter for details */
typedef void (*rplot_plottrack_fctn)(int,int,int,double, double,double,double,double,double,double);

static void*                rplot_handle    = 0;
static rplot_plottrack_fctn rplot_plottrack = 0; /*pointer to function*/

static void
loadrplotlib(void)
{
  void *handle;
  char buff[200];
  char* homedir = 0x0;
  
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

  bool root3 = false;
  
  if (root3)
   {
     handle = dlopen( "libCint.so",   RTLD_GLOBAL | RTLD_LAZY);
     if (!handle) {
       fprintf (stderr, "%s\n", dlerror());
       return;
     }
   }
  
  handle = dlopen( "libCore.so",   RTLD_GLOBAL | RTLD_LAZY);
  if (!handle) {
    fprintf (stderr, "%s\n", dlerror());
    return;
  }
  
  handle = dlopen( "libRIO.so",   RTLD_GLOBAL | RTLD_LAZY);
  if (!handle) {
    fprintf (stderr, "%s\n", dlerror());
    return;
  }

  handle = dlopen( "libMathCore.so",   RTLD_GLOBAL | RTLD_LAZY);
  if (!handle) {
    fprintf (stderr, "%s\n", dlerror());
    return;
  }
  
  handle = dlopen( "libNet.so",   RTLD_GLOBAL | RTLD_LAZY);
  if (!handle) {
    fprintf (stderr, "%s\n", dlerror());
    return;
  }
  
  handle = dlopen( "libTree.so",   RTLD_GLOBAL | RTLD_LAZY);
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


  printf("Loading $(HOME)/.madx/plugins/librplot.so\n");
  homedir = getenv("HOME");
  if (strlen(homedir) > 150)
   {
     printf("Home directory name too long: %s",homedir);
   }
  
  sprintf(buff,"%s/.madx/plugins/librplot.so",homedir);

  rplot_handle = dlopen (buff, RTLD_GLOBAL | RTLD_LAZY );
  if (!rplot_handle) {

    printf("Loading librplot.so\n");
    rplot_handle = dlopen ("librplot.so", RTLD_GLOBAL | RTLD_LAZY );
    if (!rplot_handle) {
      fprintf (stderr, "%s\n", dlerror());
      return;
    }

    return;
  }

#if 0 // not used
  warning("rplot.c: loadrplotlib()",
          "Plugin support is not enabled in this version. Unable to use rplot.");
#endif
}

static void
unloadrplotlib(void)
{
  dlclose(rplot_handle);
  rplot_handle = 0;
  rplot_plottrack = 0;
}

#endif // _PLUGIN

void
plottrack(int* particleno, int* obspoint, int* turn,
          double* x, double* xp, double* y, double* yp, double* dpOverP, double* p,
          double* length)
{
  (void)particleno, (void)obspoint, (void)turn;
  (void)x, (void)xp, (void)y, (void)yp;
  (void)dpOverP, (void)p, (void)length; 

#ifdef ROOT_PLOT

  if (MadxPlotter::GetDebug() >9)
  {
    printf("rlot.c: plottrack: %d    %d %d %10.7f %10.7f %10.7f %10.7f %10.7f %10.7f  %10.7f \n",
           *particleno, *obspoint,*turn, *x,     *xp,    *y,     *yp,    *dpOverP, *p, *length);
  }

  MadxPlotter::Instance()->Fill(*particleno, *obspoint, *turn, *x,  *xp,  *y,  *yp,  *dpOverP, *p ,*length);
#elif defined _PLUGIN

  if (rplot_plottrack == 0) return;

  rplot_plottrack(*particleno, *obspoint, *turn, *x,  *xp,  *y,  *yp,  *dpOverP, *p ,*length);

#endif
}

// not used in madx
void
plottwiss(int* obspoint,
          double* betax, double* alfax,
          double* betay, double* alfay,
          double* betaz, double* alfaz,
          double* length)
{
  (void)obspoint;
  (void)betax, (void)alfax, (void)betay, (void)alfay, (void)betaz, (void)alfaz;
  (void)length; 

#ifdef ROOT_PLOT
  if (MadxPlotter::GetDebug() >9)
  {
    printf("rlot.c: plottwiss:   %d   %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n",
           *obspoint,  *betax, *alfax, *betay, *alfay, *betaz, *alfaz, *length);
  }

  MadxPlotter::Instance()->Fill(*obspoint, *betax, *alfax, *betay, *alfay, *betaz, *alfaz, *length);
#endif
}

void
rplotfinish(void)
{
/*terminates plotter*/
#ifdef ROOT_PLOT
  MadxPlotter::Instance()->Finish();/*writes and deletes the current plot*/
#elif defined _PLUGIN

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

  rfinish();

#endif

}

void
newrplot(void)
{
/*adds new plotter*/
#ifdef ROOT_PLOT
  MadxPlotter::Instance()->NewPlot();
#elif defined _PLUGIN

  typedef void (*rvfun)();
  rvfun newrplot;
  char *error;


  if (rplot_handle == 0)
  {
    loadrplotlib();
    if (rplot_handle == 0)
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

  newrplot();
#else
  warning("rplot.c: newrplot()",
          "Plugin support is not enabled in this version. Unable to use rplot.");
#endif
}

#if 0 // not used
static void
plotter(void)
{
/*adds new plotter*/
#ifdef ROOT_PLOT
  MadxPlotter::Instance()->Plotter();
#endif
}
#endif

#if 0 // not used
static void
print(void)
{
/*adds new plotter*/
#ifdef ROOT_PLOT
  MadxPlotter::Instance()->Plotter();
#endif
}
#endif

void
rviewer(void)
{
/*adds new plotter*/

#ifdef ROOT_PLOT

  printf("Compiled with librplot.so support\n");

  MadxViewer::Instance();

#elif defined _PLUGIN

  typedef void (*rvfun)();
  rvfun fctn = 0;
  char *error;

  loadrplotlib();

  if (rplot_handle == 0) return;

  printf("Looking for runrviewer() \n");
  fctn = (rvfun)dlsym(rplot_handle, "runrviewer");
  if ((error = dlerror()) != 0)  {
    fprintf (stderr, "%s\n", error);
    return;
  }

  fctn();

  unloadrplotlib();
#else
  warning("rplot.c: rviewer()",
          "Plugin support is not enabled in this version. Unable to use rviewer.");
#endif
}

void
madxv_setknobname(int* n, const char* name)
{
  (void)n, (void)name;
  
#ifdef ROOT_PLOT
  MadxViewer::Instance()->SetKnobName(*n,name);
#elif defined _PLUGIN

  typedef void (*sknfctn)(int*,const char*);
  sknfctn fctn;
  char *error;

  if(rplot_handle == 0) return;

  fctn = (sknfctn)dlsym(rplot_handle, "madxviewer_setknobname");
  if ((error = dlerror()) != 0)
  {
    fprintf (stderr, "%s\n", error);
    return;
  }

  fctn(n,name);

#endif
}

void
madxv_setfctnname(int* n, const char* name)
{
  (void)n, (void)name;
  
#ifdef ROOT_PLOT
  MadxViewer::Instance()->SetFunctionName(*n,name);
#elif defined _PLUGIN

  typedef void (*sknfctn)(int*,const char*);
  sknfctn fctn;
  char *error;

  if(rplot_handle == 0) return;

  fctn = (sknfctn)dlsym(rplot_handle, "madxviewer_setfctnname");
  if ((error = dlerror()) != 0)
  {
    fprintf (stderr, "%s\n", error);
    return;
  }

  fctn(n,name);

#endif
}

void
madxv_setfunctionat(int* el, int* n,  const char* name)
{
  (void)el, (void)n, (void)name;
  
#ifdef ROOT_PLOT
  MadxViewer::Instance()->SetFunctionAt(*el,*n,name);
#elif defined _PLUGIN

  typedef void (*sknfctn)(int*,int*,const char*);
  sknfctn fctn;
  char *error;

  if(rplot_handle == 0) return;

  fctn = (sknfctn)dlsym(rplot_handle, "madxviewer_setfunctionat");
  if ((error = dlerror()) != 0)
  {
    fprintf (stderr, "%s\n", error);
    return;
  }

  fctn(el,n,name);

#endif
}

