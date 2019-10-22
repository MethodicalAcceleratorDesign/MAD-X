#ifdef _USE_X11

/* X include files */
#define _HPUX_SOURCE
#include <unistd.h> /* hbu for sleep */
#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mad_extrn_f.h"
#include "mad_gxx11c.h"
#include "mad_err.h"
#include "mad_mem.h"

// #define MAXCOL 256 // not used

/* declarations */
static char                   logo[] = "GXplot";
static Display*               mydisplay;
static Window                 mywindow;
static GC                     mygc[4];  /* one gc per linestyle */
static XEvent                 myevent;
// static KeySym                 mykey; // not used
static XSizeHints             myhint;
static XWMHints               mywmhint;
static XSetWindowAttributes   xswa;
static int                    i, done, myscreen, style; // , ncolors // not used
static unsigned int           wwidth, wwidth_mm, wheight, wheight_mm;
static float                  xwidth, xwidth_mm, xheight, xheight_mm;
static float                  w_lx, w_ly, wx_fact, wy_fact;
static float                  v_corrf;
static unsigned long          myforeground, mybackground, valuemask;
static Colormap               cmap;
static int                    colored;
static XColor                 color, colore; //, colors[MAXCOL];
static XWindowAttributes      theAtt;
static XGCValues              gcv;

static char solid[1]  = {0},
            dashed[2] = {3,3},
            dotted[2] = {1,3},
            dot_dashed[4] = {1,3,3,3};
static char *dash_list[] = {solid, dashed, dotted, dot_dashed};
static int dash_list_length[] = {0, 2, 2, 4};

void
wopen(int *uswid, int *ushi)
{
  int argc;
  char **argv;

  argc = 0;
  argv = NULL;

  /* initialization */
  mydisplay = XOpenDisplay((char *) getenv("DISPLAY"));
  if (!mydisplay)
  {
    fprintf(stderr, "cannot open display\n"); exit(1);
  }
  myscreen = DefaultScreen(mydisplay);

  /* default pixel values */

  mybackground = WhitePixel(mydisplay, myscreen);
  myforeground = BlackPixel(mydisplay, myscreen);

  /* default program-specified window position and size */

  myhint.x       = 100; myhint.y = 200;
  myhint.width   = *uswid; myhint.height = *ushi;
  myhint.flags   = PPosition | PSize;
  mywmhint.flags = InputHint;
  mywmhint.input = True;

  /* window creation */

  mywindow = XCreateSimpleWindow(mydisplay, DefaultRootWindow(mydisplay),
                                 myhint.x, myhint.y, myhint.width, myhint.height, 5,
                                 myforeground, mybackground);
  XSetStandardProperties(mydisplay, mywindow, logo, logo,
                         None, argv, argc, &myhint);
  XSetWMHints(mydisplay, mywindow, &mywmhint);
  xswa.backing_store = Always;
  valuemask = CWBackingStore;
  XChangeWindowAttributes(mydisplay, mywindow, valuemask, &xswa);

  /* GC creation and initialization */

  mygc[0] = XCreateGC(mydisplay, mywindow, 0, 0);
  XSetBackground(mydisplay, mygc[0], mybackground);
  XSetForeground(mydisplay, mygc[0], myforeground);
  for (i = 1;i < 4 ;i++ )
  {
    gcv.line_style = LineOnOffDash;
    mygc[i] = XCreateGC(mydisplay, mywindow, GCLineStyle, &gcv);
    XSetDashes(mydisplay, mygc[i], 0, dash_list[i],
               dash_list_length[i]);
    XSetBackground(mydisplay, mygc[i], mybackground);
    XSetForeground(mydisplay, mygc[i], myforeground);
    XSetFillStyle (mydisplay, mygc[i], FillSolid);
    XSetFillRule  (mydisplay, mygc[i], WindingRule);
  }
  XGetWindowAttributes(mydisplay, mywindow, &theAtt);
//  ncolors = theAtt.visual->map_entries; // not used
  cmap    = theAtt.colormap;
  colored = DisplayPlanes( mydisplay, myscreen ) > 1;
  if(colored) cmap = DefaultColormap( mydisplay, myscreen);

  /* input event selection */

  XSelectInput(mydisplay, mywindow,
               ButtonPressMask | KeyPressMask | ExposureMask);

  /* window mapping */

  XMapRaised(mydisplay, mywindow);
  XFlush(mydisplay);
  wwidth = DisplayWidth(mydisplay, myscreen);
  wwidth_mm = DisplayWidthMM(mydisplay, myscreen);
  wheight = DisplayHeight(mydisplay, myscreen);
  wheight_mm = DisplayHeightMM(mydisplay, myscreen);
  if (wwidth != 0 && wwidth_mm != 0 &&
      wheight != 0 && wheight_mm != 0)
  {
    xwidth     = wwidth;
    xwidth_mm  = wwidth_mm;
    xheight    = wheight;
    xheight_mm = wheight_mm;
    v_corrf    = (xheight * xwidth_mm)
      / (xwidth * xheight_mm);
    if ((v_corrf - 1.) * (v_corrf - 1.) <= 0.002)
    {
      v_corrf = 1.;
    }

  }
  sleep(3);
}

void
wclose(void)
{

  /* termination */

  for (i = 0;i < 4 ;i++ )
  {
    XFreeGC(mydisplay, mygc[i]);
  }
  XDestroyWindow(mydisplay, mywindow);
  XCloseDisplay(mydisplay);

}

void
wclrwk(int *i1,int *i2)
{
  (void)i1, (void)i2;
  /* clear workstation */
  XClearWindow(mydisplay, mywindow);

}

void
wpl(int *np, float *xp, float *yp)
{
  /* plot PolyLine of np points with coordinates (xp,yp) */
  const char *rout_name = "wpl";

  XPoint          *points;
  int             i, k;

  /* polyline drawing */
  points = mycalloc(rout_name, *np, sizeof *points);
  for (i = 0;i < *np ;i++ )
  {
    points[i].x = (xp[i] - w_lx) * wx_fact + 0.5;
    k           = (yp[i] - w_ly) * wy_fact + 0.5;
    points[i].y = myhint.height - k;
  }
  XDrawLines(mydisplay, mywindow, mygc[style], points, *np, 0);
  myfree(rout_name, points);

}

void
wfa(int *np, float *xp, float *yp)
{
  /* plot area-filled PolyLine of np points with coordinates (xp,yp) */
  const char *rout_name = "wfa";

  XPoint          *points;
  int             i, k;

  /* polyline drawing */
  points = mycalloc(rout_name, *np, sizeof *points);
  for (i = 0;i < *np ;i++ )
  {
    points[i].x = (xp[i] - w_lx) * wx_fact + 0.5;
    k              = (yp[i] - w_ly) * wy_fact + 0.5;
    points[i].y = myhint.height - k;
  }
  XFillPolygon(mydisplay, mywindow, mygc[style], points, *np,
               Nonconvex, 0);
  myfree(rout_name, points);

}

void
wswn(float *wlx, float *wxfact,
     float *wly, float *wyfact)
{
  /* store window lower world coord.s, and conversion factors */
  w_lx = *wlx; wx_fact = *wxfact;
  w_ly = *wly; wy_fact = *wyfact;

}

void
wtx(float *xp,float *yp, char *string)
{
  int k, x, y;
  x = (*xp - w_lx) * wx_fact + 0.5;
  k              = (*yp - w_ly) * wy_fact + 0.5;
  y = myhint.height - k;
  XDrawString(mydisplay, mywindow, mygc[style],
              x, y, string, strlen(string));
}

void
wwait(void)
{

  /* main event reading loop */
  done = 0;
  while (done == 0)
  {
    /* read the next event */
    XNextEvent(mydisplay, &myevent);

    switch (myevent.type)
    {

      /* process keyboard input */
      case KeyPress:
        done = 1;
        break;
    }
  }
}

void
wsetci(char *uscol)
{
  if (colored && XAllocNamedColor(mydisplay, cmap, uscol, &color, &colore))
    for(int i=0; i<4; i++)
      XSetForeground(mydisplay, mygc[i], color.pixel);
}

void
wsetls(int *ls)
{
  style = *ls;
}
/*
 * return null terminated and blank trimmed string
 */
 
void
wstring( char *s, int *l )
{
  int i, loc;

  loc = *l;
  for( i = 0; i < *l; i++ )
    while( loc > 0 && s[loc-1] == ' ' )
      loc--;

  s[loc] = '\0';
}

#else // case _USE_X11

#include <time.h>
#include "mad_extrn_f.h"
#include "mad_gxx11c.h"

#endif

void
cbyt(int* source, int* s_pos, int* target, int* t_pos, int* n)
/* inserts n_bit byte from source at position s_pos in target at t_pos.
   Attention: least significant bit is #1, positions are those of the least
   significant bit in byte */
{
  int mask = 1;
  int so = *source, tg = *target;
  so >>= (*s_pos - 1);
  mask <<= *n;
  mask -= 1;
  so &= mask;
  mask <<= (*t_pos - 1);  so <<= (*t_pos - 1);
  tg &= ~mask; tg |= so;
  *target = tg;
}

void
mydtime(int* year, int* month, int* day, int* hour, int* minute, int* sec)
{
  time_t _time;
  struct tm* tm;
  time(&_time); /* initialize timing */
  tm = localtime(&_time); /* split system time */
  *year = tm->tm_year%100;
  *month = tm->tm_mon+1;
  *day = tm->tm_mday;
  *hour = tm->tm_hour;
  *minute = tm->tm_min;
  *sec = tm->tm_sec;
}

