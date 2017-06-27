#include "madx.h"

// types and constants

#define MIN_DOUBLE 1.e-36

struct table;
struct aper_node            /* aperture limit node */
{
  char name[NAME_L];
  double n1;
  double s;
  char apertype[NAME_L];
  double aperture[4];
  double aper_tol[3];
  double deltap_twiss;
};

struct aper_e_d             /* element displacement */
{
  char name[NAME_L];        /* element name */
  int curr;                 /* # of rows */
  double tab[E_D_MAX][3];   /* the table of read values */
};

// static interface

static struct aper_e_d* true_tab;

/* struct aper_e_d* offs_tab;*/
static struct table* offs_tab;

static struct aper_node*
aperture(char *table, struct node* use_range[], struct table* tw_cp, int *tw_cnt, struct aper_node*);

static int
aper_rectellipse(double* ap1, double* ap2, double* ap3, double* ap4, int* quarterlength, double tablex[], double tabley[])
{ // build a quadrant of a polygon based on a rectangular-ellipse shape; aper_fill_quadrants() completes the polygon.
  double x, y, angle, alfa, theta, dangle;
  int i, napex;

  int debug = get_option("debug");
  if (debug)
    printf("+++ aper_rectellipse: %10.5f  %10.5f  %10.5f  %10.5f %d\n",*ap1,*ap2,*ap3,*ap4,*quarterlength);

  if ( *ap1 < MIN_DOUBLE*1.e10 || *ap2 < MIN_DOUBLE*1.e10) {
    fatal_error("Illegal negative or too small value in aper_rectellipse", "for ap1 or ap2 (rectangle)");
  }
  if ( *ap3 < MIN_DOUBLE*1.e10 || *ap4 < MIN_DOUBLE*1.e10) {
    fatal_error("Illegal negative or too small value in aper_rectellipse", "for ap3 or ap4 (ellipse)");
  }

  /* find angles where rectangle and ellipse cross */

  if ( (*ap1) >= (*ap3) ) alfa = 0. ;
  // the horizontal rectangle extent is larger than the ellipse
  // the curved part starts at angle 0
  else {
    y = (*ap3)*sqrt(1-((*ap1)*(*ap1))/((*ap3)*(*ap3)));

    if (y > (*ap2)) // the rectangle is contained within the ellipse; there is no curved part.
      alfa = atan2((*ap2),(*ap1));
    else            // the rectangle extends beyond the ellipse; there is a curved part
      alfa = atan2(y,*ap1); // this angle is not the geometrical angle
  }

  if ( (*ap2) >= (*ap4) ) theta = 0.;
  // the vertical rectangle extent is larger than the ellipse
  // the curved part extends all the way to pi/2
  else {
    x = (*ap4) * sqrt(1 - ((*ap2)*(*ap2))/((*ap4)*(*ap4)));

    if (x > (*ap1)) // the rectangle is contained within the ellipse; there is no curved part.
      theta = atan2((*ap1),(*ap2));
    else            // the rectangle extends beyond the ellipse; there is a curved part
      theta = atan2(x,*ap2); // this angle is not the geometrical angle
  }

  // at this point we know if we have a full ellipse (alfa=0, theta=0), a full rectangle (alfa+theta=pi/2) or a mixed curve.
  // we can calculate the number of apexes (napex+1) and the interval depending on the shape

  if (*quarterlength) {
    napex = 9 ;
    dangle = ((pi/2) - (alfa+theta)) / napex;
  }

  else if (fabs(alfa+theta-pi/2) < MIN_DOUBLE * 1.e10) { //rectangle, single point, zero interval
    napex = 0 ;
    dangle = 0.;
  }

  else { // general rectellipse,  1 <= napex <= 19 intervals
    napex = 1 + floor(18 * fabs(1-(alfa+theta)/(pi/2))) ;
    dangle = ((pi/2)-(alfa+theta))/napex;
  }

  /*write coordinates for first quadrant*/

  if (napex == 0) { // case of rectangle
    i=0;
    tablex[i] = *ap1;
    tabley[i] = *ap2;
    i++;
  }

  else { // there is a curved part.
    for ( i=0 ; i<=napex; i++ )	{
      angle = alfa + i*dangle; // this angle is not the geometrical angle
      tablex[i] = (*ap3)*cos(angle);
      tabley[i] = (*ap4)*sin(angle);

      // the geometrical angle alfag is such that tan(alfag) = tabley[i]/tablex[i] = ap4/ap3 * tan(alfa)

      if (i >= MAXARRAY/4) fatal_error("Memory full in aper_rectellipse", "Number of coordinates exceeds set limit");
      // should give the value of the MAXARRAY set limit that is exceeded
      // warn("aper_rectellipse: number of coordinates exceeds MAXARRAY = ", MAXARRAY)
    }
  }

  *quarterlength = i-1;
  if (debug) printf("quarterlength : %d\n", *quarterlength);

  return 0;
}

static void
aper_race(double xshift, double yshift, double r, double angle, double* x, double* y)
{
  /* this function calculates the displacement of the beam centre
     due to tolerance uncertainty for every angle */
  /* NEW VERSION of aper_race, 20feb08 BJ, potential zero-divide issues cleared */

  double angle0, angle1, angle2, alfa, gamma, theta;
  int quadrant;

  if (xshift==0 && yshift==0 && r==0) {
    *x=0; *y=0;
    return;
  }

  quadrant = angle/(pi/2) + 1;

  switch (quadrant) /*adjusting angle to first quadrant*/
    {
    case 1: /* angle = angle; */    break;
    case 2: angle = pi - angle;     break;
    case 3: angle = angle - pi;     break;
    case 4: angle = twopi - angle;  break;
    }

  if (angle == pi/2) {
    *x=0;
    *y=yshift+r;
  }

  else {
    /*finding where arc starts and ends*/
    angle0 = atan2( yshift , xshift + r );
    angle1 = atan2( r + yshift , xshift );

    /*different methods is needed, depending on angle*/
    if (angle <= angle0 + MIN_DOUBLE * 1.e10 ) {
      *x = xshift + r;
      *y = tan(angle) * (xshift+r);
    }
    else if (angle < angle1) {
      /* if this is a circle, angle2 useless */
      if (!xshift && !yshift)  angle2 = 0;
      else angle2 = atan2( yshift , xshift );

      alfa = fabs(angle-angle2);
      if (alfa < MIN_DOUBLE * 1.e10) {
	/* alfa==0 is a simpler case */
        *x = cos(angle) * (r + sqrt(xshift*xshift + yshift*yshift));
        *y = sin(angle) * (r + sqrt(xshift*xshift + yshift*yshift));
      }
      else {
	/* solving sine rule w.r.t. gamma */
        gamma = asin(sqrt(xshift*xshift + yshift*yshift)/r*sin(alfa));
	/*theta is the last corner in the triangle*/
        theta = pi - (alfa+gamma);
        *x = cos(angle) * r * sin(theta)/sin(alfa);
        *y = sin(angle) * r * sin(theta)/sin(alfa);
      }
    }
    /* upper flat part */
    else {
      *y = r + yshift;
      *x = (r+yshift) * tan(pi/2-angle);
    }
  }

  return;
}

static void
aper_adj_quad(double angle, double x, double y, double* xquad, double* yquad)
{
  int quadrant;
  quadrant = angle / (pi/2) + 1;
  switch (quadrant)
    {
    case 1: *xquad=x;  *yquad=y;  break;
    case 2: *xquad=-x; *yquad=y;  break;
    case 3: *xquad=-x; *yquad=-y; break;
    case 4: *xquad=x;  *yquad=-y; break;
    }
}

static void
aper_adj_halo_si(double ex, double ey, double betx, double bety, double bbeat,
		 double halox[], double haloy[], int halolength, double haloxsi[], double haloysi[])
{
  int j;

  for (j=0; j <= halolength + 1; j++) {
    haloxsi[j] = halox[j] * bbeat * sqrt(ex*betx);
    haloysi[j] = haloy[j] * bbeat * sqrt(ey*bety);
  }
}

static int
aper_chk_inside(double p, double q, double pipex[], double pipey[], int pipelength )
{
// winding number test for a point in a polygon
// Input:   p,q = a point,
//          pipex[], pipey[] = vertex points of the polygon with pipe[len]=pipe[0]
// Return:  wn = the winding number (=0 only when point is outside polygon)
// source : wn_PnInPoly() at http://geomalgorithms.com/a03-_inclusion.html
 int    wn = 0;    // the  winding number counter

 // loop through all edges of the polygon
 for (int i=0; i<=pipelength; i++) { // edge from V[i] to  V[i+1]
   if (pipey[i] <= q  &&  pipey[i+1]  > q) {
     // first vertex is below point; second vertex is above; upward crossing
     if ( (pipex[i+1] - pipex[i]) * (q - pipey[i]) - (p - pipex[i]) * (pipey[i+1] - pipey[i])  > 0 ) {
       // Point left of  edge
       ++wn;
       continue;
     }
   }
   if (pipey[i] > q  &&  pipey[i+1]  <= q) {
     // first vertex is above point; second vertex is below; downward crossing
     if ( (pipex[i+1] - pipex[i]) * (q - pipey[i]) - (p - pipex[i]) * (pipey[i+1] - pipey[i])  < 0) {
       // Point right of  edge
       --wn;
       continue;
     }
   }
 }

 return wn;
}

static int
aper_intersect(double a1, double b1, double c1, double a2, double b2, double c2,
		double *xm, double *ym)
{ // calculates the intersection point between two general equation lines
  // returns 0 if lines are parallel, 1 otherwise

  if (fabs(a2*b1 - a1*b2) < MIN_DOUBLE*1.e10) return 0;

  *xm = (b2*c1 - b1*c2) / (a2*b1 - a1*b2);
  *ym = (a1*c2 - a2*c1) / (a2*b1 - a1*b2);
  return 1;
}

static void
aper_linepar(double x1,double y1,double x2,double y2,double *a,double *b,double *c)
{ // general line equation, non reduced:  a*x + b*y + c = 0
  *a = y1 - y2;
  *b = x2 - x1;
  *c = (x1-x2)*y1 + (y2-y1)*x1;
  return;
}

static int
aper_on_line(double xm, double ym, double x1, double y1,
	      double x2, double y2, double dist_limit)
{
  // first check colinearity of (p1,p2) and (p1,M); return 0 if false
  if ( fabs((x2-x1)*(ym-y1) - (xm-x1)*(y2-y1)) > dist_limit ) return 0;
  // check M is between p1 and p2, with tolerance for point overlap
  if ( xm <= fmax(x1,x2)+dist_limit && xm >= fmin(x1,x2)-dist_limit &&
       ym <= fmax(y1,y2)+dist_limit && ym >= fmin(y1,y2)-dist_limit ) return 1;
  return 0; // outside of (p1,p2)
}

static void
aper_fill_quadrants(double polyx[], double polyy[], int quarterlength, int* halolength)
{/* 2013-03-21 ghislain: given the data for the upper right quadrant computed in aper_rectellipse
     and contained in the polyx and polyy tables,
     mirrors this data to the other three quadrants across the x and y axes.
     quarterlength is the length of, or number of points in, the first quadrant.*/

 /* 2015-Jul-21 ghislain: added guards to avoid duplicate points close to the axes.
     These duplicate points could cause problems in aper_calc since such a pair of
     points does not define a line segment to check intersection. This fixed a bug
     reported by R. De Maria (https://svnweb.cern.ch/trac/madx/ticket/346) where aperture
     calculations were wrong when halor == halov exactly. Changing halor by a very small
     amount was enough to get back to proper results because it pulled the pair of points
     slighlty apart.
     As implemented the check is for proximity to the axis within 1e-10 mm  which
     should be good enough. On top of this if proximity is detected the value is reset
     to exactly zero. */


  int i,j;
  int debug = get_option("debug");

  if (debug) printf("+++ aper_fill_quadrants: quarterlength = %d\n", quarterlength);

  // The counter i starts at quarterlength+1, ie the first point to be mirrored.
  i = quarterlength + 1;

  /*copying first quadrant coordinates to second quadrant*/
  for (j=quarterlength; j>=0; j--) {
    // 2015-Jul-21  16:45:36  ghislain: avoid duplicate points on axes
    if (polyx[j] < 1.e-10) {polyx[j] = 0. ; continue;} // do not mirror a point at x = 0
    polyx[i] = -polyx[j];
    polyy[i] =  polyy[j];
    i++;
  }

  /*copying first quadrant coordinates to third quadrant*/
  for (j=0; j<=quarterlength; j++) {
    // 2015-Jul-21  16:45:36  ghislain: avoid duplicate points on axes
    if (polyy[j] < 1.e-10) {polyy[j] = 0. ; continue;} // do not mirror a point at y = 0
    polyx[i] = -polyx[j];
    polyy[i] = -polyy[j];
    i++;
  }

  /*copying first quadrant coordinates to fourth quadrant*/
  for (j=quarterlength; j>=0; j--) {
    // 2015-Jul-21  16:45:36  ghislain: avoid duplicate points on axes
    if (polyx[j] < 1.e-10) {polyx[j] = 0. ; continue;} // do not mirror a point at x = 0
    polyx[i] =  polyx[j];
    polyy[i] = -polyy[j];
    i++;
  }

  /*sets the last point equal to the first, to complete the shape.
    Necessary for compatibility with aper_calc function*/
  // 2015-Jul-21  16:45:36  ghislain: avoid duplicate points on axes
  if ( polyy[0] < 1.e-10) { // do not mirror a point at y = 0
    i--;
  } else {
    polyx[i] = polyx[0];
    polyy[i] = polyy[0];
  }

  *halolength=i-1;

  if (debug) {
    for (j=0; j<=i; j++) printf("  %d  %10.5e  %10.5e \n", j, polyx[j], polyy[j]);
    printf("\n");
  }

  return;
}

static void
aper_read_twiss(const char* table, int* jslice, double* s, double* x, double* y,
		double* px, double* py,
                double* betx, double* bety, double* dx, double* dy)
{
  double_from_table_row(table, "s", jslice, s);
  double_from_table_row(table, "x", jslice, x);
  double_from_table_row(table, "y", jslice, y);
  double_from_table_row(table, "px", jslice, px);
  double_from_table_row(table, "py", jslice, py);
  double_from_table_row(table, "betx", jslice, betx);
  double_from_table_row(table, "bety", jslice, bety);
  double_from_table_row(table, "dx", jslice, dx);
  double_from_table_row(table, "dy", jslice, dy);

  return;
}

static int
aper_external_file(char *file, double tablex[], double tabley[])
{ /* receives the name of file containing coordinates. Puts coordinates into tables. */
  /* returns status -1 in case of error and number of points read otherwise */
  int i=0;
  FILE *filept;

  /* no file provided */
  if (file == NULL) return -1;

  /* file cannot be opened in read mode */
  if ((filept=fopen(file, "r")) == NULL) {
    warning("Can not open file: ", file);
    return -1;
  }

  /* build table */
  while (2==fscanf(filept, "%lf %lf", &tablex[i], &tabley[i])) {
    i++;
    if (i >= MAXARRAY) {
      // should give the value of the MAXARRAY set limit that is exceeded
      // warn("Memory full in aper_external_file; number of coordinates exceeds MAXARRAY = ", MAXARRAY)
      fatal_error("Memory full in aper_external_file. ", "Number of coordinates exceeds set limit");
    }
  }

  /* closing the shape: a last point is inserted in table
     with coordinates equal to those of the first point */
    tablex[i]=tablex[0];
    tabley[i]=tabley[0];
    fclose(filept);

    return i-1;
}

static int
aper_build_screen(char* apertype, double* ap1, double* ap2, double* ap3, double* ap4, int* pipelength, double pipex[], double pipey[])
{
  int i, quarterlength=0;

  /* 2013-03-21 -- ghislain: changed name from aper_bs;
     the same function is referenced as build_pipe in the documentation of Ivar Waarum */

  /* "var1 .. 4" represents values in the aperture array of each element  */
  /*  After they are read:                                                */
  /* *ap1 = half width rectangle                                          */
  /* *ap2 = half height rectangle                                         */
  /* *ap3 = half horizontal axis ellipse                                  */
  /* *ap4 = half vertical axis ellipse                                    */
  /*      returns 1 on success, 0 on failure          */

  // 2013-Apr-18  14:23:40  ghislain: added check for invalid values.
  // 2013-Apr-18  14:25:14  ghislain: added RECTCIRCLE type
  // 2014-Jun-27  11:14:27  ghislain: fixed bug in check of invalid values for RACETRACK

  (*ap1) = (*ap2) = (*ap3) = (*ap4) = 0;

  int debug = get_option("debug");
  if (debug)
    printf("+++ aper_build_screen; apertype = '%s' quarterlength = %d\n",apertype, quarterlength);

  if (!strcmp(apertype,"circle")) {
    *ap3 = get_aperture(current_node, "var1"); /*radius circle*/
    if ( (*ap3) <= 0. ) {
      if (debug)
	printf("+++ aper_build screen, circle parameters: %10.5f %10.5f %10.5f %10.5f  -- exiting 0\n", *ap1, *ap2, *ap3, *ap4);
      return 0;
    }
    // make a square just containing the circle
    *ap1 = *ap2 = *ap4 = *ap3;

    aper_rectellipse(ap1, ap2, ap3, ap4, &quarterlength, pipex, pipey);
    aper_fill_quadrants(pipex, pipey, quarterlength, pipelength);
    return 1;
  }

  else if (!strcmp(apertype,"ellipse")) {
    *ap3 = get_aperture(current_node, "var1"); /*half hor axis ellipse*/
    *ap4 = get_aperture(current_node, "var2"); /*half ver axis ellipse*/
    if ( (*ap3) <= 0 || (*ap4) <= 0) {
      if (debug)
	printf("+++ aper_build screen, ellipse parameters: %10.5f %10.5f %10.5f %10.5f  -- exiting 0\n", *ap1, *ap2, *ap3, *ap4);
      return 0;
    }
    // make a rectangle just containing the ellipse
    *ap1 = *ap3;   *ap2 = *ap4;

    aper_rectellipse(ap1, ap2, ap3, ap4, &quarterlength, pipex, pipey);
    aper_fill_quadrants(pipex, pipey, quarterlength, pipelength);
    return 1;
  }

  else if (!strcmp(apertype,"rectangle")) {
    *ap1 = get_aperture(current_node, "var1");      /*half width rect*/
    *ap2 = get_aperture(current_node, "var2");      /*half height rect*/
    if ( (*ap1) <= 0 || (*ap2) <= 0) {
      if (debug)
	printf("+++ aper_build screen, rectangle parameters: %10.5f %10.5f %10.5f %10.5f  -- exiting 0\n", *ap1, *ap2, *ap3, *ap4);
      return 0;
    }
    // make a circle containing the rectangle
    *ap3 = *ap4 = sqrt( (*ap1)*(*ap1) + (*ap2)*(*ap2) );

    aper_rectellipse(ap1, ap2, ap3, ap4, &quarterlength, pipex, pipey);
    aper_fill_quadrants(pipex, pipey, quarterlength, pipelength);
    return 1;
  }

  else if (!strcmp(apertype,"lhcscreen") || !strcmp(apertype, "rectcircle")) {
    // the type lhcscreen should be deprecated at some point to keep MAD agnostic...
    *ap1=get_aperture(current_node, "var1"); /*half width rect*/
    *ap2=get_aperture(current_node, "var2"); /*half height rect*/
    *ap3=get_aperture(current_node, "var3"); /*radius circle*/

    if ( (*ap1) <= 0 || (*ap2) <= 0 || (*ap3) <= 0) {
      if (debug)
	printf("+++ aper_build screen, rectcircle parameters: %10.5f %10.5f %10.5f %10.5f  -- exiting 0\n", *ap1, *ap2, *ap3, *ap4);
      return 0;
    }
    // ensure the ellipse is a circle
    *ap4 = *ap3;

    aper_rectellipse(ap1, ap2, ap3, ap4, &quarterlength, pipex, pipey);
    aper_fill_quadrants(pipex, pipey, quarterlength, pipelength);
    return 1;
  }

  else if (!strcmp(apertype,"rectellipse")) {
    *ap1=get_aperture(current_node, "var1"); /*half width rect*/
    *ap2=get_aperture(current_node, "var2"); /*half height rect*/
    *ap3=get_aperture(current_node, "var3"); /*half hor axis ellipse*/
    *ap4=get_aperture(current_node, "var4"); /*half ver axis ellipse*/

    if ( (*ap1) <= 0 || (*ap2) <= 0 || (*ap3) <= 0 || (*ap4) <= 0) {
      if (debug)
	printf("+++ aper_build screen, rectellipse parameters: %10.5f %10.5f %10.5f %10.5f  -- exiting 0\n", *ap1, *ap2, *ap3, *ap4);
      return 0;
    }
    aper_rectellipse(ap1, ap2, ap3, ap4, &quarterlength, pipex, pipey);
    aper_fill_quadrants(pipex, pipey, quarterlength, pipelength);
    return 1;
  }

  else if (!strcmp(apertype,"racetrack")) {
    *ap1=get_aperture(current_node, "var1"); /*half width extension*/
    *ap2=get_aperture(current_node, "var2"); /*half height extension*/
    *ap3=get_aperture(current_node, "var3"); /*horizontal semi-axis*/
    *ap4=get_aperture(current_node, "var4"); /*vertical semi-axis*/

    // 2014-Jun-27  11:14:27  ghislain:
    // change check from ap1 or ap2<=0  to ap1 or ap2 or ap3 < 0
    // zero horizontal or vertical explosion factors, and zero radius should be allowed.
    // 2015-Mar-09  14:43:33  ghislain: change meaning of parameters: ap1 and ap2 are now full rectangle extension,
    // 2015-Mar-10  10:18:46  ghislain: change rounded corners from circular to generalized elliptical shape
    if ( (*ap1) < 0 || (*ap2) < 0 || (*ap3) < 0 || (*ap4) < 0 || (*ap1) < (*ap3) || (*ap2) < (*ap4)) {
      if (debug)
	printf("+++ aper_build screen, racetrack parameters: %10.5f %10.5f %10.5f %10.5f  -- exiting 0\n", *ap1, *ap2, *ap3, *ap4);
      return 0;
    }
    // special call to build an ellipse first: note that we cannot invoque ap1 or ap2
    aper_rectellipse(ap3, ap4, ap3, ap4, &quarterlength, pipex, pipey);

    /* displace the quartercircle */
    for (i=0;i<=quarterlength;i++) {
      pipex[i] += (*ap1) - (*ap3);
      pipey[i] += (*ap2) - (*ap4);
    }

    aper_fill_quadrants(pipex, pipey, quarterlength, pipelength);
    return 1;
  }

  // 2015-Feb-05  18:35:01  ghislain: adding octagon aperture type
  else if (!strcmp(apertype,"octagon")) {
    *ap1=get_aperture(current_node, "var1"); /*half width rect ; >= 0 */
    *ap2=get_aperture(current_node, "var2"); /*half height rect ; >=0 */
    *ap3=get_aperture(current_node, "var3"); /*first angle ; >=0, <=pi/2 */
    *ap4=get_aperture(current_node, "var4"); /*second angle ; >=0, <=pi/2, >= *ap3*/

    if ( (*ap1) < 0 || (*ap2) < 0 || (*ap3) < 0 || (*ap4) < 0 || (*ap3) > pi/2 || (*ap4) > pi/2 || (*ap4) < (*ap3))  {
      if (debug)
	printf("+++ aper_build screen, octagon parameters: %10.5f %10.5f %10.5f %10.5f  -- exiting 0\n", *ap1, *ap2, *ap3, *ap4);
      return 0;
    }

    quarterlength = 1;
    pipex[0] = (*ap1);                      pipey[0] = (*ap1) * tan((*ap3));
    pipex[1] = (*ap2) * tan(pi/2 - (*ap4)); pipey[1] = (*ap2);

    aper_fill_quadrants(pipex, pipey, quarterlength, pipelength);
    return 1;
  }

  else if (strlen(apertype)) {
    // general case, assume the given type is a filename
    *pipelength = aper_external_file(apertype, pipex, pipey);
    *ap1 = *ap2 = *ap3 = *ap4 = 0;
    if (*pipelength > -1) return 1; else return 0;
  }

  *pipelength = -1;
  return 0;

}

static int
aper_tab_search(int cnt, struct aper_e_d* tab, char* name, int* pos)
{
  /* looks for node *name in tab[], returns 1 if found, and its pos */
  int i=-1, found=0;

  while (i < cnt && found == 0) {
    i++;
    if (strcmp(name,tab[i].name) == 0) found=1;
  }
  *pos=i;

  return found;
}

static int
aper_tab_search_tfs(struct table* t, char* name, double* row)
{
  /* looks for node *name in tab[], returns 1 if found, and its pos
     NAME                        S_IP         X_OFF        DX_OFF       DDX_OFF         Y_OFF        DY_OFF       DDY_OFF
     function value return:
     1 ok
     0  not found
     -1 column does not exist
  */

  int i=0, found=0;
  int name_pos, s_ip_pos, x_off_pos, dx_off_pos, ddx_off_pos, y_off_pos, dy_off_pos, ddy_off_pos;

  /* get column positions */
  mycpy(c_dum->c, "name");     if ((name_pos = name_list_pos(c_dum->c, t->columns)) < 0) return -1;
  mycpy(c_dum->c, "s_ip");     if ((s_ip_pos = name_list_pos(c_dum->c, t->columns)) < 0) return -1;
  mycpy(c_dum->c, "x_off");    if ((x_off_pos = name_list_pos(c_dum->c, t->columns)) < 0) return -1;
  mycpy(c_dum->c, "dx_off");   if ((dx_off_pos = name_list_pos(c_dum->c, t->columns)) < 0) return -1;
  mycpy(c_dum->c, "ddx_off");  if ((ddx_off_pos = name_list_pos(c_dum->c, t->columns)) < 0) return -1;
  mycpy(c_dum->c, "y_off");    if ((y_off_pos = name_list_pos(c_dum->c, t->columns)) < 0) return -1;
  mycpy(c_dum->c, "dy_off");   if ((dy_off_pos = name_list_pos(c_dum->c, t->columns)) < 0) return -1;
  mycpy(c_dum->c, "ddy_off");  if ((ddy_off_pos = name_list_pos(c_dum->c, t->columns)) < 0) return -1;

  while (i < t->curr && found == 0) {
    i++;
    if( !strcmp(t->s_cols[name_pos][i-1],name)) {
      row[1] = t->d_cols[s_ip_pos][i-1];
      row[2] = t->d_cols[x_off_pos][i-1];
      row[3] = t->d_cols[dx_off_pos][i-1];
      row[4] = t->d_cols[ddx_off_pos][i-1];
      row[5] = t->d_cols[y_off_pos][i-1];
      row[6] = t->d_cols[dy_off_pos][i-1];
      row[7] = t->d_cols[ddy_off_pos][i-1];
      found = 1;
    }
  }

  return found;
}

static int
aper_e_d_read(char* e_d_name, struct aper_e_d** e_d_tabp, int* cnt, char* refnode)
{
  /* Reads data for special displacements of some magnets */
  int i=1, j, k, l, e_d_flag=0, curr_e_d_max = E_D_LIST_CHUNK, new_e_d_max;
  char comment[100]="empty";
  char *strpt;
  FILE *e_d_pt;
  struct aper_e_d* e_d_tab_loc;
  struct aper_e_d* e_d_tab = *e_d_tabp;


  if (e_d_name != NULL) {
    if((e_d_pt = fopen(e_d_name,"r")) == NULL) {
      printf("\nFile does not exist: %s\n",e_d_name);
    }
    else {
      /* part for reading reference node */
      while (strncmp(comment,"reference:",10) && i != EOF) {
        /*fgets(buf, 100, e_d_pt);*/
        i = fscanf(e_d_pt, "%s", comment);
        stolower(comment);
      }

      if (i == EOF) rewind(e_d_pt);
      else {
        if (strlen(comment) != 10) {
          strpt=strchr(comment,':');
          strpt++;
          strcpy(refnode, strpt);
        }
        else i = fscanf(e_d_pt, "%s", refnode);

        stolower(refnode);
        strcat(refnode, ":1");
      }
      /* end reading reference node */

      i=0;

      while (i != EOF) {
        i=fscanf(e_d_pt, "%s", e_d_tab[*cnt].name);
        /*next while-loop treats comments*/
        while ( e_d_tab[*cnt].name[0] == '!' && i != EOF) {
          if (fgets(comment, 100, e_d_pt) == NULL) { /* discard */ ; }
          i=fscanf(e_d_pt, "%s", e_d_tab[*cnt].name);
        }

        stolower(e_d_tab[*cnt].name);

        if (i != EOF) {
          strcat(e_d_tab[*cnt].name, ":1");

          k=0; j=3;
          while (j == 3 && k < E_D_MAX) {
            j=fscanf(e_d_pt, "%lf %lf %lf", &e_d_tab[*cnt].tab[k][0],
                     &e_d_tab[*cnt].tab[k][1],
                     &e_d_tab[*cnt].tab[k][2]);
            k++;

            if (e_d_tab[*cnt].curr == E_D_MAX) printf("\nToo many points of x,y displacement...\n");
          }

          e_d_tab[*cnt].curr=k-2;

          ++*cnt;

          if (*cnt == curr_e_d_max) { /* grow e_d array */
            /* printf("\nToo many special elements...(less than %d expected)\n", E_D_MAX); */
            new_e_d_max = curr_e_d_max + E_D_LIST_CHUNK;
            printf("\ngrowing e_d_max array to %d\n", new_e_d_max);
            e_d_tab_loc = mycalloc("Aperture", new_e_d_max, sizeof *e_d_tab_loc);

            for( l=0 ; l < curr_e_d_max; l++)
              e_d_tab_loc[l] = e_d_tab[l];

            myfree("Aperture", e_d_tab);

            e_d_tab = e_d_tab_loc;
            curr_e_d_max = new_e_d_max;
          }

          i=j;
        }
      } /* while !EOF */

      printf("\nUsing extra displacements from file \"%s\"\n",e_d_name);
      e_d_flag=1; fclose(e_d_pt);
      --*cnt;
    }
  }

  *e_d_tabp = e_d_tab;

  return e_d_flag;
}

static struct table*
aper_e_d_read_tfs(char* e_d_name, int* cnt, char* refnode)
{
  /* Reads displacement data in tfs format */
  struct table* t = NULL;
  struct char_p_array* tcpa = NULL;
  struct name_list* tnl = NULL;
  int i, k, error = 0;
  short  sk;
  char *cc, *tmp, *name;

  int tempcount = 0;

  (void)cnt;

  if (e_d_name == NULL) return NULL;

  printf("\n Reading offsets from tfs \"%s\"\n",e_d_name);

  if ((tab_file = fopen(e_d_name, "r")) == NULL) {
    warning("cannot open file:", e_d_name); return NULL;
  }

  const char* sep=" \"\n";

  while (fgets(aux_buff->c, aux_buff->max, tab_file))  {
    tempcount++;

    cc = strtok(aux_buff->c, sep);

    if (*cc == '@') {
      if ((tmp = strtok(NULL, sep)) != NULL
	  && strcmp(tmp, "REFERENCE") == 0) { /* search for reference node */
	if ((name = strtok(NULL, sep)) != NULL) { /* skip format */
	  if ((name = strtok(NULL, sep)) != NULL) {
	    strcpy(refnode, name);
	    stolower(refnode);
	    strcat(refnode, ":1");
	  }
	}
      }
    }

    else if (*cc == '*' && tnl == NULL) { /* search for column names and register them*/
      tnl = new_name_list("table_names", 20);
      while ((tmp = strtok(NULL, sep)) != NULL)
	add_to_name_list(permbuff(stolower(tmp)), 0, tnl);
    }

    else if (*cc == '$' && tcpa == NULL) { /* search for formats */
      if (tnl == NULL) {
	warning("formats before names","skipped"); return NULL;
      }
      tcpa = new_char_p_array(20);
      while ((tmp = strtok(NULL, sep)) != NULL) {
	if (tcpa->curr == tcpa->max) grow_char_p_array(tcpa);
	if (strcmp(tmp, "%s") == 0)       tnl->inform[tcpa->curr] = 3;
	else if (strcmp(tmp, "%hd") == 0) tnl->inform[tcpa->curr] = 1;
	else if (strcmp(tmp, "%d") == 0)  tnl->inform[tcpa->curr] = 1;
	else                              tnl->inform[tcpa->curr] = 2;
	tcpa->p[tcpa->curr++] = permbuff(tmp);
      }
    }

    else {
      if(t == NULL) {
	if (tcpa == NULL) {
	  warning("TFS table without formats,","skipped"); error = 1;
	}
	else if (tnl == NULL) {
	  warning("TFS table without column names,","skipped"); error = 1;
	}
	else if (tnl->curr == 0) {
	  warning("TFS table: empty column name list,","skipped");
	  error = 1;
	}
	else if (tnl->curr != tcpa->curr) {
	  warning("TFS table: number of names and formats differ,", "skipped");
	  error = 1;
	}
	if (error) {
	  delete_name_list(tnl); return NULL;
	}
	if(e_d_name != NULL) {
	  t = new_table(e_d_name, "OFFSETS",    500, tnl);
	} else {
	  t = new_table(e_d_name, "OFFSETS",    500, tnl);
	}
	t->curr = 0;
      }

      for (i = 0; i < tnl->curr; i++) {
	if (t->curr == t->max) grow_table(t);
	tmp = tcpa->p[i];
	if (strcmp(tmp,"%s") == 0)  {
	  char buf[strlen(cc)+3]; // reserve enough space for strcat
	  t->s_cols[i][t->curr] = tmpbuff(stolower(strcat(strcpy(buf, cc), ":1")));
	}
	else if (strcmp(tmp,"%d") == 0 ) {
	  sscanf(cc, tmp, &k); t->d_cols[i][t->curr] = k;
	}
	else if (strcmp(tmp,"%hd") == 0 ) {
	  sscanf(cc, tmp, &sk); t->d_cols[i][t->curr] = sk;
	}
	else sscanf(cc, tmp, &t->d_cols[i][t->curr]);
	if (i+1 < tnl->curr) {
	  if ((cc =strtok(NULL, sep)) == NULL) {
	    warning("incomplete table line starting with:", aux_buff->c);
	    return NULL;
	  }
	}
      }
      t->curr++;
    }
  }

  fclose(tab_file);
  t->origin = 1;
  /*  next line commented : avoid memory error at 2nd APERTURE command */
  /*  when the offset file has the same same, BJ 8apr2009 */
  /*  add_to_table_list(t, table_register);*/
  return t;
}

static void
aper_header(struct table* aper_t, struct aper_node *lim)
  /* puts beam and aperture parameters at start of the aperture table */
{
  int i, h_length = 26; // not used, nint=1;
  double vtmp[4]; // not used, deltap_twiss;
  char tmp[NAME_L], name[NAME_L], *stmp;

  strncpy(name, lim->name, sizeof name);

  /* =================================================================*/
  /* ATTENTION: if you add header lines, augment h_length accordingly */
  /* =================================================================*/


  /* many modif to make the header being standard; BJ 25feb2008 */

  if (aper_t == NULL) return;
  stmp = command_par_string("pipefile", this_cmd->clone);
  if (stmp) h_length++;
  stmp = command_par_string("halofile", this_cmd->clone);
  if (stmp) h_length += 1; else h_length += 4;

//  printf("\nheader %d \n",h_length);

  /* beam properties */
  if (aper_t->header == NULL)  aper_t->header = new_char_p_array(h_length);
  strncpy(tmp, current_sequ->name, sizeof tmp);
  table_add_header(aper_t, "@ SEQUENCE         %%%02ds \"%s\"", strlen(tmp),stoupper(tmp));
  i = get_string("beam", "particle", tmp);
  table_add_header(aper_t, "@ PARTICLE         %%%02ds \"%s\"",i,stoupper(tmp));
  table_add_header(aper_t, "@ MASS             %%le  %F", get_value("beam", "mass"));
  table_add_header(aper_t, "@ ENERGY           %%le  %F", get_value("beam", "energy"));
  table_add_header(aper_t, "@ PC               %%le  %F", get_value("beam", "pc"));
  table_add_header(aper_t, "@ GAMMA            %%le  %F", get_value("beam", "gamma"));
  // 2015-Mar-03  12:05:58  ghislain: added
  table_add_header(aper_t, "@ BETA             %%le  %F", get_value("beam", "beta"));
  // end addition
  // 2015-Mar-11  15:30:15  ghislain: changed to get emittances from BEAM command, not from input parameters.
  table_add_header(aper_t, "@ EXN              %%le  %G", get_value("beam", "exn"));
  table_add_header(aper_t, "@ EYN              %%le  %G", get_value("beam", "eyn"));
  table_add_header(aper_t, "@ DQF              %%le  %F", command_par_value("dqf", this_cmd->clone));
  table_add_header(aper_t, "@ BETAQFX          %%le  %F", command_par_value("betaqfx", this_cmd->clone));
  table_add_header(aper_t, "@ PARAS_DX         %%le    %g", command_par_value("dparx", this_cmd->clone));
  table_add_header(aper_t, "@ PARAS_DY         %%le    %g", command_par_value("dpary", this_cmd->clone));
  table_add_header(aper_t, "@ DP_BUCKET_SIZE   %%le  %F", command_par_value("dp", this_cmd->clone));

  // LD: summary table is corrupted by embedded twiss, use saved value
  // double_from_table_row("summ","deltap",&nint,&deltap_twiss);
  table_add_header(aper_t, "@ TWISS_DELTAP     %%le  %F", lim->deltap_twiss);

  table_add_header(aper_t, "@ CO_RADIUS        %%le  %F", command_par_value("cor", this_cmd->clone));
  table_add_header(aper_t, "@ BETA_BEATING     %%le  %F", command_par_value("bbeat", this_cmd->clone));
  table_add_header(aper_t, "@ NB_OF_ANGLES     %%d       %g", command_par_value("nco", this_cmd->clone)*4);

  /* if a filename with halo coordinates is given, need not show halo */
  stmp = command_par_string("halofile", this_cmd->clone);
  if (stmp)
  {
    strncpy(tmp, stmp, sizeof tmp);
    table_add_header(aper_t, "@ HALOFILE         %%%02ds \"%s\"",strlen(tmp),stoupper(tmp));
  }
  else
  {
    i = command_par_vector("halo", this_cmd->clone, vtmp);
    table_add_header(aper_t, "@ HALO_PRIM        %%le       %g",vtmp[0]);
    table_add_header(aper_t, "@ HALO_R           %%le       %g",vtmp[1]);
    table_add_header(aper_t, "@ HALO_H           %%le       %g",vtmp[2]);
    table_add_header(aper_t, "@ HALO_V           %%le       %g",vtmp[3]);
  }
  /* show filename with pipe coordinates if given */
  stmp = command_par_string("pipefile", this_cmd->clone);
  if (stmp)
  {
    strncpy(tmp, stmp, sizeof tmp);
    table_add_header(aper_t, "@ PIPEFILE         %%%02ds \"%s\"", strlen(tmp), stoupper(tmp) );
  }

  /* 2013-Nov-18  13:49:45  ghislain: changed n1min and at_element strings to uppercase in output to TFS */
  table_add_header(aper_t, "@ N1MIN            %%le   %g", lim->n1);
  set_value("beam","n1min",&lim->n1);
  table_add_header(aper_t, "@ AT_ELEMENT       %%%02ds  \"%s\"", strlen(name), stoupper(name));
}

static void
aper_surv(double init[], int nint)
{
  struct in_cmd* aper_survey;
  struct name_list* asnl;
  int aspos;

  /* Constructs artificial survey command, the result is the  */
  /* table 'survey' which can be accessed from all functions. */
  /* init[0] = x0, init[1] = y0, init[2] = z0,                */
  /* init[3] = theta0, init[4] = phi0, init[5] = psi0         */

  aper_survey = new_in_cmd(10);
  aper_survey->type = 0;
  aper_survey->clone = aper_survey->cmd_def = clone_command(find_command("survey",defined_commands));
  asnl = aper_survey->cmd_def->par_names;
  aspos = name_list_pos("table", asnl);
  aper_survey->cmd_def->par->parameters[aspos]->string = permbuff("survey");
  aper_survey->cmd_def->par_names->inform[aspos] = 1;

  aspos = name_list_pos("x0", asnl);
  aper_survey->cmd_def->par->parameters[aspos]->double_value = init[0];
  aper_survey->cmd_def->par_names->inform[aspos] = 1;

  aspos = name_list_pos("y0", asnl);
  aper_survey->cmd_def->par->parameters[aspos]->double_value = init[1];
  aper_survey->cmd_def->par_names->inform[aspos] = 1;

  aspos = name_list_pos("z0", asnl);
  aper_survey->cmd_def->par->parameters[aspos]->double_value = init[2];
  aper_survey->cmd_def->par_names->inform[aspos] = 1;

  aspos = name_list_pos("theta0", asnl);
  aper_survey->cmd_def->par->parameters[aspos]->double_value = init[3];
  aper_survey->cmd_def->par_names->inform[aspos] = 1;

  aspos = name_list_pos("phi0", asnl);
  aper_survey->cmd_def->par->parameters[aspos]->double_value = init[4];
  aper_survey->cmd_def->par_names->inform[aspos] = 1;

  aspos = name_list_pos("psi0", asnl);
  aper_survey->cmd_def->par->parameters[aspos]->double_value = init[5];
  aper_survey->cmd_def->par_names->inform[aspos] = 1;

/* frs: suppressing the survey file created by the internal survey command */
  aspos = name_list_pos("file", asnl);
  aper_survey->cmd_def->par->parameters[aspos]->string = NULL;
  aper_survey->cmd_def->par_names->inform[aspos] = 0;

  current_survey=aper_survey->clone;
  pro_survey(aper_survey);

  double_from_table_row("survey","x",&nint, &init[0]);
  double_from_table_row("survey","y",&nint, &init[1]);
  double_from_table_row("survey","z",&nint, &init[2]);
  double_from_table_row("survey","theta",&nint, &init[3]);
  double_from_table_row("survey","phi",&nint, &init[4]);
  double_from_table_row("survey","psi",&nint, &init[5]);
}

static void
aper_trim_ws(char* string, int len)
{
  /* Replaces the first ws or : in a string with a '\0', */
  /* thus translating a FORTRAN-like attribute string to */
  /* C compatibility, or washes the ':1' from node names */
  int c=0;
  while (string[c] && string[c]!=' ' && c<len-1) c++;
  string[c]=0;
  if (c<len-1) string[c+1]=' '; /*adds a ws to avoid two \0 in a row*/
}

static void
aper_write_table(char* name, double* n1, double* n1x_m, double* n1y_m,
		 double* rtol, double* xtol, double* ytol, double* xoffset, double* yoffset,
		 char* apertype,double* ap1,double* ap2,double* ap3,double* ap4,
		 double* on_ap, double* on_elem, double* spec,double* s,
		 double* x, double* y, double* px, double* py,
		 double* betx, double* bety,double* dx, double* dy, char *table)
{
  string_to_table_curr(table, "name", name);
  double_to_table_curr(table, "n1", n1);
  double_to_table_curr(table, "n1x_m", n1x_m);
  double_to_table_curr(table, "n1y_m", n1y_m);
  double_to_table_curr(table, "rtol", rtol);
  double_to_table_curr(table, "xtol", xtol);
  double_to_table_curr(table, "ytol", ytol);
  double_to_table_curr(table, "xoffset", xoffset);
  double_to_table_curr(table, "yoffset", yoffset);
  string_to_table_curr(table, "apertype", apertype);
  double_to_table_curr(table, "aper_1", ap1);
  double_to_table_curr(table, "aper_2", ap2);
  double_to_table_curr(table, "aper_3", ap3);
  double_to_table_curr(table, "aper_4", ap4);
  double_to_table_curr(table, "on_ap", on_ap);
  double_to_table_curr(table, "on_elem", on_elem);
  double_to_table_curr(table, "spec", spec);
  double_to_table_curr(table, "s", s);
  double_to_table_curr(table, "x", x);
  double_to_table_curr(table, "y", y);
  double_to_table_curr(table, "px", px);
  double_to_table_curr(table, "py", py);
  double_to_table_curr(table, "betx", betx);
  double_to_table_curr(table, "bety", bety);
  double_to_table_curr(table, "dx", dx);
  double_to_table_curr(table, "dy", dy);

  augment_count(table);
}


static double
aper_calc(double p, double q, double* minhl,
	  double halox[], double haloy[], int halolength, double haloxadj[],double haloyadj[],
	  double pipex[], double pipey[], int pipelength, double notsimple)
{ // 2015-Jul-30  ghislain: partially rewritten from original to limit number of calculations
  // in loop to find ratios.
  int i=0, j=0, c=0;
  double dist_limit = 1.e-10;
  double a1, b1, c1, a2, b2, c2;
  double xm, ym, h, l;

  double tmphalox[MAXARRAY], tmphaloy[MAXARRAY];

  double bbxmax, bbymax, bbxmin, bbymin; // rectangular bounding box for aperture

  for (c=0; c <= halolength + 1; c++) {
    haloxadj[c] = halox[c] + p;
    haloyadj[c] = haloy[c] + q;
  }

  /* if halo centre is outside of beam pipe : zero aperture! */
  if (!aper_chk_inside(p, q, pipex, pipey, pipelength)) {
    *minhl = 0;
    return -1;
  }

  /* determine a rectangular bounding box for the beampipe */
  bbxmin=pipex[0]; bbymin=pipey[0]; bbxmax=pipex[0]; bbymax=pipey[0];
  for (j=1; j<=pipelength; j++) {
    bbxmin = fmin(bbxmin,pipex[j]); bbxmax = fmax(bbxmax,pipex[j]);
    bbymin = fmin(bbymin,pipey[j]); bbymax = fmax(bbymax,pipey[j]);
  }
  bbxmax *= 1.01; bbymax *= 1.01; bbxmin *= 1.01; bbymin *= 1.01; // 1% clearance


  /* not simply connex beampipes... adding apexes to the halo */
  c = 0;
  if (notsimple) {
    for (j=0; j<=halolength; j++) {
      tmphalox[c] = haloxadj[j];
      tmphaloy[c] = haloyadj[j];
      c++;

      /* Find general line equation for line joining two consecutive halo apex points*/
      aper_linepar(haloxadj[j], haloyadj[j], haloxadj[j+1], haloyadj[j+1], &a2, &b2, &c2);

      for (i=0; i<=pipelength; i++) {
	/* Find general line equation for line joining halo centre to pipe apex point*/
	aper_linepar(p, q, pipex[i], pipey[i], &a1, &b1, &c1);

	/* find intersection coordinates; cycle if lines are parallel */
	if (0 == aper_intersect(a1, b1, c1, a2, b2, c2, &xm, &ym)) continue; // parallel lines
	/* intersection point must be on halo line segment */
	if (0 == aper_on_line(xm, ym, haloxadj[j], haloyadj[j], haloxadj[j+1], haloyadj[j+1], dist_limit)) continue;
	/* pipe apex must be between halo center and intersection point */
	if (0 == aper_on_line(pipex[i], pipey[i], p, q, xm, ym, dist_limit)) continue;

	tmphalox[c] = xm;
	tmphaloy[c] = ym;
	c++;
      }
    }

    halolength = c - 1;
    for (j=0; j <= halolength; j++) {
      haloxadj[j] = tmphalox[j];
      haloyadj[j] = tmphaloy[j];
    }
  } // notsimple

  /*Calculate smallest ratio:*/
  i=0; // start from first pipeline apex

  for (j=0; j <= halolength; j++) { // loop over all halo apexes

    /* Find general line equation for line joining halo centre to halo apex point*/
    aper_linepar(p, q, haloxadj[j], haloyadj[j], &a1, &b1, &c1);

    for (;;) { // infinite loop
      /* find general equation of current pipeline segement */
      aper_linepar(pipex[i], pipey[i], pipex[i+1], pipey[i+1], &a2, &b2, &c2);

      if ( // not parallel lines
	  0 != aper_intersect(a1, b1, c1, a2, b2, c2, &xm, &ym)
	  // intersection point is inside the bounding box
	  &&  ( xm < bbxmax && xm > bbxmin && ym < bbymax && ym > bbymin)
	  // intersection point is inside pipe line segment
	  && 0 != aper_on_line(xm, ym, pipex[i], pipey[i], pipex[i+1], pipey[i+1], dist_limit)
	  // halo center is not inside segement defined by halo apex and intersection point
	  && 0 != aper_on_line(p,q,haloxadj[j],haloyadj[j],xm,ym,dist_limit) )
	      break; // valid point is found

      if (++i == pipelength + 1) i = 0; // cycle through the pipeline
    }

    h = sqrt((xm-p)*(xm-p) + (ym-q)*(ym-q));
    l = sqrt((haloxadj[j]-p)*(haloxadj[j]-p) + (haloyadj[j]-q)*(haloyadj[j]-q));
    if (h/l < *minhl) *minhl = h/l;
  }

  return 0;
}


// public interface

double
get_apertol(struct node* node, const char* par)
  /* returns aper_tol parameter 'i' where i is integer at the end of par;
     e.g. aptol_1 gives i = 1 etc. (count starts at 1) */
{
  int i, k, n = strlen(par);
  double val = zero, vec[100];
  for (i = 0; i < n; i++)  if (isdigit(par[i])) break;
  if (i == n) return val;
  sscanf(&par[i], "%d", &k); k--;
  if ((n = element_vector(node->p_elem, "aper_tol", vec)) > k)  val = vec[k];
  return val;
}

double
get_aperture(struct node* node, const char* par)
  /* returns aperture parameter 'i' where i is integer at the end of par;
     e.g. aper_1 gives i = 1 etc. (count starts at 1) */
{
  int i, k, n = strlen(par);
  double val = zero, vec[100];
  for (i = 0; i < n; i++)  if (isdigit(par[i])) break;
  if (i == n) return val;
  sscanf(&par[i], "%d", &k); k--;
  if ((n = element_vector(node->p_elem, "aperture", vec)) > k)  val = vec[k];
  return val;
}

void
pro_aperture(struct in_cmd* cmd)
{
  struct node *use_range[2];
  struct table* tw_cp;
  char *file, *range, tw_name[NAME_L], table_ap[NAME_L]="aperture", *table=table_ap;
  int tw_cnt, rows;
  double interval;

  embedded_twiss_cmd = cmd;

  /* check for valid sequence, beam and Twiss table */
  if (current_sequ == NULL || current_sequ->ex_start == NULL || sequence_length(current_sequ) == 0)
    fatal_error("Aperture module - no active sequence", "");

  if (attach_beam(current_sequ) == 0)
    fatal_error("Aperture module - sequence without beam:", current_sequ->name);

  if (!sequ_check_valid_twiss(current_sequ)) {
    warning("Aperture module - no valid TWISS table present", "Aperture command ignored");
    return;
  }

  range = command_par_string("range", this_cmd->clone);
  if (get_ex_range(range, current_sequ, use_range) == 0) {
    warning("Illegal range.","Aperture command ignored");
    return;
  }

  current_node = use_range[0];

  /* navigate to starting point in Twiss table */
  tw_cp = current_sequ->tw_table;
  tw_cnt = 1; /* table starts at 1 */

  while (1) {
    string_from_table_row(tw_cp->name, "name", &tw_cnt, tw_name);
    aper_trim_ws(tw_name, NAME_L);

    if (!strcmp(tw_name,current_node->name)) break;

    if (++tw_cnt > tw_cp->curr) {
      warning("Aperture command ignored: could not find range start in Twiss table:", current_node->name);
      return;
    }
  }

  /* approximate # of needed rows in aperture table */
  interval = command_par_value("interval", this_cmd->clone);
  rows = current_sequ->n_nodes + 2 * (sequence_length(current_sequ)/interval);

  /* make empty aperture table */
  aperture_table = make_table(table, table, ap_table_cols, ap_table_types, rows);
  aperture_table->dynamic = 1;
  add_to_table_list(aperture_table, table_register);

  /* calculate apertures and fill table */
  struct aper_node limit_node = { "none", -1, -1, "none", {-1,-1,-1,-1}, {-1,-1,-1}, 0.0 };
  struct aper_node *limit_pt = &limit_node;

  limit_pt = aperture(table, use_range, tw_cp, &tw_cnt, limit_pt);

  if (limit_pt->n1 != -1) {
    printf("\n\nAPERTURE LIMIT: %s, n1: %g, at: %g\n\n", limit_pt->name, limit_pt->n1, limit_pt->s);

    aper_header(aperture_table, limit_pt);

    file = command_par_string("file", this_cmd->clone);
    if (file != NULL) out_table(table, aperture_table, file);

    // 2015-Jul-31  11:41:59  ghislain: aperture twiss file for output of twiss table ! not needed
    //if (strcmp(aptwfile,"dummy")) out_table(tw_cp->name, tw_cp, aptwfile);
  }
  else warning("Could not run aperture command.","Aperture command ignored");

  /* set pointer to updated Twiss table */
  current_sequ->tw_table=tw_cp;
}

static struct aper_node*
aperture(char *table, struct node* use_range[], struct table* tw_cp, int *tw_cnt, struct aper_node *lim_pt)
{
  int stop=0, nint=1, jslice=1, first, ap=1;
  int true_flag, true_node=0, offs_node=0, do_survey=0;
  int truepos=0, true_cnt=0, offs_cnt=0;
  int halo_q_length=1, halolength, pipelength, namelen=NAME_L, ntol;
  double surv_init[6]={0, 0, 0, 0, 0, 0};
  double surv_x=zero, surv_y=zero;
  double xa=0, xb=0, xc=0, ya=0, yb=0, yc=0;
  double on_ap=1, on_elem=0;
  double ex, ey;
  double dqf, betaqfx, dp, dparx, dpary;
  double cor, bbeat, nco, halo[4], interval, spec, notsimple;
  double s=0, x=0, y=0, px=0, py=0, betx=0, bety=0, dx=0, dy=0, ratio, n1, length; // nr not used
  double xeff=0,yeff=0;
  double n1x_m, n1y_m;
  double s_start, s_curr, s_end;
  double node_s=-1, node_n1=-1;
  double aper_tol[3], ap1, ap2, ap3, ap4;
  double dispx, dispy, tolx, toly;
  double dispxadj=0, dispyadj=0, coxadj, coyadj, tolxadj=0, tolyadj=0;
  double angle, dangle, deltax, deltay;
  double xshift, yshift, r;
  double halox[MAXARRAY], haloy[MAXARRAY], haloxsi[MAXARRAY], haloysi[MAXARRAY];
  double haloxadj[MAXARRAY], haloyadj[MAXARRAY];
  double pipex[MAXARRAY], pipey[MAXARRAY];
  double parxd,paryd;
  char *halofile, *truefile, *offsfile;
  char refnode[NAME_L]="";
  char *cmd_refnode;
  char apertype[NAME_L];
  char name[NAME_L];
  char tol_err_mess[80] = "";
  int code;

  // 2014-Sep-18  17:19:52  ghislain: attempt to read offset values from element attributes...
  double aper_offset[2], xoffset, yoffset;
  int noffset;

  struct node* rng_glob[2];
  // struct aper_node limit_node = {"none", -1, -1, "none", {-1,-1,-1,-1}, {-1,-1,-1}};
  // struct aper_node* lim_pt = &limit_node;

  int is_zero_len;

  int debug = get_option("debug");

  int rc_deltap=0;

  true_tab = mycalloc("Aperture", E_D_LIST_CHUNK, sizeof *true_tab);
  /* offs_tab = mycalloc("Aperture",E_D_LIST_CHUNK,sizeof *offs_tab); */

  printf("\nProcessing apertures from %s to %s...\n",use_range[0]->name,use_range[1]->name);

  /* read command parameters */
  halofile = command_par_string("halofile", this_cmd->clone);

  /* removed IW 240205 */
  /*  pipefile = command_par_string("pipefile", this_cmd->clone); */

  // 2015-Mar-03  17:25:49  ghislain: get geometric emittances from BEAM command
  ex = get_value("beam","ex");
  ey = get_value("beam","ey");

  dqf = command_par_value("dqf", this_cmd->clone);
  betaqfx = command_par_value("betaqfx", this_cmd->clone);
  dp = command_par_value("dp", this_cmd->clone);
  dparx = command_par_value("dparx", this_cmd->clone);
  dpary = command_par_value("dpary", this_cmd->clone);
  cor = command_par_value("cor", this_cmd->clone);
  bbeat = command_par_value("bbeat", this_cmd->clone);
  nco = command_par_value("nco", this_cmd->clone);
  command_par_vector("halo", this_cmd->clone, halo); // nhalopar = // not used
  interval = command_par_value("interval", this_cmd->clone);
  spec = command_par_value("spec", this_cmd->clone);
  notsimple = command_par_value("notsimple", this_cmd->clone);
  truefile = command_par_string("trueprofile", this_cmd->clone);
  offsfile = command_par_string("offsetelem", this_cmd->clone);

  cmd_refnode = command_par_string("refnode", this_cmd->clone);

  /* calculate delta angle */
  dangle = twopi/(nco*4);

  // 2013-Nov-13  16:20:29  ghislain:
  // attempt to extract maximum parameters from twiss summary table
  // double_from_table_row("summ","dxmax",&nint,&dqf);
  // printf ("+++++++ dqf from TWISS %12.6g\n",dqf);
  // TODO:  add some error checking.

  /* fetch deltap as set by user in the former TWISS command */
  /* will be used below for displacement associated to parasitic dipersion */

  if ( (rc_deltap = double_from_table_row("summ","deltap",&nint,&lim_pt->deltap_twiss)) != 0) {
    printf ("+++++++ deltap from TWISS could not be read; assume 0.0d0\n");
  }
  else printf ("+++++++ deltap from TWISS %12.6g\n",lim_pt->deltap_twiss);


  /* check if trueprofile and offsetelem files exist */
  true_flag = aper_e_d_read(truefile, &true_tab, &true_cnt, refnode);
  /* offs_flag = aper_e_d_read(offsfile, &offs_tab, &offs_cnt, refnode);*/
  offs_tab = aper_e_d_read_tfs(offsfile, &offs_cnt, refnode);

  if (cmd_refnode != NULL) {
      strcpy(refnode, cmd_refnode);
      strcat(refnode, ":1");
  }
  if (strcmp(refnode,"")) printf("\nreference node: %s\n",refnode);

  /* build halo polygon based on input ratio values or coordinates */
  if ((halolength = aper_external_file(halofile, halox, haloy)) > -1) {
    printf("\n Reading halo from file \"%s\"\n",halofile);
  }
  else if (aper_rectellipse(&halo[2], &halo[3], &halo[1], &halo[1], &halo_q_length, halox, haloy)) {
    warning("Not valid parameters for halo. ", "Unable to make polygon.");

    /* IA */
    myfree("Aperture",true_tab);
    myfree("Aperture",offs_tab);

    return lim_pt;
  }
  else aper_fill_quadrants(halox, haloy, halo_q_length, &halolength);

  /* check for externally given pipe polygon */
  /* changed this recently, IW 240205 */
  /*  pipelength = aper_external_file(pipefile, pipex, pipey);
      if ( pipelength > -1) ext_pipe=1; */

  /* get initial twiss parameters, from start of first element in range */
  aper_read_twiss(tw_cp->name, tw_cnt, &s_end, &x, &y, &px, &py, &betx, &bety, &dx, &dy);
  // LD: shift further results by one step (?) and finish outside the table
  //  (*tw_cnt)++;
  aper_adj_halo_si(ex, ey, betx, bety, bbeat, halox, haloy, halolength, haloxsi, haloysi);

  /* calculate initial normal+parasitic disp. */
  /* modified 27feb08 BJ */
  parxd = dparx * sqrt(betx/betaqfx) * dqf;
  paryd = dpary * sqrt(bety/betaqfx) * dqf;

  /* Initialize n1 limit value */
  lim_pt->n1=999999;

  int n = 0;

  while (!stop) { // loop over nodes
    ++n;

    strcpy(name,current_node->name);
    aper_trim_ws(name, NAME_L);

    is_zero_len = 0;

    /* the first node in a sequence can not be sliced, hence: */
    if (current_sequ->range_start == current_node) first=1; else first=0;

    double_from_table_row(tw_cp->name, "s", tw_cnt, &s_end);
    length  = node_value("l");
    s_start = s_end - length;
    s_curr  = s_start;

    node_string("apertype", apertype, &namelen);
    aper_trim_ws(apertype, NAME_L);

    if (!strncmp("drift",name,5)) on_elem=-999999;
    else on_elem=1;

    if ( (offs_tab != NULL) && (strcmp(refnode, name) == 0)) do_survey=1; // current name is refnode; switch survey on.

    if (debug) printf("\nname: %s, ref: %s, do_survey: %d, true_flag: %d\n",name,refnode,do_survey,true_flag);

    // 2015-Mar-19  09:07:37  ghislain:
    code = node_value("mad8_type");
    if (code == 20)
      warning("Found deprecated ECOLLIMATOR element;"," Should be replaced by COLLIMATOR");
    if (code == 21)
      warning("Found deprecated RCOLLIMATOR element;"," Should be replaced by COLLIMATOR");

    /* read data for tol displacement of halo */
    get_node_vector("aper_tol",&ntol,aper_tol);
    if (ntol == 3) {
      r = aper_tol[0];
      xshift = aper_tol[1];
      yshift = aper_tol[2];
    }
    else r = xshift = yshift = 0.;
    //if (debug) printf("\nname: %s, x-shift: %f, y-shift: %f\n",name,xshift, yshift);

    xoffset=yoffset=0;
    // 2014-Sep-18  17:19:52  ghislain: attempt to read offset values from element attributes...
    /* read data for aper_offset */
    get_node_vector("aper_offset", &noffset, aper_offset);
    if (noffset == 2) {
      xoffset = aper_offset[0];
      yoffset = aper_offset[1];
    }
    else xoffset=yoffset=0;
    if (debug) printf("name: %s, x-offset: %f, y-offset: %f\n", name, xoffset, yoffset);
    xoffset=yoffset=0; // 2015-Nov-20  15:43:45  ghislain: until clarified in MAD-X meeting

    /*read aperture data and make polygon tables for beam pipe*/
    /* IW 250205 */
    /*  if (ext_pipe == 0) */
    ap = aper_build_screen(apertype, &ap1, &ap2, &ap3, &ap4, &pipelength, pipex, pipey);

    if (ap == 0 || first == 1) {
      /* if no pipe can be built, the n1 is set to inf and Twiss parms read for reference*/
      n1=999999; n1x_m=999999; n1y_m=999999; on_ap=-999999; nint=1;

      aper_read_twiss(tw_cp->name, tw_cnt, &s_end, &x, &y, &px, &py, &betx, &bety, &dx, &dy);

      aper_write_table(name, &n1, &n1x_m, &n1y_m, &r, &xshift, &yshift, &xoffset, &yoffset,
		       apertype, &ap1, &ap2, &ap3, &ap4, &on_ap, &on_elem, &spec,
                       &s_end, &x, &y, &px, &py, &betx, &bety, &dx, &dy, table);
      on_ap=1;

      double_to_table_row(tw_cp->name, "n1", tw_cnt, &n1);
      (*tw_cnt)++;

      /* calc disp and adj halo to have ready for next node */
      /* modified 27feb08 BJ */
      parxd = dparx * sqrt(betx/betaqfx) * dqf;
      paryd = dpary * sqrt(bety/betaqfx) * dqf;

      aper_adj_halo_si(ex, ey, betx, bety, bbeat, halox, haloy, halolength, haloxsi, haloysi);

     /*do survey to have ready init for next node */
      if (do_survey) {
        rng_glob[0] = current_sequ->range_start;
        rng_glob[1] = current_sequ->range_end;
        current_sequ->range_start = current_sequ->range_end = current_node;
        aper_surv(surv_init, nint);
        double_from_table_row("survey","x",&nint, &surv_x);
        double_from_table_row("survey","y",&nint, &surv_y);
        current_sequ->range_start = rng_glob[0];
        current_sequ->range_end = rng_glob[1];
      }
    }    /* end loop 'if no pipe ' */

    else {
      node_n1   = 999999;
      true_node = 0;
      offs_node = 0;

      /* calculate the number of slices per node */
      if (true_flag == 0)
        nint = length / interval;
      else {
        true_node = aper_tab_search(true_cnt, true_tab, name, &truepos);
        if (true_node) nint = true_tab[truepos].curr;
        else nint = length / interval;
      }

      if (debug) printf("name: %s, nint: %d, length: %f\n", name, nint, length);

      if (!nint) nint = 1;

      /* do not interpolate 0-length elements*/
      if (fabs(length) < MIN_DOUBLE ) is_zero_len = 1;

      /* slice the node, call survey if necessary, make twiss for slices*/
      interpolate_node(&nint);

      /* do survey */
      if (do_survey) {
        double offs_row[8] = { 0 };

        aper_surv(surv_init, nint);

        offs_node = aper_tab_search_tfs(offs_tab, name, offs_row);
        if (offs_node) {
          xa=offs_row[4];
          xb=offs_row[3];
          xc=offs_row[2];
          ya=offs_row[7];
          yb=offs_row[6];
          yc=offs_row[5];
	  if (debug) {
	    printf("\nusing offset:");
	    printf("\n xa: %f, xb: %f, xc: %f \n ya: %f, yb: %f, yc: %f\n", xa, xb, xc, ya, yb, yc);
	  }
        }
      }

      embedded_twiss();

      /* Treat each slice, for all angles */
      for (jslice=0; jslice <= nint; jslice++) {
        ratio=999999;

        if (jslice == 0) {
          // parameters from previous node will be used
          s_curr += 1.e-12;     /*to get correct plot at start of elements*/
          s = 0;                /*used to calculate survey adjustments */
        }
	else {
          aper_read_twiss("embedded_twiss_table", &jslice, &s, &x, &y, &px, &py,
			  &betx, &bety, &dx, &dy);

	  if(debug) printf("embedded twiss for slice %d: s= %f betx= %f bety= %f dx= %f dy= %f\n",
			   jslice, s, betx, bety, dx, dy);

          s_curr = s_start + s;
          aper_adj_halo_si(ex, ey, betx, bety, bbeat, halox, haloy, halolength, haloxsi, haloysi);

          /* calculate normal+parasitic disp.*/
          /* modified 27feb08 BJ */
          parxd = dparx * sqrt(betx/betaqfx) * dqf;
          paryd = dpary * sqrt(bety/betaqfx) * dqf;

          if (do_survey) {
            double_from_table_row("survey","x",&jslice,&surv_x);
            double_from_table_row("survey","y",&jslice,&surv_y);
          }
        }

        /* BJ 3 APR 2009 : introduced xeff and yeff in order to avoid
           interferences with x and y as used for the first slice
	         (re-use from end of former node) */
        xeff = x;
        yeff = y;

	/* offset adjustment */
	if (noffset == 2) {
	  xeff -= xoffset;
	  yeff -= yoffset;
	}

        /* survey adjustments */
        if (offs_node) {
          xeff += surv_x - (xa*s*s + xb*s + xc);
          yeff += surv_y - (ya*s*s + yb*s + yc);
        }

        /* discrete adjustments */
        if (true_node) {
          xeff += true_tab[truepos].tab[jslice][1];
          yeff += true_tab[truepos].tab[jslice][2];
        }

	if (debug) printf("\n adjustments xeff: %f, yeff: %f\n",xeff,yeff);

        for (angle=0; angle<twopi; angle+=dangle) {
          /* new 27feb08 BJ */
          dispx = bbeat * (fabs(dx)*dp + parxd*(fabs(lim_pt->deltap_twiss)+dp) );
          dispy = bbeat * (fabs(dy)*dp + paryd*(fabs(lim_pt->deltap_twiss)+dp) );

          /*adjust dispersion to worst-case for quadrant*/
          aper_adj_quad(angle, dispx, dispy, &dispxadj, &dispyadj);

          /*calculate displacement co+tol for each angle*/
          coxadj = cor * cos(angle); coyadj = cor * sin(angle);

          /* Error check added 20feb08 BJ */
          if ( xshift < 0 || yshift < 0 || r < 0 ) {
            sprintf(tol_err_mess,"In element : %s\n",name);
            fatal_error("Illegal negative tolerance",tol_err_mess);
          }

          aper_race(xshift, yshift, r, angle, &tolx, &toly);
          aper_adj_quad(angle, tolx, toly, &tolxadj, &tolyadj);

          /* add all displacements */
          deltax = coxadj + tolxadj + xeff + dispxadj;
          deltay = coyadj + tolyadj + yeff + dispyadj;

          /* send beta adjusted halo and its displacement to aperture calculation */
          aper_calc(deltax, deltay, &ratio, haloxsi, haloysi, halolength, haloxadj, haloyadj,
                    pipex, pipey, pipelength, notsimple);

	  if (debug) printf("\n Angle: %f deltax: %f deltay: %f minratio: %f\n", angle, deltax, deltay, ratio);
        }

        //nr = ratio * halo[1];
        //n1 = nr / (halo[1]/halo[0]); /* ratio r/n = 1.4 */
	n1 = ratio * halo[0]; // 2015-Jul-30  17:23:26  ghislain: replaced above two lines.


	if (debug) printf("\n Found ratio: %f n1: %f \n",ratio,n1);

        n1x_m = n1 * bbeat * sqrt(betx*ex);
        n1y_m = n1 * bbeat * sqrt(bety*ey);

	/* Change below, BJ 23oct2008                              */
	/* test block 'if (n1 < node_n1)' included in test block   */

        if (is_zero_len == 0 || jslice == 1) {
          aper_write_table(name, &n1, &n1x_m, &n1y_m, &r, &xshift, &yshift, &xoffset, &yoffset,
			   apertype, &ap1, &ap2, &ap3, &ap4, &on_ap, &on_elem, &spec, &s_curr,
                           &xeff, &yeff, &px, &py, &betx, &bety, &dx, &dy, table);

	  /* save node minimum n1 */
          if (n1 < node_n1) {
              node_n1 = n1;
              node_s = s_curr;
          }
	} // if is_zero_len
      } // for jslice

      reset_interpolation();

      /* insert minimum node value into Twiss table */
      double_to_table_row(tw_cp->name, "n1", tw_cnt, &node_n1);
      (*tw_cnt)++;

      /* save range minimum n1 */
      if (node_n1 < lim_pt->n1) {
        strcpy(lim_pt->name,name);
        lim_pt->n1 = node_n1;
        lim_pt->s = node_s;
        strcpy(lim_pt->apertype,apertype);
        lim_pt->aperture[0] = ap1;
        lim_pt->aperture[1] = ap2;
        lim_pt->aperture[2] = ap3;
        lim_pt->aperture[3] = ap4;
        lim_pt->aper_tol[0] = r;
        lim_pt->aper_tol[1] = xshift;
        lim_pt->aper_tol[2] = yshift;
      }
    } // if !(ap == 0 || first == 1)

    if (!strcmp(current_node->name,use_range[1]->name)) stop=1;
    if (!advance_node()) stop=1;
  } // while !stop

  // if an offset table was provided and the do_survey flag is still not set,
  // the reference node was certainly not found within the range given; hence no offset could be treated.
  if ( (offs_tab != NULL) && (do_survey == 0)) {
    warning("Offset reference node was not found in the active range", "Offsets were not used.");
    printf("\nOffset reference node: %s active range: %s / %s\n",
	    refnode, use_range[0]->name, use_range[1]->name);
  }

  myfree("Aperture",true_tab);
  if (offs_tab != NULL) myfree("Aperture",offs_tab);

  return lim_pt;
}

