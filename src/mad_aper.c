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
    // 2013-Apr-16  15:38:43  ghislain: following formula is correct but can be rewritten
    // y=sqrt((*ap3)*(*ap3)-(*ap1)*(*ap1));
    // alfa=atan2(y,*ap1);
    // 2013-Apr-23  12:15:51 ghislain --
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
    // 2013-Apr-16  15:41:13  ghislain: following formula is correct but can be rewritten
    // x=sqrt(((*ap3)*(*ap3)) * (1 - ((*ap2)*(*ap2)) / ((*ap4)*(*ap4))));
    // y=sqrt((*ap3)*(*ap3)-x*x);
    // theta=atan2(x,y);
    // 2013-Apr-23  12:15:51 ghislain --
    x = (*ap4)*sqrt(1-((*ap2)*(*ap2))/((*ap4)*(*ap4)));
    
    if (x > (*ap1)) // the rectangle is contained within the ellipse; there is no curved part.
      theta=atan2((*ap1),(*ap2)); 
    else            // the rectangle extends beyond the ellipse; there is a curved part
      theta = atan2(x,*ap2); // this angle is not the geometrical angle
  }

  // at this point we know if we have a full ellipse (alfa=0, theta=0), a full rectangle (alfa+theta=pi/2) or a mixed curve.
  // we can calculate the number of apexes (napex+1) and the interval depending on the shape
  
  if (*quarterlength) { 
    napex=9 ; 
    dangle=((pi/2)-(alfa+theta))/napex; 
  } 

  else if (fabs(alfa+theta-pi/2) < MIN_DOUBLE * 1.e10) { //rectangle, single point, zero interval  
    napex=0 ; 
    dangle=0.; 
  }
  // following case is useless and covered by the general rectellipse. 
  //  else if (fabs(alfa+theta) < MIN_DOUBLE * 1.e10) // ellipse, 20 points, 19 intervals            
  //  { napex=19 ; 
  //    dangle=(pi/2)/napex; } 

  else { // general rectellipse,  1 <= napex <= 19 intervals
    napex = 1 + floor(18 * fabs(1-(alfa+theta)/(pi/2))) ; 
    dangle=((pi/2)-(alfa+theta))/napex; 
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
      tablex[i]=(*ap3)*cos(angle);
      tabley[i]=(*ap4)*sin(angle);
      
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
aper_adj_quad(double angle, double x, double y, double* xquad, double* yquad)
{
  int quadrant;
  quadrant=angle/(pi/2)+1;
  switch (quadrant) 
    {
    case 1: *xquad=x; *yquad=y; break;
    case 2: *xquad=-x; *yquad=y; break;
    case 3: *xquad=-x; *yquad=-y; break;
    case 4: *xquad=x; *yquad=-y; break;
    }
}

static void
aper_adj_halo_si(double ex, double ey, double betx, double bety, double bbeat, double halox[], double haloy[], int halolength, double haloxsi[], double haloysi[])
{
  int j;

  for (j=0;j<=halolength+1;j++) {
    haloxsi[j]=halox[j]*bbeat*sqrt(ex*betx);
    haloysi[j]=haloy[j]*bbeat*sqrt(ey*bety);
  }
}

static double
aper_online(double xm, double ym, double startx, double starty,
            double endx, double endy, double dist_limit)

/* BJ 13.02.2009.
   - added check |x-e| < dist_limit
   - removed useless calculations of sqrt
   - made consistent use of dist_limit and min_double */
/* BJ 13.02.2009
   check if point x = (xm,ym) is in the segment [s,e] with
   s = (startx,starty) and e = (endx,endy) by computing 
   cosfi = (x-s).(x-e) / |x-s||x-e|. cosfi = -1 : x is in
   first check if |x-s| and |x-e| are not too small.If yes for one of them : in
   if OK , the zero divide check must be superfluous. But keep it anyway.
*/

{
  double cosfi=1 , aaa, sss, eee;

  sss = sqrt((xm-startx)*(xm-startx)+ (ym-starty)*(ym-starty));
  eee = sqrt((xm-endx)*(xm-endx)    + (ym-endy)*(ym-endy));

  if ( sss <= dist_limit)      cosfi=-1;
  else if ( eee <= dist_limit) cosfi=-1;
  else {
    aaa = sss * eee ;
    if ( aaa < MIN_DOUBLE )
	fatal_error("Attempt to zero divide ", "In aper_online");
    cosfi=  ((xm-startx)*(xm-endx)+(ym-starty)*(ym-endy)) / aaa;
  }

  if (cosfi <= -1+dist_limit) cosfi=-1;
  
  return cosfi;
}

static void
aper_race(double xshift, double yshift, double r, double angle, double* x, double* y)
{
  /* NEW VERSION of aper_race, 20feb08 BJ, potential zero-divide issues cleared */

  double angle0, angle1, angle2, alfa, gamma, theta;
  int quadrant;

  /* this function calculates the displacement of the beam centre
     due to tolerance uncertainty for every angle */

  if (xshift==0 && yshift==0 && r==0) {
    *x=0; *y=0;
    return;
  }
  
  quadrant=angle/(pi/2)+1;
  
  switch (quadrant) /*adjusting angle to first quadrant*/
    {
    case 1: angle=angle; break;
    case 2: angle=pi-angle; break;
    case 3: angle=angle-pi; break;
    case 4: angle=twopi-angle; break;
    }

  if (angle==pi/2) {
    *x=0;
    *y=yshift+r;
  }

  else {
    /*finding where arc starts and ends*/
    angle0=atan2( yshift , xshift+r );
    angle1=atan2( r+yshift , xshift );
    
    /*different methods is needed, depending on angle*/
    if (angle <= angle0 + MIN_DOUBLE * 1.e10 ) {
      *x=xshift+r;
      *y=tan(angle)*(xshift+r);
    }
    else if (angle<angle1) {
      /* if this is a circle, angle2 useless */
      if (!xshift && !yshift)  angle2 = 0;
      else angle2 = atan2( yshift , xshift );

      alfa = fabs(angle-angle2);
      if (alfa < MIN_DOUBLE * 1.e10) {
	/* alfa==0 is a simpler case */
        *x=cos(angle)*(r+sqrt(xshift*xshift+yshift*yshift));
        *y=sin(angle)*(r+sqrt(xshift*xshift+yshift*yshift));
      }
      else {
	/* solving sine rule w.r.t. gamma */
        gamma=asin(sqrt(xshift*xshift+yshift*yshift)/r*sin(alfa));
	/*theta is the last corner in the triangle*/
        theta=pi-(alfa+gamma);
        *x=cos(angle)*r*sin(theta)/sin(alfa);
        *y=sin(angle)*r*sin(theta)/sin(alfa);
      }
    }
    /* upper flat part */
    else {
      *y=r+yshift;
      *x=(r+yshift)*tan(pi/2-angle);
    }
  }

  return;
}

static int
aper_chk_inside(double p, double q, double pipex[], double pipey[], double dist_limit, int pipelength)
{
  int i;
  double n12, salfa, calfa, alfa=0;

  /*checks first whether p,q is exact on a pipe coordinate*/
  for (i=0;i<=pipelength;i++) {
    if (-1 == aper_online(p, q, pipex[i], pipey[i], pipex[i+1], pipey[i+1], dist_limit))
      return 0; 
  }

  /*calculates and adds up angle from centre between all coordinates*/
  for (i=0;i<=pipelength;i++) {
    n12=sqrt(((pipex[i]-p)*(pipex[i]-p) + (pipey[i]-q)*(pipey[i]-q))
             * ((pipex[i+1]-p)*(pipex[i+1]-p) + (pipey[i+1]-q)*(pipey[i+1]-q)));

    salfa=((pipex[i]-p)*(pipey[i+1]-q) - (pipey[i]-q)*(pipex[i+1]-p))/n12;

    calfa=((pipex[i]-p)*(pipex[i+1]-p) + (pipey[i]-q)*(pipey[i+1]-q))/n12;
    
    alfa += atan2(salfa, calfa);
  }

  /*returns yes to main if total angle is at least twopi*/
  if (sqrt(alfa*alfa)>=(twopi-dist_limit)) return 1;

  return 0;
}

static void
aper_intersect(double a1, double b1, double a2, double b2, double x1, double y1, double x2, double y2,
               int ver1, int ver2, double *xm, double *ym)
{
  (void)y1;
  
  if (ver1 && ver2 && x1==x2) {
    *xm=x2;
    *ym=y2;
  }

  else if (ver1) {
    *xm=x1;
    *ym=a2*x1+b2;
  }

  else if (ver2) {
    *xm=x2;
    *ym=a1*x2+b1;
  }

  else {
    *xm=(b1-b2)/(a2-a1);
    *ym=a1*(*xm)+b1;
  }

  return;
}

static int
aper_linepar(double x1,double y1,double x2,double y2,double *a,double *b)
{
  if ( fabs(x1-x2)< MIN_DOUBLE) { 
    // line is vertical
    *a=0;
    *b=y1-(*a)*x1;
    return 1;
  } 
  else {
    // general line equation
    *a=(y1-y2)/(x1-x2);
    *b=y1-(*a)*x1;
    return 0;
  }
}

static void
aper_fill_quadrants(double polyx[], double polyy[], int quarterlength, int* halolength)
{/* 2013-03-21 ghislain: given the data for the upper right quadrant computed in aper_rectellipse 
     and contained in the polyx and polyy tables,
     mirrors this data to the other three quadrants across the x and y axes.
     quarterlength is the length of, or number of points in, the first quadrant.*/

  int i,j;
  int debug = get_option("debug");

  if (debug) printf("+++ aper_fill_quadrants: quarterlength = %d\n", quarterlength);
  
  // The counter i starts at quarterlength+1, ie the first point to be mirrored.
  i=quarterlength+1;

  /*copying first quadrant coordinates to second quadrant*/
  for (j=quarterlength; j>=0; j--) {
    polyx[i]= -polyx[j];
    polyy[i]=  polyy[j];
    i++;
  }

  /*copying first quadrant coordinates to third quadrant*/
  for (j=0; j<=quarterlength; j++) {
    polyx[i]= -polyx[j];
    polyy[i]= -polyy[j];
    i++;
  }

  /*copying first quadrant coordinates to fourth quadrant*/
  for (j=quarterlength; j>=0; j--) {
    polyx[i]=  polyx[j];
    polyy[i]= -polyy[j];
    i++;
  }

  /*sets the last point equal to the first, to complete the shape.
    Necessary for compatibility with aper_calc function*/
  polyx[i]=polyx[0];
  polyy[i]=polyy[0];

  *halolength=i-1;
  
  if (debug) {
    for (j=0;j<=i;j++) printf("  %d  %10.5f  %10.5f \n", j, polyx[j], polyy[j]);
    printf("\n");
  }

  return;
}

static void
aper_read_twiss(char* table, int* jslice, double* s, double* x, double* y,
                double* betx, double* bety, double* dx, double* dy)
{
  double_from_table_row(table, "s", jslice, s);
  double_from_table_row(table, "x", jslice, x);
  double_from_table_row(table, "y", jslice, y);
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

  /* 2013-03-21 -- ghislain: changed name from aper_bs; the same function is referenced as build_pipe in the documentation of Ivar Waarum */

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

  (*ap1)=(*ap2)=(*ap3)=(*ap4)=0;

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

  else if (!strcmp(apertype,"marguerite")) {
    printf("\nApertype %s not supported.", apertype);
    return 0;
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
    *ap1=get_aperture(current_node, "var1"); /*half width rect*/
    *ap2=get_aperture(current_node, "var2"); /*half height rect*/
    *ap3=get_aperture(current_node, "var3"); /*radius circle*/
    
    *ap4 = *ap3; // curved part is a circle

    // 2014-Jun-27  11:14:27  ghislain: 
    // change check from ap1 or ap2<=0  to ap1 or ap2 or ap3 < 0
    // zero horizontal or vertical explosion factors, and zero radius should be allowed.
    if ( (*ap1) < 0 || (*ap2) < 0 || (*ap3) < 0 ) { 
      if (debug) 
	printf("+++ aper_build screen, racetrack parameters: %10.5f %10.5f %10.5f %10.5f  -- exiting 0\n", *ap1, *ap2, *ap3, *ap4); 
      return 0;
    }
    // special call to build a circle first: note that we cannot invoque ap1 or ap2
    aper_rectellipse(ap3, ap3, ap3, ap4, &quarterlength, pipex, pipey);
      
    /* displace the quartercircle */
    for (i=0;i<=quarterlength;i++) {
      pipex[i] += (*ap1); 
      pipey[i] += (*ap2); 
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
    pipex[0] = (*ap1); pipey[0] = (*ap1)*tan((*ap3));
    pipex[1] = (*ap2)*tan(pi/2 - (*ap4)); pipey[1] = (*ap2);

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
          fgets(comment, 100, e_d_pt);
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
  int i, h_length = 25; // not used, nint=1;
  double dtmp, vtmp[4]; // not used, deltap_twiss;
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
  sprintf(c_dum->c, v_format("@ SEQUENCE         %%%02ds \"%s\""),strlen(tmp),stoupper(tmp));
  aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);
  i = get_string("beam", "particle", tmp);
  sprintf(c_dum->c, v_format("@ PARTICLE         %%%02ds \"%s\""),i,stoupper(tmp));
  aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = get_value("beam", "mass");
  sprintf(c_dum->c, v_format("@ MASS             %%le  %F"), dtmp);
  aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = get_value("beam", "energy");
  sprintf(c_dum->c, v_format("@ ENERGY           %%le  %F"), dtmp);
  aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = get_value("beam", "pc");
  sprintf(c_dum->c, v_format("@ PC               %%le  %F"), dtmp);
  aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = get_value("beam", "gamma");
  sprintf(c_dum->c, v_format("@ GAMMA            %%le  %F"), dtmp);
  aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);


  /* aperture command properties */

  /* 2013-Nov-14  15:45:23  ghislain: The global parameters that have a default in the 
     dictionary or can be input from other commands like BEAM, should be obtained adequately, not 
     from the cmd input.
 */
  dtmp = command_par_value("exn", this_cmd->clone);
  sprintf(c_dum->c, v_format("@ EXN              %%le  %F"), dtmp);
  aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = command_par_value("eyn", this_cmd->clone);
  sprintf(c_dum->c, v_format("@ EYN              %%le  %F"), dtmp);
  aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = command_par_value("dqf", this_cmd->clone);
  sprintf(c_dum->c, v_format("@ DQF              %%le  %F"), dtmp);
  aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = command_par_value("betaqfx", this_cmd->clone);
  sprintf(c_dum->c, v_format("@ BETAQFX          %%le  %F"), dtmp);
  aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = command_par_value("dparx", this_cmd->clone);
  sprintf(c_dum->c, v_format("@ PARAS_DX         %%le       %g"), dtmp);
  aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = command_par_value("dpary", this_cmd->clone);
  sprintf(c_dum->c, v_format("@ PARAS_DY         %%le       %g"), dtmp);
  aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = command_par_value("dp", this_cmd->clone);
  sprintf(c_dum->c, v_format("@ DP_BUCKET_SIZE   %%le  %F"), dtmp);
  aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);

  // LD: summary table is corrupted by embedded twiss, use saved value
  // double_from_table_row("summ","deltap",&nint,&deltap_twiss);
  sprintf(c_dum->c, v_format("@ TWISS_DELTAP     %%le  %F"), lim->deltap_twiss);
  aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);

  dtmp = command_par_value("cor", this_cmd->clone);
  sprintf(c_dum->c, v_format("@ CO_RADIUS        %%le  %F"), dtmp);
  aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = command_par_value("bbeat", this_cmd->clone);
  sprintf(c_dum->c, v_format("@ BETA_BEATING     %%le  %F"), dtmp);
  aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);
  dtmp = command_par_value("nco", this_cmd->clone);
  sprintf(c_dum->c, v_format("@ NB_OF_ANGLES     %%d   %g"), dtmp*4);
  aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);

  /* if a filename with halo coordinates is given, need not show halo */
  stmp = command_par_string("halofile", this_cmd->clone);
  if (stmp)
  {
    strncpy(tmp, stmp, sizeof tmp);
    sprintf(c_dum->c, v_format("@ HALOFILE         %%%02ds \"%s\""),strlen(tmp),stoupper(tmp));
    aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);
  }
  else
  {
    i = command_par_vector("halo", this_cmd->clone, vtmp);
    sprintf(c_dum->c, v_format("@ HALO_PRIM        %%le       %g"),vtmp[0]);
    aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);
    sprintf(c_dum->c, v_format("@ HALO_R           %%le       %g"),vtmp[1]);
    aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);
    sprintf(c_dum->c, v_format("@ HALO_H           %%le       %g"),vtmp[2]);
    aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);
    sprintf(c_dum->c, v_format("@ HALO_V           %%le       %g"),vtmp[3]);
    aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);
  }
  /* show filename with pipe coordinates if given */
  stmp = command_par_string("pipefile", this_cmd->clone);
  if (stmp)
  {
    strncpy(tmp, stmp, sizeof tmp);
    sprintf(c_dum->c, v_format("@ PIPEFILE         %%%02ds \"%s\""), strlen(tmp), stoupper(tmp) );
    aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);
  }

  /* 2013-Nov-18  13:49:45  ghislain: changed n1min and at_element strings to uppercase in output to TFS */
  sprintf(c_dum->c, v_format("@ N1MIN            %%le   %g"), lim->n1);
  aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);
  set_value("beam","n1min",&lim->n1);

  sprintf(c_dum->c, v_format("@ AT_ELEMENT       %%%02ds  \"%s\""), strlen(name), stoupper(name) );
  aper_t->header->p[aper_t->header->curr++] = tmpbuff(c_dum->c);
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
  aper_survey->cmd_def->par->parameters[aspos]->string = "survey";
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
  int c=0;

  /* Replaces the first ws or : in a string with a '\0', */
  /* thus translating a FORTRAN-like attribute string to */
  /* C compatibility, or washes the ':1' from node names */

  while (string[c] && string[c]!=' ' && c<len-1) c++;

  string[c]=0;
  if (c<len-1) string[c+1]=' '; /*adds a ws to avoid two \0 in a row*/
}

static void
aper_write_table(char* name, double* n1, double* n1x_m, double* n1y_m,
                  double* rtol, double* xtol, double* ytol,
                  char* apertype,double* ap1,double* ap2,double* ap3,double* ap4,
                  double* on_ap, double* on_elem, double* spec,double* s,
                  double* x, double* y, double* betx, double* bety,double* dx, double* dy,
                  char *table)
{
  string_to_table_curr(table, "name", name);
  double_to_table_curr(table, "n1", n1);
  double_to_table_curr(table, "n1x_m", n1x_m);
  double_to_table_curr(table, "n1y_m", n1y_m);
  double_to_table_curr(table, "rtol", rtol);
  double_to_table_curr(table, "xtol", xtol);
  double_to_table_curr(table, "ytol", ytol);
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
  double_to_table_curr(table, "betx", betx);
  double_to_table_curr(table, "bety", bety);
  double_to_table_curr(table, "dx", dx);
  double_to_table_curr(table, "dy", dy);

  augment_count(table);
}

static double
aper_calc(double p, double q, double* minhl, double halox[], double haloy[],
          int halolength,double haloxadj[],double haloyadj[],
          double newhalox[], double newhaloy[], double pipex[], double pipey[],
          int pipelength, double notsimple)
{
  int i=0, j=0, c=0, ver1, ver2;
  double dist_limit=0.0000000001; // 1.e-10
  double a1, b1, a2, b2, xm, ym, h, l;

  for (c=0; c<=halolength+1; c++) {
    haloxadj[c]=halox[c]+p;
    haloyadj[c]=haloy[c]+q;
  }

  c=0;

  /* if halo centre is inside beam pipe, calculate smallest H/L ratio */
  if (aper_chk_inside(p, q, pipex, pipey, dist_limit, pipelength)) {
    if (notsimple) {
      /* Add extra apexes first:*/
      for (j=0; j<=halolength; j++) {
        newhalox[c]=haloxadj[j];
        newhaloy[c]=haloyadj[j];
        c++;

        for (i=0;i<=pipelength;i++) {
          /*Find a and b parameters for line*/
          ver1=aper_linepar(p, q, pipex[i], pipey[i], &a1, &b1);
          ver2=aper_linepar(haloxadj[j], haloyadj[j], haloxadj[j+1], haloyadj[j+1], &a2, &b2);

          /*find meeting coordinates for infinitely long lines*/
          aper_intersect(a1, b1, a2, b2, pipex[i], pipey[i], haloxadj[j], haloyadj[j], ver1, ver2, &xm, &ym);

          /*eliminate intersection points not between line limits*/
          if (-1 == aper_online(xm, ym, haloxadj[j], haloyadj[j], haloxadj[j+1], haloyadj[j+1], dist_limit)) { /*halo line*/
            if (-1 != aper_online(p, q, pipex[i], pipey[i], xm, ym, dist_limit)) { /*test line*/
              newhalox[c]=xm;
              newhaloy[c]=ym;
              c++;
            }
          }
        }
      }

      halolength=c-1;
      for (j=0; j<=halolength; j++) {
        haloxadj[j]=newhalox[j];
        haloyadj[j]=newhaloy[j];
      }
    }

    /*Calculate smallest ratio:*/
    for (i=0; i<=pipelength; i++) {
      for (j=0; j<=halolength; j++) {
        /*Find a and b parameters for line*/
        ver1 = aper_linepar(p, q, haloxadj[j], haloyadj[j], &a1, &b1);
        ver2 = aper_linepar(pipex[i], pipey[i], pipex[i+1], pipey[i+1], &a2, &b2);

        /*find meeting coordinates for infinitely long lines*/
        aper_intersect(a1, b1, a2, b2, haloxadj[j], haloyadj[j], pipex[i], pipey[i], ver1, ver2, &xm, &ym);

        /*eliminate intersection points not between line limits*/
        if (-1 == aper_online(xm, ym, pipex[i], pipey[i], pipex[i+1], pipey[i+1], dist_limit)) { /*pipe line*/
          if (-1 != aper_online(p, q, haloxadj[j], haloyadj[j], xm, ym, dist_limit)) {  /*test line*/
            h=sqrt((xm-p)*(xm-p)+(ym-q)*(ym-q));
            l=sqrt((haloxadj[j]-p)*(haloxadj[j]-p) + (haloyadj[j]-q)*(haloyadj[j]-q));
            if (h/l < *minhl) *minhl=h/l;            
          }
        }
      }
    }
  }
  else { /*if halo centre is outside of beam pipe*/    
    *minhl=0;
    return -1;
  }

  return 0;
}

// public interface

double
get_apertol(struct node* node, char* par)
  /* returns aper_tol parameter 'i' where i is integer at the end of par;
     e.g. aptol_1 gives i = 1 etc. (count starts at 1) */
{
  int i, k, n = strlen(par);
  double val = zero, vec[100];
  for (i = 0; i < n; i++)  if(isdigit(par[i])) break;
  if (i == n) return val;
  sscanf(&par[i], "%d", &k); k--;
  if ((n = element_vector(node->p_elem, "aper_tol", vec)) > k)  val = vec[k];
  return val;
}

double 
get_aperture(struct node* node, char* par)
  /* returns aperture parameter 'i' where i is integer at the end of par;
     e.g. aper_1 gives i = 1 etc. (count starts at 1) */
{
  int i, k, n = strlen(par);
  double val = zero, vec[100];
  for (i = 0; i < n; i++)  if(isdigit(par[i])) break;
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
  tw_cp=current_sequ->tw_table;
  tw_cnt=1; /* table starts at 1 */

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
  aperture_table=make_table(table, table, ap_table_cols, ap_table_types, rows);
  aperture_table->dynamic=1;
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
    
    if (strcmp(aptwfile,"dummy")) out_table(tw_cp->name, tw_cp, aptwfile);
  }
  else warning("Could not run aperture command.","Aperture command ignored");

  /* set pointer to updated Twiss table */
  current_sequ->tw_table=tw_cp;
}

static struct aper_node*
aperture(char *table, struct node* use_range[], struct table* tw_cp, int *tw_cnt, struct aper_node *lim_pt)
{
  int stop=0, nint=1, jslice=1, first, ap=1; // , err not used
  int true_flag, true_node=0, offs_node=0, do_survey=0;
  int truepos=0, true_cnt=0, offs_cnt=0;
  int halo_q_length=1, halolength, pipelength, namelen=NAME_L, ntol; // nhalopar, not used
  double surv_init[6]={0, 0, 0, 0, 0, 0};
  double surv_x=zero, surv_y=zero;
  double xa=0, xb=0, xc=0, ya=0, yb=0, yc=0;
  double on_ap=1, on_elem=0;
  double mass, energy, exn, eyn, dqf, betaqfx, dp, dparx, dpary;
  double cor, bbeat, nco, halo[4], interval, spec, ex, ey, notsimple;
  double s=0, x=0, y=0, betx=0, bety=0, dx=0, dy=0, ratio, n1, nr, length;
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
  double haloxadj[MAXARRAY], haloyadj[MAXARRAY], newhalox[MAXARRAY], newhaloy[MAXARRAY];
  double pipex[MAXARRAY], pipey[MAXARRAY];
  double parxd,paryd;
  char *halofile, *truefile, *offsfile;
  char refnode[NAME_L]="";
  char *cmd_refnode;
  char apertype[NAME_L];
  char name[NAME_L];
  char tol_err_mess[80] = "";

  
  // 2014-Sep-18  17:19:52  ghislain: attempt to read offset values from element attributes...
  //double aper_offset[2], xoffset, yoffset;
  //int noffset;

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
 
  exn = command_par_value("exn", this_cmd->clone);
  eyn = command_par_value("eyn", this_cmd->clone);
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

  mass = get_value("beam", "mass");
  energy = get_value("beam", "energy");

  /* 2013-Nov-13  16:20:29  ghislain: attempt to extract relevant parameters 
     from BEAM command instead of internal parameters. 
     This works but needs a bit more thoughts, in conjuction with fetching other parameters 
     from Twiss table

  // exn = get_value("beam", "exn");
  // eyn = get_value("beam", "eyn");

  // and attempt to extract maximum parameters from twiss summary table
  // double_from_table_row("summ","dxmax",&nint,&dqf);
  // printf ("+++++++ dqf from TWISS %12.6g\n",dqf);
  // TODO:  add some error checking.

   end of ghislain: attempt... */

  /* fetch deltap as set by user in the former TWISS command */
  /* will be used below for displacement associated to parasitic dipersion */

  if ( (rc_deltap = double_from_table_row("summ","deltap",&nint,&lim_pt->deltap_twiss)) != 0) {
    printf ("+++++++ deltap from TWISS could not be read; assume 0.0d0\n");
  }
  else printf ("+++++++ deltap from TWISS %12.6g\n",lim_pt->deltap_twiss);

  /* calculate emittance and delta angle */
  /* mad_beam.c says :     ex = exn / (4 * beta * gamma); */
  /* Warning: 1- MAD uses a different definition for emittance
              2- This assumes beta = 1, explicitly ultra-relativistic particles */
  ex = mass*exn/energy; ey = mass*eyn/energy;
  dangle = twopi/(nco*4);

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
  if ((halolength = aper_external_file(halofile, halox, haloy)) > -1)
    ;
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
  aper_read_twiss(tw_cp->name, tw_cnt, &s_end, &x, &y, &betx, &bety, &dx, &dy);
  // LD: shift further results by one step (?) and finish outside the table
  //  (*tw_cnt)++;
  aper_adj_halo_si(ex, ey, betx, bety, bbeat, halox, haloy, halolength, haloxsi, haloysi);

  /* calculate initial normal+parasitic disp. */
  /* modified 27feb08 BJ */
  parxd = dparx*sqrt(betx/betaqfx)*dqf;
  paryd = dpary*sqrt(bety/betaqfx)*dqf;

  /* Initialize n1 limit value */
  lim_pt->n1=999999;

  int n = 0;

  while (!stop) {
    ++n;

    strcpy(name,current_node->name);
    aper_trim_ws(name, NAME_L);

    is_zero_len = 0;

    /* the first node in a sequence can not be sliced, hence: */
    if (current_sequ->range_start == current_node) first=1; else first=0;

    length=node_value("l");
    double_from_table_row(tw_cp->name, "s", tw_cnt, &s_end);
    s_start=s_end-length;
    s_curr=s_start;

    node_string("apertype", apertype, &namelen);
    aper_trim_ws(apertype, NAME_L);

    if (!strncmp("drift",name,5)) on_elem=-999999;
    else on_elem=1;

    if ( (offs_tab != NULL) && (strcmp(refnode, name) == 0)) do_survey=1; // current name is refnode; switch survey on.

    if (debug) printf("\nname: %s, ref: %s, do_survey: %d, true_flag: %d\n",name,refnode,do_survey,true_flag);

    /* read data for tol displacement of halo */
    get_node_vector("aper_tol",&ntol,aper_tol);
    if (ntol == 3) {
      r = aper_tol[0];
      xshift = aper_tol[1];
      yshift = aper_tol[2];
    }
    else r=xshift=yshift=0;
    //if (debug) printf("\nname: %s, x-shift: %f, y-shift: %f\n",name,xshift, yshift);

    // 2014-Sep-18  17:19:52  ghislain: attempt to read offset values from element attributes...
    /* read data for aper_offset */
    //get_node_vector("aper_offset",&noffset,aper_offset);
    //if (noffset == 2) {
    //  xoffset = aper_offset[0];
    //  yoffset = aper_offset[1];
    //}
    //else xoffset=yoffset=0;
    //if (debug) printf("\nname: %s, x-offset: %f, y-offset: %f\n",name,xoffset, yoffset);

    /*read aperture data and make polygon tables for beam pipe*/
    /* IW 250205 */
    /*  if (ext_pipe == 0) */
    ap=aper_build_screen(apertype, &ap1, &ap2, &ap3, &ap4, &pipelength, pipex, pipey);

    if (ap == 0 || first == 1) {
      /* if no pipe can be built, the n1 is set to inf and Twiss parms read for reference*/
      n1=999999; n1x_m=999999; n1y_m=999999; on_ap=-999999; nint=1;

      aper_read_twiss(tw_cp->name, tw_cnt, &s_end, &x, &y, &betx, &bety, &dx, &dy);

      aper_write_table(name, &n1, &n1x_m, &n1y_m, &r, &xshift, &yshift, apertype,
                       &ap1, &ap2, &ap3, &ap4, &on_ap, &on_elem, &spec,
                       &s_end, &x, &y, &betx, &bety, &dx, &dy, table);
      on_ap=1;

      double_to_table_row(tw_cp->name, "n1", tw_cnt, &n1);
      (*tw_cnt)++;

      /* calc disp and adj halo to have ready for next node */
      /* modified 27feb08 BJ */
      parxd = dparx*sqrt(betx/betaqfx)*dqf;
      paryd = dpary*sqrt(bety/betaqfx)*dqf;

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
      node_n1=999999;
      true_node=0;
      offs_node=0;

      /* calculate the number of slices per node */
      if (true_flag == 0)
        nint=length/interval;
      else {
        true_node=aper_tab_search(true_cnt, true_tab, name, &truepos);
        if (true_node) nint=true_tab[truepos].curr;
        else nint=length/interval;
      }

      if (debug) printf("name: %s, nint: %d, length: %f\n", name, nint, length);

      if (!nint) nint=1;

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
          s=0;                  /*used to calculate survey adjustments */
        } 
	else {
          aper_read_twiss("embedded_twiss_table", &jslice, &s, &x, &y, &betx, &bety, &dx, &dy);
	  
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

        for (angle=0; angle < twopi; angle += dangle) {
          /* new 27feb08 BJ */
          dispx = bbeat*(fabs(dx)*dp + parxd*(fabs(lim_pt->deltap_twiss)+dp) );
          dispy = bbeat*(fabs(dy)*dp + paryd*(fabs(lim_pt->deltap_twiss)+dp) );

          /*adjust dispersion to worst-case for quadrant*/
          aper_adj_quad(angle, dispx, dispy, &dispxadj, &dispyadj);

          /*calculate displacement co+tol for each angle*/
          coxadj=cor*cos(angle); coyadj=cor*sin(angle);

          /* Error check added 20feb08 BJ */
          if ( xshift < 0 || yshift < 0 || r < 0 ) {
            sprintf(tol_err_mess,"In element : %s\n",name);
            fatal_error("Illegal negative tolerance",tol_err_mess);
          }

          aper_race(xshift,yshift,r,angle,&tolx,&toly);
          aper_adj_quad(angle, tolx, toly, &tolxadj, &tolyadj);

          /* add all displacements */
          deltax = coxadj + tolxadj + xeff + dispxadj;
          deltay = coyadj + tolyadj + yeff + dispyadj;

          /* send beta adjusted halo and its displacement to aperture calculation */
          aper_calc(deltax,deltay,&ratio,haloxsi,haloysi,
                    halolength,haloxadj,haloyadj,newhalox,newhaloy,
                    pipex,pipey,pipelength,notsimple);
        }

        nr=ratio*halo[1];
        n1=nr/(halo[1]/halo[0]); /* ratio r/n = 1.4 */

        n1x_m=n1*bbeat*sqrt(betx*ex);
        n1y_m=n1*bbeat*sqrt(bety*ey);

	/* Change below, BJ 23oct2008                              */
	/* test block 'if (n1 < node_n1)' included in test block   */

        if (is_zero_len == 0 || jslice == 1) {
          aper_write_table(name, &n1, &n1x_m, &n1y_m, &r, &xshift, &yshift, apertype,
                           &ap1, &ap2, &ap3, &ap4, &on_ap, &on_elem, &spec, &s_curr,
                           &xeff, &yeff, &betx, &bety, &dx, &dy, table);

	  /* save node minimum n1 */
          if (n1 < node_n1) {
              node_n1=n1;
              node_s=s_curr;
          }
	} // if is_zero_len
      } // for jslice

      reset_interpolation(&nint);

      /* insert minimum node value into Twiss table */
      double_to_table_row(tw_cp->name, "n1", tw_cnt, &node_n1);
      (*tw_cnt)++;

      /* save range minimum n1 */
      if (node_n1 < lim_pt->n1) {
        strcpy(lim_pt->name,name);
        lim_pt->n1=node_n1;
        lim_pt->s=node_s;
        strcpy(lim_pt->apertype,apertype);
        lim_pt->aperture[0]=ap1;
        lim_pt->aperture[1]=ap2;
        lim_pt->aperture[2]=ap3;
        lim_pt->aperture[3]=ap4;
        lim_pt->aper_tol[0]=r;
        lim_pt->aper_tol[1]=xshift;
        lim_pt->aper_tol[2]=yshift;
      }
    } // if !(ap == 0 || first == 1)

    if (!strcmp(current_node->name,use_range[1]->name)) stop=1;
    if (!advance_node()) stop=1;
  } // while !stop

  // 2013-Sep-17  16:45:34  ghislain: 
  // if an offset table was provided and the do_survey flag is still not set, 
  // the reference node was certainly not found within the range given; hence no offset could be treated.
  if ( (offs_tab != NULL) && (do_survey == 0)) {
    printf("\nWarning: Offset reference node %s was not found in the active range %s to %s.",
	    refnode,use_range[0]->name,use_range[1]->name); 
    printf("\nOffsets were not used.");
    printf("\nThe active range must contain the reference node for offsets to be taken into account.\n"); 
  }  

  myfree("Aperture",true_tab);
  if (offs_tab != NULL) myfree("Aperture",offs_tab);

  return lim_pt;
}

