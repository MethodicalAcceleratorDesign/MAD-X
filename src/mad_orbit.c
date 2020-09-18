#include "madx.h"

// private types

struct val_mic {
  double before[2];
  double after[2];
};

struct id_mic {
  int id_ttb;
  int enable;
  struct val_mic val;
  struct node* p_node;
  struct id_mic *next;
  struct id_mic *previous;
};

struct id_mic2 {
  int id_ttb[2];
  int enable;
  struct val_mic val;
  struct node* p_node;
  struct node* p_node_s1;
  struct node* p_node_s2;
  struct id_mic2 *next;
  struct id_mic2 *previous;
};

struct orb_cor {
  double qx0;
  double qy0;
  double units;
  struct id_mic *cor_table;
  struct id_mic *mon_table;
};

struct orb_cor2 {
  double qx0;
  double qy0;
  double units;
  struct id_mic2 *cor_table;
  struct id_mic2 *mon_table;
};

// forward declarations

static void correct_correct(struct in_cmd*);
static void correct_usemonitor(struct in_cmd*);
static void correct_usekick(struct in_cmd*);
static void correct_putorbit(struct in_cmd*);
static void correct_getorbit(struct in_cmd*); /* empty */
static void correct_option(struct in_cmd*);
static void correct_readcorr(struct in_cmd*);
static void correct_setcorr(struct in_cmd*);
// static void correct_prtcorr(struct in_cmd*);

static void correct_correct1(struct in_cmd*);
static int pro_correct_getactive(int ip, int *nm, int *nx, int *nc,
				 double *corvec, double *monvec, char *conm);
static void pro_correct_write_results(double *monvec, double *resvec, double *corvec,
				      int *nx, int *nc, int *nm, int imon, int icor, int ip,
				      int resout);
static void pro_correct_fill_mon_table(int ip, char *name, double old, double new_);
static void pro_correct_fill_corr_table(int ip, char *name, double old, double new_);
static void pro_correct_make_mon_table(void);
static void pro_correct_make_corr_table(void);
static double* pro_correct_response_line(int ip, int nc, int nm);
static double* pro_correct_response_ring(int ip, int nc, int nm);
static int pro_correct_filter(int iplane, double sigcut);
static void pro_correct_write_cocu_table(void);
static void pro_correct_prtwiss(void);
static int pro_correct_getcorrs(struct in_cmd*);
static int pro_correct_getorbit_ext(struct in_cmd*);
static int pro_correct_getorbit(struct in_cmd*);
static int pro_correct_gettables(int iplane, struct in_cmd*);
static int pro_correct_getcommands(struct in_cmd*);
// static void    pro_correct_option(struct in_cmd*);

static void correct_correct2(struct in_cmd*);
static void pro_correct2_fill_mon_table(int ip, char *name, double old, double new_);
static void pro_correct2_fill_corr_table(int b, int ip, char *name, double old, double new_);
static void pro_correct2_make_mon_table(void);
static void pro_correct2_make_corr_table(void);
static double* pro_correct2_response_ring(int ip, int nc, int nm);
static int pro_correct2_getactive(int ip, int *nm, int *nx, int *nc,
				  double *corvec, double *monvec, char *conm);
static int pro_correct2_getcorrs(struct in_cmd*);
static int pro_correct2_getorbit(struct in_cmd*);
static int pro_correct2_gettables(int iplane, struct in_cmd*);
static void pro_correct2_write_results(double *monvec, double *resvec, double *corvec,
				       int *nx, int *nc, int *nm, int imon, int icor, int ip,
				       int resout);

static void fill_orbit_table(struct table* t_out, struct table* t_in);

// private interface

static double caverage(double *r, int m) {
  // 2014-Jan-07  17:33:59  ghislain: added function
  double xave   = { 0.0 };
  int i;

  for (i = 0; i < m; i++) {
    xave = xave + r[i];
  }
  xave = xave / m;

  return (xave);
}

static double cstddev(double *r, int m) {
  // 2014-Jan-07  17:33:59  ghislain: changed function name from crms to cstddev, and variable names
  double xave   = { 0.0 };
  double xstdev = { 0.0 };
  int i;

  for (i = 0; i < m; i++) {
    xave = xave + r[i];
  }
  xave = xave / m;
  for (i = 0; i < m; i++) {
    xstdev = xstdev + (xave - r[i]) * (xave - r[i]);
  }
  xstdev = sqrt(xstdev / m);

  return (xstdev);
}

static double crms(double *r, int m) {
  // 2014-Jan-07  17:33:59  ghislain: added function
  double xrms = { 0.0 };
  int i;

  for (i = 0; i < m; i++) {
    xrms = xrms +  r[i]*r[i];
  }
  xrms = sqrt(xrms / m);

  return (xrms);
}

static double cptp(double *r, int m) {
  // 2014-Jan-07  17:33:59  ghislain: changed function name from cprp (typo ?)
  double xhi = { -9999. };
  double xlo = { 9999. };
  double xptp = { 0.0 };
  int i;

  for (i = 0; i < m; i++) {
    if (r[i] < xlo) xlo = r[i];
    if (r[i] > xhi) xhi = r[i];
  }
  xptp = xhi - xlo;

  return (xptp);
}

static double copk(double *r, int m) {
  double xpk = { -9999. };
  int i;

  for (i = 0; i < m; i++) {
    if (fabs(r[i]) > xpk) xpk = fabs(r[i]);
  }

  return (xpk);
}

static int c_micit(double *dmat, char *conm, double *monvec, double *corvec,
       double *resvec, int *nx, double rms, int imon, int icor, int niter) {
  const char *rout_name = "c_micit";
  int *ny;
  int ifail;
  double *ax, *cinx, *xinx, *resx;
  double *rho, *ptop, *rmss, *xrms, *xptp, *xiter;

  /* allocate auxiliary vectors used by correction algorithms */
  ny    = mycalloc_atomic("c_micit_ny"   , icor , sizeof *ny);
  ax    = mycalloc_atomic("c_micit_ax"   , imon*icor, sizeof *ax);
  cinx  = mycalloc_atomic("c_micit_cinx" , icor , sizeof *cinx);
  xinx  = mycalloc_atomic("c_micit_xinx" , imon , sizeof *xinx);
  resx  = mycalloc_atomic("c_micit_resx" , imon , sizeof *resx);
  rho   = mycalloc_atomic("c_micit_rho"  , 3*icor , sizeof *rho);
  ptop  = mycalloc_atomic("c_micit_ptop" , icor , sizeof *ptop);
  rmss  = mycalloc_atomic("c_micit_rmss" , icor , sizeof *rmss);
  xrms  = mycalloc_atomic("c_micit_xrms" , icor , sizeof *xrms);
  xptp  = mycalloc_atomic("c_micit_xptp" , icor , sizeof *xptp);
  xiter = mycalloc_atomic("c_micit_xiter", icor , sizeof *xiter);

  micit_(dmat, conm, monvec, corvec, resvec, nx, &rms, &imon, &icor, &niter,
   ny, ax, cinx, xinx, resx, rho, ptop, rmss, xrms, xptp, xiter, &ifail);

  myfree(rout_name, ny);
  myfree(rout_name, ax);
  myfree(rout_name, cinx);
  myfree(rout_name, xinx);
  myfree(rout_name, resx);
  myfree(rout_name, rho);
  myfree(rout_name, ptop);
  myfree(rout_name, rmss);
  myfree(rout_name, xrms);
  myfree(rout_name, xptp);
  myfree(rout_name, xiter);

  return (ifail);
}

static void c_haveit(double *dmat, double *monvec, double *corvec,
         double *resvec, int *nx, int imon, int icor) {
  const char *rout_name = "c_haveit";
  double *cb, *xmeas, *xres, *y, *z, *xd;

  cb    = mycalloc_atomic("c_haveit_cb" ,   icor , sizeof *cb);
  xmeas = mycalloc_atomic("c_haveit_xmeas", imon , sizeof *xmeas);
  xres  = mycalloc_atomic("c_haveit_xres" , imon , sizeof *xres);
  y     = mycalloc_atomic("c_haveit_y" ,    icor*imon , sizeof *y);
  z     = mycalloc_atomic("c_haveit_z" ,    icor*icor , sizeof *z);
  xd    = mycalloc_atomic("c_haveit_xd" ,   icor , sizeof *xd);

  haveit_(dmat, monvec, corvec, resvec, nx, &imon, &icor, cb, xmeas, xres, y, z, xd);

  myfree(rout_name, cb);
  myfree(rout_name, xmeas);
  myfree(rout_name, xres);
  myfree(rout_name, y);
  myfree(rout_name, z);
  myfree(rout_name, xd);
}

static int c_svddec(double *dmat, int imon, int icor, int *sing,
		     double *sngcut, double *sngval) {
  const char *rout_name = "c_svddec";
  int flag;
  int debug = get_option("debug");

  double *s, *u, *v;
  double *ws, *wv;
  int *sw;

  s  = mycalloc_atomic("c_svddec_s" , icor*imon, sizeof *s);
  u  = mycalloc_atomic("c_svddec_u" , icor*imon, sizeof *u);
  v  = mycalloc_atomic("c_svddec_v" , icor*imon, sizeof *v);
  ws = mycalloc_atomic("c_svddec_ws", icor , sizeof *ws);
  wv = mycalloc_atomic("c_svddec_wv", icor , sizeof *wv);
  sw = mycalloc_atomic("c_svddec_sw", icor , sizeof *sw);

  svddec_(dmat, s, u, v, ws, wv, sw, sngcut, sngval,
	  &imon, &icor, &flag, sing, &debug);

  myfree(rout_name, s);
  myfree(rout_name, u);
  myfree(rout_name, v);
  myfree(rout_name, ws);
  myfree(rout_name, wv);
  myfree(rout_name, sw);

  return flag;
}

static int c_svdcorr(double *dmat, double *xin, double *cor, double *res,
		      int *nx, int imon, int icor) {
  const char *rout_name = "c_svdcorr";
  int flag;
  int debug = get_option("debug");

  double *s, *u, *v, *w, *ut, *vt, *wt;
  double *xp, *wv, *ws;
  int *sw;

  s  = mycalloc_atomic("c_svdcorr_s" , imon*icor, sizeof *s);
  u  = mycalloc_atomic("c_svdcorr_u" , imon*icor, sizeof *u);
  v  = mycalloc_atomic("c_svdcorr_v" , icor*icor, sizeof *v);
  w  = mycalloc_atomic("c_svdcorr_w" , icor*icor, sizeof *w);
  ut = mycalloc_atomic("c_svdcorr_ut", icor*imon, sizeof *ut);
  vt = mycalloc_atomic("c_svdcorr_vt", icor*icor, sizeof *vt);
  wt = mycalloc_atomic("c_svdcorr_wt", icor*icor, sizeof *wt);

  xp = mycalloc_atomic("c_svdcorr_xp", imon, sizeof *xp);
  ws = mycalloc_atomic("c_svdcorr_xp", icor, sizeof *ws);
  wv = mycalloc_atomic("c_svdcorr_xp", icor, sizeof *wv);

  sw = mycalloc_atomic("c_svdcorr_sw", icor, sizeof *sw);

  svdcorr_(dmat, s, u, v, w, ut, vt, wt, xin, cor, res, xp, ws,
	   wv, sw, nx, &imon, &icor, &flag, &debug);

  myfree(rout_name, s);
  myfree(rout_name, u);
  myfree(rout_name, v);
  myfree(rout_name, w);
  myfree(rout_name, ut);
  myfree(rout_name, vt);
  myfree(rout_name, wt);
  myfree(rout_name, sw);
  myfree(rout_name, xp);
  myfree(rout_name, ws);
  myfree(rout_name, wv);

  return flag;
}

static void fill_orbit_table(struct table* t_out, struct table* t_in)
/* fills a table with orbit values at monitor positions */
{
  int i, j, pos;
  t_out->curr = 0;
  for (i = 0; i < t_in->curr; i++) {
    if (strstr(t_in->s_cols[1][i], "monitor")) {
      for (j = 0; j < t_out->num_cols; j++) {
  if ((pos = name_list_pos(t_out->columns->names[j], t_in->columns)) > -1) {
    if (t_out->columns->inform[j] < 3)
      t_out->d_cols[j][t_out->curr] = t_in->d_cols[pos][i];
    else
      t_out->s_cols[j][t_out->curr] = tmpbuff( t_in->s_cols[pos][i]);
  } else {
    if (t_out->columns->inform[j] < 3)
      t_out->d_cols[j][t_out->curr] = zero;
    else
      t_out->s_cols[j][t_out->curr] = tmpbuff(blank);
  }
      }
      t_out->curr++;
    }
  }
}

static void correct_setcorr(struct in_cmd* cmd) {
  // 2013-Dec-04  08:56:16  ghislain: The command SETCORR is not documented.

  /* read the correctors from named table and stores
     them in the nodes of the sequence at
     "->chkick" and "->cvkick". Subsequent Twiss will
     use them correctly.
     ===> Must be preceded by a call to "read_table"
     ===> (unless table exists in memory !)

     ===> Watch out, does not yet take care of existing corrector
     ===> settings already present in sequence

     ===> Uses table with name specified in parameter: table=
  */

  int i, ix;

  struct node *ndexe;
  struct node *nextnode;

  char name[NAME_L];
  char slname[NAME_L];

  char nname[NAME_L];
  char slnname[NAME_L];

  const char* namtab;

  double xnew, ynew;

  /* set up pointers to current sequence for later use */
  struct sequence* mysequ = current_sequ;
  nextnode = mysequ->ex_start;
  ndexe = mysequ->ex_end;

  /* printf("Pointers: %d %d %d\n",mysequ,nextnode,ndexe); */

  if ((namtab = command_par_string("table", cmd->clone)) != NULL ) {
    printf("Want to use named table: %s\n", namtab);
    if (table_exists(namtab)) {
      printf("The table ==> %s <=== was found \n", namtab);
    } else {
      /* fatal_error("Corrector table requested, but not existing:",namtab); */
      /* exit(-77); */
      printf("No such corrector table in memory: %s\n", namtab);
    }
  } else {
    if (get_option("debug")) {
      printf("No table name requested\n");
      printf("Use default name\n");
    }
    namtab = "corr";
  }

  i = 1;
  ix = 0;

  while (ix == 0) {

    // 2013-Dec-04  10:48:48  ghislain: changed ix= into ix+=
    ix += string_from_table_row(namtab, "name", &i, name);
    ix += double_from_table_row(namtab, "px.correction", &i, &xnew);
    ix += double_from_table_row(namtab, "py.correction", &i, &ynew);

    if (ix == 0) { // 2013-Dec-02  18:19:55  ghislain: what was that supposed to test ?
                   // Before: just the last return value above, Now: any of the above return values ?
      stolower(name);
      strcpy(slname, strip(name));
      supp_tb(slname);

      /* printf("corrs: %s %d %e %e %e %e\n",name,ix,xold,yold,xnew,ynew); */
      nextnode = mysequ->ex_start;
      while (nextnode != ndexe) {
	stolower(name);
	strcpy(slname, strip(name));
	supp_tb(slname);

	strcpy(nname, nextnode->name);
	stolower(nname);
	strcpy(slnname, strip(nname));
	supp_tb(slnname);

	if (strcmp(slname, slnname) == 0) {
	  nextnode->chkick += xnew;
	  nextnode->cvkick += ynew;
	  nextnode = ndexe;
	} else {
	  nextnode = nextnode->next;
	}
      }
    }
    i++;
  }
  return;
}

static void correct_readcorr(struct in_cmd* cmd) {
  // 2013-Dec-04  08:56:16  ghislain: The command READCORR is not documented.

  /* read the correctors from table "corr" and stores
     them in the nodes of the sequence at
     "->chkick" and "->cvkick". Subsequent Twiss will
     use them correctly.
     ===> Must be preceded by a call to "read_table"

     ===> Watch out, does not yet take care of existing corrector
     ===> settings already present in sequence
     2013-Dec-04  10:43:32  ghislain: I think this is done already (see code)

     ===> Always uses table with name "corr", will change ...
  */

  int i, ix;

  struct node *ndexe;
  struct node *nextnode;

  char name[NAME_L];
  char lname[NAME_L];
  char slname[NAME_L];
  char* uslname;

  char nname[NAME_L];
  char lnname[NAME_L];
  char slnname[NAME_L];
  char* uslnname;

  double xnew, ynew;

  /* set up pointers to current sequence for later use */
  struct sequence* mysequ = current_sequ;
  nextnode = mysequ->ex_start;
  ndexe = mysequ->ex_end;

  /* printf("Pointers: %d %d %d\n",mysequ,nextnode,ndexe); */

  (void) cmd;
  i = 1;
  ix = 0;

  while (ix == 0) {
    // loop over all elements in table corr (indexed with i)

    // 2013-Dec-04  10:48:48  ghislain: changed ix= into ix+=
    ix += string_from_table_row("corr", "name", &i, name);
    ix += double_from_table_row("corr", "px.correction", &i, &xnew);
    ix += double_from_table_row("corr", "py.correction", &i, &ynew);

    if (ix == 0) { // 2013-Dec-02  18:19:55  ghislain: what was that supposed to test ?
                   // Before: just the last return value above, Now: any of the above return values ?
      nextnode = mysequ->ex_start;
      while (nextnode != ndexe) {
	// loop over all elements in the sequence; look for the corrector "name"
	strcpy(lname, name);
	stolower(lname);
	strcpy(slname, strip(lname));
	uslname = supp_tb(slname);

	strcpy(nname, nextnode->name);
	strcpy(lnname, nname);
	stolower(lnname);
	strcpy(slnname, strip(lnname));
	uslnname = supp_tb(slnname);

	if (strcmp(uslname, uslnname) == 0) {
	  // found element "nname" matches the corrector "name"
	  // add the kicks read from corr table for corrector "name" to the
	  // already existing kick values for element "nname"
	  nextnode->chkick += xnew;
	  nextnode->cvkick += ynew;
	  // and call it quit for this corrector by placing ourselves at end of sequence
	  // which will break the loop over sequence nodes.
	  nextnode = ndexe;
	} else
	  {
	    // the "nname" element does not match the "name" corector; continue along the sequence.
	    nextnode = nextnode->next;
	  }
      }
    }
    i++;
  }
  return;
}

#if 0
static void correct_prtcorr(struct in_cmd* cmd)
{
  // 2013-Dec-04  11:06:35  ghislain: inserted routine received from W Herr.
  /*
    static char atm[6][4] = {"hmon","vmon","moni","hkic","vkic","kick"};
  */

  // int i, ix, j;

  char   *corrfil;
  char   corrfil1[100];

  struct table *seqcorr_table = NULL;

  struct node *ndexe;
  struct node *nextnode;

  strcpy(corrfil1,"\0");

  if(seqcorr_table == NULL) {
    seqcorr_table = make_table("seqcorr", "seqcorr", corr_table_cols, corr_table_types, 15000);
    //seqcorr_table = make_table("seqcorr", "seqcorr", seqcorr_table_cols, seqcorr_table_types, 15000);
    add_to_table_list(seqcorr_table, table_register);
  }

  /* set up pointers to current sequence for later use */
  struct sequence* mysequ = current_sequ;
  nextnode = mysequ->ex_start;
  ndexe = mysequ->ex_end;

  /* printf("Pointers: %d %d %d\n",mysequ,nextnode,ndexe); */

  nextnode = mysequ->ex_start;
  while (nextnode != ndexe) {
    if((strncmp(atc[0],nextnode->base_name,4) == 0) ||
       (strncmp(atc[1],nextnode->base_name,4) == 0) ||
       (strncmp(atc[2],nextnode->base_name,4) == 0)) {
      string_to_table_curr("seqcorr","name",nextnode->name);
      //string_to_table_curr("seqcorr","type",nextnode->base_name);
      //double_to_table_curr("seqcorr","s",&nextnode->position);
      //double_to_table_curr("seqcorr","hkick",&nextnode->chkick);
      double_to_table_curr("seqcorr","px.correction",&nextnode->chkick);
      //double_to_table_curr("seqcorr","vkick",&nextnode->cvkick);
      double_to_table_curr("seqcorr","py.correction",&nextnode->cvkick);
      augment_count("seqcorr");
    }
    nextnode = nextnode->next;
  }

  if ((corrfil = command_par_string("file",cmd->clone)) != NULL) {
    strcpy(corrfil1,corrfil);
  } else {
    strcpy(corrfil1,"seqcorr.out");
  }

  out_table("seqcorr",seqcorr_table,corrfil1);
  return;
}
#endif

static void correct_correct2(struct in_cmd* cmd)
/* Steering routine for orbit corrections of two beams */
{
  const char *rout_name = "correct_correct2";

  /*
    struct name_list* spos = sequences->list;
    struct table *twb1;
    struct table *twb2;
    int idrop;
    int pos;
    struct timeb tp;
    int sflag, svdflg;
    double  sigcut;
  */

  int ix, im, ip;
  int i, j, nnnseq;
  int imon, icor;
  int ncorr, nmon;
  int niter;
  int resout;

  int twism;
  int ifail;
  double rms;
  double rrms;
  double tmp1, tmp2, tmp3, tmp4;
  char *clist, *mlist;           /* file names for monitor and corrector output */
  char clist1[100], clist2[100]; /* file names for corrector output ring 1 and ring 2 */
  double *dmat = { NULL };       /* response matrix, double precision */
  double *corvec, *monvec;       /* vectors to hold measured orbit and correctors */
  double *resvec;                /* vector to hold corrected orbit */
  char *conm;                    /* vector to hold corrector names (for MICADO) */
  //  int     *sing;  // not used /* array to store pointer to singular correctors */
  int *nm, *nx, *nc;
  struct id_mic2 *c;
  //    struct id_mic2   *m;

  int debug = get_option("debug");

  /* If only Twiss summary is required prepare and write it */
  if ((twism = command_par_value("twissum", cmd->clone)) > 0) {

    if (!ftdata && !(ftdata = fopen("twiss.sum" , "w")))
      fatal_error("Cannot open file twiss.sum with write access", ", MAD-X terminates ");

    j = 1;
    if ((nnnseq = get_variable("n")) == 0)   nnnseq = twism;

    double_from_table_row("summ", "xcomax", &j, &tmp1);
    double_from_table_row("summ", "xcorms", &j, &tmp2);
    double_from_table_row("summ", "ycomax", &j, &tmp3);
    double_from_table_row("summ", "ycorms", &j, &tmp4);
    fprintf(ftdata, " T: %d %e %e %e %e\n", nnnseq, tmp1, tmp2, tmp3, tmp4);
    printf("TWISSUM: Data from twiss summary written to twiss.summ; aborting correction\n");
    fflush(ftdata);
    return; // abort the correction here
  }

  strcpy(clist1, "\0");
  strcpy(clist2, "\0");

  printf("for two beams orbit corrections ...\n");

  ip = pro_correct_getcommands(cmd);
  im = pro_correct2_gettables(ip, cmd);
  ncorr = im % 30000;
  nmon = im / 30000;
  printf("%d monitors and %d correctors found in input\n", nmon, ncorr);

  if (nmon == 0) {
    printf("No monitor found in input, no correction done\n");
    return;
  }

  if (ncorr == 0) {
    printf("No corrector found in input, no correction done\n");
    return;
  }

  /* Prepare file descriptors for the output */
  resout = command_par_value("resout", cmd->clone);
  if (resout > 0) {
    if (!fddata && !(fddata = fopen("corr.out" , "w")))
      fatal_error("Cannot open file corr.out with write access", ", MAD-X terminates ");
    if (!fcdata && !(fcdata = fopen("stren.out", "w")))
      fatal_error("Cannot open file stren.out with write access", ", MAD-X terminates ");
    if (!fgdata && !(fgdata = fopen("plot.orb", "w")))
      fatal_error("Cannot open file plot.orb with write access", ", MAD-X terminates ");
  }

  /* allocate vectors used by correction algorithms */
  nx =     mycalloc_atomic("correct_correct2_nx",     ncorr,    sizeof *nx);
  nc =     mycalloc_atomic("correct_correct2_nc",     ncorr,    sizeof *nc);
  nm =     mycalloc_atomic("correct_correct2_nm",     nmon,     sizeof *nm);
  //sing = mycalloc_atomic("correct_correct2_sing",   ncorr*2,  sizeof *sing); // not used
  corvec = mycalloc_atomic("correct_correct2_corvec", ncorr,    sizeof *corvec);
  monvec = mycalloc_atomic("correct_correct2_monvec", nmon,     sizeof *monvec);
  resvec = mycalloc_atomic("correct_correct2_resvec", nmon,     sizeof *resvec);
  conm =   mycalloc_atomic("correct_correct2_conm",   ncorr*16, sizeof *conm);

  /* get original settings of correctors from input Twiss-table */
  pro_correct2_getcorrs(cmd);

  /* get input orbit, default is from input Twiss-table */
  pro_correct2_getorbit(cmd);

  /* find and prepare enabled correctors and monitors, may be repeated */
  ix = pro_correct2_getactive(ip, nm, nx, nc, corvec, monvec, conm);
  icor = ix % 30000;
  imon = ix / 30000;
  printf("%d monitors and %d correctors enabled\n", imon, icor);

  if (debug) {
    for (i = 0; i < icor; i++)
      printf("C: %d %d \n", nx[i], nc[i]);
    for (i = 0; i < imon; i++)
      printf("M: %d %e \n", nm[i], monvec[i]);
  }

  // RING is only valid case
  if (strcmp("ring", command_par_string("flag", cmd->clone)) == 0) {
    if (dmat != NULL )
      myfree(rout_name, dmat);
    /* icor and imon used to set up correct matrix size !! */
    dmat = pro_correct2_response_ring(ip, icor, imon);
  } else {
    fatal_error("Invalid machine type other than RING with option TWOBEAM", ", MAD-X terminates ");
  }

  /* MICADO correction, get desired number of correctors from command */
  corrl = command_par_value("corrlim", cmd->clone);
  set_variable("corrlim", &corrl);

  if (strcmp("micado", command_par_string("mode", cmd->clone)) == 0) {
    printf("enter MICADO correction ...\n");

    niter = command_par_value("ncorr", cmd->clone);
    if (niter == 0) {
      printf("Requested %d correctors (\?\?\?) set to %d\n", niter, icor);
      niter = icor;
    } else if (niter < 0) {
      printf("Requested %d correctors (\?\?\?) set to 0\n", niter);
      niter = 0;
    } else if (niter > icor) {
      printf("Fewer correctors available than requested by ncorr\n");
      printf("you want %d,  you get %d\n", niter, icor);
      printf("ncorr reset to %d\n", icor);
      niter = icor;
    }

    rms = 1000.0 * command_par_value("error", cmd->clone);

    ifail = c_micit(dmat, conm, monvec, corvec, resvec, nx, rms, imon, icor, niter);
    if (debug) printf("Back from micado %d\n", ifail);

    if (ifail != 0) {
      printf("MICADO correction completed with error code %d\n\n", ifail);
      warning("MICADO back with error", ", no correction done");
    }

    rrms = crms(monvec, imon);
    printf("RMS before %e\n", rrms);
    rrms = crms(resvec, imon);
    printf("RMS after  %e\n", rrms);

    if (fgdata != NULL ) {
      for (i = 0; i < nmon; i++) {
	fprintf(fgdata, "%e %e \n", monvec[i], resvec[i]);
      }
      fflush(fgdata);
    }

    c = correct_orbit12->cor_table;
    for (i = 0; i < icor; i++) {
      printf("%s %e\n", c[nc[i]].p_node->name, corvec[nx[i] - 1]);
    }
    printf("\n");

    if (ifail != 0) {
      printf("MICADO correction completed with error code %d\n\n", ifail);
      warning("MICADO back with error", ", no correction done");
    }
    if (ifail == 0) {
      pro_correct2_write_results(monvec, resvec, corvec, nx, nc, nm, imon,
				 icor, ip, resout);
    }
  }

  /* write corrector output to tfs table */
  if ((clist = command_par_string("clist", cmd->clone)) != NULL ) {
    strcat(clist1, clist); strcat(clist1, "_1");
    strcat(clist2, clist); strcat(clist2, "_2");
    out_table("corr1", corr_table1, clist1);
    out_table("corr2", corr_table2, clist2);
  }

  /* write monitor output to tfs table */
  if ((mlist = command_par_string("mlist", cmd->clone)) != NULL ) {
    out_table("mon", mon_table, mlist);
  }

  /* Clean up at the end of the module */
  myfree(rout_name, nm);
  myfree(rout_name, dmat);
  myfree(rout_name, nx);
  myfree(rout_name, nc);
  myfree(rout_name, corvec);
  myfree(rout_name, monvec);
  myfree(rout_name, resvec);
  myfree(rout_name, conm);
}

static int pro_correct2_gettables(int iplane, struct in_cmd* cmd) {
  const char *rout_name = "pro_correct2_gettables";

  struct id_mic2 *cor_l1, *cor_l2;
  struct id_mic2 *mon_l1, *mon_l2;
  struct id_mic2 *cor_l12, *mon_l12;
  struct id_mic2 *prt;

  // struct table *ttb; // not used

  struct table *b1 = NULL;
  struct table *b2 = NULL;

  char* orbtab1;
  char* orbtab2;

  int ebl1, ebl2;

  int j, k;

  int cntm1 = 0;
  int cntc1 = 0;
  int cntm2 = 0;
  int cntc2 = 0;
  int cntm12 = 0;
  int cntc12 = 0;

  double ounits;

  /*
    static char atm[6][4] = {"hmon","vmon","moni","hkic","vkic","kick"};
  */

  /* Get access to tables, for orbit and model the default is twiss_table */
  if ((orbtab1 = command_par_string("beam1tab", cmd->clone)) != NULL ) {
    printf("Want to use orbit from: %s\n", orbtab1);
    if (!(b1 = find_table(orbtab1))) {
      fatal_error("Beam 1 ORBIT table requested, but not provided:", orbtab1);
    }
  } else {
    // Jun 25, 2013 2:41:26 PM ghislain : FIXME - ??? empty else statement
  }

  if ((orbtab2 = command_par_string("beam2tab", cmd->clone)) != NULL ) {
    printf("Want to use orbit from: %s\n", orbtab2);
    if (!(b2 = find_table(orbtab2))) {
      fatal_error("Beam 2 ORBIT table requested, but not provided:", orbtab2);
    }
  } else {
    // Jun 25, 2013 2:41:01 PM ghislain : FIXME - ??? empty else statement
  }

  /* store as globals for later use */
  if ((b1 != NULL) && (b2 != NULL)){
    twiss_table_beam1 = b1;
    twiss_table_beam2 = b2;
  } else {
    fatal_error("Beam 1 and 2 orbit tables not found:",orbtab1);
  }

  /* reserve space for orbit correction structures */
  if (correct_orbit12 == NULL )
    correct_orbit12 = mycalloc("pro_correct2_gettables", 1, sizeof *correct_orbit12);

  if (correct_orbit12->cor_table != NULL )
    myfree(rout_name, correct_orbit12->cor_table);
  if (correct_orbit12->mon_table != NULL )
    myfree(rout_name, correct_orbit12->mon_table);

  correct_orbit12->cor_table = mycalloc("pro_correct2_gettables_cor",30200, sizeof *correct_orbit12->cor_table);
  correct_orbit12->mon_table = mycalloc("pro_correct2_gettables_mon",30200, sizeof *correct_orbit12->mon_table);

  /* orbit table available, get units, if defined */
  if ((ounits = command_par_value("units", cmd->clone)) > 0)
    correct_orbit12->units = ounits;
  else
    correct_orbit12->units = 1.0;

  // ttb = model_table; // not used
  /* no more need, we have b1 and b2 as pointers .. */

  correct_orbit12->mon_table->previous = NULL;
  correct_orbit12->mon_table->next = NULL;
  correct_orbit12->cor_table->previous = NULL;
  correct_orbit12->cor_table->next = NULL;

  mon_l1 = correct_orbit12->mon_table;
  cor_l1 = correct_orbit12->cor_table;

  for (j = 0; j < b1->curr; j++) {
    if ((strncmp(atm[iplane - 1], b1->p_nodes[j]->base_name, 4) == 0)
	|| (strncmp(atm[2], b1->p_nodes[j]->base_name, 4) == 0)) {
      /*    printf("1m: %s %ld\n", b1->p_nodes[j]->name, strstr(".b2", b1->p_nodes[j]->name)); */
      if (strstr(b1->p_nodes[j]->name, ".b1") != NULL ) {
	mon_l1->id_ttb[0] = j;
	mon_l1->id_ttb[1] = -1;
	mon_l1->enable = b1->p_nodes[j]->enable;
	mon_l1->p_node = b1->p_nodes[j];
	mon_l1->next = mon_l1;
	mon_l1->next++;
	mon_l1++;
	cntm1++;
      } else {
	/*      printf("Removed: %s\n",b1->p_nodes[j]->name); */
      }
    }

    if ((strncmp(atc[iplane - 1], b1->p_nodes[j]->base_name, 4) == 0)
	|| (strncmp(atc[2], b1->p_nodes[j]->base_name, 4) == 0)) {
      /*    printf("1c: %s %ld\n", b1->p_nodes[j]->name, b1->p_nodes[j]->name); */
      if (strstr(b1->p_nodes[j]->name, ".b1") != NULL ) {
	cor_l1->id_ttb[0] = j;
	cor_l1->id_ttb[1] = -1;
	cor_l1->enable = b1->p_nodes[j]->enable;
	cor_l1->p_node = b1->p_nodes[j];
	cor_l1->p_node_s1 = b1->p_nodes[j];
	cor_l1->p_node_s2 = NULL;

	if (command_par_value("corzero", cmd->clone) > 0) {
	  if (iplane == 1)
	    cor_l1->p_node_s1->chkick = 0.0;
	  if (iplane == 2)
	    cor_l1->p_node_s1->cvkick = 0.0;
	}

	cor_l1->next = cor_l1;
	cor_l1->next++;
	cor_l1++;
	cntc1++;
      } else {
	/* printf("Removed: %s\n",b1->p_nodes[j]->name); */
      }
    }
  }

  mon_l2 = mon_l1;
  cor_l2 = cor_l1;

  for (j = 0; j < b2->curr; j++) {
    if ((strncmp(atm[iplane - 1], b2->p_nodes[j]->base_name, 4) == 0)
	|| (strncmp(atm[2], b2->p_nodes[j]->base_name, 4) == 0)) {
      /* printf("2m: %s %ld\n", b2->p_nodes[j]->name, b2->p_nodes[j]->name); */
      if (strstr(b2->p_nodes[j]->name, ".b2") != NULL ) {
	mon_l2->id_ttb[0] = -1;
	mon_l2->id_ttb[1] = j;
	mon_l2->enable = b2->p_nodes[j]->enable;
	mon_l2->p_node = b2->p_nodes[j];
	mon_l2->next = mon_l2;
	mon_l2->next++;
	mon_l2++;
	cntm2++;
      } else {
	/* printf("Removed: %s\n",b2->p_nodes[j]->name); */
      }
    }

    if ((strncmp(atc[iplane - 1], b2->p_nodes[j]->base_name, 4) == 0)
  || (strncmp(atc[2], b2->p_nodes[j]->base_name, 4) == 0)) {
      /*    printf("2c: %s %ld\n", b2->p_nodes[j]->name, b2->p_nodes[j]->name); */
      if (strstr(b2->p_nodes[j]->name, ".b2") != NULL ) {
	cor_l2->id_ttb[0] = -1;
	cor_l2->id_ttb[1] = j;
	cor_l2->enable = b2->p_nodes[j]->enable;
	cor_l2->p_node = b2->p_nodes[j];
	cor_l2->p_node_s2 = b2->p_nodes[j];
	cor_l2->p_node_s1 = NULL;

	if (command_par_value("corzero", cmd->clone) > 0) {
	  if (iplane == 1)
	    cor_l2->p_node_s2->chkick = 0.0;
	  if (iplane == 2)
	    cor_l2->p_node_s2->cvkick = 0.0;
	}

	cor_l2->next = cor_l2;
	cor_l2->next++;
	cor_l2++;
	cntc2++;
      } else {
	/* printf("Removed: %s\n",b2->p_nodes[j]->name); */
      }
    }
  }

  mon_l12 = mon_l2;
  cor_l12 = cor_l2;
  for (j = 0; j < b1->curr; j++) {
    if ((strncmp(atm[iplane - 1], b1->p_nodes[j]->base_name, 4) == 0)
	|| (strncmp(atm[2], b1->p_nodes[j]->base_name, 4) == 0)) {
      /*    printf("12m: %s \n", b1->p_nodes[j]->name); */
      if ((strstr(b1->p_nodes[j]->name, ".b1") == NULL )&&
	  (strstr(b1->p_nodes[j]->name,".b2") == NULL)){
	mon_l12->id_ttb[0] = j;
	for (k=0; k < b2->curr; k++) {
	  if(strcmp(b2->p_nodes[k]->name,b1->p_nodes[j]->name) == 0) {
	    mon_l12->id_ttb[1] = k;
	  }
	}
	mon_l12->enable = b1->p_nodes[j]->enable;
	mon_l12->p_node = b1->p_nodes[j];
	mon_l12->next = mon_l12;
	mon_l12->next++; mon_l12++;
	cntm12++;
      } else {
	/* printf("Removed: %s\n",b1->p_nodes[j]->name); */
      }
    }

    if((strncmp(atc[iplane-1],b1->p_nodes[j]->base_name,4) == 0) ||
       (strncmp(atc[2], b1->p_nodes[j]->base_name,4) == 0)) {
      /* printf("12c: %s \n", b1->p_nodes[j]->name);     */
      if((strstr(b1->p_nodes[j]->name,".b1") == NULL) &&
	 (strstr(b1->p_nodes[j]->name,".b2") == NULL)) {
	cor_l12->id_ttb[0] = j;
	for (k=0; k < b2->curr; k++) {
	  if(strcmp(b2->p_nodes[k]->name,b1->p_nodes[j]->name) == 0) {
	    cor_l12->id_ttb[1] = k;
	  }
	}
	cor_l12->p_node = b1->p_nodes[j];
	cor_l12->p_node_s1 = b1->p_nodes[cor_l12->id_ttb[0]];
	cor_l12->p_node_s2 = b2->p_nodes[cor_l12->id_ttb[1]];
	ebl1 = b1->p_nodes[cor_l12->id_ttb[0]]->enable;
	ebl2 = b2->p_nodes[cor_l12->id_ttb[1]]->enable;
	cor_l12->enable = ebl1*ebl2;

	if (command_par_value("corzero", cmd->clone) > 0) {
	  if(iplane == 1) cor_l12->p_node_s1->chkick = 0.0;
	  if(iplane == 2) cor_l12->p_node_s1->cvkick = 0.0;
	  if(iplane == 1) cor_l12->p_node_s2->chkick = 0.0;
	  if(iplane == 2) cor_l12->p_node_s2->cvkick = 0.0;
	}

	cor_l12->next = cor_l12;
	cor_l12->next++; cor_l12++;
	cntc12++;
      } else {
	/* printf("Removed: %s\n",b1->p_nodes[j]->name); */
      }
    }
  }

  /* terminate linked list   */
  mon_l12--;
  mon_l12->next = NULL;
  cor_l12--;
  cor_l12->next = NULL;

  printf("mons and corrs (beam 1)   : %ld %ld\n", (long int) cntm1,  (long int) cntc1);
  printf("mons and corrs (beam 2)   : %ld %ld\n", (long int) cntm2,  (long int) cntc2);
  printf("mons and corrs (beam 1+2) : %ld %ld\n", (long int) cntm12, (long int) cntc12);

  if (get_option("debug")) {
    prt = correct_orbit12->mon_table;
    while (prt != NULL ) {
      printf("Monitors beam12: %s %ld %ld\n", prt->p_node->name,
	     (long int) prt->id_ttb[0], (long int) prt->id_ttb[1]);
      prt = prt->next;
    }

    prt = correct_orbit12->cor_table;
    while (prt != NULL ) {
      printf("Correctors beam12: %s %ld %ld\n", prt->p_node->name,
	     (long int) prt->id_ttb[0], (long int) prt->id_ttb[1]);
      prt = prt->next;
    }
  }

  if (corr_table1 == NULL ) {
    corr_table1 = make_table("corr1", "corr1", corr_table_cols, corr_table_types, 15000);
    add_to_table_list(corr_table1, table_register);
  }

  if (corr_table2 == NULL ) {
    corr_table2 = make_table("corr2", "corr2", corr_table_cols, corr_table_types, 15000);
    add_to_table_list(corr_table2, table_register);
  }

  pro_correct2_make_corr_table();

  if (mon_table == NULL ) {
    mon_table = make_table("mon", "mon", mon_table_cols, mon_table_types, 15000);
    add_to_table_list(mon_table, table_register);
    pro_correct2_make_mon_table();
  }

  // 2013-Jun-24  12:21:49  ghislain:
  // following is a kludge to return a single value but has to be decoded on other side.
  if( cntc1+cntc2+cntc12 >= 30000)
    fatal_error("Found more than 30000 correctors; decoding in mad_orbit.c will fail",
    "Please report this issue to MAD developpers (mad@cern.ch)");
  return 30000 * (cntm1 + cntm2 + cntm12) + cntc1 + cntc2 + cntc12;
}

static int pro_correct2_getorbit(struct in_cmd* cmd) {
  struct id_mic2 *m; /* access to tables for monitors and correctors */
  double **da1;
  double **da2;
  double xlimit;

  char strx[40];
  char stry[40];

  int posx, posy, pospx, pospy;

  int debug = get_option("debug");

  da1 = twiss_table_beam1->d_cols;
  da2 = twiss_table_beam2->d_cols;

  m = correct_orbit12->mon_table;

  strcpy(strx, "x");
  strcpy(stry, "y");

  if ((posx = name_list_pos(strx, twiss_table_beam1->columns)) < 0)
    fatal_error("orbit x not found in input table", ", MAD-X terminates ");

  if ((posy = name_list_pos(stry, twiss_table_beam1->columns)) < 0)
    fatal_error("orbit y not found in input table", ", MAD-X terminates ");

  if (debug) {
    if ((pospx = name_list_pos("px", twiss_table_beam1->columns)) < 0)
      warning("orbit px not found in input table", ", MAD-X continues ");

    if ((pospy = name_list_pos("py", twiss_table_beam1->columns)) < 0)
      warning("orbit py not found in input table", ", MAD-X continues ");

    printf("====c1===>  %d %d %d %d \n", posx, posy, pospx, pospy);
  }

  while (m) {

    /* If correction to target orbit, subtract the wanted orbit ... */
    if (m->id_ttb[0] > 0) {
      m->val.before[0] = m->p_node->other_bv * da1[9] [m->id_ttb[0]] * 1000;
      m->val.before[1] = m->p_node->other_bv * da1[11][m->id_ttb[0]] * 1000;
    } else if (m->id_ttb[1] > 0) {
      m->val.before[0] = m->p_node->other_bv * da2[9] [m->id_ttb[1]] * 1000;
      m->val.before[1] = m->p_node->other_bv * da2[11][m->id_ttb[1]] * 1000;
    } else {
      fatal_error("Unforeseen case in pro_correct2_getorbit", ", MAD-X terminates ");
    }

    if (par_present("monon", cmd->clone)) {
      xlimit = command_par_value("monon", cmd->clone);
      if (frndm() > xlimit) {
	m->enable = 0;
	printf("Monitor %s disabled\n", m->p_node->name);
      }
    }

    if (debug) {
      printf("m-list: %d %d %s %s\n", m->id_ttb[0], m->id_ttb[1], m->p_node->name, m->p_node->base_name);
      printf("initial reading: %e %e\n\n", m->val.before[0], m->val.before[1]);
    }

    m = m->next;
  }
  return (0);
}

static int pro_correct2_getcorrs(struct in_cmd* cmd) {
  struct id_mic2 *c; /* access to tables for monitors and correctors */

  (void) cmd;

  int debug = get_option("debug");

  c = correct_orbit12->cor_table;
  while (c) {
    if (c->id_ttb[0] > 0) {
      c->val.before[0] = c->p_node_s1->chkick * 1000.;
      c->val.before[1] = c->p_node_s1->cvkick * 1000.;
    } else if (c->id_ttb[1] > 0) {
      c->val.before[0] = c->p_node_s2->chkick * 1000.;
      c->val.before[1] = c->p_node_s2->cvkick * 1000.;
    }

    if (debug) {
      printf("c-list: %d %d %s %s\n", c->id_ttb[0], c->id_ttb[1], c->p_node->name, c->p_node->base_name);
      printf("initial strengths: %e %e\n", c->val.before[0], c->val.before[1]);
    }

    c = c->next;
  }

  return (0);
}

static int pro_correct2_getactive(int ip, int *nm, int *nx, int *nc,
          double *corvec, double *monvec, char *conm) {
  int imon, icor;
  int imona, icora;
  struct id_mic2 *m, *c;

  int debug = get_option("debug");

  m = correct_orbit12->mon_table;
  imon = 0;
  imona = 0;

  while (m) {
    if (debug) {
      printf("from list: %d %d %s %s\n", m->id_ttb[0], m->id_ttb[1], m->p_node->name, m->p_node->base_name);
      printf("orbit readings: %d %f %f\n", ip, m->val.before[0], m->val.before[1]);
    }
    if (m->enable == 1) {
      monvec[imon] = m->val.before[ip - 1];
      nm[imon] = imona;
      imon++;
    }
    imona++;
    m = m->next;
  };

  c = correct_orbit12->cor_table;
  icor = 0;
  icora = 0;

  while (c) {
    if (debug) {
      printf("from list: %d %d %d %s %s\n", c->enable, c->id_ttb[0], c->id_ttb[1], c->p_node->name, c->p_node->base_name);
      printf("kicker readings: %f %f\n", c->val.before[0], c->val.before[1]);
    }
    if (c->enable == 1) {
      corvec[icor] = c->val.before[ip - 1];
      nx[icor] = icora;
      nc[icor] = icora;
      strcpy(conm, c->p_node->name);
      conm += 16;
      /* printf("nc: %d %d \n",icor,nc[icor]); */
      icor++;
    }
    icora++;
    c = c->next;
  };

  // 2013-Jun-24  12:21:49  ghislain:
  // following is a kludge to return a single value but has to be decoded on other side.
  if(icor >= 30000)
    fatal_error("Found more than 30000 correctors; decoding in mad_orbit.c will fail",
    "Please report this issue to MAD developpers (mad@cern.ch)");

  return (30000 * imon + icor);
}

static double* pro_correct2_response_ring(int ip, int nc, int nm) {
  int ic, im;
  struct id_mic2 *m, *c; /* access to tables for monitors and correctors */

  double **da1;
  double **da2;
  double bx_c, by_c, pix_c, piy_c;
  double bx_m, by_m, pix_m, piy_m;
  double qx0, qy0;
  double respx, respy;
  double *dmat;
  int *imat;
  int mp;
  int i_zero, i_one;
  int icb;
  int i, j;

  int debug = get_option("debug");

  ic = 0;
  im = 0;
  i_zero = 0;
  i_one = 1;

  da1 = twiss_table_beam1->d_cols;
  da2 = twiss_table_beam2->d_cols;

  dmat = mycalloc_atomic("pro_correct2_response_ring", nc*nm, sizeof *dmat);
  imat = mycalloc_atomic("pro_correct2_response_ring", nc*nm, sizeof *imat);

  /* initialize imat: */
  for (i = 0; i < nc; i++) {
    for (j = 0; j < nm; j++) {
      setupi_(&i_zero, imat, &j, &i, &nm, &nc);
    }
  }

  c = correct_orbit12->cor_table;
  ic = 0;

  while (c) {
    if (debug)
      printf("corrector flag: %d\n", c->enable);

    if (c->enable == 1) {

      for (icb = 0; icb < 2; icb++) { // correcting for two beams respectively 0 and 1
	if (c->id_ttb[icb] > 0) {

	  if (icb == 0) { // beam1
	    correct_orbit12->qx0 = da1[5][twiss_table_beam1->curr - 1];
	    correct_orbit12->qy0 = da1[8][twiss_table_beam1->curr - 1];
	    qx0 = correct_orbit12->qx0;
	    qy0 = correct_orbit12->qy0;
	    if (c->id_ttb[icb] > 0) {
	      bx_c = da1[3][c->id_ttb[icb]];
	      by_c = da1[6][c->id_ttb[icb]];
	      pix_c = da1[5][c->id_ttb[icb]];
	      piy_c = da1[8][c->id_ttb[icb]];
	    } else {
	      bx_c = 0.0;
	      by_c = 0.0;
	      pix_c = 0.0;
	      piy_c = 0.0;
	    }
	  } else { // beam2
	    correct_orbit12->qx0 = da2[5][twiss_table_beam2->curr - 1];
	    correct_orbit12->qy0 = da2[8][twiss_table_beam2->curr - 1];
	    qx0 = correct_orbit12->qx0;
	    qy0 = correct_orbit12->qy0;
	    if (c->id_ttb[icb] > 0) {
	      bx_c = da2[3][c->id_ttb[icb]];
	      by_c = da2[6][c->id_ttb[icb]];
	      pix_c = da2[5][c->id_ttb[icb]];
	      piy_c = da2[8][c->id_ttb[icb]];
	    } else {
	      bx_c = 0.0;
	      by_c = 0.0;
	      pix_c = 0.0;
	      piy_c = 0.0;
	    }
	  }

	  m = correct_orbit12->mon_table;
	  im = 0;

	  while (m) {
	    if (debug)
	      printf("monitor flag: %d\n", m->enable);
	    if (m->enable == 1) {
	      if ((m->id_ttb[icb] > 0) && (c->id_ttb[icb] > 0)) {
		if (m->id_ttb[icb] > 0) {
		  if (icb == 0) { // beam1
		    mp = m->id_ttb[icb];
		    bx_m = da1[3][mp];
		    by_m = da1[6][mp];
		    pix_m = da1[5][mp];
		    piy_m = da1[8][mp];
		  } else { // beam2
		    mp = m->id_ttb[icb];
		    bx_m = da2[3][mp];
		    by_m = da2[6][mp];
		    pix_m = da2[5][mp];
		    piy_m = da2[8][mp];
		  }
		} else {
		  bx_m = 0.0;
		  by_m = 0.0;
		  pix_m = 0.0;
		  piy_m = 0.0;
		}

		respx = 0.0;
		respy = 0.0;

		/*  print Twiss parameters ... */
		if (debug) {
		  printf("%s %d %e %e %e %e -- %s %e %e %e %e\n",
			 c->p_node->name, icb, bx_c, by_c,
			 pix_c, piy_c, m->p_node->name, bx_m,
			 by_m, pix_m, piy_m);
		}

		if (ip == 1) { // x plane
		  respx = sqrt(bx_m * bx_c) / (2.0 * sin(pi * qx0) ) * cos((fabs(pix_m - pix_c) * twopi) - qx0 * pi);
		  setup_(&respx, dmat, &im, &ic, &nm, &nc);
		} else if (ip == 2) { // y plane
		  respy = sqrt(by_m * by_c) / (2.0 * sin(pi * qy0)) * cos((fabs(piy_m - piy_c) * twopi) - qy0 * pi);
		  setup_(&respy, dmat, &im, &ic, &nm, &nc);
		}

		if ((fabs(respy) > 0.000006) || (fabs(respx) > 0.000006)) {
		  if (debug) printf("true %d %d\n", ic, im);
		  setupi_(&i_one, imat, &im, &ic, &nm, &nc);
		} else {
		  if (debug) printf("false \n");
		  setupi_(&i_zero, imat, &im, &ic, &nm, &nc);
		}

		if (debug) printf("Response:  %d %d %e %e %e \n", ic, im, respx, respy, fabs(respy));

	      }
	      im++;
	    }
	    m = m->next;
	  }
	}
      }
      ic++;
    }
    c = c->next;
  }

  if (debug) {
    printf("\n");
    primat_(imat, &nm, &nc);
    printf("\n");
    prdmat_(dmat, &nm, &nc);
    printf("\n");
    printf("\n");
  }

  myfree("pro_correct2_response_ring", imat);

  return dmat;
}

static void pro_correct2_write_results(double *monvec, double *resvec, double *corvec,
				       int *nx, int *nc, int *nm, int imon, int icor, int ip,
				       int resout)
{
  /*                                              */
  /* Writes a summary of the correction           */
  /* Writes correctors strengths into sequences   */
  /* Fills TFS tables for correctors and monitors */
  /* Fills the 'stren.out' output                 */
  /* Makes various prints on request              */
  /*                                              */
  int i;
  int rst;
  double corrm;
  struct id_mic2 *m, *c; /* access to tables for monitors and correctors */

  m = correct_orbit12->mon_table;
  c = correct_orbit12->cor_table;

  if (fddata != NULL ) {
    rst = get_variable("n");
    fprintf(fddata, "%d %d %e %e %e %e %e %e\n", ip, rst,
      cptp(monvec, imon), cptp(resvec, imon), crms(monvec, imon),
      crms(resvec, imon), copk(monvec, imon), copk(resvec, imon));
  }

  if (print_correct_opt > 0) {
    printf("CORRECTION SUMMARY:   \n\n");
    printf("                   average [mm]  std.dev. [mm]      RMS [mm]        peak-to-peak [mm]\n\n");
    printf("before correction: %f        %f          %f        %f \n",
	   caverage(monvec,imon), cstddev(monvec,imon), crms(monvec, imon), cptp(monvec, imon));
    printf("after correction:  %f        %f          %f        %f \n\n\n",
	   caverage(resvec,imon), cstddev(resvec,imon), crms(resvec, imon), cptp(resvec, imon));
  }

  if (print_correct_opt > 1) {
    printf("Monitor:  Before:     After:    Difference:\n");
    printf("           (mm)        (mm)         (mm)   \n");
  }

  for (i = 0; i < imon; i++) {
    if (print_correct_opt > 1)
      printf("%s   %-4.3f     %-4.3f     %-4.3f\n",
	     m[nm[i]].p_node->name, monvec[i], resvec[i], resvec[i] - monvec[i]);

    m[nm[i]].val.after[ip - 1] = resvec[i];
    pro_correct2_fill_mon_table(ip, m[nm[i]].p_node->name, monvec[i], resvec[i]);
  }

  corrm = copk(corvec, icor);

  if (corrm > corrl) {
    printf("Max strength: %e should be less than corrector strength limit: %e\n", corrm, corrl);
    warning("maximum corrector strength larger than limit.","");
  } else {
    printf("Max strength: %e is below corrector strength limit: %e\n", corrm, corrl);
  }

  set_variable("corrmax", &corrm);

  if (print_correct_opt > 1) {
    printf("Max strength: %e\n", copk(corvec, icor));
    printf("Corrector:  Before:     After:    Difference:\n");
    printf("             (mrad)     (mrad)       (mrad)  \n");
  }

  for (i = 0; i < icor; i++) { /* loop over all correctors */

    c[nc[i]].val.after[ip - 1] = corvec[nx[i] - 1];
    if (print_correct_opt > 1) {
      printf("%s %-3.6f %-3.6f %-3.6f\n", c[nc[i]].p_node->name,
       c[nc[i]].val.before[ip - 1],
       corvec[nx[i] - 1] + c[nc[i]].val.before[ip - 1],
       corvec[nx[i] - 1]);
    }

    if (ip == 1) {
      /* Fill horizontal corrections for beam 1  */
      if (c[nc[i]].id_ttb[0] > 0) {
	c[nc[i]].p_node_s1->chkick += c[nc[i]].p_node_s1->other_bv * 0.001 * corvec[nx[i] - 1];
	pro_correct2_fill_corr_table(0, ip, c[nc[i]].p_node->name, c[nc[i]].val.before[ip - 1] * 0.001,
				     c[nc[i]].p_node_s1->chkick);
	/* ??? c[nc[i]].p_node_s1->other_bv*0.001*corvec[nx[i]-1]); */
	if (fcdata != NULL ) {
	  fprintf(fcdata, "%s->hkick = %e; \t! [1] %d\n", strip(c[nc[i]].p_node->name),
		  c[nc[i]].p_node_s1->other_bv * 0.001 * corvec[nx[i] - 1], resout);
	}
      }
      /* Fill horizontal corrections for beam 2  */
      if (c[nc[i]].id_ttb[1] > 0) {
	c[nc[i]].p_node_s2->chkick += 0.001 * corvec[nx[i] - 1];
	pro_correct2_fill_corr_table(1, ip, c[nc[i]].p_node->name, c[nc[i]].val.before[ip - 1] * 0.001,
				     c[nc[i]].p_node_s2->chkick);
	/* ??? c[nc[i]].p_node_s2->other_bv*0.001*corvec[nx[i]-1]); */
	if (fcdata != NULL ) {
	  fprintf(fcdata, "%s->hkick = %e; \t! [2] %d\n", strip(c[nc[i]].p_node->name),
		  0.001 * corvec[nx[i] - 1], resout);
	}
      }

    } else if (ip == 2) {
      /* Fill vertical corrections for beam 1  */
      if (c[nc[i]].id_ttb[0] > 0) {
	c[nc[i]].p_node_s1->cvkick += c[nc[i]].p_node_s1->other_bv * 0.001 * corvec[nx[i] - 1];
	pro_correct2_fill_corr_table(0, ip, c[nc[i]].p_node->name, c[nc[i]].val.before[ip - 1] * 0.001,
				     c[nc[i]].p_node_s1->cvkick);
	/* ??? c[nc[i]].p_node_s1->other_bv*0.001*corvec[nx[i]-1]); */
	if (fcdata != NULL ) {
	  fprintf(fcdata, "%s->vkick = %e; \t! [1] %d\n", strip(c[nc[i]].p_node->name),
		  c[nc[i]].p_node_s1->other_bv * 0.001 * corvec[nx[i] - 1], resout);
	}
      }
      if (c[nc[i]].id_ttb[1] > 0) {
	/* Fill vertical corrections for beam 2  */
	c[nc[i]].p_node_s2->cvkick += 0.001 * corvec[nx[i] - 1];
	pro_correct2_fill_corr_table(1, ip, c[nc[i]].p_node->name, c[nc[i]].val.before[ip - 1] * 0.001,
				     c[nc[i]].p_node_s2->cvkick);
	/* ??? c[nc[i]].p_node_s2->other_bv*0.001*corvec[nx[i]-1]); */
	if (fcdata != NULL ) {
	  fprintf(fcdata, "%s->vkick = %e; \t! [2] %d\n", strip(c[nc[i]].p_node->name),
		  0.001 * corvec[nx[i] - 1], resout);
	}
      }
    }

  } /* end loop over correctors */

  if (fcdata != NULL ) fflush(fcdata);
  if (fddata != NULL ) fflush(fddata);

}

static void correct_correct1(struct in_cmd* cmd)
/* Steering routine for orbit corrections of one beam */
{
  const char *rout_name = "correct_correct";
  int ix, im, ip, idrop;
  int j, nnnseq;
  int imon, icor;
  int ncorr, nmon;
  int niter;
  int resout;

  int twism;
  int ifail, sflag;
  double rms;
  double sngcut, sngval;
  double tmp1, tmp2, tmp3, tmp4;
  double sigcut;           /* number of sigmas (normalized) for filter cut */
  char *clist, *mlist;     /* file names for monitor and corrector output */
  double *dmat = { NULL }; /* response matrix, double precision */
  double *monvec, *corvec; /* vectors to hold measured orbit and correctors */
  double *resvec;          /* vector to hold corrected orbit (result) */
  char *conm;              /* vector to hold corrector names (for MICADO) */
  int *sing;               /* array to store pointer to singular correctors */
  int *nm, *nx, *nc;
  struct id_mic *corl;

  int debug = get_option("debug");

  /* If only Twiss summary is required prepare and write it */
  if ((twism = command_par_value("twissum", cmd->clone)) > 0) {
    if (ftdata == NULL ) {
      if ((ftdata = fopen("twiss.summ", "w")) == NULL )
	fatal_error("Cannot open file twiss.summ with write access", ", MAD-X terminates ");
    }
    j = 1;
    if ((nnnseq = get_variable("n")) == 0)
      nnnseq = twism;

    double_from_table_row("summ", "xcomax", &j, &tmp1);
    double_from_table_row("summ", "xcorms", &j, &tmp2);
    double_from_table_row("summ", "ycomax", &j, &tmp3);
    double_from_table_row("summ", "ycorms", &j, &tmp4);
    fprintf(ftdata, " T: %d %e %e %e %e\n", nnnseq, tmp1, tmp2, tmp3, tmp4);
    printf("TWISSUM: Data from twiss summary written to twiss.summ; aborting correction\n");
    fflush(ftdata);
    return; // abort the correction here
  }

  ip = pro_correct_getcommands(cmd);
  im = pro_correct_gettables(ip, cmd);
  ncorr = im % 30000;
  nmon = im / 30000;
  printf("%d monitors and %d correctors found in input\n", nmon, ncorr);

  if (nmon == 0) {
    printf("No monitor found in input, no correction done\n");
    return;
  }

  if (ncorr == 0) {
    printf("No corrector found in input, no correction done\n");
    return;
  }

  /* Prepare file descriptors for the output */
  resout = command_par_value("resout", cmd->clone);
  if (resout > 0) {
    if (!fddata && !(fddata = fopen("corr.out" , "w")))
      fatal_error("Cannot open file corr.out with write access", ", MAD-X terminates ");
    if (!fcdata && !(fcdata = fopen("stren.out", "w")))
      fatal_error("Cannot open file stren.out with write access", ", MAD-X terminates ");
  }

  /* allocate vectors used by correction algorithms */
  nx =     mycalloc("correct_correct_nx",    ncorr,   sizeof(int));
  nc =     mycalloc("correct_correct_nc",    ncorr,   sizeof(int));
  nm =     mycalloc("correct_correct_nm",    nmon,    sizeof(int));
  sing =   mycalloc("correct_correct_sing",  ncorr*2, sizeof(int));
  corvec = mycalloc("correct_correct_corvec",ncorr,   sizeof(double));
  monvec = mycalloc("correct_correct_monvec",nmon,    sizeof(double));
  resvec = mycalloc("correct_correct_resvec",nmon,    sizeof(double));
  conm =   mycalloc("correct_correct_conm",  ncorr*16,sizeof(char));

  /* get original settings of correctors from input Twiss-table */
  pro_correct_getcorrs(cmd);

  /* get input orbit, default is from input Twiss-table */
  /* if flag "extern" is true: can be from external table */
  if (command_par_value("extern", cmd->clone))
    pro_correct_getorbit_ext(cmd);
  else
    pro_correct_getorbit(cmd);

  /* find and prepare enabled correctors and monitors, may be repeated */
  ix = pro_correct_getactive(ip, nm, nx, nc, corvec, monvec, conm);
  icor = ix % 30000;
  imon = ix / 30000;
  printf("%d monitors and %d correctors enabled\n", imon, icor);

  /* normalized cut on beam position, if requested */
  if ((sigcut = command_par_value("moncut", cmd->clone)) > 0) {
    idrop = pro_correct_filter(ip, sigcut);
    printf("Disabled %d monitors with %-2.2f sigma cut\n", idrop, sigcut);
    ix = pro_correct_getactive(ip, nm, nx, nc, corvec, monvec, conm);
    icor = ix % 30000;
    imon = ix / 30000;
    printf("After filter of %-2.2f sigma:\n", sigcut);
    printf("%d monitors and %d correctors enabled\n", imon, icor);
  }

  /* set up response matrix for ring or line */
  corl = correct_orbit->cor_table;

  // RING
  if (strcmp("ring", command_par_string("flag", cmd->clone)) == 0) {
    if (dmat != NULL ) myfree(rout_name, dmat);
    /* icor and imon used to set up correct matrix size !! */
    dmat = pro_correct_response_ring(ip, icor, imon);

    // SVD CONDITIONING
    if (command_par_value("cond", cmd->clone) == 1) {
      sngcut = command_par_value("sngcut", cmd->clone);
      sngval = command_par_value("sngval", cmd->clone);
      printf("SVD conditioning requested ...\n");
      if (debug) printf("Conditioning parameters: %e %e\n", sngcut, sngval);

      sflag = c_svddec(dmat, imon, icor, sing, &sngcut, &sngval);
      printf("Initially found %d singular values\n", sflag);

      for (ix = 0; ix < sflag; ix++) {
	corl[nx[sing[2 * ix + 0]]].enable = 0;
	if (debug)
	  printf("Removed:   %d %s\n", nx[sing[2 * ix + 0]], corl[nx[sing[2 * ix + 0]]].p_node->name);
      }

      /* find and prepare enabled correctors and monitors, may be repeated */
      ix = pro_correct_getactive(ip, nm, nx, nc, corvec, monvec, conm);
      icor = ix % 30000;
      imon = ix / 30000;
      printf("After SVD conditioning:             \n");
      printf("%d monitors and %d correctors enabled\n\n", imon, icor);

      if (dmat != NULL ) myfree(rout_name, dmat);

      /* icor and imon used to set up correct matrix size !! */
      dmat = pro_correct_response_ring(ip, icor, imon);
      sflag = c_svddec(dmat, imon, icor, sing, &sngcut, &sngval);
      printf("Finally found %d singular values\n", sflag);
    } //end SVD Conditioning

  } // END RING

  // LINE
  else if (strcmp("line", command_par_string("flag", cmd->clone)) == 0) {
    if (dmat != NULL ) myfree(rout_name, dmat);
    printf("make response for line\n");
    dmat = pro_correct_response_line(ip, icor, imon);

    // SVD CONDITIONING
    if (command_par_value("cond", cmd->clone) == 1) {
      sngcut = command_par_value("sngcut", cmd->clone);
      sngval = command_par_value("sngval", cmd->clone);
      printf("SVD conditioning requested ...\n");
      if (debug)
	printf("Conditioning parameters: %e %e\n", sngcut, sngval);

      sflag = c_svddec(dmat, imon, icor, sing, &sngcut, &sngval);
      printf("Initially found %d singular values\n", sflag);

      for (ix = 0; ix < sflag; ix++) {
	corl[nx[sing[2 * ix + 0]]].enable = 0;
	if (debug)
	  printf("Removed:   %d %s\n", nx[sing[2 * ix + 0]], corl[nx[sing[2 * ix + 0]]].p_node->name);
      }

      ix = pro_correct_getactive(ip, nm, nx, nc, corvec, monvec, conm);
      icor = ix % 30000;
      imon = ix / 30000;
      printf("After SVD conditioning:             \n");
      printf("%d monitors and %d correctors enabled\n\n", imon, icor);

      if (dmat != NULL ) myfree(rout_name, dmat);

      /* icor and imon used to set up correct matrix size !! */
      dmat = pro_correct_response_ring(ip, icor, imon);
      sflag = c_svddec(dmat, imon, icor, sing, &sngcut, &sngval);
      printf("Finally found %d singular values\n", sflag);
    } //end SVD Conditioning
  } //END LINE

  else { // neither ring nor line
    fatal_error("Invalid machine type in CORRECT command", ", MAD-X terminates ");
  }

  if (debug) {
    pro_correct_prtwiss();
    pro_correct_write_cocu_table();
  }

  corrl = command_par_value("corrlim", cmd->clone);
  set_variable("corrlim", &corrl);

  /* Switch block between LSQ, SVD and MICADO correction methods... */

  /* LSQ correction, use all available correctors */
  if (strcmp("lsq", command_par_string("mode", cmd->clone)) == 0) {
    c_haveit(dmat, monvec, corvec, resvec, nx, imon, icor);
    pro_correct_write_results(monvec, resvec, corvec, nx, nc, nm, imon, icor, ip, resout);
  }

  /* SVD correction, use all available correctors */
  else if (strcmp("svd", command_par_string("mode", cmd->clone)) == 0) {
    sflag = c_svdcorr(dmat, monvec, corvec, resvec, nx, imon, icor);
    pro_correct_write_results(monvec, resvec, corvec, nx, nc, nm, imon, icor, ip, resout);
  }

  /* MICADO correction, get desired number of correctors from command */
  else if (strcmp("micado", command_par_string("mode", cmd->clone)) == 0) {
    printf("enter MICADO correction ...\n");

    niter = command_par_value("ncorr", cmd->clone);
    if (niter == 0) {
      printf("Requested %d correctors (\?\?\?) set to %d\n", niter, icor);
      niter = icor;
    } else if (niter < 0) {
      printf("Requested %d correctors (\?\?\?) set to 0\n", niter);
      niter = 0;
    } else if (niter > icor) {
      printf("Fewer correctors available than requested by ncorr\n");
      printf("you want %d,  you get %d\n", niter, icor);
      printf("ncorr reset to %d\n", icor);
      niter = icor;
    }

    /* specifies the maximum RMS value of the orbit to be reached by the correction algorithm.*/
    rms = 1000.0 * command_par_value("error", cmd->clone);

    ifail = c_micit(dmat, conm, monvec, corvec, resvec, nx, rms, imon, icor, niter);

    if (ifail == 0)
      pro_correct_write_results(monvec, resvec, corvec, nx, nc, nm, imon, icor, ip, resout);
    else {
      printf("MICADO correction completed with error code %d\n\n", ifail);
      warning("MICADO back with error", ", no correction done");
    }
  }

  else { // neither LSQ, nor SVD, nor MICADO correction type
    fatal_error("Invalid correction mode in CORRECT command", ", MAD-X terminates ");
  }

  /* write corrector output to tfs table */
  if ((clist = command_par_string("clist", cmd->clone)) != NULL )
    out_table("corr", corr_table, clist);

  /* write monitor output to tfs table */
  if ((mlist = command_par_string("mlist", cmd->clone)) != NULL )
    out_table("mon", mon_table, mlist);

  /* Clean up at the end of the module */
  myfree(rout_name, nm);
  myfree(rout_name, dmat);
  myfree(rout_name, nx);
  myfree(rout_name, nc);
  myfree(rout_name, corvec);
  myfree(rout_name, monvec);
  myfree(rout_name, resvec);
  myfree(rout_name, conm);
  return;
}

static void correct_correct(struct in_cmd* cmd)
/* Steering routine for orbit corrections */
{
  /*
    const char *rout_name = "correct_correct";
  */
  char *orbtab1, *orbtab2;

  if (command_par_value("tworing", cmd->clone)) {
    printf("Want to correct orbit for two rings\n");

    // the following tests only whether a parameter was supplied;
    // the validity of the parameter is tested in pro_correct2_gettables
    if ((orbtab1 = command_par_string("beam1tab", cmd->clone)) == NULL )
      fatal_error("Two beam correction requested but no table supplied for beam 1", orbtab1);
    if ((orbtab2 = command_par_string("beam2tab", cmd->clone)) == NULL )
      fatal_error("Two beam correction requested but no table supplied for beam 2", orbtab2);

    printf("Want to use orbits from: %s and : %s\n", orbtab1, orbtab2);

    correct_correct2(cmd);

  } else {
    printf("Want to correct orbit of a single ring\n");

    if ((orbtab1 = command_par_string("beam1tab", cmd->clone)) != NULL )
      warning("Single beam correction requested but beam 1 table supplied; specified table ignored:", orbtab1);

    if ((orbtab2 = command_par_string("beam2tab", cmd->clone)) != NULL )
      warning("Single beam correction requested but beam 2 table supplied; specified table ignored:", orbtab2);

    correct_correct1(cmd);
  }
}

#if 0 // not used...
static void pro_correct_option(struct in_cmd* cmd)
{
  int i;
  int val, seed;

  int debug = get_option("debug");

  if (debug) {
    fprintf(prt_file, "in coption routine\n");
    for(i=0;i<cmd->tok_list->curr;i++) {
      fprintf(prt_file, "command(s): %s\n",cmd->tok_list->p[i]);
    }
  }

  if (par_present("seed", cmd->clone))
    {
    seed = command_par_value("seed", cmd->clone);
    init55(seed);
    }

  val = command_par_value("print", cmd->clone);

  if(val == 0) {
    if (debug) fprintf(prt_file, "print option not set\n");
    print_correct_opt = 0;
  } else {
    if (debug) fprintf(prt_file, "print option set\n");
    print_correct_opt = val;
  }
}
#endif

static int pro_correct_getcommands(struct in_cmd* cmd) {

  int iplane = 1;
  char plane[20];

  if (get_option("debug"))
    printf("enter CORRECT module\n");

  if (current_sequ == NULL || current_sequ->ex_start == NULL ) {
    warning("CORRECT, but no active sequence:", "ignored");
    return (-1);
  }

  strcpy(plane, command_par_string("plane", cmd->clone));
  if      (strcmp("x", plane) == 0) iplane = 1;
  else if (strcmp("y", plane) == 0) iplane = 2;
  else if (strcmp("h", plane) == 0) iplane = 1;
  else if (strcmp("v", plane) == 0) iplane = 2;
  else {
    printf("No valid plane specified, x plane used \n");
    iplane = 1;
  }

  return iplane;
}

static int pro_correct_gettables(int iplane, struct in_cmd* cmd) {

  const char *rout_name = "pro_correct_gettables";

  struct id_mic *cor_l;
  struct id_mic *mon_l;

  struct table *ttb;

  char* orbtab;
  char* tartab;
  char* modtab;

  int j;

  int cntm = 0;
  int cntc = 0;

  double ounits;

  /*
    static char atm[6][4] = {"hmon","vmon","moni","hkic","vkic","kick"};
  */

  int corzero;

  int debug = get_option("debug");

  /* Get access to tables, for orbit and model the default is twiss_table */

  if ((orbtab = command_par_string("orbit", cmd->clone)) != NULL ) {
    printf("Want to use orbit from: %s\n", orbtab);
    if (!(orbin_table = find_table(orbtab))) {
      fatal_error("ORBIT table for correction requested, but not provided:", orbtab);
    }
  } else { // the orbit table is the twiss table
    if ((orbin_table = twiss_table) == NULL ) {
      fatal_error("ORBIT cannot be obtained from non-existing TWISS table",
		  "You MUST run TWISS before trying to correct the orbit");
    } else {
      if (debug)
	printf("orbit from TWISS table at address: %p\n", (void*) twiss_table);
    }
  }

  if ((tartab = command_par_string("target", cmd->clone)) != NULL ) {
    printf("Want to use target orbit from: %s\n", tartab);
    if (!(target_table = find_table(tartab))) {
      fatal_error("TARGET table for correction requested, but not provided:", tartab);
    }
  } else {
    if (debug)
      printf("No target orbit requested\n");
  }

  if ((modtab = command_par_string("model", cmd->clone)) != NULL ) {
    printf("Want to use model orbit from: %s\n", modtab);
    if (!(model_table = find_table(modtab))) {
      fatal_error("MODEL table for correction requested, but not provided:", modtab);
    }
  } else {
    if ((model_table = twiss_table) == NULL ) {
      fatal_error("MODEL cannot be obtained from non-existing TWISS table",
		  "You MUST run TWISS before trying to correct the orbit");
    } else {
      if (debug)
	printf("model from TWISS table at address: %p\n", (void*) twiss_table);
    }
  }

  if (debug)
    printf( "The orbit, twiss, target and model tables are at addresses: %p %p %p %p\n",
	    (void*) orbin_table, (void*) twiss_table, (void*) target_table, (void*) model_table);

  if (correct_orbit == NULL )
    correct_orbit = mycalloc("pro_correct_gettables", 1, sizeof *correct_orbit);

  if (debug) printf("-0-\n");
  // if(corr_table == NULL) {
  corr_table = make_table("corr", "corr", corr_table_cols, corr_table_types, 15000);
  add_to_table_list(corr_table, table_register);
  pro_correct_make_corr_table();
  // }

  if (debug) printf("-1-\n");
  // if(mon_table == NULL) {
  mon_table = make_table("mon", "mon", mon_table_cols, mon_table_types, 15000);
  add_to_table_list(mon_table, table_register);
  pro_correct_make_mon_table();
  // }

  if (debug) printf("-2-\n");

  if (correct_orbit->cor_table != NULL ) myfree(rout_name, correct_orbit->cor_table);
  if (correct_orbit->mon_table != NULL ) myfree(rout_name, correct_orbit->mon_table);
  correct_orbit->cor_table = mycalloc("pro_correct_gettables_cor", 30200, sizeof *correct_orbit->cor_table);
  correct_orbit->mon_table = mycalloc("pro_correct_gettables_mon", 30200, sizeof *correct_orbit->mon_table);

  /* orbit table available, get units, if defined */
  //  2013-Jun-24  12:06:03  ghislain: FIXME - units option of CORRECT command is only partially documented!!!
  if ((ounits = command_par_value("units", cmd->clone)) > 0)
    correct_orbit->units = ounits;
  else correct_orbit->units = 1.0;

  ttb = model_table;
  correct_orbit->mon_table->previous = NULL;
  correct_orbit->mon_table->next = NULL;
  correct_orbit->cor_table->previous = NULL;
  correct_orbit->cor_table->next = NULL;

  if (debug) printf("-3-\n");

  mon_l = correct_orbit->mon_table;
  cor_l = correct_orbit->cor_table;

  corzero = command_par_value("corzero", cmd->clone);

  // go through the model table and build chained lists of monitors and correctors
  for (j = 0; j < ttb->curr; j++) {
    if (!ttb->p_nodes[j]->base_name) continue;
    if ((strncmp(atm[iplane - 1], ttb->p_nodes[j]->base_name, 4) == 0)
	|| (strncmp(atm[2], ttb->p_nodes[j]->base_name, 4) == 0)) {
      mon_l->id_ttb = j;
      mon_l->enable = ttb->p_nodes[j]->enable;
      mon_l->p_node = ttb->p_nodes[j];
      mon_l->next = mon_l;
      mon_l->next++;
      mon_l++;
      cntm++;
    }
    if ((strncmp(atc[iplane - 1], ttb->p_nodes[j]->base_name, 4) == 0)
	|| (strncmp(atc[2], ttb->p_nodes[j]->base_name, 4) == 0)) {
      cor_l->id_ttb = j;
      cor_l->enable = ttb->p_nodes[j]->enable;
      cor_l->p_node = ttb->p_nodes[j];

      if (corzero > 0) {
        if (iplane == 1) cor_l->p_node->chkick = 0.0;
        if (iplane == 2) cor_l->p_node->cvkick = 0.0;
      }

      cor_l->next = cor_l;
      cor_l->next++;
      cor_l++;
      cntc++;
    }
  }

  if (debug) printf("-4-\n");

  mon_l--;
  mon_l->next = NULL;
  cor_l--;
  cor_l->next = NULL;

  if (debug) printf("done: %d %d\n", cntm, cntc);

  // 2013-Jun-24  12:21:49  ghislain:
  // following is a kludge to return a single value but has to be decoded on other side.
  if(cntc >= 30000)
    fatal_error("Found more than 30000 correctors; decoding in mad_orbit.c will fail",
		"Please report this issue to MAD developpers (mad@cern.ch)");
  return (30000 * cntm + cntc);
}

static int pro_correct_getorbit(struct in_cmd* cmd) {
  struct id_mic *m; /* access to tables for monitors and correctors */
  struct table *ttb;
  struct table *tar = NULL;
  double **da1;
  double **da2 = NULL;
  double xlimit;
  double rx, ry, dpsi;

  char strx[40];
  char stry[40];

  int posx, posy, pospx, pospy;
  int tosx = -1;
  int tosy = -1;
  int tospx, tospy;

  int debug = get_option("debug");

  ttb = orbin_table;
  da1 = ttb->d_cols;

  if (target_table != NULL ) {
    tar = target_table;
    da2 = tar->d_cols;
  }

  m = correct_orbit->mon_table;

  strcpy(strx, "x");
  strcpy(stry, "y");

  if ((posx = name_list_pos(strx, ttb->columns)) < 0)
    fatal_error("orbit x not found in input table", ", MAD-X terminates ");

  if ((posy = name_list_pos(stry, ttb->columns)) < 0)
    fatal_error("orbit y not found in input table", ", MAD-X terminates ");

  if (debug) {
    if ((pospx = name_list_pos("px", ttb->columns)) < 0)
      fatal_error("orbit px not found in input table", ", MAD-X terminates ");

    if ((pospy = name_list_pos("py", ttb->columns)) < 0)
      fatal_error("orbit py not found in input table", ", MAD-X terminates ");

    printf("====c1===>  %d %d %d %d \n", posx, posy, pospx, pospy);
  }

  if (command_par_string("target", cmd->clone) != NULL ) {
    if ((tosx = name_list_pos("x", tar->columns)) < 0)
      fatal_error("target orbit x not found in table", ", MAD-X terminates ");

    if ((tosy = name_list_pos("y", tar->columns)) < 0)
      fatal_error("target orbit y not found in table", ", MAD-X terminates ");

    if (debug) {
      if ((tospx = name_list_pos("px", tar->columns)) < 0)
  fatal_error("target orbit px not found in table", ", MAD-X terminates ");

      if ((tospy = name_list_pos("py", tar->columns)) < 0)
  fatal_error("target orbit px not found in table", ", MAD-X terminates ");

      printf("====c1===>  %d %d %d %d \n", tosx, tosy, tospx, tospy);
    }
  }

  while (m) {

    if (command_par_string("target", cmd->clone) != NULL ) {
      /* If correction to target orbit, subtract the target orbit to correct only difference... */
      m->val.before[0] = (da1[posx][m->id_ttb] - da2[tosx][m->id_ttb]) * 1000. * correct_orbit->units;
      m->val.before[1] = (da1[posy][m->id_ttb] - da2[tosy][m->id_ttb]) * 1000. * correct_orbit->units;
    } else {
      m->val.before[0] = da1[posx][m->id_ttb] * 1000. * correct_orbit->units;
      m->val.before[1] = da1[posy][m->id_ttb] * 1000. * correct_orbit->units;
    }

    /* monon=xlimit determines the fraction of available monitors */
    if (par_present("monon", cmd->clone)) {
      xlimit = command_par_value("monon", cmd->clone);
      if (frndm() > xlimit) {
	m->enable = 0;
	printf("Monitor %s disabled\n", m->p_node->name);
      }
    }

    /* scaling error should come first, monitor alignment not scaled ... */
    if (par_present("monscale", cmd->clone)) {
      if ((command_par_value("monscale", cmd->clone)) == 1) {
	if (m->p_node->p_al_err != NULL ) {
	  if (debug) {
	    printf("m-list: %d %s %s\n", m->id_ttb, m->p_node->name, m->p_node->base_name);
	    printf("scales: %e %e\n", m->p_node->p_al_err->a[12], m->p_node->p_al_err->a[13]);
	  }
	  m->val.before[0] = m->val.before[0] * (1.0 + m->p_node->p_al_err->a[12]);
	  m->val.before[1] = m->val.before[1] * (1.0 + m->p_node->p_al_err->a[13]);
	}
      }
    }

    /* monitor misalignment after all other reading manipulations ! */
    if (par_present("monerror", cmd->clone)) {
      if ((command_par_value("monerror", cmd->clone)) == 1) {
	if (m->p_node->p_al_err != NULL ) {
	  if (debug) {
	    printf("m-list: %d %s %s\n", m->id_ttb, m->p_node->name, m->p_node->base_name);
	    printf("errors: %e %e \n", m->p_node->p_al_err->a[6], m->p_node->p_al_err->a[7]);
	  }
	  dpsi = m->p_node->p_al_err->a[5];
	  rx = m->val.before[0];
	  ry = m->val.before[1];
	  if (debug) printf("\nA: %e %e %e\n", m->val.before[0], m->val.before[1], dpsi);
	  m->val.before[0] =  rx * cos(dpsi) + ry * sin(dpsi);
	  m->val.before[1] = -rx * sin(dpsi) + ry * cos(dpsi);
	  if (debug) printf("B: %e %e %e\n", m->val.before[0], m->val.before[1], dpsi);
	  m->val.before[0] += m->p_node->p_al_err->a[6] * 1000.;
	  m->val.before[1] += m->p_node->p_al_err->a[7] * 1000.;
	  if (debug) printf("C: %e %e %e\n", m->val.before[0], m->val.before[1], dpsi);
	}
      }
    }

    m = m->next;
  }

  return (0);
}

static int pro_correct_getorbit_ext(struct in_cmd* cmd) {
  int j;

  struct id_mic *m; /* access to tables for monitors and correctors */
  struct table *ttb;
  struct table *tar = NULL;
  double **da1;
  double **da2 = NULL;
  double xlimit;
  char name[NAME_L];
  char l1name[NAME_L];
  char l2name[NAME_L];
  char l3name[NAME_L];
  char l4name[NAME_L];
  double rx, ry, dpsi;

  char *nam_col;
  char *x_col;
  char *y_col;

  char strx[40];
  char stry[40];
  char strn[40];

  int posx, posy, pospx, pospy;
  int tosx = -1;
  int tosy = -1;
  int tospx, tospy;

  int yok;

  int jjx, jjy, jj;

  int debug = get_option("debug");

  ttb = orbin_table;
  da1 = ttb->d_cols;

  if (target_table != NULL ) {
    tar = target_table;
    da2 = tar->d_cols;
  }

  m = correct_orbit->mon_table;

  if ((x_col = command_par_string("x_col", cmd->clone)) != NULL ) {
    printf("X orbit in column: %s\n", x_col);
    strcpy(strx, x_col);
  } else {
    strcpy(strx, "x");
  }

  if ((y_col = command_par_string("y_col", cmd->clone)) != NULL ) {
    printf("Y orbit in column: %s\n", y_col);
    strcpy(stry, y_col);
  } else {
    strcpy(stry, "y");
  }

  if ((nam_col = command_par_string("name_col", cmd->clone)) != NULL ) {
    printf("names in column: %s\n", nam_col);
    // strcpy(strn, "name"); 2015-Apr-22  16:24:34  ghislain: bug fix, but never used
    strcpy(strn, nam_col);
  } else {
    strcpy(strn, "name");
  }

  if ((posx = name_list_pos(strx, ttb->columns)) < 0)
    fatal_error("orbit x not found in input table", ", MAD-X terminates ");

  if ((posy = name_list_pos(stry, ttb->columns)) < 0)
    fatal_error("orbit y not found in input table", ", MAD-X terminates ");

  if (debug) {
    if ((pospx = name_list_pos("px", ttb->columns)) < 0)
      warning("orbit px not found in input table", ", MAD-X continues ");

    if ((pospy = name_list_pos("py", ttb->columns)) < 0)
      warning("orbit py not found in input table", ", MAD-X continues ");

    printf("====c1===>  %d %d %d %d \n", posx, posy, pospx, pospy);
  }

  if (command_par_string("target", cmd->clone) != NULL ) {
    if ((tosx = name_list_pos("x", tar->columns)) < 0)
      fatal_error("target orbit x not found in table", ", MAD-X terminates ");

    if ((tosy = name_list_pos("y", tar->columns)) < 0)
      fatal_error("target orbit y not found in table", ", MAD-X terminates ");

    if (debug) {
      if ((tospx = name_list_pos("px", tar->columns)) < 0)
       warning("target orbit px not found in table", ", MAD-X continues ");

      if ((tospy = name_list_pos("py", tar->columns)) < 0)
       warning("target orbit py not found in table", ", MAD-X continues ");

      printf("====c1===>  %d %d %d %d \n", tosx, tosy, tospx, tospy);
    }
  }

  if (debug) {
    printf("Number in table: %d\n", ttb->curr);
    // 2013-Jun-24  13:42:16  ghislain: ????
    // for (j=1; j < (ttb->curr)+1; j++) {
    //   string_from_table_row(ttb->name, "name", &j, name);
    // }
  }

  jj = 0;

  while (m) { // loop over the monitors
    strcpy(l1name, m->p_node->name);
    stolower(l1name);
    strcpy(l2name, strip(l1name));
    supp_tb(l2name);

    if (debug)
      printf("monitor name: %s\n", l2name);

    jjx = -1;
    jjy = -1;
    jj++;
    yok = 0;

    for (j = 1; j <= ttb->curr; j++) {
      string_from_table_row(ttb->name, "name", &j, name);
      strcpy(l3name, name);
      stolower(l3name);
      strcpy(l4name, strip(l3name));
      supp_tb(l4name);
      if (strlen(l4name) == strlen(l2name)) {
        if (strncmp(l4name, l2name, strlen(l2name)) == 0) {
          jjx = j - 1;
          jjy = jj - 1;
          yok = 1;
          if (debug)
            printf("monitor names found: %s %s %d\n", l2name, l4name, yok);
        }
      }
    }





    if (debug)
      printf("jjx, jjy, yok : %d %d %d\n", jjx, jjy, yok);

    if ((jjy >= 0) && (yok == 1)) { // the monitor was found

      if (command_par_string("target", cmd->clone) != NULL ) { // target exists
	// If correction to target orbit, subtract the target orbit, then correct to zero ...

        if (debug) {
          printf("x ==> %d %d %e %e\n", jjx, m->id_ttb, da1[posx][jjx], da2[tosx][jjy]);
          printf("y ==> %e %e\n",                       da1[posy][jjx], da2[tosy][jjy]);
        }

        m->val.before[0] = (da1[posx][jjx] - da2[tosx][jjy]) * 1000. * correct_orbit->units;
        m->val.before[1] = (da1[posy][jjx] - da2[tosy][jjy]) * 1000. * correct_orbit->units;

        if (debug)
          printf("bxy ==> %s %d %e %e\n", m->p_node->name, jjx, m->val.before[0], m->val.before[1]);

      } else { // no target was given

        if (debug) {
          printf("x ==> %e \n", da1[posx][jjx]);
          printf("y ==> %e \n", da1[posy][jjx]);
        }

        m->val.before[0] = da1[posx][jjx] * 1000. * correct_orbit->units;
        m->val.before[1] = da1[posy][jjx] * 1000. * correct_orbit->units;

        if (debug)
          printf("bxy ==> %s %d %e %e\n", m->p_node->name, jjx, m->val.before[0], m->val.before[1]);
      }

      /* monon=xlimit determines the fraction of available monitors */
      if (par_present("monon", cmd->clone)) {
        xlimit = command_par_value("monon", cmd->clone);
        if (frndm() > xlimit) {
          m->enable = 0;
          printf("Monitor %s disabled\n", m->p_node->name);
        }
      }

      /* scaling error should come first, monitor alignment not scaled ... */
      if (par_present("monscale", cmd->clone)) {
        if ((command_par_value("monscale", cmd->clone)) == 1) {
          if (m->p_node->p_al_err != NULL ) {

            if (debug) {
              printf("m-list: %d %s %s\n", m->id_ttb, m->p_node->name, m->p_node->base_name);
              printf("scales: %e %e\n", m->p_node->p_al_err->a[12], m->p_node->p_al_err->a[13]);
            }

            m->val.before[0] = m->val.before[0] * (1.0 + m->p_node->p_al_err->a[12]);
            m->val.before[1] = m->val.before[1] * (1.0 + m->p_node->p_al_err->a[13]);
          }
        }
      }

      /* monitor misalignment after all other reading manipulations ! */
      if (par_present("monerror", cmd->clone)) {
        if ((command_par_value("monerror", cmd->clone)) == 1) {
          if (m->p_node->p_al_err != NULL ) {

            if (debug) {
              printf("m-list: %d %s %s\n", m->id_ttb, m->p_node->name, m->p_node->base_name);
              printf("errors: %e %e \n", m->p_node->p_al_err->a[6], m->p_node->p_al_err->a[7]);
            }

            dpsi = m->p_node->p_al_err->a[5];
            rx = m->val.before[0];
            ry = m->val.before[1];

            if (debug)
              printf("\nA: %e %e %e\n", m->val.before[0], m->val.before[1], dpsi);

            m->val.before[0] = rx * cos(dpsi) + ry * sin(dpsi);
            m->val.before[1] = -rx * sin(dpsi) + ry * cos(dpsi);

            if (debug)
              printf("B: %e %e %e\n", m->val.before[0], m->val.before[1], dpsi);

            m->val.before[0] += m->p_node->p_al_err->a[6] * 1000.;
            m->val.before[1] += m->p_node->p_al_err->a[7] * 1000.;

            if (debug)
              printf("C: %e %e %e\n", m->val.before[0], m->val.before[1], dpsi);
          }
        }
      }
    } // end of treatment of the monitor that was found to exist
    else {
      m->enable = 0; /* Only enable monitors found in input */
    }

    m = m->next;
  }

  return 0;
}

static int pro_correct_getcorrs(struct in_cmd* cmd) {
  struct id_mic *c; /* access to tables for monitors and correctors */

  (void) cmd;

  int debug = get_option("debug");

  c = correct_orbit->cor_table;
  while (c) {
    c->val.before[0] = c->p_node->chkick * 1000.;
    c->val.before[1] = c->p_node->cvkick * 1000.;

    if (debug) {
      printf("c-list: %d %s %s\n", c->id_ttb, c->p_node->name, c->p_node->base_name);
      printf("initial strengths: %e %e\n", c->val.before[0], c->val.before[1]);
    }

    c = c->next;
  }

  return (0);
}

static void pro_correct_prtwiss(void) {
  int i, j;
  int pr_cols;
  struct table *ttb;
  double **da1;

  ttb = model_table;

  printf(" %d %d\n", ttb->curr, ttb->num_cols);
  for (i = 0; i < ttb->curr; i++) {
    printf(" %s %s\n", ttb->s_cols[0][i], ttb->s_cols[1][i]);
  }

  da1 = ttb->d_cols;
  for (j = 0; j < ttb->curr; j++) {
    printf("\n\n");
    printf("from table: %s \n", ttb->node_nm->p[j]);
    printf("from node:  %s \n", ttb->p_nodes[j]->name);
    printf(" %s %s\n", ttb->s_cols[0][j], ttb->s_cols[1][j]);

    pr_cols = ttb->num_cols;
    pr_cols = 19; /* print only for 20 columns */
    for (i = 0; i < pr_cols; i++) {
      if (&da1[i][0] != NULL ) {
	printf("%-8s %f\n", twiss_table_cols[i], da1[i][j]);
      }
    }
  }
  return;
}

static void pro_correct_write_cocu_table(void) {
  int i, j;
  int pr_cols;
  int cp[13] = { 1, 0, 2, 9, 11, 3, 6, 4, 7, 5, 8, 15, 17 };
  struct table *ttb;
  double **da1;
  FILE *fp1;

  fp1 = fopen("cocu_in.opt", "w");
  ttb = model_table;

  pr_cols = ttb->num_cols;
  pr_cols = 13; /* print only for 19 columns */
  fprintf(fp1, "*");
  for (i = 0; i < pr_cols; i++) {
    fprintf(fp1, "%-8s ", twiss_table_cols[cp[i]]);
  }

  da1 = ttb->d_cols;
  for (j = 0; j < ttb->curr; j++) {
    fprintf(fp1, "\n%s %s ", ttb->s_cols[1][j], ttb->s_cols[0][j]);
    for (i = 2; i < pr_cols; i++) {
      if (&da1[cp[i]][0] != NULL ) {
	fprintf(fp1, " %f", da1[cp[i]][j]);
      }
    }
  }
  return;
}

static int pro_correct_filter(int iplane, double sigcut) {
  int im, ip, icnt; // ic, no used
  struct id_mic *m; /* access to tables for monitors */

  struct table *ttb;
  static char pl[2] = "xy";
  double **da1;
  double bx_m = -9999.;
  double xsig;
  double xmea;
  double xn;

  int debug = get_option("debug");

  ttb = model_table;
  da1 = ttb->d_cols;
  im = 0;
  icnt = 0;
  ip = iplane - 1;

  printf("A (normalized) cut of %-2.2f is requested\n", sigcut);

  m = correct_orbit->mon_table;
  xmea = 0.0;
  while (m) {
    if (debug)
      printf("monitor flag: %d\n", m->enable);

    if (m->enable == 1) {
      if (ip == 0)
	bx_m = da1[3][m->id_ttb];
      else if (ip == 1)
	bx_m = da1[6][m->id_ttb];
      else {
	// 2013-Dec-06  17:01:03  ghislain: failure case
      }

      xn = m->val.before[ip] / sqrt(bx_m);
      xmea += xn;
      if (debug) {
	printf("==> %s %-4.3f %-4.3f \n", m->p_node->name, bx_m, m->val.before[ip]);
	printf("==> %-4.3f \n", xn);
      }
      im++;
    }
    m = m->next;
  }

  xmea = xmea / im;
  if (debug)
    printf("Mean values: %-4.3f \n", xmea);

  m = correct_orbit->mon_table;
  im = 0;
  xsig = 0.0;
  while (m) {
    if (m->enable == 1) {
      if (ip == 0)
	bx_m = da1[3][m->id_ttb];
      else if (ip == 1)
	bx_m = da1[6][m->id_ttb];
      else {
	// 2013-Dec-06  17:02:41  ghislain: failure case
      }
      xn = m->val.before[ip] / sqrt(bx_m);
      xsig += (xmea - xn) * (xmea - xn);
      im++;
    }
    m = m->next;
  }

  xsig = sqrt(xsig / im);
  if (debug)
    printf("Sigma values: %-4.3f \n", xsig);

  m = correct_orbit->mon_table;
  while (m) {
    if (m->enable == 1) {
      if (ip == 0)
	bx_m = da1[3][m->id_ttb];
      else if (ip == 1)
	bx_m = da1[6][m->id_ttb];
      else {
	// 2013-Dec-06  17:03:47  ghislain: failure case
      }
      xn = (m->val.before[ip] / sqrt(bx_m)) - xmea;
      if (fabs(xn) > (sigcut * xsig)) {
	printf( "disabled %s %c = %-4.3f (%-4.3f), limit is %-2.2f*%-4.3f\n",
		m->p_node->name, pl[ip], xn, m->val.before[ip], sigcut, xsig);
	m->enable = 0;
	icnt++;
      }
    }
    m = m->next;
  }

  return (icnt);
}

static double*
pro_correct_response_ring(int ip, int nc, int nm) {
  int ic, im;
  struct id_mic *m, *c; /* access to tables for monitors and correctors */

  struct table *ttb;
  double **da1;
  double bx_c, by_c, pix_c, piy_c;
  double bx_m, by_m, pix_m, piy_m;
  double qx0, qy0;
  double respx, respy;
  double *dmat;

  int debug = get_option("debug");

  ttb = model_table;
  da1 = ttb->d_cols;
  ic = 0;
  im = 0;

  dmat = mycalloc_atomic("pro_correct_response_ring", nc*nm, sizeof *dmat);

  correct_orbit->qx0 = da1[5][ttb->curr - 1];
  correct_orbit->qy0 = da1[8][ttb->curr - 1];
  qx0 = correct_orbit->qx0;
  qy0 = correct_orbit->qy0;

  c = correct_orbit->cor_table;
  ic = 0;

  while (c) {
    if (debug)
      printf("corrector flag: %d\n", c->enable);

    if (c->enable == 1) {
      bx_c = da1[3][c->id_ttb];
      by_c = da1[6][c->id_ttb];
      pix_c = da1[5][c->id_ttb];
      piy_c = da1[8][c->id_ttb];
      m = correct_orbit->mon_table;
      im = 0;
      while (m) {
        if (debug)
        printf("monitor flag: %d\n", m->enable);
	if (m->enable == 1) {
	  bx_m = da1[3][m->id_ttb];
	  by_m = da1[6][m->id_ttb];
	  pix_m = da1[5][m->id_ttb];
	  piy_m = da1[8][m->id_ttb];
	  respx = 0.0;
	  respy = 0.0;
	  if (ip == 1) {
	    respx = sqrt(bx_m * bx_c) / (2.0 * sin(pi * qx0)) * cos((fabs(pix_m - pix_c) * twopi) - qx0 * pi);
	    setup_(&respx, dmat, &im, &ic, &nm, &nc);
	  } else if (ip == 2) {
	    respy = sqrt(by_m * by_c) / (2.0 * sin(pi * qy0)) * cos((fabs(piy_m - piy_c) * twopi) - qy0 * pi);
	    setup_(&respy, dmat, &im, &ic, &nm, &nc);
	  }
	  im++;
	}
	m = m->next;
      };
      ic++;
    }
    c = c->next;
  }

  return dmat;
}

static double*
pro_correct_response_line(int ip, int nc, int nm) {
  int ic, im;
  struct id_mic *m, *c; /* access to tables for monitors and correctors */

  struct table *ttb;
  double **da1;
  double bx_c, by_c, pix_c, piy_c;
  double bx_m, by_m, pix_m, piy_m;
  double respx, respy;
  double *dmat;

  ttb = model_table;
  da1 = ttb->d_cols;
  ic = 0;
  im = 0;

  dmat = mycalloc_atomic("pro_correct_response_ring", nc*nm, sizeof *dmat);

  correct_orbit->qx0 = da1[5][ttb->curr - 1];
  correct_orbit->qy0 = da1[8][ttb->curr - 1];

  c = correct_orbit->cor_table;
  ic = 0;
  while (c) {
    if (c->enable == 1) {
      bx_c = da1[3][c->id_ttb];
      by_c = da1[6][c->id_ttb];
      pix_c = da1[5][c->id_ttb];
      piy_c = da1[8][c->id_ttb];
      m = correct_orbit->mon_table;
      im = 0;
      while (m) {
	if (m->enable == 1) {
	  bx_m = da1[3][m->id_ttb];
	  by_m = da1[6][m->id_ttb];
	  pix_m = da1[5][m->id_ttb];
	  piy_m = da1[8][m->id_ttb];
	  respx = 0.0;
	  respy = 0.0;
	  if (ip == 1) {
	    if (pix_m > pix_c) {
	      respx = sqrt(bx_m * bx_c) * sin((pix_m - pix_c) * twopi);
	    } else {
	      respx = 0.0;
	    }
	    setup_(&respx, dmat, &im, &ic, &nm, &nc);
	  } else if (ip == 2) {
	    if (piy_m > piy_c) {
	      respy = sqrt(by_m * by_c) * sin((piy_m - piy_c) * twopi);
	    } else {
	      respy = 0.0;
	    }
	    setup_(&respy, dmat, &im, &ic, &nm, &nc);
	  }
	  im++;
	}
	m = m->next;
      }
      ic++;
    }
    c = c->next;
  }

  return dmat;
}

static void pro_correct_make_corr_table(void) {
  struct table *ttb;
  int j;
  /*
    static char atm[5][4] = {"hmon","vmon","hkic","vkic","kick"};
  */
  /*
    ttb = orbin_table;
  */

  ttb = model_table;

  if(model_table == 0x0)
   {
     fatal_error("pro_correct_make_corr_table "," Model table does not exist");
   }
  /*   else
       {
       printf("pro_correct_make_corr_table model table OK \n");
       printf("pro_correct_make_corr_table model table name %s\n",model_table->name);
       printf("pro_correct_make_corr_table model table type %s\n",model_table->type);
       printf("pro_correct_make_corr_table model table curr %d\n",model_table->curr);

       }
  */
  for (j = 0; j < ttb->curr; j++)
    {
      /*
	 printf("pro_correct_make_corr_table node %d \n",j);
	 printf("                         node_addr %#x\n",ttb->p_nodes[j]);
	 printf("                         base_name %#x\n",ttb->p_nodes[j]->base_name);

      */
      if (!ttb->p_nodes[j]->base_name)
	{
	  continue;
	}

      if ( (strncmp(atc[0], ttb->p_nodes[j]->base_name, 4) == 0)
	   || (strncmp(atc[1], ttb->p_nodes[j]->base_name, 4) == 0)
	   || (strncmp(atc[2], ttb->p_nodes[j]->base_name, 4) == 0)) {
	string_to_table_curr("corr", "name", ttb->p_nodes[j]->name);
	augment_count("corr");
      }
    }
}

static void pro_correct2_make_corr_table(void) {
  struct id_mic2 *ttb;
  /*
    static char atm[5][4] = {"hmon","vmon","hkic","vkic","kick"};
  */

  ttb = correct_orbit12->cor_table;

  while (ttb != NULL ) {
    if ((strncmp(atc[0], ttb->p_node->base_name, 4) == 0)
	|| (strncmp(atc[1], ttb->p_node->base_name, 4) == 0)
	|| (strncmp(atc[2], ttb->p_node->base_name, 4) == 0)) {
      if (ttb->id_ttb[0] > 0) {
        string_to_table_curr("corr1", "name", ttb->p_node->name);
        augment_count("corr1");
      }

      if (ttb->id_ttb[1] > 0) {
        string_to_table_curr("corr2", "name", ttb->p_node->name);
        augment_count("corr2");
      }
    }
    ttb = ttb->next;
  }
}

static void pro_correct_make_mon_table(void) {
  struct table *ttb;
  int j;
  /*
    static char atm[3][4] = {"hmon","vmon","moni"};
  */

  ttb = model_table;

  for (j = 0; j < ttb->curr; j++) {
    if (!ttb->p_nodes[j]->base_name) continue;
    if ((strncmp(atm[0], ttb->p_nodes[j]->base_name, 4) == 0)
     || (strncmp(atm[1], ttb->p_nodes[j]->base_name, 4) == 0)
     || (strncmp(atm[2], ttb->p_nodes[j]->base_name, 4) == 0)) {
      string_to_table_curr("mon", "name", ttb->p_nodes[j]->name);
      augment_count("mon");
    }
  }
}

static void pro_correct2_make_mon_table(void) {
  struct id_mic2 *ttb;
  /*
    static char atm[3][4] = {"hmon","vmon","moni"};
  */

  ttb = correct_orbit12->mon_table;

  while (ttb != NULL ) {
    if ((strncmp(atm[0], ttb->p_node->base_name, 4) == 0)
     || (strncmp(atm[1], ttb->p_node->base_name, 4) == 0)
     || (strncmp(atm[2], ttb->p_node->base_name, 4) == 0)) {
      string_to_table_curr("mon", "name", ttb->p_node->name);
      augment_count("mon");
    }
    ttb = ttb->next;
  }
}

static void pro_correct_fill_corr_table(int ip, char *name, double old, double new) {
  struct table *cor;

  int j;

  cor = corr_table;

  for (j = 0; j < cor->curr; j++) {
    if (strcmp(name, cor->s_cols[0][j]) == 0) {
      cor->d_cols[ip][j] = old;
      cor->d_cols[ip + 2][j] = new;
    }
  }
}

static void pro_correct2_fill_corr_table(int b, int ip, char *name, double old, double new) {
  struct table *cor = NULL;

  int j;

  if ((b != 1) && (b != 0)) {
    char buf[64];
    sprintf(buf, "%d", b);
    fatal_error("Invalid beam requested:", buf);
  }

  if (b == 0) cor = corr_table1;
  if (b == 1) cor = corr_table2;

  for (j = 0; j < cor->curr; j++) {
    if (strcmp(name, cor->s_cols[0][j]) == 0) {
      cor->d_cols[ip][j] = old;
      cor->d_cols[ip + 2][j] = new;
    }
  }
}

static void pro_correct_fill_mon_table(int ip, char *name, double old, double new) {
  struct table *mon;

  int j;

  mon = mon_table;

  for (j = 0; j < mon->curr; j++) {
    if (strcmp(name, mon->s_cols[0][j]) == 0) {
      mon->d_cols[ip][j] = old * 0.001;
      mon->d_cols[ip + 2][j] = new * 0.001;
    }
  }
}

static void pro_correct2_fill_mon_table(int ip, char *name, double old, double new) {
  struct table *mon;

  int j;

  mon = mon_table;

  for (j = 0; j < mon->curr; j++) {
    if (strcmp(name, mon->s_cols[0][j]) == 0) {
      mon->d_cols[ip][j] = old * 0.001;
      mon->d_cols[ip + 2][j] = new * 0.001;
    }
  }
}

static void pro_correct_write_results(double *monvec, double *resvec, double *corvec,
				      int *nx, int *nc, int *nm, int imon, int icor, int ip,
				      int resout)
/*                                              */
/* Writes a summary of the correction           */
/* Writes correctors strengths into sequences   */
/* Fills TFS tables for correctors and monitors */
/* Fills the 'stren.out' output                 */
/* Makes various prints on request              */
/*                                              */
{
  int i;
  int rst;
  double corrm;
  struct id_mic *m, *c; /* access to tables for monitors and correctors */

  m = correct_orbit->mon_table;
  c = correct_orbit->cor_table;

  if (fddata != NULL ) {
    rst = get_variable("n");
    fprintf(fddata, "%d %d %e %e %e %e %e %e\n", ip, rst,
      cptp(monvec, imon), cptp(resvec, imon), crms(monvec, imon),
      crms(resvec, imon), copk(monvec, imon), copk(resvec, imon));
  }

  // 2014-Jan-07  17:33:59  ghislain: change of output format to print RMS and STDDEV values
  if (print_correct_opt > 0) {
    printf("CORRECTION SUMMARY:   \n\n");
    printf("                   average [mm]   std.dev. [mm]      RMS [mm]        peak-to-peak [mm]\n\n");
    printf("before correction: %f        %f          %f        %f \n",
     caverage(monvec,imon), cstddev(monvec,imon), crms(monvec, imon), cptp(monvec, imon));
    printf("after correction:  %f        %f          %f        %f \n\n\n",
     caverage(resvec,imon), cstddev(resvec,imon), crms(resvec, imon), cptp(resvec, imon));
  }

  if (print_correct_opt > 1) {
    printf("Monitor:  Before:     After:    Difference:\n");
    printf("           (mm)        (mm)         (mm)   \n");
  }

  for (i = 0; i < imon; i++) {
    if (print_correct_opt > 1) {
      printf("%s   %-4.3f     %-4.3f     %-4.3f\n", m[nm[i]].p_node->name,
       monvec[i], resvec[i], resvec[i] - monvec[i]);
    }
    m[nm[i]].val.after[ip - 1] = resvec[i];
    pro_correct_fill_mon_table(ip, m[nm[i]].p_node->name, monvec[i], resvec[i]);
  }

  corrm = copk(corvec, icor);

  if (corrm > corrl) {
    printf("Max strength: %e should be less than corrector strength limit: %e\n", corrm, corrl);
    warning("maximum corrector strength larger than limit","");
  } else {
    printf("Max strength: %e is below corrector strength limit: %e\n", corrm, corrl);
  }

  set_variable("corrmax", &corrm);
  if (print_correct_opt > 1) {
    printf("Max strength: %e\n", copk(corvec, icor));
    printf("Corrector:  Before:     After:    Difference:\n");
    printf("             (mrad)     (mrad)       (mrad)  \n");
  }

  if (fcdata != NULL ) fprintf(fcdata, "\n! RESOUT = %d\n", resout);

  for (i = 0; i < icor; i++) { /* loop over all correctors */

    if (print_correct_opt > 1) {
      printf("%s %-3.6f %-3.6f %-3.6f\n", c[nc[i]].p_node->name,
       c[nc[i]].val.before[ip - 1],
       corvec[nx[i] - 1] + c[nc[i]].val.before[ip - 1],
       corvec[nx[i] - 1]);
    }

    c[nc[i]].val.after[ip - 1] = corvec[nx[i] - 1];
    if (ip == 1) {
      c[nc[i]].p_node->chkick += c[nc[i]].p_node->other_bv * 0.001 * corvec[nx[i] - 1];
      pro_correct_fill_corr_table(ip, c[nc[i]].p_node->name,
                                       c[nc[i]].val.before[ip - 1] * 0.001,
                                       c[nc[i]].p_node->chkick);
      /*                          c[nc[i]].p_node->other_bv*0.001*corvec[nx[i]-1]); */

      if (fcdata != NULL ) {
	fprintf(fcdata, "%s->hkick = %e; \t! %d\n", strip(c[nc[i]].p_node->name),
		c[nc[i]].p_node->other_bv * 0.001 * corvec[nx[i] - 1], resout);
      }


    } else if (ip == 2) {
      c[nc[i]].p_node->cvkick += c[nc[i]].p_node->other_bv * 0.001 * corvec[nx[i] - 1];
      pro_correct_fill_corr_table(ip, c[nc[i]].p_node->name,
                                      c[nc[i]].val.before[ip - 1] * 0.001,
                                      c[nc[i]].p_node->cvkick);
      /*                          c[nc[i]].p_node->other_bv*0.001*corvec[nx[i]-1]); */
      if (fcdata != NULL ) {
	fprintf(fcdata, "%s->vkick = %e; \t! %d\n", strip(c[nc[i]].p_node->name),
		c[nc[i]].p_node->other_bv * 0.001 * corvec[nx[i] - 1], resout);
      }
    }
  } /* end of loop over correctors */

  if (fcdata != NULL ) fflush(fcdata);
  if (fddata != NULL ) fflush(fddata);

}

static int pro_correct_getactive(int ip, int *nm, int *nx, int *nc,
         double *corvec, double *monvec, char *conm) {
  int imon, icor;
  int imona, icora;
  struct id_mic *m, *c;

  int debug = get_option("debug");

  m = correct_orbit->mon_table;
  imon = 0;
  imona = 0;
  while (m) {
    if (debug) {
      printf("from list: %d %d %s %s ", m->enable, m->id_ttb, m->p_node->name, m->p_node->base_name);
      printf("\t\t orbit readings: %d %f %f\n", ip, m->val.before[0], m->val.before[1]);
    }
    if (m->enable == 1) {
      monvec[imon] = m->val.before[ip - 1];
      nm[imon] = imona;
      imon++;
    }
    imona++;
    m = m->next;
  }

  c = correct_orbit->cor_table;
  icor = 0;
  icora = 0;
  while (c) {
    if (debug) {
      printf("from list: %d %d %s %s ", c->enable, c->id_ttb, c->p_node->name, c->p_node->base_name);
      printf("\t\t kicker readings: %f %f\n", c->val.before[0], c->val.before[1]);
    }
    if (c->enable == 1) {
      corvec[icor] = c->val.before[ip - 1];
      nx[icor] = icora;
      nc[icor] = icora;
      strcpy(conm, c->p_node->name);
      conm += 16;
      icor++;
    }
    icora++;
    c = c->next;
  }

  // 2013-Jun-24  12:21:49  ghislain:
  // following is a kludge to return a single value but has to be decoded on other side.
  if(icor >= 30000)
    fatal_error("Found more than 30000 correctors; decoding in mad_orbit.c will fail",
    "Please report this issue to MAD developpers (mad@cern.ch)");
  return (30000 * imon + icor);
}

static void correct_option(struct in_cmd* cmd) {
  int i;

  int debug = get_option("debug");

  if (debug) {
    fprintf(prt_file, "in coption routine\n");
    for (i = 0; i < cmd->tok_list->curr; i++) {
      fprintf(prt_file, "command(s): %s\n", cmd->tok_list->p[i]);
    }
  }

  if (par_present("seed", cmd->clone)) {
      int seed = command_par_value("seed", cmd->clone);
      init55(seed);
  }

  print_correct_opt = command_par_value("print", cmd->clone);

  if (debug) {
    if (print_correct_opt == 0)
      fprintf(prt_file, "print option not set\n");
    else
      fprintf(prt_file, "print option set\n");
  }

}

static void correct_getorbit(struct in_cmd* cmd) {
  (void) cmd;
}

static void correct_putorbit(struct in_cmd* cmd) {
// Jun 25, 2013 3:33:03 PM ghislain : DOC - this option is documented as deprecated but still alive
  int i;
  struct name_list* nl;
  char* filename = command_par_string("file", cmd->clone);
  char* table_name;

  current_twiss = clone_command(find_command("twiss", defined_commands));

  nl = current_twiss->par_names;
  for (i = 0; i < nl->curr; i++)
    nl->inform[i] = 0;

  pro_twiss();

  table_name = permbuff("orbit");
  orbit_table = make_table(table_name, "orbit", orbit_table_cols, orbit_table_types, current_sequ->n_nodes);
  add_to_table_list(orbit_table, table_register);
  fill_orbit_table(orbit_table, orbin_table);
  out_table("orbit", orbit_table, filename);

  current_twiss = delete_command(current_twiss);
}

static void correct_usekick(struct in_cmd* cmd) {
  char temp[12];
  int count = set_enable("kicker", cmd);
  sprintf(temp, "%d", count);
  put_info(temp, "corrector(s) affected");
}

static void correct_usemonitor(struct in_cmd* cmd) {
  char temp[12];
  int count = set_enable("monitor", cmd);
  sprintf(temp, "%d", count);
  put_info(temp, "monitor(s) affected");
}

// public interface

void store_orbit(struct command* comm, double* orbit) {
  if (par_present("x",  comm)) orbit[0] = command_par_value("x",  comm);
  if (par_present("px", comm)) orbit[1] = command_par_value("px", comm);
  if (par_present("y",  comm)) orbit[2] = command_par_value("y",  comm);
  if (par_present("py", comm)) orbit[3] = command_par_value("py", comm);
  if (par_present("t",  comm)) orbit[4] = command_par_value("t",  comm);
  if (par_present("pt", comm)) orbit[5] = command_par_value("pt", comm);
}

void pro_correct(struct in_cmd* cmd) {
  if      (strcmp(cmd->tok_list->p[0], "correct") == 0)    correct_correct(cmd);
  else if (strcmp(cmd->tok_list->p[0], "usekick") == 0)    correct_usekick(cmd);
  else if (strcmp(cmd->tok_list->p[0], "usemonitor") == 0) correct_usemonitor(cmd);
  else if (strcmp(cmd->tok_list->p[0], "getorbit") == 0)   correct_getorbit(cmd); // FIXME obsolete command; should be flagged and not call anything...
  else if (strcmp(cmd->tok_list->p[0], "putorbit") == 0)   correct_putorbit(cmd); // FIXME obsolete command; should be flagged and not call anything...
  else if (strcmp(cmd->tok_list->p[0], "readmytable") == 0) read_table(cmd);
  else if (strcmp(cmd->tok_list->p[0], "readcorr") == 0)   correct_readcorr(cmd); //FIXME not documented
  else if (strcmp(cmd->tok_list->p[0], "setcorr") == 0)    correct_setcorr(cmd); // FIXME not documented
  // else if (strcmp(cmd->tok_list->p[0], "prtcorr") == 0)    correct_prtcorr(cmd); // FIXME not documented
  else if (strcmp(cmd->tok_list->p[0], "coption") == 0)    correct_option(cmd);
}

void f_ctof(int *j, char *string, int *nel) {
  long i, flg = 0;

  for (i = 0; i < *nel; i++) {
    if (flg == 1) {
      string[i] = ' ';
      continue;
    }

    if (string[i] == '\0') {
      string[i] = ' ';
      flg = 1;
      continue;
    }
  }
  *j = i;
}

