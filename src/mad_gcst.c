#include "mad_def.h"
#include "mad_gcst.h"

/* Temporary file: global constants
   these constants will be split over their respective modules...
   and should be const pointers to constant values...
*/

// madx version and date from makefile

#define mkstr(a)  mkstr_(a)
#define mkstr_(a) #a

const char * const version_name   = mkstr(_VERSION);
const char * const version_arch   = sizeof(void*) == 4 ? "32" : sizeof(void*) == 8 ? "64" : "??";
const char * const version_ostype = mkstr(_VERSION_OSTYPE);
const char * const version_date   = mkstr(_VERSION_DATE);
const int          version_num    = _VERSION_NUM;

#undef mkstr
#undef mkstr_

// madx constants

const char* const functs[] = {
  "dummyfunction", "abs", "sqrt", "exp", "log", "log10",
  "sin", "cos", "tan", "asin", "acos",
  "atan", "sinh", "cosh", "tanh",
  "ranf", "gauss", "tgauss", "table", "exist",
  "floor", "ceil", "round", "frac",
  "erf", "erfc", "sinc",
  ""}; /* keep "" ! */

const char* const op_string = "-+*/^";
const char* const file_string = "file"; /* to avoid local in routine alias */
const char* const vrai = "true";        /* to avoid local in routine alias */
const char* const faux = "false";       /* to avoid local in routine alias */
const int n_match = 17;                 /* # of match token lists in cmd_match_base */
const int s_match[] =              /* position of first token of command below */
{0, 1, 4, 8, 13, 17, 22, 25, 29, 32, 36, 39, 43, 45, 48, 50, 52, 56};

const int t_match[] = /* order in which the commands are matched */
{0, 1, 16, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

const char* const cmd_match_base[] =
{ /*  0 */ "@cmd",
  /*  1 */ "@name", ":", "@cmd",
  /*  2 */ "int", "const", "@name", "=",
  /*  3 */ "int", "const", "@name", ":", "=",
  /*  4 */ "real", "const", "@name", "=",
  /*  5 */ "real", "const", "@name", ":", "=",
  /*  6 */ "int", "@name", "=",
  /*  7 */ "int", "@name", ":", "=",
  /*  8 */ "real", "@name", "=",
  /*  9 */ "real", "@name", ":", "=",
  /* 10 */ "const", "@name", "=",
  /* 11 */ "const", "@name", ":", "=",
  /* 12 */ "@name", "=",
  /* 13 */ "@name", ":", "=",
  /* 14 */ "@name", ":",
  /* 15 */ "@name", "@name",
  /* 16 */ "shared", "@name", ":", "@cmd" };

/* aperture types and # of parameters, needed for twiss table */

const char* const aperture_types[] =
{
  "circle", "ellipse", "rectangle", "lhcscreen", "rectcircle",
  "rectellipse", "racetrack", "octagon",
  " "  /* blank terminates */
};

const int aperture_npar[] =
{
  1, 2, 2, 3, 3,
  4, 4, 4
};

/* table descriptors: type 1 = int, type 2 = double, type 3 = string;
   internally, however, int are stored as double */

const int ap_table_types[] =
{
  3, 2, 2, 2, 3,
  2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2
};

const char* const ap_table_cols[] =
{
  "name", "n1", "n1x_m", "n1y_m", "apertype",
  "aper_1", "aper_2", "aper_3", "aper_4",
  "rtol", "xtol", "ytol", "xoffset", "yoffset",
  "s", "betx", "bety", "dx", "dy", "x", "y", "px", "py",
  "on_ap", "on_elem", "spec", "x_pos_hit", "y_pos_hit",
  " "  /* blank terminates */
};

const int survey_table_types[] =
{
  3, 3, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2,
  1, 1, 2,
  /*== dealt with the new property v_pos as for mech_sep */
  2,3
  /*==*/
};

const char* const survey_table_cols[] =
{
  "name", "keyword", "s", "l", "angle", "x",
  "y", "z", "theta", "phi", "psi", "globaltilt", "tilt",
  "slot_id", "assembly_id", "mech_sep",
  /*== dealt with the new property v_pos asc for mech_sep */
  "v_pos", "comments",
  /*==*/
  " "  /* blank terminates */
};



const int efield_table_types[] =
{
  3, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2,
  /* AL: RF-multipolar error */
  2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2
};

const char* const efield_table_cols[] =
{
  "name",
  "k0l", "k0sl", "k1l", "k1sl",
  "k2l", "k2sl", "k3l", "k3sl", "k4l",
  "k4sl", "k5l", "k5sl", "k6l", "k6sl",
  "k7l", "k7sl", "k8l", "k8sl", "k9l",
  "k9sl", "k10l", "k10sl", "k11l", "k11sl",
  "k12l", "k12sl", "k13l", "k13sl", "k14l",
  "k14sl", "k15l", "k15sl", "k16l", "k16sl",
  "k17l", "k17sl", "k18l", "k18sl", "k19l",
  "k19sl", "k20l", "k20sl",
  "dx", "dy", "ds", "dphi", "dtheta",
  "dpsi", "mrex", "mrey", "mredx", "mredy",
  "arex", "arey", "mscalx", "mscaly",
  /* AL: RF-multipolar errors */
  "rfm_freq", "rfm_harmon", "rfm_lag",
  "p0l", "p0sl", "p1l", "p1sl",
  "p2l", "p2sl", "p3l", "p3sl", "p4l",
  "p4sl", "p5l", "p5sl", "p6l", "p6sl",
  "p7l", "p7sl", "p8l", "p8sl", "p9l",
  "p9sl", "p10l", "p10sl", "p11l", "p11sl",
  "p12l", "p12sl", "p13l", "p13sl", "p14l",
  "p14sl", "p15l", "p15sl", "p16l", "p16sl",
  "p17l", "p17sl", "p18l", "p18sl", "p19l",
  "p19sl", "p20l", "p20sl",
  " "  /* blank terminates */
};


const char* const sxf_table_names[] =
{
  "l","angle", "k0","k0s","k1","k1s",
  "e1","e2","k2","k2s","h1",
  "h2","hgap","fint","k3","k3s",
  "lrad","knl","ksl","ks","volt",
  "lag","harmon","betrf","pg",
  "shunt","tfill","eloss","ex","ey",
  "hkick","vkick","xsize","ysize","sigx",
  "sigy","xma","yma","charge", "freq",
  " " /* blank terminates */
};

const int twiss_opt_end = 33; /* last column filled by twiss module */
const int twiss_mult_end = 78; /* last multipole column filled
                            by complete_twiss_table */
const int twiss_fill_end = 102; /* last standard column filled
                            by complete_twiss_table */
/*== increased twiss_fill_end from 96 to 97 to accomodate for v_pos */

/* warning: modify routine complete_twiss_table in case of changes */
const int twiss_table_types[] =
{
  3, 3, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2,
  2, 2, 2, 2,
  2, 2, 2, 2,
  2, 2, 2, 2,
  2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2,
  2, 2, 2, 2,
  1, 1, 2,
  2, /* v_pos */
  2, 2, 2, 2, 2,
  2, 3, 3,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2,
  2, 2, 2,
  2, 2, 2,
  2, 2, 2,
  2, 2, 2,
  2, 2, 2,
  2, 2, 2,
  2, 2, 2,
  2, 2, 2,
  /* delta_p dependency terms */
  2,2,2, /* beta11p, beta12p, beta13p */
  2,2,2, /* beta21p, beta22p, beta23p  */
  2,2,2, /* beta31p, beta32p, beta33p  */
  2,2,2, /* alfa11p, alfa12p, alfa13p */
  2,2,2, /* alfa21p, alfa22p, alfa23p */
  2,2,2, /* alfa31p, alfa32p, alfa33p */
  2,2,2, /* gama11p, gama12p, gama13p */
  2,2,2, /* gama21p, gama22p, gama23p */
  2,2,2, /* gama31p, gama32p, gama33p */
  /* end of delta_p dependency terms */
  2, 2, 2, 2,
  /* derivatives of dispersion w.r.t. delta_p */
  2, 2, 2, 2,
  2, 2, 2, 2, /* second order derivatives */
  2, 2, 2, 2, /* third order derivatives */
  /* end of dispersion derivatives w.r.t. delta_p */
  2, 2, 2, /* mu1, mu2, mu3 */
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2
};

const char* const twiss_table_cols[] =
{
  "name", "keyword", "s", "betx", "alfx",
  "mux", "bety", "alfy", "muy", "x",
  "px", "y", "py", "t", "pt",
  "dx", "dpx", "dy", "dpy", "wx",
  "phix", "dmux", "wy", "phiy", "dmuy",
  "ddx", "ddpx", "ddy", "ddpy", "r11",
  "r12", "r21", "r22", "energy", "l",
  "angle", "k0l", "k0sl", "k1l", "k1sl",
  "k2l", "k2sl", "k3l", "k3sl", "k4l",
  "k4sl", "k5l", "k5sl", "k6l", "k6sl",
  "k7l", "k7sl", "k8l", "k8sl", "k9l",
  "k9sl", "k10l", "k10sl", "k11l", "k11sl",
  "k12l", "k12sl", "k13l", "k13sl",
  "k14l", "k14sl", "k15l", "k15sl",
  "k16l", "k16sl", "k17l", "k17sl",
  "k18l", "k18sl", "k19l", "k19sl",
  "k20l", "k20sl", "ksi", "hkick",
  "vkick", "tilt", "e1", "e2", "h1",
  "h2", "hgap", "fint", "fintx",
  "volt", "lag", "freq", "harmon",
  "slot_id","assembly_id","mech_sep",
  /*== dealt with the new property v_pos as for mech_sep */
  "v_pos",
  "bbcharge","xma", "yma", "sigx", "sigy",
  /*==*/
  "lrad","parent","comments",
  "re11", "re12", "re13", "re14", "re15", "re16",
  "re21", "re22", "re23", "re24", "re25", "re26",
  "re31", "re32", "re33", "re34", "re35", "re36",
  "re41", "re42", "re43", "re44", "re45", "re46",
  "re51", "re52", "re53", "re54", "re55", "re56",
  "re61", "re62", "re63", "re64", "re65", "re66",
  "kmax", "kmin", "calib", "polarity", "alfa",
  "beta11", "beta12", "beta13",
  "beta21", "beta22", "beta23",
  "beta31", "beta32", "beta33",
  "alfa11", "alfa12", "alfa13",
  "alfa21", "alfa22", "alfa23",
  "alfa31", "alfa32", "alfa33",
  "gama11", "gama12", "gama13",
  "gama21", "gama22", "gama23",
  "gama31", "gama32", "gama33",
  /* delta_p dependency: derivatives of the above Twiss parameters */
  "beta11p","beta12p","beta13p",
  "beta21p","beta22p","beta23p",
  "beta31p","beta32p","beta33p",
  "alfa11p", "alfa12p","alfa13p",
  "alfa21p", "alfa22p","alfa23p",
  "alfa31p", "alfa32p","alfa33p",
  "gama11p", "gama12p","gama13p",
  "gama21p", "gama22p","gama23p",
  "gama31p", "gama32p","gama33p",
  /* end of delta_p dependency */
  "disp1", "disp2", "disp3","disp4",
  /* derivatives of dispersion w.r.t. delta_p */
  "disp1p", "disp2p", "disp3p", "disp4p",
  "disp1p2", "disp2p2", "disp3p2", "disp4p2", /* second order derivatives */
  "disp1p3", "disp2p3", "disp3p3", "disp4p3", /* third order derivatives */
  /* end of dispersion derivatives w.r.t. delta_p */
  "mu1", "mu2", "mu3",
  /* IT  sigma matrix */
  "sig11", "sig12", "sig13", "sig14", "sig15", "sig16",
  "sig21", "sig22", "sig23", "sig24", "sig25", "sig26",
  "sig31", "sig32", "sig33", "sig34", "sig35", "sig36",
  "sig41", "sig42", "sig43", "sig44", "sig45", "sig46",
  "sig51", "sig52", "sig53", "sig54", "sig55", "sig56",
  "sig61", "sig62", "sig63", "sig64", "sig65", "sig66",
  /*  "eign11", "eign12", "eign13", "eign14", "eign15", "eign16",
      "eign21", "eign22", "eign23", "eign24", "eign25", "eign26",
      "eign31", "eign32", "eign33", "eign34", "eign35", "eign36",
      "eign41", "eign42", "eign43", "eign44", "eign45", "eign46",
      "eign51", "eign52", "eign53", "eign54", "eign55", "eign56",
      "eign61", "eign62", "eign63", "eign64", "eign65", "eign66",*/
  "n1",
  " "  /* blank terminates */
};

const int twiss_sector_table_types[] = {
  3, 2,
  2, 2, 2, 2, 2, 2,
  /* 36 elements for the R-matrix */
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  /* 216 elements for the T-matrix */
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2
};

const char* const twiss_sector_table_cols[] = {
  "name", "pos",
  "k1", "k2", "k3", "k4", "k5", "k6",
  "r11", "r21", "r31", "r41", "r51", "r61",
  "r12", "r22", "r32", "r42", "r52", "r62",
  "r13", "r23", "r33", "r43", "r53", "r63",
  "r14", "r24", "r34", "r44", "r54", "r64",
  "r15", "r25", "r35", "r45", "r55", "r65",
  "r16", "r26", "r36", "r46", "r56", "r66",
  "t111", "t211", "t311", "t411", "t511", "t611",
  "t121", "t221", "t321", "t421", "t521", "t621",
  "t131", "t231", "t331", "t431", "t531", "t631",
  "t141", "t241", "t341", "t441", "t541", "t641",
  "t151", "t251", "t351", "t451", "t551", "t651",
  "t161", "t261", "t361", "t461", "t561", "t661",
  "t112", "t212", "t312", "t412", "t512", "t612",
  "t122", "t222", "t322", "t422", "t522", "t622",
  "t132", "t232", "t332", "t432", "t532", "t632",
  "t142", "t242", "t342", "t442", "t542", "t642",
  "t152", "t252", "t352", "t452", "t552", "t652",
  "t162", "t262", "t362", "t462", "t562", "t662",
  "t113", "t213", "t313", "t413", "t513", "t613",
  "t123", "t223", "t323", "t423", "t523", "t623",
  "t133", "t233", "t333", "t433", "t533", "t633",
  "t143", "t243", "t343", "t443", "t543", "t643",
  "t153", "t253", "t353", "t453", "t553", "t653",
  "t163", "t263", "t363", "t463", "t563", "t663",
  "t114", "t214", "t314", "t414", "t514", "t614",
  "t124", "t224", "t324", "t424", "t524", "t624",
  "t134", "t234", "t334", "t434", "t534", "t634",
  "t144", "t244", "t344", "t444", "t544", "t644",
  "t154", "t254", "t354", "t454", "t554", "t654",
  "t164", "t264", "t364", "t464", "t564", "t664",
  "t115", "t215", "t315", "t415", "t515", "t615",
  "t125", "t225", "t325", "t425", "t525", "t625",
  "t135", "t235", "t335", "t435", "t535", "t635",
  "t145", "t245", "t345", "t445", "t545", "t645",
  "t155", "t255", "t355", "t455", "t555", "t655",
  "t165", "t265", "t365", "t465", "t565", "t665",
  "t116", "t216", "t316", "t416", "t516", "t616",
  "t126", "t226", "t326", "t426", "t526", "t626",
  "t136", "t236", "t336", "t436", "t536", "t636",
  "t146", "t246", "t346", "t446", "t546", "t646",
  "t156", "t256", "t356", "t456", "t556", "t656",
  "t166", "t266", "t366", "t466", "t566", "t666",
  " " /* blank terminates */
};


const int ptc_twiss_summary_table_types[] =
  {
    2, 2, 2, 2, 2, 2, 2,/* "length", "alpha_c", "alpha_c_p", "alpha_c_p2", "alpha_c_p3", "eta_c", "gamma_tr", */
    2, 2, 2, 2, 2,	    /* "q1", "q2", "dq1", "dq2", "qs", */
    2, 2, 	    /* "beta_x_min","beta_x_max", */
    2, 2, 	    /* "beta_y_min","beta_y_max", */
    2, 2,	    /* "beta11min","beta11max", */
    2, 2, 	    /* "beta12min","beta12max", */
    2, 2, 	    /* "beta13min","beta13max", */
    2, 2, 	    /* "beta21min","beta21max", */
    2, 2, 	    /* "beta22min","beta22max", */
    2, 2, 	    /* "beta23min","beta23max", */
    2, 2, 	    /* "beta31min","beta31max", */
    2, 2, 	    /* "beta32min","beta32max", */
    2, 2, 	    /* "beta33min","beta33max", */
    2, 2, 	    /* "disp1min", "disp1max", */
    2, 2, 	    /* "disp2min", "disp2max", */
    2, 2, 	    /* "disp3min", "disp3max", */
    2, 2, 	    /* "disp4min", "disp4max", */
    2,	    /* "deltap", */
    2,2,2,	    /* "orbit_x","orbit_px","orbit_y", */
    2,2,2,	    /* "orbit_py","orbit_pt","orbit_-cT", */
    2,2,2,2,2,2,	    /* "xcorms","ycorms","pxcorms","pycorms","tcorms","ptcorms", */
    2,2,2,2,2,2,	    /* "xcomax","ycomax","pxcomax","pycomax","tcomax","ptcomax", */
    2,2,2,2,2,2	    /* "xcomin","ycomin","pxcomin","pycomin","tcomin","ptcomin", */
  };
const char* const ptc_twiss_summary_table_cols[] = {
  "length", "alpha_c", "alpha_c_p", "alpha_c_p2", "alpha_c_p3", "eta_c", "gamma_tr",
  "q1", "q2", "dq1", "dq2", "qs",
  "beta_x_min","beta_x_max",
  "beta_y_min","beta_y_max",
  "beta11min","beta11max",
  "beta12min","beta12max",
  "beta13min","beta13max",
  "beta21min","beta21max",
  "beta22min","beta22max",
  "beta23min","beta23max",
  "beta31min","beta31max",
  "beta32min","beta32max",
  "beta33min","beta33max",
  "disp1min","disp1max",
  "disp2min","disp2max",
  "disp3min","disp3max",
  "disp4min","disp4max",
  "deltap",
  "orbit_x","orbit_px","orbit_y",
  "orbit_py","orbit_pt","orbit_t",
  "xcorms","ycorms","pxcorms","pycorms","tcorms","ptcorms",
  "xcomax","ycomax","pxcomax","pycomax","tcomax","ptcomax",
  "xcomin","ycomin","pxcomin","pycomin","tcomin","ptcomin",
  " " /* blank terminates */
};

const int ibs_table_types[] =
{
  3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2
};

const char* const ibs_table_cols[] =
{
  "name", "s", "dels", "tli", "txi", "tyi", "betx","alfx","dx","dpx", "bety","alfy","dy","dpy",
  " "  /* blank terminates */
};

const int bb6d_ixy_types[]=
{
  2,2,2,2,2,2,2
};

const char* const bb6d_ixy_cols[]=
{
  "turn","n_macro_surv","n_for_i","ex_rms","ey_rms","sigma_p","sigma_z",
  " "  /* blank terminates */
};

const int map_tab_types[]=
{
  3,2,1,1,1,1,1,1,1,1,1
};

const char* const map_tab_cols[]=
{
  "name","coef","n_vector","nv","order","nx","nxp","ny","nyp","ndeltap","nt",
  " "  /* blank terminates */
};

const int normal_res_types[] =
{
  3, 1, 1, 1, 1, 2
};

const char* const normal_res_cols[] =
{
  "name", "order1", "order2", "order3", "order4", "value",
  " "  /* blank terminates */
};

const int sodd_detune_5_types[] =
{
  1, 1, 2, 1, 1
};

const char* const sodd_detune_5_cols[] =
{
  "multipoleorder", "plane", "detuning", "h_inv_order", "v_inv_order",
  " "  /* blank terminates */
};

const int sodd_distort1_8_types[] =
{
  2, 2, 2, 2, 2, 2, 2, 2
};

const char* const sodd_distort1_8_cols[] =
{
  "multipoleorder", "cosine", "sine", "amplitude", "j", "k", "l", "m",
  " "  /* blank terminates */
};

const int sodd_distort1_11_types[] =
{
  1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1
};

const char* const sodd_distort1_11_cols[] =
{
  "multipoleorder", "location", "resonance", "position[m]", "cosine", "sine", "amplitude", "j", "k", "l", "m",
  " "  /* blank terminates */
};

const int sodd_distort2_9_types[] =
{
  1, 1, 2, 2, 2, 1, 1, 1, 1
};

const char* const sodd_distort2_9_cols[] =
{
  "multipoleorder1", "multipoleorder2", "cosine", "sine", "amplitude", "j", "k", "l", "m",
  " "  /* blank terminates */
};

const int touschek_table_types[] =
{
  3, 2, 2, 2, 2
};

const char* const touschek_table_cols[] =
{
  "name", "s", "tli", "tliw", "tlitot",
  " "  /* blank terminates */
};

const int mon_table_types[] =
{
  3, 2, 2, 2, 2
};

const char* const mon_table_cols[] =
{
  "name", "x.old", "y.old", "x", "y",
  " "  /* blank terminates */
};

const int corr_table_types[] =
{
  3, 2, 2, 2, 2
};

const char* const corr_table_cols[] =
{
  "name", "px.old", "py.old", "px.correction", "py.correction",
  " "  /* blank terminates */
};

const int orbit_table_types[] =
{
  3, 2, 2, 1,
};

const char* const orbit_table_cols[] =
{
  "name", "x", "y", "status",
  " "  /* blank terminates */
};

const int special_comm_cnt[] =
{
  3, 5, 7, 6, 5, 4,
  0
};

const char* const special_comm_desc[] = /* ">?" = skip from start including char. at ? */
{
  "if(", "else{", "elseif(", "while(", ">:macro", ">:line",
  " "  /* blank terminates , line must remain last */
};

const int summ_table_types[] =
{
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 
  //2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2,
  2, //for nflips
};

const char* const summ_table_cols[] =
{
  "length", "orbit5", "alfa", "gammatr", "q1",
  "dq1", "betxmax", "dxmax", "dxrms", "xcomax",
  "xcorms", "q2", "dq2", "betymax", "dymax",
  "dyrms", "ycomax", "ycorms", "deltap",
  //"synch_1","synch_2","synch_3","synch_4","synch_5",
  "synch_1","synch_2","synch_3","synch_4","synch_5",
  "synch_6","synch_8",
   "nflips", //for nflips
  " "  /* blank terminates */
};

// different layout for SUMM table columns
// char* summ_table_cols[] =
// {
//   "length", "orbit5", "alfa", "gammatr", "deltap",
//   "q1", "dq1", "betxmax", "dxmax", "dxrms", "xcomax", "xcorms",
//   "q2", "dq2", "betymax", "dymax", "dyrms", "ycomax", "ycorms",
//   "synch_1","synch_2","synch_3","synch_4","synch_5",
//   " "  /* blank terminates */
// };

const int trackone_table_types[] =
{
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2
};

const char* const trackone_table_cols[] =
{
  "number", "turn", "x", "px", "y", "py", "t", "pt", "s", "e",
  " "  /* blank terminates */
};

const int track_table_types[] =
{
  1, 1, 2, 2, 2, 2, 2, 2, 2, 2
};

const char* const track_table_cols[] =
{
  "number", "turn", "x", "px", "y", "py", "t", "pt", "s", "e",
  " "  /* blank terminates */
};

const char* const dist_table_cols[] =
{
  "number", "x", "px", "y", "py", "t", "pt",
  " "  /* blank terminates */
};

const int dist_table_types[] =
{
  1, 2, 2, 2, 2, 2, 2
};


const int track_table_cols_len = sizeof track_table_cols / sizeof track_table_cols[0];

const int tracksumm_table_types[] =
{
  1, 1, 2, 2, 2, 2, 2, 2, 2, 2
};

const char* const tracksumm_table_cols[] =
{
  "number", "turn", "x", "px", "y", "py", "t", "pt", "s", "e",
  " "  /* blank terminates */
};

const int ptcnodetrack_table_types[] =
{  1,        3,      1,         1,           1,      2,       2,   2,   2,    2,   2,    2,   2,    2 };

const char* const ptcnodetrack_table_cols[] =
{"number", "name", "elnumber","trnumber" , "turn","s_slice", "s", "x", "px", "y", "py", "t", "pt", "s",
 " "  /* blank terminates */
};

const int trackloss_table_types[] =
{
  1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3
};

const char* const trackloss_table_cols[] =
{
  "number", "turn", "x", "px", "y", "py", "t", "pt", "s", "e", "element",
  " "  /* blank terminates */
};

const int dynap_table_types[] =
{
  2,2,2,2,
  2,2,2,2,2,2,
  2,2,2,2,2
};

const char* const dynap_table_cols[] =
{
  "dktrturns", "xend", "pxend", "yend",
  "pyend", "tend", "ptend", "wxmin", "wxmax", "wymin", "wymax",
  //"wxymin", "wxymax", "smear", "yapunov",
  "wxymin", "wxymax", "smear", "lyapunov",
  " "  /* blank terminates */
};

const int dynaptune_table_types[] =
{
  2,2,2,2,2
};

const char* const dynaptune_table_cols[] =
{
  "x", "y", "tunx", "tuny", "dtune",
  " "  /* blank terminates */
};

/* Definition of "select_ptc_normal" parameters for "ptc_normal" */
const char* const names[]=
{
  "dx","dpx","dy","dpy","q1","q2","dq1","dq2","anhx","anhy","haml","gnfu","eign"," "
};

const char* const atm[] =
{
 "hmon","vmon","moni"," "
};

const char* const atc[] =
{
 "hkic","vkic","kick"," "
};

const char* atc_type = 0;
int   atc_flag = 0;

const char* atm_type = 0;
int   atm_flag = 0;


const char* const nonlin_table_cols[] =
{
  "name", "nickname", "basevariable", "value",
  "order", "order_x", "order_px","order_y",
           "order_py","order_pt","order_t", " "
};

const int nonlin_table_types[] =
{
  3, 3, 3, 2,
  1, 1, 1, 1,
     1, 1, 1
};

