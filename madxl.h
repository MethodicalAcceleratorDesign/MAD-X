/* defined constants for word lengths etc. */
#define ALIGN_MAX 14        /* alignment error array length */
#define FIELD_MAX 42        /* field error array length */
#define SEQ_DUMP_LEVEL 0    /* chooses amount of dumped output */
#define NAME_L 24           /* internal name length */
#define FNAME_L 240         /* for file names */
#define AUX_LG 10000        /* for all sorts of ancillary buffers */
#define INVALID 1.e20       /* used for erroneous value requests */
#define MAX_ITEM  1000      /* initial # of items in tok_list etc. */
#define MAX_D_ITEM 30000    /* storage for doubles */
#define MAX_LINE 20000      /* max. input line length */
#define MAX_LOOP 100        /* max. count for (possibly circular) calls */
#define MAX_COND 100        /* max. nesting level for "if" and "while" */
#define CHAR_BUFF_SIZE 100000 /* size of each dynamic char_buff member */
#define IN_BUFF_SIZE 500000 /* size of buffer for command groups */
#define LINE_FILL 70        /* max. line length -2 for "save" output */
#define LINE_F_MAD8 70      /* the same, for mad-8 format */
#define ELEM_OUT 1000       /* character buffer for writing one element */
#define MAX_RAND 1000000000 /* for random generator */
#define NR_RAND 55          /* for random generator */
#define NJ_RAND 24          /* for random generator */
#define ND_RAND 21          /* for random generator */
#define MATCH_WORK 10       /* no. of work spaces in matching */
#define TABLE_KEY 16        /* max. table keyword is TABLE_KEY-1 */
#define TRACK_ROWS 1200     /* initial length of track tables */

char* const functs[] = {"dummyfunction", "abs", "sqrt", "exp", "log", "log10",
                        "sin", "cos", "tan", "asin", "acos",
                        "atan", "sinh", "cosh", "tanh", "ranf",
                        "gauss", "tgauss", "table",
                        ""}; /* keep "" ! */

const char op_string[] = "-+*/^";
char file_string[] = "file"; /* to avoid local in routine alias */

const int n_match = 17; /* # of match token lists in cmd_match_base */
const int s_match[] = /* position of first token of command below */
{0, 1, 4, 8, 13, 17, 22, 25, 29, 32, 36, 39, 43, 45, 48, 50, 52, 56};

const int t_match[] = /* order in which the commands are matched */
{0, 1, 16, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
const char cmd_match_base[][TABLE_KEY] = 
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
  /* 16 */ "shared", "@name", ":", "@cmd"};

/* table descriptors: type 1 = int, type 2 = double, type 3 = string;
   internally, however, int are stored as double */

int survey_table_types[] =
{
3, 2, 2, 2, 2,
2, 2, 2, 2, 2
};

char survey_table_cols[][TABLE_KEY] =
{
"name", "s", "l", "angle", "x",
"y", "z", "theta", "phi", "psi",
" "  /* blank terminates */
};

int efield_table_types[] =
{
3, 2, 2, 2, 2,
2, 2, 2, 2, 2,
2, 2, 2, 2, 2,
2, 2, 2, 2, 2,
2, 2, 2
};

char efield_table_cols[][TABLE_KEY] =
{
"name", 
"k0l", "k0sl", "k1l", "k1sl", 
"k2l", "k2sl", "k3l", "k3sl", "k4l", 
"k4sl", "k5l", "k5sl", "k6l", "k6sl", 
"k7l", "k7sl", "k8l", "k8sl", "k9l", 
"k9sl", "k10l", "k10sl",
" "  /* blank terminates */
};

int twiss_opt_end = 33; /* last column filled by twiss module */
int twiss_fill_end = 60; /* last column filled by twiss_complete */
int twiss_table_types[] =
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
2, 3,
2, 2, 2, 2, 2,
2, 2, 2, 2, 2,
2, 2, 2, 2, 2,
2, 2, 2, 2, 2,
2, 2, 2, 2, 2,
2, 2, 2, 2, 2,
2, 2, 2, 2, 2,
2
};

char twiss_table_cols[][TABLE_KEY] =
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
"k9sl", "k10l", "k10sl", "ks", "hkick", 
"vkick", "parent",
"re11", "re12", "re13", "re14", "re15",
"re16", "re21", "re22", "re23", "re24",
"re25", "re26", "re31", "re32", "re33",
"re34", "re35", "re36", "re41", "re42",
"re43", "re44", "re45", "re46", "re51",
"re52", "re53", "re54", "re55", "re56", 
"re61", "re62", "re63", "re64", "re65",
"re66", 
" "  /* blank terminates */
};

int ibs_table_types[] =
{
  3, 2, 2, 2, 2, 2
};

char ibs_table_cols[][TABLE_KEY] =
{
  "name", "s", "dels", "tli", "txi", "tyi",
" "  /* blank terminates */
};

int mon_table_types[] =
{
3, 2, 2, 2, 2
};

char mon_table_cols[][TABLE_KEY] =
{
  "name", "x(old)", "y(old)", "x(new)", "y(new)",
" "  /* blank terminates */
};

int corr_table_types[] =
{
3, 2, 2, 2, 2
};

char corr_table_cols[][TABLE_KEY] =
{
  "name", "px(old)", "py(old)", "px(new)", "py(new)",
" "  /* blank terminates */
};

char inverted_mag_forces[][TABLE_KEY] =
{
"angle", "k0", "k0s", "kick", "hkick", "vkick", 
"k1", "k1s", "k2", "k2s", "k3", "k3s", "k4", "k4s",
" "  /* blank terminates */
};

int orbit_table_types[] =
{
  3, 2, 2, 1,
};

char orbit_table_cols[][TABLE_KEY] =
{
  "name", "x", "y", "status",
" "  /* blank terminates */
};

int special_comm_cnt[] =
{
  3, 5, 7, 6, 5, 4,
0
};

char special_comm_desc[][TABLE_KEY] = /* ">?" = skip from start including char. at ? */
{
  "if(", "else{", "elseif(", "while(", ">:macro", ">:line",
" "  /* blank terminates , line must remain last */
};

int summ_table_types[] =
{
2, 2, 2, 2, 2,
2, 2, 2, 2, 2,
2, 2, 2, 2, 2,
2, 2, 2, 2
};

char summ_table_cols[][TABLE_KEY] =
{
"length", "orbit5", "alfa", "gammatr", "q1",
"dq1", "betxmax", "dxmax", "dxrms", "xcomax",
"xcorms", "q2", "dq2", "betymax", "dymax", 
"dyrms", "ycomax", "ycorms", "deltap",
" "  /* blank terminates */
};

int trackone_table_types[] =
{
1, 1, 2, 2, 2, 2,
2, 2
};

char trackone_table_cols[][TABLE_KEY] =
{
"turn", "number", "x", "px", "y", "py",
"t", "pt",
" "  /* blank terminates */
};

int track_table_types[] =
{
1, 2, 2, 2, 2,
2, 2
};

char track_table_cols[][TABLE_KEY] =
{
"turn", "x", "px", "y", "py",
"t", "pt",
" "  /* blank terminates */
};

int tracksumm_table_types[] =
{
1, 1, 2, 2, 2, 2,
2, 2
};

char tracksumm_table_cols[][TABLE_KEY] =
{
"number", "turn", "x", "px", "y", "py",
"t", "pt",
" "  /* blank terminates */
};

int dynap_table_types[] =
{
2,2,2,2,2,
2,2,2,2,2,
2,2,2,2,2
};

char dynap_table_cols[][TABLE_KEY] =
{
"dynapfrac", "dktrturns", "xend", "pxend", "yend",
"pyend", "tend", "wxmin", "wxmax", "wymin", "wymax",
"wxymin", "wxymax", "smear", "yapunov",
" "  /* blank terminates */
};

int dynaptune_table_types[] =
{
  2,2,2,2,2
};

char dynaptune_table_cols[][TABLE_KEY] =
{
"x", "y", "tunx", "tuny", "dtune",
" "  /* blank terminates */
};


