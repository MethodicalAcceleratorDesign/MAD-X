#include "madx.h"
#include "mad_6track_name_mangler.h"

/*---------------------------------------------------------------------*
 *                                                                      *
 *                           CERN                                       *
 *                                                                      *
 *     European Organization for Nuclear Research                       *
 *                                                                      *
 *     Program name: c6t: MAD-Sixtrack Format Converter                 *
 *                                                                      *
 *     Author and contact:   Hans GROTE                                 *
 *                           SL Division                                *
 *                           CERN                                       *
 *                           CH-1211 GENEVA 23                          *
 *                           SWITZERLAND                                *
 *                      Tel. [041] (022) 767 49 61                      *
 *                           Hans.Grote@cern.ch                         *
 *                                                                      *
 *     Converted to MAD-X by Mark HAYES                                 *
 *     Followed up by Frank Schmidt                                     *
 *                                                                      *
 *     Copyright  CERN,  Geneva  2000  -  Copyright  and  any   other   *
 *     appropriate  legal  protection  of  this  computer program and   *
 *     associated documentation reserved  in  all  countries  of  the   *
 *     world.                                                           *
 *                                                                      *
 *     Organizations collaborating with CERN may receive this program   *
 *     and documentation freely and without charge.                     *
 *                                                                      *
 *     CERN undertakes no obligation  for  the  maintenance  of  this   *
 *     program,  nor responsibility for its correctness,  and accepts   *
 *     no liability whatsoever resulting from its use.                  *
 *                                                                      *
 *     Program  and documentation are provided solely for the use  of   *
 *     the organization to which they are distributed.                  *
 *                                                                      *
 *     This program  may  not  be  copied  or  otherwise  distributed   *
 *     without  permission. This message must be retained on this and   *
 *     any other authorized copies.                                     *
 *                                                                      *
 *     The material cannot be sold. CERN should be  given  credit  in   *
 *     all references.                                                  *
 *                                                                      *
 *---------------------------------------------------------------------*/

/* 17.08.2004 - FS fix print-out of special f34 file needed as input file
   for the sodd program. */
/* 15/03/2004 - FS fixing faulty variable passing to "create_aperture" */
/* 01/07/2003 - FS added the "arbitrary matrix" element */
/* 21/03/2003 - FS fixed segmentation fault which was due to a faulty
   free-ing of object that had already been freed before */
/* 10/07/2002 - MH fixed missing mcdo bug, caused by recursion up
   element tree to unexpanded double_array */
/* 20/06/2002 - MH fixed double declarations and memory leaks because the
   original c6t was only ment to be run once - but not this one! */
/* 19/06/2002 - MH found last 'bug' in rhic sequence... due to micron
   length quadrupole and rounding errors caused by them */
/* 29/04/2002 - MH&HG made it copy all the collimators across */
/* 23/04/2002 - MH changed ref_delta so that =0 for quads and higher
   this directly effects fort.3 */
/* 14/04/2002 - MH changed calloc to mycalloc for HG error checking */
/* extract Sixtrack input files from DOOM */
/* question 1: change BEAM defaults to those of MAD8 */
/* question 2: att_lcavity ? */
/* #define _call_tree_ */

/* JMJ, 7/11/2002 commenting out the following
   to see if it helps  for Visual Fortran ....
   already in madxn.c

   #include <string.h>
   #include <stdio.h>
   #include <stdlib.h>
   #include <sys/types.h>
   #include <ctype.h>
   #include <math.h>
   #include <time.h>

   and I moved
   #include "c6t.h"
   to madxn.c

*/

#define BASE_TYPES 100    /* maximum no. of element types allowed */
#define EL_COUNT   100    /* initial length of element type list */

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502
#endif

// types

/* MADX name and internal codes    :
   circle=1=CR, rectangle=2=RE, ellipse=3=EL, rectcircle=lhcscreen=4=RC,
   rectellipse=5=RL, racetrack=6=RT, octagon=7=OC */
struct aper_struct {
    int apply;
    char name[255];
    char style[255];
    double value[8]; // 2015-Jul-31 ghislain: adding more parameters to aperture
};

struct c6t_element
{
    char name[48],           /* name including occurrence count if > 1 */
      org_name[48],            /* original name */
      base_name[48];           /* basic type */
    struct c6t_element* previous;
    struct c6t_element* next;
    struct c6t_element* equiv;   /* pointer to first identical element */
    int flag;                /* treatment flag_1 or _2 or _3 (see type_info) */
    int force;               /* magnet flag (flag_5) (see type_info) */
    int c_drift;             /* treatment flag_4 (see type_info) */
    int split;               /* treatment flag_6 (see type_info) */
    int n_values;            /* length of value */
    int w_flag;              /* 0 if not, 1 if written on fort.2 */
    int out_1;               /* output parameter 1, fort.2 */
    int na_err;              /* current no. of alignment errors */
    int nf_err;              /* current no. of field errors */
    int nc_pos;              /* component count, only multipoles */
    int npole_sign;          /* sign inversion flag for even (created) npoles */
    int keep_in;             /* if not 0, do not yank */
    int mult_order;          /* error reference comp., only multipoles */
    int f3_flag;             /* for multipole def. on fc.3 */
    int occ_cnt;             /* occurrence count */
    int twtab_row;           /* row number in twiss table */
    double position;         /* s position in sequence [m] */
    double rad_length;       /* radiation length of multipoles [m] */
    double ref_radius;       /* reference radius for multipole errors [m] */
    double ref_delta;        /* reference delta for multipole errors */
    double out_2;            /* output parameter 2, fort.2 */
    double out_3;            /* output parameter 3, fort.2 */
    double out_4;            /* output parameter 4, fort.2 */
    double out_5;            /* output parameter 5, fort.2 */
    double out_6;            /* output parameter 6, fort.2 */
    double out_7;            /* output parameter 7, fort.2 */
    double* value;           /* element strength etc. values */
    struct object* p_al_err; /* pointer to alignment error object */
    struct object* p_fd_err; /* pointer to field error object */
    struct object* p_ph_err; /* pointer to field phase error array AL: */
    double rfm_freq;         /* frequency of the rf-multipole fields  AL: */
    int tilt_err;            /* allow write_f8 to dump tilt as well */
    int do_not_free;         /* avoid free crash */
};

struct c6t_el_list /* contains list of element pointers */
{
    int   max,                /* max. pointer array size */
          curr;               /* current occupation */
    char  base_name[48];
    struct c6t_element** elem; /* element pointer list */
};

struct block
{
    char    name[48];
    double  length;
    int     flag;              /* if 0 take element, else block */
    struct c6t_element* first;
    struct c6t_element* last;
    struct block* previous;
    struct block* next;
    struct block* equiv;
    struct c6t_el_list* elements;
};

struct li_list /* contains list of list pointers */
{
    int curr;               /* current occupation */
    struct c6t_el_list* member[BASE_TYPES]; /* list pointer list */
};

struct type_info /* info about types */
{
    char name[48];    /* base_type */
    /* flag meanings - 0: skip, 1: linear, >1: non-linear,
       2: convert to multipole (temporarily), 3: cavity
       4: make 2 if in explicit list, else skip
       5: only split */
    int flag_1,  /* for length = 0 */
        flag_2,  /* for length > 0, normal */
        flag_3,  /* for length > 0, skew */
        flag_4,  /* if > 0: make drift, print warning when encountered */
        flag_5,  /* if > 0: magnet (for k0n * l) */
        flag_6;  /* if length > 0: 0 = no split
                    1 = split; if flag_2(_3) = 1: two identical + zero m.
                    if flag_2(_3) = 2: two drift + full m. */
};

/* (some) constants and structure definitions for DOOM */

/*
  element definition: (F: Fortran, C: C/C++)
  d.p. array
  word / e_type = 1            = 2                               = 3
  C    F
  0    1    l      [m]          l [m]                             l [m]
  1    2    rhoinv [1/m]        volt [MV]                          kick
  2    3    e1                  ex [MV/m]                           .
  3    4    e2                  ey [MV/m]                           .
  4    5    h1                  freq [MHz]                          .
  5    6    h2                  lag [2 Pi]                          .
  6    7    tilt                tilt                               kick
  7    8    ks                  betrf                       C:7-42 F:8-43: rm
  8    9    hgap [m]            pg {MW]                 C:43-258 F:44-259: tm
  9   10    fint [Tm]           shunt [MOhm/m]
  10   11    angle = K_0*l       tfill [micro sec]
  11   12    lrad                harmon
  12   13    k0 or k0*l (l=0)    xsize (coll.) or xma (beam-beam) or x (mon.)
  13   14    k0s or k0s*l        ysize (coll.) or yma (beam_beam) or y (mon.)
  14   15    k1 or k1*l (l=0)    sigx
  15   16    k1s or k1s*l        sigy
  16   17    k2 or k2*l          charge
  17   18    k2s  etc.           npart (# particles in opposite beam)

  int array: as d.p. array, containing expression flag:
  ex_flag = 1   value
  ex_flag > 1   expression

  name array:as d.p. array, pointers to parameter names if ex_flag > 0
*/

/*
  parameter definition:
  int  array:
  1    exflag            1 if value, > 1 if expression
  d.p. array:
  1    value             (always)
  char array:
  string            expression as read if exflag > 1
*/

struct object
{
  char key[48];          /* d.b. key */
/* The order of the first 11 variables below is FIXED */
  int ma_time,            /* start of control part =
                             major time at creation or last modification */
      mi_time,            /* minor time at creation or last modification */
      l_int,              /* length of integer array */
      l_dble,             /* length of double array */
      l_char,             /* length of string */
      l_obj,              /* length of object pointer array */
      c_int,              /* occupation of integer array */
      c_dble,             /* occupation of double array */
      c_char,             /* occupation of string */
      c_obj;              /* occupation of object and names pointer array */

  char par_name[24],      /* parent name */
       base_name[24],     /* basic type name (e.g. QUADRUPOLE, DRIFT,..) */
       obj_type[24];      /* object type such as ELEMENT, TWISS_SUMMARY etc. */

  int*    a_int;          /* integer array */
  double* a_dble;         /* d.p. array */
  char*   a_char;         /* string */
  struct object* parent;  /* pointer to parent object */
  struct object** p_obj;  /* object pointer array */
  char**  names;          /* name pointers into a_char */
};

/* already defined as 42 in fulll.h */
/*#define FIELD_MAX 40*/        /* field error array length */
#define KEY_LENGTH 100      /* from DOOM */
#define MM_KEEP 2           /* no. of element name starts to keep */
#define MULTI_MAX 24        /* element array length for multipoles */
#define NT34 5              /* no. of element types in special fort.34 */
#define LINES_MAX 3         /* structure output line max. names */
#define SEQ_DUMP_LEVEL 0    /* chooses amount of dumped output */

static void add_c6t_drifts(void);
static void add_split_list(struct c6t_element*);
static void add_to_ellist(struct c6t_element*);
static void app_factor(double, double*, int);
// static void arr_print(double*, int); // not used
static void assign_att(void);
static void att_aperture(struct c6t_element*);
static void att_beambeam(struct c6t_element*);
static void att_colli(struct c6t_element*);
static void att_decapole(struct c6t_element*);
static void att_drift(struct c6t_element*);
static int f34_values(struct c6t_element*, int*, double*);
static void att_hkicker(struct c6t_element*);
static void att_kicker(struct c6t_element*);
static void att_lcavity(struct c6t_element*);
static void att_marker(struct c6t_element*);
static void att_matrix(struct c6t_element*);
static void att_multipole(struct c6t_element*);
static void att_octupole(struct c6t_element*);
static void att_quadrupole(struct c6t_element*);
static void att_rbend(struct c6t_element*);
static void att_rfcavity(struct c6t_element*);
static void att_crabcavity(struct c6t_element*);
static void att_dipedge(struct c6t_element*);
static void att_solenoid(struct c6t_element*);
static void att_hacdipole(struct c6t_element*);
static void att_vacdipole(struct c6t_element*);
static void att_sbend(struct c6t_element*);
static void att_sextupole(struct c6t_element*);
static void att_vkicker(struct c6t_element*);
static void att_rfmultipole(struct c6t_element*);
static void att_xrotation(struct c6t_element*);
static void att_yrotation(struct c6t_element*);
static void att_srotation(struct c6t_element*);
static void att_sixmarker(struct c6t_element* el);
static void att_wire(struct c6t_element* el);
static void att_undefined(struct c6t_element*);
static void clean_c6t_element(struct c6t_element*);
static struct c6t_element* create_aperture(const char* ,const char* ,double, double, double, double, double, double, double,
             struct double_array*);
static void concat_drifts(void);
static void conv_elem(void);
static void c6t_finish(void);
static void c6t_init(void);
static struct c6t_element* convert_madx_to_c6t(struct node*,  int ncombined);
// static void dump_c6t_element(struct c6t_element*); // not used
// static void dump_c6t_sequ(int); // not used
// static void dump_types(int); // not used
static void equiv_elem(void);
static struct block* get_block_equiv(struct block*);
static void get_args(struct in_cmd*);
static void get_error_refs(struct c6t_element*);
static int get_flag(struct c6t_element*, struct type_info*);
// static struct c6t_element* get_from_ellist(char*, char*); // not used
static void get_multi_refs(void);
static int get_next_name(char*, char);
// static void gnu_file(struct c6t_element*); // not used
static void grow_ellist(struct c6t_el_list*);
static int ident_el(struct c6t_element*, struct c6t_element*);
static int ident_zero(struct c6t_element*);
static int in_keep_list(struct c6t_element*);
static void invert_normal(int, double*);
// static void invert_skew(int, double*); // not used
static void link_behind(struct c6t_element*, struct c6t_element*);
static void link_c6t_in_front(struct c6t_element*, struct c6t_element*);
static struct c6t_element* make_c6t_element(struct node*, int ncombined);
static struct object* make_obj(const char*, int, int, int, int);
static void make_multipole(struct c6t_element*);
static void mod_errors(void);
static void mod_lcavity(struct c6t_element*);
static void mod_multipole(struct c6t_element*);
static void mod_octupole(struct c6t_element*);
static void mod_quadrupole(struct c6t_element*);
static void mod_rbend(struct c6t_element*);
static void mod_rfcavity(struct c6t_element*);
static void mod_crabcavity(struct c6t_element*);
// static void mod_dipedge(struct c6t_element*); // not defined
// static void mod_solenoid(struct c6t_element*); // not defined
// static void mod_hacdipole(struct c6t_element*); // not defined
// static void mod_vacdipole(struct c6t_element*); // not defined
static void mod_sextupole(struct c6t_element*);
static void multi_loop(void);
static struct c6t_element* new_c6t_element(int, const char*, const char*);
static struct block* new_block(void);
static void post_multipoles(void);
static double power_of(double, int);
static void pre_multipole(struct c6t_element*);
static void pro_elem(struct node*, int ncombined);
static void process_c6t(void);
static void read_sequ(void);
static void remove_from_ellist(struct c6t_element*);
static void replace_c6t(struct c6t_element*, struct c6t_element*);
static void split(void);
static void split_kicker(struct c6t_element*);
static void split_other(struct c6t_element*);
static void split_special(struct c6t_element*);
static void supp_elem(void);
static void supp_small_comp(struct c6t_element*);
static void treat_split(struct c6t_element*);
static void yank(struct c6t_element*);
static void write_all_el(void);
static void write_blocks(void);
static void write_c6t_element(struct c6t_element*);
static void write_f16_errors(void);
static void write_f8_errors(void);
static void write_f3_aper(void);
static void write_f3_aux(void);
static void write_f3_matrix(void);
static void write_f3_wire(void);
static void write_f3_entry(const char*, struct c6t_element*);
static void write_f3_mult(struct c6t_element*);
static void write_f34_special(void);
static void write_struct(void);
static void setup_output_string(void);
static void write_f3_rfmultipoles(struct c6t_element*);
static int my_table_row(struct table*, char*);
static void write_f3_sixmarker(struct element* el);
/* routines used from makethin.c */
static struct li_list types;

static struct block   *first_block; //, *last_block; not used
static struct block*   prev_block;
static struct block*   current_block = NULL;

static int virgin_c6t = 1;

static struct c6t_element *first_in_sequ, *last_in_sequ_org; // *last_in_sequ, // not used
static struct c6t_element* prev_element;
static struct c6t_element* current_element = NULL;
// static struct c6t_element* debug_element = NULL; // not used
static struct c6t_el_list* split_list = NULL;
static struct aper_struct tag_aperture;

static struct object *p_err_zero;  /* pointer to error object with all zeroes */

static int last_row = 0;

static char el_info[][60] = /* see type_info definition */
/*           l=0 l>0,normal l>0,skew ->drift make_k*l split */
{"aperture     2       2       2       0       0       0",
 "beambeam     2       2       2       0       0       0",
 "beamint      0       1       1       1       0       0",
 "collimator   2       1       1       0       0       0",
 "drift        0       1       1       0       0       0",
 "decapole     2       2       2       0       1       2",
 "ecollimator  2       1       1       0       0       0",
 "elseparator  0       1       1       1       0       0",
 "hkicker      5       5       5       1       0       3",
 "hmonitor     0       1       1       1       0       0",
 "instrument   0       1       1       1       0       0",
 "placeholder  0       1       1       1       0       0",
 "kicker       6       6       6       1       0       3",
 "tkicker      6       6       6       1       0       3",
 "lcavity      3       3       3       0       0       2",
 "marker       4       0       0       0       0       0",
 "matrix       2       2       2       0       0       0",
 "monitor      0       1       1       1       0       0",
 "multipole    2       2       2       0       0       0",
 "octupole     2       2       2       0       1       2",
 "quadrupole   2       1       2       0       1       1",
 "rbend        2       1       1       0       1       1",
 "rcollimator  2       1       1       0       0       0",
 "rfcavity     3       3       3       0       0       2",
 "sbend        2       1       1       0       1       1",
 "sextupole    2       2       2       0       1       2",
 "vkicker      5       5       5       1       0       3",
 "vmonitor     0       1       1       1       0       0",
 "crabcavity   3       3       3       0       0       2",
 "dipedge      2       2       2       0       0       0",
 "solenoid     2       2       2       0       0       0",
 "hacdipole    3       3       3       0       0       2",
 "vacdipole    3       3       3       0       0       2",
 "rfmultipole  2       0       2       0       0       2",
 "rfdipole     2       0       2       0       0       2",
 "rfquadrupole 2       0       2       0       0       2",
 "rfsextupole  2       0       2       0       0       2",
 "rfoctupole   2       0       2       0       0       2",
 "xrotation    2       2       2       0       0       0",
 "yrotation    2       2       2       0       0       0",
 "srotation    2       2       2       0       0       0",
 "sixmarker    2       2       2       0       0       2",
 "wire         2       2       2       0       0       2"
};

/* no. of valid element types */
enum { N_TYPES = sizeof el_info / sizeof *el_info };

static struct type_info* t_info[N_TYPES];

static char mpole_names[][16] = {"dipole", "quadrupole", "sextupole",
                          "octupole", "decapole", "multipole"};
static char acro_list[20];   /* list for name starts */
static int acro_cnt[20];    /* counters for name starts */
static char tmp_name[KEY_LENGTH];
static char name_format[70];
static char name_format_short[6];
static char name_format_error[62];
static char name_format_aper[61];
static char name_format_3[40];
static char name_format_4[40];
static char name_format_5[40];
static char name_format_6[60];
static int general_rf_req = 50299 ;
//static char name_format[80]; /*This is used by fprint to determin the length of the names"*/

static int
  block_count = 0,     /* current block count for naming */
  elem_cnt = 0,        /* element count */
  acro_occ = 0,        /* acro list occupation */
  align_cnt = 0,       /* element with align errors count */
  field_cnt = 0,       /* element with field errors count */
  f3_cnt = 0,          /* f3 write flag */
  f3aux_cnt = 0,       /* f3aux write flag */
  f8_cnt = 0,          /* f8 write count */
  f16_cnt = 0,         /* f16 write count */
  f34_cnt = 0,         /* f34 write count */
  special_flag = 1,    /* produce special output file from twiss */
  cavall_flag = 0,     /* if 0 dump all cavities into first */
  markall_flag = 0,    /* if 1 dump all markers into first */
  long_names_flag = 0,    /* if 1 dump all markers into first */
  multicol_flag = 0,   /* if 1 dump multi-column STRUCTURE block */
  aperture_flag = 0,   /* if 1 insert apertures into structure */
//  radius_flag = 0, // not used    /* change the default reference radius */
  split_flag = 0,      /* if 1 keep zero multipoles after split */
  mangle_flag = 0,     /* if 1 truncate to 14 chars and mangle names */
  mult_auto_off = 1,   /* if 1 code does not process zero value
                          multipoles;
                          if 0 process up to order max_mult_ord */
  max_mult_ord = 20,   /* Process up to this order for mult_auto_off = 0 */
  multi_type = -1,     /* is set to multipole type if any found */
  cavity_count = 0;    /* count cavities in output */

static double
  sequ_length,         /* length of  sequence */
  sequ_start,
  sequ_end,
  total_voltage = 0,
  harmon = 0,
//  freq = 0, // not used
  error_matrix[FIELD_MAX],
  tmp_buff[FIELD_MAX];

static const double ten   = 10;
static const double c1p3 = 1.e3;
static const double eps_6 = 1.e-6;
static const double eps_9 = 1.e-9;
static const double eps_12 = 1.e-12;
static double ref_def = 0.017;
static double six_version = 0;

static FILE *f2, *f3, *f3aux, *f3aper, *f8, *f16, *f34;

// private functions

static void
add_c6t_drifts(void)
{
  int af;
  struct c6t_element *d1;
  double pos = sequ_start, dl, el2;
  char c[24];
  current_element = first_in_sequ;
  while (current_element != NULL)
  {
    el2 = current_element->value[0] / two;
    dl = current_element->position - el2 - pos;
    if (dl + eps_9 < zero)
    {
      printf(
        "+=+=+= c6t fatal - negative drift in front of %s, length %f\n",
        current_element->name, dl);
      exit(1);
    }
    else if (dl > eps_9) // create an appropriate drift
    {
      af = get_next_name(c, 'd');
      d1 = new_c6t_element(1, c, "drift");
      d1->value[0] = dl; d1->flag = 1;
      link_c6t_in_front(d1, current_element);
      d1->position = pos + dl / two;
      if (af != 0)  add_to_ellist(d1);
    }
    pos = current_element->position + el2;
    current_element = current_element->next;
  }
}

static void
add_split_list(struct c6t_element* el)
{
  int i;
  const char *rout_name = "c6t:add_split_list";
  if (split_list == NULL)
  {
    split_list = mycalloc(rout_name,1, sizeof *split_list);
    split_list->elem = mycalloc(rout_name,EL_COUNT, sizeof *split_list->elem);
    split_list->max = EL_COUNT;
  }
  else if (split_list->curr == split_list->max) grow_ellist(split_list);
  for (i = 0; i < split_list->curr; i++) if (split_list->elem[i] == el) return;
  split_list->elem[split_list->curr++] = el;
}

static void
add_to_ellist( /* adds element to correct object list */
  struct c6t_element* p_elem)
{
  int j;
  const char *rout_name = "c6t:add_to_ellist";

#ifdef _call_tree_
  puts("+++++++ add_to_ellist");
#endif
  for (j = 0; j < types.curr; j++)
  {
    if (strcmp(types.member[j]->base_name, p_elem->base_name) == 0)
    {
      if (types.member[j]->curr == types.member[j]->max)
        grow_ellist(types.member[j]);
      types.member[j]->elem[types.member[j]->curr++] = p_elem;
      return;
    }
  }
  /* type list does not exist - create */
  if (types.curr == BASE_TYPES)
  {
    printf("+++ fatal - %s overruns type buffer of %d types\n",
           p_elem->base_name, types.curr);
    exit(1);
  }
  types.member[types.curr] = mycalloc(rout_name,1,sizeof *types.member[0]);
  types.member[types.curr]->elem = mycalloc(rout_name,EL_COUNT, sizeof *types.member[0]->elem);
  types.member[types.curr]->elem[types.member[types.curr]->curr++] = p_elem;
  types.member[types.curr]->max = EL_COUNT;
  strcpy(types.member[types.curr]->base_name, p_elem->base_name);
  types.curr++;
}

static void
app_factor(double fact, double* array, int count)
{
  int i;
  for (i = 0; i < count; i++) array[i] *= fact;
}

#if 0 // not used
static void
arr_print(double array[], int occ)
{
  int i;
  for (i = 0; i < occ; i++)
  {
    printf(" %12.4e", array[i]); if ((i+1)%5 == 0) printf("\n");
  }
  printf("\n");
}
#endif

static void
assign_att(void)
{
  struct c6t_element *el;
  int i, j;

  for (i = 0; i < types.curr; i++)  /* loop over base types */
  {
    for (j = 0; j < types.member[i]->curr; j++) /* loop over el. in type */
    {
      el = types.member[i]->elem[j];
      if (el->flag > 0 && el->equiv == el)  /* all others ignored */
      {
        if      (strcmp(el->base_name, "aperture") == 0) att_aperture(el);
        else if (strcmp(el->base_name, "beambeam") == 0) att_beambeam(el);
        else if (strcmp(el->base_name, "collimator") == 0) att_colli(el);
        else if (strcmp(el->base_name, "decapole") == 0) att_decapole(el);
        else if (strcmp(el->base_name, "drift") == 0) att_drift(el);
        else if (strcmp(el->base_name, "ecollimator") == 0) att_colli(el);
        else if (strcmp(el->base_name, "hkicker") == 0) att_hkicker(el);
        else if (strcmp(el->base_name, "kicker") == 0) att_kicker(el);
        else if (strcmp(el->base_name, "tkicker") == 0) att_kicker(el);
        else if (strcmp(el->base_name, "lcavity") == 0) att_lcavity(el);
        else if (strcmp(el->base_name, "marker") == 0) att_marker(el);
        else if (strcmp(el->base_name, "matrix") == 0) att_matrix(el);
        else if (strcmp(el->base_name, "multipole") == 0) att_multipole(el);
        else if (strcmp(el->base_name, "octupole") == 0) att_octupole(el);
        else if (strcmp(el->base_name, "quadrupole") == 0) att_quadrupole(el);
        else if (strcmp(el->base_name, "rbend") == 0) att_rbend(el);
        else if (strcmp(el->base_name, "rcollimator") == 0) att_colli(el);
        else if (strcmp(el->base_name, "rfcavity") == 0) att_rfcavity(el);
        else if (strcmp(el->base_name, "crabcavity") == 0) att_crabcavity(el);
        else if (strcmp(el->base_name, "dipedge") == 0) att_dipedge(el);
        else if (strcmp(el->base_name, "solenoid") == 0) att_solenoid(el);
        else if (strcmp(el->base_name, "hacdipole") == 0) att_hacdipole(el);
        else if (strcmp(el->base_name, "vacdipole") == 0) att_vacdipole(el);
        else if (strcmp(el->base_name, "sbend") == 0) att_sbend(el);
        else if (strcmp(el->base_name, "sextupole") == 0) att_sextupole(el);
        else if (strcmp(el->base_name, "vkicker") == 0) att_vkicker(el);
        else if (strcmp(el->base_name, "rfmultipole") == 0) att_rfmultipole(el);
        else if (strcmp(el->base_name, "xrotation") == 0) att_xrotation(el);
        else if (strcmp(el->base_name, "yrotation") == 0) att_yrotation(el);
        else if (strcmp(el->base_name, "srotation") == 0) att_srotation(el);
        else if (strcmp(el->base_name, "sixmarker") == 0) att_sixmarker(el);
        else if (strcmp(el->base_name, "wire") == 0) att_wire(el);
        else att_undefined(el);
      }
    }
  }
}

static void
att_aperture(struct c6t_element* el)
{
  el->out_1 = 0;
  el->out_2 = 0.0;
  el->out_3 = 0.0;
  el->out_4 = 0.0;
}

static void
att_beambeam(struct c6t_element* el)
{

  double beamx=zero,beamy=zero;
  if (double_from_table_row("twiss","x",&(el->twtab_row),&beamx) != 0 ||
      double_from_table_row("twiss","y",&(el->twtab_row),&beamy) != 0)
  {
    warning("c6t: beambeam element not found in twiss table","");
  }
  el->out_1 = 20;
  el->out_2 = - c1p3*(el->value[12] - beamx);
  el->out_3 = - c1p3*(el->value[13] - beamy);
  el->out_4 = el->value[16];
  el->out_5 = pow(c1p3*el->value[14], 2);
  el->out_6 = pow(c1p3*el->value[15], 2);
  el->out_7 = 0;

}

static void
att_colli(struct c6t_element* el)
/* collimator, ecollimator, rcollimator - make drift, do not concatenate */
{
  el->out_1 = 0; el->out_4 = el->value[0];
}

static void
att_decapole(struct c6t_element* el)
{
  if (el->value[20] != zero)
  {
    el->out_1 = 5; el->out_2 = -el->value[20]/24;
  }
  else if (el->value[21] != zero)
  {
    el->out_1 = -5; el->out_2 = el->value[21]/24;
  }
  else el->out_1 = 0;
}

static void
att_drift(struct c6t_element* el)
{
  el->out_4 = el->value[0];
}

static void
att_hkicker(struct c6t_element* el)
{
  el->out_1 = 1; el->out_2 = el->value[12];
}

static void
att_kicker(struct c6t_element* el)
{
  (void)el;
}

static void
att_lcavity(struct c6t_element* el)
{
  double lag = -el->value[5];
  el->out_1 = 12;
  el->out_2 = cavall_flag == 0 ? total_voltage : el->value[1];
  el->out_3 = 0; /* ??? harmon = p_beam->a_dble[41] / el->value[11]; */
  printf("harmon: %e\n", harmon);
  if (lag < -0.5) lag +=1.;
  else if (lag > 0.5) lag -=1.;
  el->out_4 = 360. * lag;
}

static void
att_marker(struct c6t_element* el)
{
  (void)el;
}

static void
att_matrix(struct c6t_element* el)
{
  el->out_1 = 22;
  el->out_2 = 0;
  el->out_3 = 0;
  el->out_4 = el->value[0];
}

static void
att_multipole(struct c6t_element* el)
{
  el->out_1 = 11;
  if (el->nc_pos == 0)
  {
    el->out_2 = el->out_3 = 1;
  }
  else
  {
    el->out_3 = el->rad_length;
    if (el->nc_pos == 12)
    {
      el->out_2 = -el->value[12]; el->out_4 = -1;
    }
    else if (el->nc_pos == 13)
    {
      el->out_2 = el->value[13]; el->out_4 = -2;
    }
  }
}

static void
att_octupole(struct c6t_element* el)
{
  if (el->value[18] != zero)
  {
    el->out_1 = 4; el->out_2 = -el->value[18]/6;
  }
  else if (el->value[19] != zero)
  {
    el->out_1 = -4; el->out_2 = el->value[19]/6;
  }
  else el->out_1 = 0;
}

static void
att_quadrupole(struct c6t_element* el)
{
  el->out_4 = el->value[0];
  if (el->value[14] != zero)
  {
    el->out_1 = 2;
    if (el->value[0] == zero) el->out_2 = -el->value[14];
    else                      el->out_3 = -el->value[14];
  }
  else if (el->value[15] != zero)
  {
    el->out_1 = -2; el->out_2 = el->value[15];
  }
  else el->out_1 = 0;
}

static void
att_rbend(struct c6t_element* el)
{
  el->out_4 = el->value[0];
  if (el->value[12] != zero)
  {
    el->out_2 = -el->value[1];
    if (el->value[14] == zero)  el->out_1 = 1;
    else
    {
      el->out_1 = 6;
      el->out_3 = -el->value[14];
    }
  }
  else if (el->value[13] != zero)
  {
    el->out_2 = el->value[1];
    if (el->value[15] == zero)  el->out_1 = 4;
    else
    {
      el->out_1 = 4;
      el->out_3 = el->value[15];
    }
  }
  else el->out_1 = 0;
}

static void
att_rfcavity(struct c6t_element* el)
{
  double lag = 0.5 - el->value[5];
  el->out_1 = 12;
  if (cavall_flag == 0)
  {
    el->out_2 = total_voltage;
    strcpy(el->name, "CAV");
  }
  else el->out_2 = el->value[1];
  el->out_3 = harmon = el->value[11];
  if (lag < -0.5) lag +=1.;
  else if (lag > 0.5) lag -=1.;
  el->out_4 = 360. * lag;
}

static void
att_crabcavity(struct c6t_element* el)
{
  double lag = -el->value[5];
  double tilt = el->value[12];
  if (fabs(tilt - M_PI/2)<eps_9)
    el->out_1 = -23;
  else
    el->out_1 = 23;
  /* AL: Discussions with RdM lead to the conclusions that a CrabCavity shouldn't be considered as a regular accelerating cavity */
  /* if (cavall_flag == 0) */
  /* { */
  /*   el->out_2 = total_voltage; */
  /*   strcpy(el->name, "CAV"); */
  /* } */
  /* else */
  el->out_2 = el->value[1];
  el->out_3 = el->value[4]; // freq = // not used
  if (lag < -0.5) lag +=1.;
  else if (lag > 0.5) lag -=1.;
  el->out_4 = 2 * M_PI * lag;
}

static void
att_dipedge(struct c6t_element* el)
{
  double corr;
  corr = 2*el->value[1]*el->value[8]*el->value[9];
  if (el->value[1] != zero && (el->value[2] != zero || corr != zero))
  {
    el->out_1 = 24;
    el->out_2 =  el->value[1]*tan(el->value[2]);
    el->out_3 = -el->value[1]*tan(el->value[2]-corr/cos(el->value[2])*
                                  (one+sin(el->value[2])*sin(el->value[2])));
  }
  else
  {
    el->out_1 = 0;
    el->out_2 = 0;
    el->out_3 = 0;
  }
  el->out_4 = 0;
}

static void
att_xrotation(struct c6t_element* el)
{
  el->out_1 = 43;
  el->out_2 = el->value[1] ; 
}
static void
att_yrotation(struct c6t_element* el)
{
  el->out_1 = 44;
  el->out_2 = el->value[1] ;
}
static void
att_srotation(struct c6t_element* el)
{
  el->out_1 = 45;
  el->out_2 = el->value[1] ;
}

static void
att_sixmarker(struct c6t_element* el)
{
  el->out_1 = el->value[1];
  el->out_2 = el->value[2];
  el->out_3 = el->value[3];
  el->out_4 = el->value[4];
  el->out_5 = el->value[5];
  el->out_6 = el->value[6];
  el->out_7 = el->value[7];

}
static void 
att_wire(struct c6t_element* el)
{
  el->out_1 = 15;

}
static void
att_solenoid(struct c6t_element* el)
{
  el->out_1 = 25;
  el->out_2 = el->value[2];
  el->out_3 = el->value[3];
  el->out_4 = el->value[0];
}

static void
att_hacdipole(struct c6t_element* el)
{
  el->out_1 = 16;
  el->out_2 = el->value[2];
  el->out_3 = el->value[3];
  el->out_4 = el->value[4];
}

static void
att_vacdipole(struct c6t_element* el)
{
  el->out_1 = -16;
  el->out_2 = el->value[2];
  el->out_3 = el->value[3];
  el->out_4 = el->value[4];
}

static void
att_sbend(struct c6t_element* el)
{
  el->out_4 = el->value[0];
  if (el->value[12] != zero)
  {
    el->out_2 = -el->value[1];
    if (el->value[14] == zero)  el->out_1 = 3;
    else
    {
      el->out_1 = 6;
      el->out_3 = -el->value[14];
    }
  }
  else if (el->value[13] != zero)
  {
    el->out_1 = 5;
    el->out_2 = el->value[1];
    el->out_3 = el->value[15];
  }
  else el->out_1 = 0;
}

static void
att_sextupole(struct c6t_element* el)
{
  if (el->value[16] != zero)
  {
    el->out_1 = 3; el->out_2 = -el->value[16]/two;
  }
  else if (el->value[17] != zero)
  {
    el->out_1 = -3; el->out_2 = el->value[17]/two;
  }
  else el->out_1 = 0;
}

static void
att_vkicker(struct c6t_element* el)
{
  el->out_1 = -1; el->out_2 = el->value[13];
}

static void
att_undefined(struct c6t_element* el)
{
  el->out_4 = el->value[0];
}

static void
att_rfmultipole(struct c6t_element* el)
{
  (void)el;
}

static void
block_it(void)
{
  struct c6t_element* el;
  const char *rout_name = "c6t:block_it";

  current_element = first_in_sequ;
  while ((el = current_element) != NULL)
  {
    current_block = new_block();
    current_block->previous = prev_block;
    current_block->next = NULL;
    if (prev_block == NULL) first_block = current_block;
    else                    prev_block->next = current_block;
    current_block->elements = mycalloc(rout_name,1,sizeof *current_block->elements);
    current_block->elements->elem = mycalloc(rout_name,EL_COUNT, sizeof *current_block->elements->elem);
    current_block->elements->max = EL_COUNT;
    current_block->first = el;
    current_block->length = el->equiv->value[0];
    current_block->elements->elem[0] = el;
    current_block->elements->curr = 1;
    if (el->flag < 2)
    {
      while (el->next != NULL && el->next->flag < 2)
      {
        el = el->next;
        current_block->length += el->equiv->value[0];
        if (current_block->elements->curr == current_block->elements->max)
          grow_ellist(current_block->elements);
        current_block->elements->elem[current_block->elements->curr++]
          = el;
      }
      current_element = el;
    }
    current_block->last = current_element;
    if (current_block->first == current_block->last &&
        current_block->last->flag >= 2)  current_block->flag = 0;
    else current_block->flag = 1;
    current_block->equiv = get_block_equiv(current_block);
    current_element = current_element->next;
    prev_block = current_block;
  }
  // last_block = current_block; // not used
}

static void
concat_drifts(void)
{
  struct c6t_element *d1, *temp, *nk;
  int flag, cnt;
  double suml, pos;
  char c[24];
  current_element = first_in_sequ;
  while (current_element != NULL)
  {
    cnt = 0;
    suml = current_element->value[0];
    pos = current_element->position - suml / two;
    if (strcmp(current_element->base_name, "drift") == 0)
    {
      temp = current_element->next;
      while (temp != NULL && strcmp(temp->base_name, "drift") == 0)
      {
        suml += temp->value[0]; cnt++;
        temp = temp->next;
      }
    }
    if (cnt > 0) /* actually concatenated something */
    {
      flag = get_next_name(c, 'd');
      d1 = new_c6t_element(1, c, "drift");  d1->flag = 1;
      d1->value[0] = suml; d1->position = pos + suml / two;
      if (flag != 0) add_to_ellist(d1);
      temp = current_element->next;
      while (temp != NULL && strcmp(temp->base_name, "drift") == 0)
      {
        nk = temp->next;
        yank(temp);
        temp = nk;
      }
      if (current_element == first_in_sequ) first_in_sequ = d1;
      replace_c6t(current_element, d1); current_element = d1;
    }
    current_element = current_element->next;
  }
}

static void
conv_elem(void)
{
  int i, j, nup;
  struct type_info* type = NULL;
  struct c6t_element* el;

  for (i = 0; i < types.curr; i++)  /* loop over base types */
  {
    for (j = 0; j < N_TYPES; j++)
    {

      if (strcmp(types.member[i]->base_name, t_info[j]->name) == 0)
      {

        type = t_info[j]; break;
      }
    }
    if (type == NULL)
    {
      printf("+=+=+= c6t fatal - type %s not defined\n",
             types.member[i]->base_name);
      exit(1);
    }
    nup = types.member[i]->curr;
    for (j = 0; j < nup; j++) /* loop over el. in type */
    {
      el = types.member[i]->elem[j];
      if (type->flag_4 > 1)
        printf("+++ warning - treated as drift: %s\n", el->name);
      el->flag = get_flag(el, type);
      if (el->flag > 0)  /* all others ignored */
      {
        if (el->value[0] < eps_9)
        {
          el->value[0] = zero;
          if (el->flag == 1)  el->flag = 0;
        }
        if (el->flag > 0)
        {
          el->c_drift = type->flag_4;
          el->force = type->flag_5;
          el->split = type->flag_6;
          if (el->split > 0) add_split_list(el);
        }
      }
    }
  }
}

static void
clean_c6t_element(struct c6t_element* cleanme)
{
  int i;
  for(i=0; i<cleanme->n_values; i++) { cleanme->value[i]=0; }
}

static struct c6t_element*
create_aperture(const char* name, const char* type, double ap1, double ap2, double ap3, double ap4,
    double offx, double offy, double tilt, struct double_array* p_al_err)
{
  struct c6t_element* aper_element;
  aper_element = new_c6t_element(9,name,"aperture");
  clean_c6t_element(aper_element);
  strcpy(aper_element->org_name,name);

  // type of element is coded by integer...
  if      (strcmp(type,"CR")==0) aper_element->value[1] = 1;
  else if (strcmp(type,"RE")==0) aper_element->value[1] = 2;
  else if (strcmp(type,"EL")==0) aper_element->value[1] = 3;
  else if (strcmp(type,"RC")==0) aper_element->value[1] = 4;
  else if (strcmp(type,"RL")==0) aper_element->value[1] = 5;
  else if (strcmp(type,"RT")==0) aper_element->value[1] = 6;
  else if (strcmp(type,"OC")==0) aper_element->value[1] = 7;
  else aper_element->value[1] = 0;

  aper_element->value[0] = 0.0; // zero length ?
  aper_element->value[2] = ap1 * 1e3; // sixtrack units are mm
  aper_element->value[3] = ap2 * 1e3;
  aper_element->value[4] = ap3 * 1e3;
  aper_element->value[5] = ap4 * 1e3;
  aper_element->value[6] = offx * 1e3;
  aper_element->value[7] = offy * 1e3;
  aper_element->value[8] = tilt / M_PI * 180 ; // sixtrack units are degrees

  if (aper_element->value[1] == 7) // Octagon
  {
    aper_element->value[4] = ap3 / M_PI * 180; // sixtrack units are degrees
    aper_element->value[5] = ap4 / M_PI * 180;
  }

  aper_element->keep_in=1;
  /* alignment errors of aperture are to be copied to alignment errors of element */
  if (p_al_err && p_al_err->curr>11)
  {
    align_cnt++;
    aper_element->na_err = p_al_err->curr;
    aper_element->p_al_err = make_obj("ALDUM",0,ALIGN_MAX,0,0);
    aper_element->p_al_err->c_dble = p_al_err->curr;
    aper_element->p_al_err->a_dble[0] = p_al_err->a[10];
    aper_element->p_al_err->a_dble[1] = p_al_err->a[11];
  }
  return aper_element;
}

static struct c6t_element*
convert_madx_to_c6t(struct node* p, int ncombined)
{
  struct command_parameter *kn_param = NULL,*ks_param = NULL,*aper_param = NULL;
  struct command_parameter *pn_param = NULL,*ps_param = NULL;
  struct c6t_element* c6t_elem = NULL;
  char t_name[255]; /* minimum NAME_MANGLER_BASE + NAME_MANGLER_SUFFIX = 16 characters + 1 */
  int i,j;
  char* cp;
  int index=-1;
  if (mangle_flag)
    NameMangler_mangle(p->name, t_name);
  else
    strcpy(t_name, p->name);

  if ((cp = strchr(t_name, ':')) != NULL) *cp = '\0';
  if ((strcmp(p->base_name,"rbend") == 0)      ||
      (strcmp(p->base_name,"sbend") == 0)      ||
      (strcmp(p->base_name,"quadrupole") == 0) ||
      (strcmp(p->base_name,"sextupole") == 0)  ||
      (strcmp(p->base_name,"octupole") == 0)   ||
      (strcmp(p->base_name,"vkicker") == 0)    ||
      (strcmp(p->base_name,"hkicker") == 0)    ||
      (strcmp(p->base_name,"tkicker") == 0)    ||
      (strcmp(p->base_name,"kicker") == 0))
  {
    c6t_elem = new_c6t_element(20,t_name,p->base_name);
    clean_c6t_element(c6t_elem);
    strcpy(c6t_elem->org_name,t_name);

    c6t_elem->value[0] = el_par_value_recurse("l",p->p_elem);
    c6t_elem->value[1] = el_par_value("rhoinv",p->p_elem);  /* not intrinsic */
    c6t_elem->value[2] = el_par_value_recurse("e1",p->p_elem);
    c6t_elem->value[3] = el_par_value_recurse("e2",p->p_elem);
    c6t_elem->value[4] = el_par_value_recurse("h1",p->p_elem);
    c6t_elem->value[5] = el_par_value_recurse("h2",p->p_elem);
    c6t_elem->value[6] = el_par_value_recurse("tilt",p->p_elem);
    c6t_elem->value[7] = el_par_value_recurse("k0s",p->p_elem);
    c6t_elem->value[8] = el_par_value_recurse("hgap",p->p_elem);
    c6t_elem->value[9] = el_par_value_recurse("fint",p->p_elem);
    c6t_elem->value[10] = el_par_value_recurse("angle",p->p_elem);
    c6t_elem->value[11] = el_par_value_recurse("lrad",p->p_elem);
    c6t_elem->value[12] = el_par_value_recurse("k0",p->p_elem);
    if (c6t_elem->value[12] == zero && c6t_elem->value[0] > zero)
    {
      c6t_elem->value[12] = c6t_elem->value[10]/c6t_elem->value[0];
    }
    c6t_elem->value[13] = el_par_value_recurse("k0s",p->p_elem);
    c6t_elem->value[14] = el_par_value_recurse("k1",p->p_elem);
    c6t_elem->value[15] = el_par_value_recurse("k1s",p->p_elem);
    c6t_elem->value[16] = el_par_value_recurse("k2",p->p_elem);
    c6t_elem->value[17] = el_par_value_recurse("k2s",p->p_elem);
    c6t_elem->value[18] = el_par_value_recurse("k3",p->p_elem);
    c6t_elem->value[19] = el_par_value_recurse("k3s",p->p_elem);
  }
  else if ((strcmp(p->base_name,"multipole") == 0))
  {
    int maxkn=0, maxks=0;
    /*      if ((kn_param = return_param_recurse("knl",p->p_elem))) maxkn=kn_param->double_array->curr; */
    /*      if ((ks_param = return_param_recurse("ksl",p->p_elem))) maxks=ks_param->double_array->curr; */
    if ((index = name_list_pos("knl",p->p_elem->def->par_names))>-1)
    {
      kn_param = p->p_elem->def->par->parameters[index];
      maxkn=kn_param->double_array->curr;
    }
    if ((index = name_list_pos("ksl",p->p_elem->def->par_names))>-1)
    {
      ks_param = p->p_elem->def->par->parameters[index];
      maxks=ks_param->double_array->curr;
    }
    if (maxkn > maxks) {j=maxkn;} else {j=maxks;}
    i=j*2+12+1;
    c6t_elem = new_c6t_element(i,t_name,p->base_name);
    clean_c6t_element(c6t_elem);
    strcpy(c6t_elem->org_name,t_name);
    c6t_elem->value[0] = el_par_value_recurse("l",p->p_elem);
    c6t_elem->value[6] = el_par_value_recurse("tilt",p->p_elem);
    c6t_elem->value[10] = el_par_value_recurse("angle",p->p_elem);
    c6t_elem->value[11] = el_par_value_recurse("lrad",p->p_elem);
    for (i=0; i<j; i++)
    {
      if (i<maxkn) c6t_elem->value[i*2+12] = kn_param->double_array->a[i];
      if (i<maxks) c6t_elem->value[i*2+13] = ks_param->double_array->a[i];
    }
  }
  else if ((strcmp(p->base_name,"matrix") == 0))
  {
    c6t_elem = new_c6t_element(43,t_name,p->base_name);
    clean_c6t_element(c6t_elem);
    strcpy(c6t_elem->org_name,t_name);
    c6t_elem->value[0] = el_par_value_recurse("l",p->p_elem);
    c6t_elem->value[1] = el_par_value_recurse("kick1",p->p_elem);
    c6t_elem->value[2] = el_par_value_recurse("kick2",p->p_elem);
    c6t_elem->value[3] = el_par_value_recurse("kick3",p->p_elem);
    c6t_elem->value[4] = el_par_value_recurse("kick4",p->p_elem);
    c6t_elem->value[5] = el_par_value_recurse("kick5",p->p_elem);
    c6t_elem->value[6] = el_par_value_recurse("kick6",p->p_elem);
    c6t_elem->value[7] = el_par_value_recurse("rm11",p->p_elem);
    c6t_elem->value[8] = el_par_value_recurse("rm12",p->p_elem);
    c6t_elem->value[9] = el_par_value_recurse("rm13",p->p_elem);
    c6t_elem->value[10] = el_par_value_recurse("rm14",p->p_elem);
    c6t_elem->value[11] = el_par_value_recurse("rm15",p->p_elem);
    c6t_elem->value[12] = el_par_value_recurse("rm16",p->p_elem);
    c6t_elem->value[13] = el_par_value_recurse("rm21",p->p_elem);
    c6t_elem->value[14] = el_par_value_recurse("rm22",p->p_elem);
    c6t_elem->value[15] = el_par_value_recurse("rm23",p->p_elem);
    c6t_elem->value[16] = el_par_value_recurse("rm24",p->p_elem);
    c6t_elem->value[17] = el_par_value_recurse("rm25",p->p_elem);
    c6t_elem->value[18] = el_par_value_recurse("rm26",p->p_elem);
    c6t_elem->value[19] = el_par_value_recurse("rm31",p->p_elem);
    c6t_elem->value[20] = el_par_value_recurse("rm32",p->p_elem);
    c6t_elem->value[21] = el_par_value_recurse("rm33",p->p_elem);
    c6t_elem->value[22] = el_par_value_recurse("rm34",p->p_elem);
    c6t_elem->value[23] = el_par_value_recurse("rm35",p->p_elem);
    c6t_elem->value[24] = el_par_value_recurse("rm36",p->p_elem);
    c6t_elem->value[25] = el_par_value_recurse("rm41",p->p_elem);
    c6t_elem->value[26] = el_par_value_recurse("rm42",p->p_elem);
    c6t_elem->value[27] = el_par_value_recurse("rm43",p->p_elem);
    c6t_elem->value[28] = el_par_value_recurse("rm44",p->p_elem);
    c6t_elem->value[29] = el_par_value_recurse("rm45",p->p_elem);
    c6t_elem->value[30] = el_par_value_recurse("rm46",p->p_elem);
    c6t_elem->value[31] = el_par_value_recurse("rm51",p->p_elem);
    c6t_elem->value[32] = el_par_value_recurse("rm52",p->p_elem);
    c6t_elem->value[33] = el_par_value_recurse("rm53",p->p_elem);
    c6t_elem->value[34] = el_par_value_recurse("rm54",p->p_elem);
    c6t_elem->value[35] = el_par_value_recurse("rm55",p->p_elem);
    c6t_elem->value[36] = el_par_value_recurse("rm56",p->p_elem);
    c6t_elem->value[37] = el_par_value_recurse("rm61",p->p_elem);
    c6t_elem->value[38] = el_par_value_recurse("rm62",p->p_elem);
    c6t_elem->value[39] = el_par_value_recurse("rm63",p->p_elem);
    c6t_elem->value[40] = el_par_value_recurse("rm64",p->p_elem);
    c6t_elem->value[41] = el_par_value_recurse("rm65",p->p_elem);
    c6t_elem->value[42] = el_par_value_recurse("rm66",p->p_elem);
  }
  else if ((strcmp(p->base_name,"rfcavity") == 0))
  {
    c6t_elem = new_c6t_element(11,t_name,p->base_name);
    clean_c6t_element(c6t_elem);
    strcpy(c6t_elem->org_name,t_name);
    c6t_elem->value[0] = el_par_value_recurse("l",p->p_elem);
    c6t_elem->value[1] = el_par_value_recurse("volt",p->p_elem);
    c6t_elem->value[4] = el_par_value_recurse("freq",p->p_elem);
    c6t_elem->value[5] = el_par_value_recurse("lag",p->p_elem);
    c6t_elem->value[7] = el_par_value_recurse("betrf",p->p_elem);
    c6t_elem->value[8] = el_par_value_recurse("pg",p->p_elem);
    c6t_elem->value[9] = el_par_value_recurse("shunt",p->p_elem);
    c6t_elem->value[10] = el_par_value_recurse("tfill",p->p_elem);
    c6t_elem->value[11] = el_par_value_recurse("harmon",p->p_elem);
  }
  else if ((strcmp(p->base_name,"crabcavity") == 0))
  {
    c6t_elem = new_c6t_element(12,t_name,p->base_name);
    clean_c6t_element(c6t_elem);
    strcpy(c6t_elem->org_name,t_name);
    c6t_elem->value[0] = el_par_value_recurse("l",p->p_elem);
    c6t_elem->value[1] = el_par_value_recurse("volt",p->p_elem);
    c6t_elem->value[4] = el_par_value_recurse("freq",p->p_elem);
    c6t_elem->value[5] = el_par_value_recurse("lag",p->p_elem);
    c6t_elem->value[7] = el_par_value_recurse("betrf",p->p_elem);
    c6t_elem->value[8] = el_par_value_recurse("pg",p->p_elem);
    c6t_elem->value[9] = el_par_value_recurse("shunt",p->p_elem);
    c6t_elem->value[10] = el_par_value_recurse("tfill",p->p_elem);
    c6t_elem->value[11] = el_par_value_recurse("harmon",p->p_elem);
    c6t_elem->value[12] = el_par_value_recurse("tilt",p->p_elem);
  }
  else if ((strcmp(p->base_name,"dipedge") == 0))
  {
    c6t_elem = new_c6t_element(11,t_name,p->base_name);
    clean_c6t_element(c6t_elem);
    strcpy(c6t_elem->org_name,t_name);
    c6t_elem->value[1] = el_par_value("h",p->p_elem);
    c6t_elem->value[2] = el_par_value_recurse("e1",p->p_elem);
    c6t_elem->value[8] = el_par_value_recurse("hgap",p->p_elem);
    c6t_elem->value[9] = el_par_value_recurse("fint",p->p_elem);
  }
  else if ((strcmp(p->base_name,"solenoid") == 0))
  {
    c6t_elem = new_c6t_element(11,t_name,p->base_name);
    clean_c6t_element(c6t_elem);
    strcpy(c6t_elem->org_name,t_name);
    c6t_elem->value[0] = el_par_value_recurse("l",p->p_elem);
    c6t_elem->value[2] = el_par_value_recurse("ks",p->p_elem);
    c6t_elem->value[3] = el_par_value_recurse("ksi",p->p_elem);
  }
  else if ((strcmp(p->base_name,"hacdipole") == 0))
  {
    c6t_elem = new_c6t_element(11,t_name,p->base_name);
    clean_c6t_element(c6t_elem);
    strcpy(c6t_elem->org_name,t_name);
    c6t_elem->value[2] = el_par_value_recurse("volt",p->p_elem);
    c6t_elem->value[3] = el_par_value_recurse("freq",p->p_elem);
    c6t_elem->value[4] = el_par_value_recurse("lag",p->p_elem);
  }
  else if ((strcmp(p->base_name,"vacdipole") == 0))
  {
    c6t_elem = new_c6t_element(11,t_name,p->base_name);
    clean_c6t_element(c6t_elem);
    strcpy(c6t_elem->org_name,t_name);
    c6t_elem->value[2] = el_par_value_recurse("volt",p->p_elem);
    c6t_elem->value[3] = el_par_value_recurse("freq",p->p_elem);
    c6t_elem->value[4] = el_par_value_recurse("lag",p->p_elem);
  }
  else if ((strcmp(p->base_name,"marker") == 0)   ||
           (strcmp(p->base_name,"instrument") == 0)    ||
           (strcmp(p->base_name,"placeholder") == 0)    ||
           (strcmp(p->base_name,"hmonitor") == 0) ||
           (strcmp(p->base_name,"vmonitor") == 0) ||
           (strcmp(p->base_name,"monitor") == 0))
  {
    c6t_elem = new_c6t_element(0,t_name,p->base_name);
    clean_c6t_element(c6t_elem);
    strcpy(c6t_elem->org_name,t_name);
    c6t_elem->value[0] = el_par_value_recurse("l",p->p_elem);
  }
  else if ((strcmp(p->base_name,"collimator") == 0))
  {
    c6t_elem = new_c6t_element(13,t_name,p->base_name);
    clean_c6t_element(c6t_elem);
    strcpy(c6t_elem->org_name,t_name);
    c6t_elem->value[0] = el_par_value_recurse("l",p->p_elem);
  }
  else if ((strcmp(p->base_name,"ecollimator") == 0) ||
     (strcmp(p->base_name,"rcollimator") == 0))
  {
    c6t_elem = new_c6t_element(13,t_name,p->base_name);
    clean_c6t_element(c6t_elem);
    strcpy(c6t_elem->org_name,t_name);
    c6t_elem->value[0] = el_par_value_recurse("l",p->p_elem);
    c6t_elem->value[12] = el_par_value_recurse("xsize",p->p_elem);
    c6t_elem->value[13] = el_par_value_recurse("ysize",p->p_elem);
  }
  else if((strcmp(p->base_name,"beambeam") == 0  ))
  {
    c6t_elem = new_c6t_element(17,t_name,p->base_name);
    clean_c6t_element(c6t_elem);
    strcpy(c6t_elem->org_name,t_name);
    c6t_elem->value[0] = el_par_value_recurse("l",p->p_elem);
    c6t_elem->value[12] = el_par_value_recurse("xma",p->p_elem);
    c6t_elem->value[13] = el_par_value_recurse("yma",p->p_elem);
    c6t_elem->value[14] = el_par_value_recurse("sigx",p->p_elem);
    c6t_elem->value[15] = el_par_value_recurse("sigy",p->p_elem);
    c6t_elem->value[16] = el_par_value_recurse("charge",p->p_elem);
    c6t_elem->value[17] = 0.0; /* npart */
  }
  else if((strcmp(p->base_name,"elseparator") == 0  ))
  {
    c6t_elem = new_c6t_element(3,t_name,p->base_name);
    clean_c6t_element(c6t_elem);
    strcpy(c6t_elem->org_name,t_name);
    c6t_elem->value[0] = el_par_value_recurse("l",p->p_elem);
    c6t_elem->value[2] = el_par_value_recurse("ex",p->p_elem);
    c6t_elem->value[3] = el_par_value_recurse("ey",p->p_elem);
  }
  else if((strcmp(p->base_name,"xrotation") == 0  || (strcmp(p->base_name,"yrotation") == 0) || (strcmp(p->base_name,"srotation") == 0)))
  {
    c6t_elem = new_c6t_element(3,t_name,p->base_name);
    clean_c6t_element(c6t_elem);
    strcpy(c6t_elem->org_name,t_name);
    c6t_elem->value[1] = el_par_value_recurse("angle",p->p_elem);
  }
  else if(strcmp(p->base_name,"sixmarker") == 0){
    c6t_elem = new_c6t_element(20,t_name,p->base_name);

    clean_c6t_element(c6t_elem);
    strcpy(c6t_elem->org_name,t_name);
    c6t_elem->value[1] = el_par_value_recurse("eltype",p->p_elem);
    int nl;
    double attrtemp [20];
    nl = element_vector(p->p_elem, "attr", attrtemp);
    write_f3_sixmarker(p->p_elem);

    for(int i=0; i< nl; i++){
        c6t_elem->value[i+2] = attrtemp[i];
    }

  }
  else if(strcmp(p->base_name,"wire") == 0){
    char snum[11];
    int ill_l, xma_l, yma_l, l_phy_l, l_int_l;
    strcat(t_name, "w_");
    sprintf(snum, "%d", ncombined);
    strcat(t_name, snum);
    c6t_elem = new_c6t_element(15,t_name,p->base_name);

    double current[20];
    double xma [20];
    double yma [20];     
    double l_int [20];
    double l_phy [20];
    ill_l   = element_vector(p->p_elem, "current", current);
    xma_l   = element_vector(p->p_elem, "xma", xma);
    yma_l   = element_vector(p->p_elem, "yma", yma);
    l_int_l = element_vector(p->p_elem, "l_int", l_int);
    l_phy_l = element_vector(p->p_elem, "l_phy", l_phy);
    if (ill_l==xma_l && ill_l==yma_l && ill_l==l_int_l && ill_l==l_phy_l){
      c6t_elem->value[1] = el_par_value("closed_orbit", p->p_elem);
      c6t_elem->value[2] = current[ncombined];
      c6t_elem->value[3] = l_int[ncombined];
      c6t_elem->value[4] = l_phy[ncombined];
      c6t_elem->value[5] = xma[ncombined]*1000;
      c6t_elem->value[6] = yma[ncombined]*1000;
      c6t_elem->value[7] = 0;
      c6t_elem->value[8] = 0;
    }
    else{
      mad_error("The length of the xma, yma and current is different for element :",p->base_name);
    }
  }

  else if (strcmp(p->base_name,"drift") == 0)
  {
    c6t_elem = new_c6t_element(0,t_name,p->base_name);
    clean_c6t_element(c6t_elem);
    strcpy(c6t_elem->org_name,t_name);
    c6t_elem->value[0] = el_par_value_recurse("l",p->p_elem);
  }
  else if (strcmp(p->base_name,"rfmultipole") == 0)
  {
    int maxkn=0, maxks=0, maxpn=0, maxps=0, mmult=0;
    if ((index = name_list_pos("knl",p->p_elem->def->par_names))>-1)
    {
      kn_param = p->p_elem->def->par->parameters[index];
      maxkn=kn_param->double_array->curr;
    }
    if ((index = name_list_pos("ksl",p->p_elem->def->par_names))>-1)
    {
      ks_param = p->p_elem->def->par->parameters[index];
      maxks=ks_param->double_array->curr;
    }
    if ((index = name_list_pos("pnl",p->p_elem->def->par_names))>-1)
    {
      pn_param = p->p_elem->def->par->parameters[index];
      maxpn=kn_param->double_array->curr;
    }
    if ((index = name_list_pos("psl",p->p_elem->def->par_names))>-1)
    {
      ps_param = p->p_elem->def->par->parameters[index];
      maxps=ks_param->double_array->curr;
    }

    if(six_version < general_rf_req)
    {
      if (maxkn>3 || maxks>3) {
        printf("warning while converting rfmultipole: components beyond octupole are ignored\n");
      }

      c6t_elem = new_c6t_element(19,t_name,p->base_name);
      clean_c6t_element(c6t_elem);
      strcpy(c6t_elem->org_name,t_name);
      // rf & general params
      c6t_elem->value[0] = el_par_value_recurse("l",p->p_elem);
      c6t_elem->value[1] = el_par_value_recurse("volt",p->p_elem);
      c6t_elem->value[2] = el_par_value_recurse("tilt",p->p_elem);
      c6t_elem->value[3] = el_par_value_recurse("freq",p->p_elem);
      // normal components
      c6t_elem->value[4] = maxkn>0?(kn_param->double_array->a[0]):0.0;
      c6t_elem->value[5] = maxkn>1?(kn_param->double_array->a[1]):0.0;
      c6t_elem->value[6] = maxkn>2?(kn_param->double_array->a[2]):0.0;
      c6t_elem->value[7] = maxkn>3?(kn_param->double_array->a[3]):0.0;
      c6t_elem->value[8] = maxpn>0?(pn_param->double_array->a[0]):0.0;
      c6t_elem->value[9] = maxpn>1?(pn_param->double_array->a[1]):0.0;
      c6t_elem->value[10] = maxpn>2?(pn_param->double_array->a[2]):0.0;
      c6t_elem->value[11] = maxpn>3?(pn_param->double_array->a[3]):0.0;
      // skew component
      c6t_elem->value[18] = maxks>0?(ks_param->double_array->a[0]):0.0;
      c6t_elem->value[12] = maxks>1?(ks_param->double_array->a[1]):0.0;
      c6t_elem->value[13] = maxks>2?(ks_param->double_array->a[2]):0.0;
      c6t_elem->value[14] = maxks>3?(ks_param->double_array->a[3]):0.0;
      c6t_elem->value[19] = maxps>0?(ps_param->double_array->a[0]):0.0;
      c6t_elem->value[15] = maxps>1?(ps_param->double_array->a[1]):0.0;
      c6t_elem->value[16] = maxps>2?(ps_param->double_array->a[2]):0.0;
      c6t_elem->value[17] = maxps>3?(ps_param->double_array->a[3]):0.0;
        
    }
    else 
    { //This is the new RF-multipoles
      if(maxkn >=maxks){
        mmult = maxkn;
      }
      else{
        mmult = maxks;
      }

      c6t_elem = new_c6t_element(11+mmult*4,t_name,p->base_name);
      clean_c6t_element(c6t_elem);
      strcpy(c6t_elem->org_name,t_name);

      c6t_elem->value[0] = el_par_value_recurse("l",p->p_elem);
      c6t_elem->value[1] = el_par_value_recurse("volt",p->p_elem);
      c6t_elem->value[2] = el_par_value_recurse("freq",p->p_elem);
      c6t_elem->value[3] = mmult;
      c6t_elem->value[6] = el_par_value_recurse("tilt",p->p_elem);


      for (int i=0; i<mmult; i++){
        c6t_elem->value[7+i*4]  = maxkn>i?(kn_param->double_array->a[i]):0.0;
        c6t_elem->value[8+i*4]  = maxpn>i?(pn_param->double_array->a[i]):0.0;
        c6t_elem->value[9+i*4]  = maxks>i?(ks_param->double_array->a[i]):0.0;
        c6t_elem->value[10+i*4] = maxps>i?(ps_param->double_array->a[i]):0.0;
      }
    }
  }
  else
  {
    printf("Element not convertible! name= %s, basename = %s\n",p->name,p->base_name);
  }


  if (c6t_elem)
  {
    for (j = 0; j < c6t_elem->n_values; j++)
      if (fabs(c6t_elem->value[j]) < eps_12)
        c6t_elem->value[j] = 0.0;

    /* check to see if this has an aperture assigned, check for aperture flag */
    if ((aperture_flag)
        && (aper_param = return_param_recurse("apertype", p->p_elem)))
    {
      tag_aperture.apply=1;
      strcpy(tag_aperture.style,aper_param->string);
      strcpy(tag_aperture.name,t_name);
      strcat(tag_aperture.name,"_AP");
      for (i=0; i<8; i++) tag_aperture.value[i] = 0. ;

      if ((aper_param = return_param_recurse("aperture", p->p_elem)))
      {
        if (aper_param->expr_list != NULL)
          update_vector(aper_param->expr_list, aper_param->double_array);
  j = 4;
  if (aper_param->double_array->curr < 4) j = aper_param->double_array->curr;
        for(i=1; i<=j; i++) tag_aperture.value[i] = aper_param->double_array->a[i-1];
      }

      if ((aper_param = return_param_recurse("aper_offset", p->p_elem)))
      {
        if (aper_param->expr_list != NULL)
          update_vector(aper_param->expr_list, aper_param->double_array);
  j = 2;
  if (aper_param->double_array->curr < 2) j = aper_param->double_array->curr;
        for(i=1; i<=j; i++) tag_aperture.value[i+4] = aper_param->double_array->a[i-1];
      }

      if ((aper_param = return_param_recurse("aper_tilt", p->p_elem)))
      {
        tag_aperture.value[7] = el_par_value("aper_tilt", p->p_elem);
      }
    }

    /* name used has to be without occ_cnt as this is added
       (only 1) in tab_name_code
       !!! fixed FS 17.08.2004 !!! */
    c6t_elem->twtab_row = my_table_row(current_sequ->tw_table,t_name);
  }

  return c6t_elem;
}

#if 0 // not used
static void
dump_c6t_element(struct c6t_element* el)
{
  int j;
  char pname[24] = "NULL", nname[24] = "NULL";
  if (el->previous != NULL) strcpy(pname, el->previous->name);
  if (el->next != NULL) strcpy(nname, el->next->name);
  printf("name: %s  base: %s  position: %f\n", el->name, el->base_name,
         el->position);
  printf("  names of - previous: %s  next: %s  equiv: %s\n",
         pname, nname, el->equiv->name);
  printf("  flag: %d  force: %d  split: %d keep: %d values: %d f_errors: %d\n",
         el->flag, el->force, el->split, el->keep_in,el->n_values, el->nf_err);
  for (j = 0; j < el->n_values; j++)
  {
    printf("%e ", el->value[j]);
    if ((j+1)%5 == 0 || j+1 == el->n_values)  printf("\n");
  }
  if (el->nf_err)  puts("field errors");
  for (j = 0; j < el->nf_err; j++)
  {
    printf("%e ", el->p_fd_err->a_dble[j]);
    if ((j+1)%5 == 0 || j+1 == el->nf_err)  printf("\n");
  }
  if (el->na_err)  puts("alignment errors");
  for (j = 0; j < el->na_err; j++)
  {
    printf("%e ", el->p_al_err->a_dble[j]);
    if ((j+1)%5 == 0 || j+1 == el->na_err)  printf("\n");
  }
}
#endif

#if 0 // not used
static void
dump_c6t_sequ(int level)
{
  double suml = zero;
  puts("+++++++++ dump sequence +++++++++");
  current_element = first_in_sequ;
  while (current_element != NULL)
  {
    suml += current_element->value[0];
    if (level > 2)  dump_c6t_element(current_element);
    else if (level > 1)  gnu_file(current_element);
    else if (level > 0 && strcmp(current_element->base_name, "drift") != 0)
      printf("%s: %s at = %f\n", current_element->name,
             current_element->equiv->name, current_element->position);
    current_element = current_element->next;
  }
  printf("=== sum of element length: %f\n", suml);
}
#endif

#if 0 // not used
static void
dump_types(int flag)
{
  int i, j;

  puts("+++++++++ dump types +++++++++");
  for (i = 0; i < types.curr; i++)
  {
    puts(types.member[i]->base_name);
    for (j = 0; j < types.member[i]->curr; j++)
      printf("       %s  %f\n", types.member[i]->elem[j]->name,
             types.member[i]->elem[j]->value[0]);
    if (flag > 0) dump_c6t_element(types.member[i]->elem[j]);
  }
}
#endif

static void
equiv_elem(void)
{
  int i, j, k;
  struct c6t_element *el, *eln;

  for (i = 0; i < types.curr; i++)  /* loop over base types */
  {
    for (j = 0; j < types.member[i]->curr; j++) /* loop over el. in type */
    {
      el = types.member[i]->elem[j];
      if (el->flag > 0)  /* all others ignored */
      {
        if (el->equiv == el /* not yet equivalenced */
            && strcmp(el->base_name,"marker") != 0 /* do not touch markers and wires*/
            && strcmp(el->base_name,"wire") != 0)
        {
          for (k = j+1; k < types.member[i]->curr; k++)
          {
            eln = types.member[i]->elem[k];
            if (eln->flag > 0
                && eln->equiv == eln
                && ident_el(el, eln) == 0
                && eln->nf_err == el->nf_err
                && strcmp(eln->base_name,"beambeam") != 0
                && strcmp(eln->base_name,"marker") != 0
                && strstr(eln->base_name,"colli") == NULL)
              eln->equiv = el;
          }
        }
      }
    }
  }
}

static int
f34_values(struct c6t_element* el, int* flags, double* values)
{
  int i, j, np, nd, cnt = 0;
  double pow, tmp[FIELD_MAX];
  for (i = 0; i < FIELD_MAX; i++)
  {
    tmp[i] = zero;
    j = i + 12;
    if (j < el->n_values && el->value[j] != zero)
    {
      if (el->value[0] != zero) tmp[i] += el->value[0] * el->value[j];
      else tmp[i] += el->value[j];
    }
    if (i < el->nf_err && el->p_fd_err->a_dble[i] != zero)
      tmp[i] += el->p_fd_err->a_dble[i];
  }
  for (i = 3; i < FIELD_MAX; i++)
  {
    if (tmp[i] != zero)
    {
      np = i / 2 + 1;
      nd = 1; for (j = 2; j < np; j++)  nd *= j;
      pow = nd;
      pow = power_of(ten, 6-3*np) / pow;
      if (i%2 == 0)
      {
        flags[cnt] = np; if (el->npole_sign) pow = -pow;
      }
      else           flags[cnt] = -np;
      values[cnt++] = pow * tmp[i];
    }
  }
  return cnt;
}

static struct block*
get_block_equiv(struct block* current)
{
  struct block* p = first_block;
  int i, k;
  while (p != current)
  {
    if (current->elements->curr == p->elements->curr)
    {
      k = 0;
      for (i = 0; i < current->elements->curr; i++)
      {
        if (strcmp(current->elements->elem[i]->equiv->name,
                   p->elements->elem[i]->equiv->name) == 0) k++;
      }
      if (k == current->elements->curr)  return p;
    }
    p = p->next;
  }
  return p;
}

static void
get_args(struct in_cmd* my_cmd)
{
  int tmp_max_mult_ord;
  double tmp_ref_def;
  if ((aperture_flag = command_par_value("aperture", my_cmd->clone)))
    put_info("c6t - aperture flag selected","");
  if ((cavall_flag = command_par_value("cavall", my_cmd->clone)))
    put_info("c6t - cavall flag selected","");
  if ((markall_flag = command_par_value("markall", my_cmd->clone)))
    put_info("c6t - markall flag selected","");
  if ((mult_auto_off = command_par_value("mult_auto_off", my_cmd->clone)))
    put_info("c6t - mult_auto_off flag selected","");
  if ((split_flag = command_par_value("split", my_cmd->clone)))
    put_info("c6t - split flag selected","");
  if ((mangle_flag = command_par_value("mangle", my_cmd->clone)))
    put_info("c6t - mangle flag selected","");
  if ((long_names_flag = command_par_value("long_names", my_cmd->clone)))
    put_info("c6t - long names flag selected","");
  if ((multicol_flag = command_par_value("multicol", my_cmd->clone)))
    put_info("c6t - multicol flag selected","");
  if ((six_version = command_par_value("six_version", my_cmd->clone)))
    printf("SixTrack Version of (or higher is assumed): %f\n",six_version);
  if (mult_auto_off &&
     (tmp_max_mult_ord = command_par_value("max_mult_ord", my_cmd->clone))>0)
  {
    max_mult_ord = tmp_max_mult_ord;
    printf("max_mult_ord set to : %d\n",max_mult_ord);
  }
  if ((tmp_ref_def = command_par_value("radius", my_cmd->clone))>0.)
  {
//    radius_flag = 1; // not_used
    ref_def = tmp_ref_def;
    printf("Reference radius set to : %f\n",ref_def);
  }

}

static void
get_error_refs(struct c6t_element* el)
{
  int i;
  double tmp;
  i = 12 + 2 * el->mult_order;
  if (i+1 < el->n_values)
  {
    tmp = fabs(el->value[i]) > fabs(el->value[i+1]) ?
      fabs(el->value[i]) : fabs(el->value[i+1]);
  }
  else if(i < el->n_values) tmp = fabs(el->value[i]);
  else tmp = 1;
  if (tmp == zero) tmp = 1;
  if (el->mult_order==0)
  {
    el->ref_delta = c1p3 * tmp * power_of(el->ref_radius, el->mult_order);
  }
  else
  {
    el->ref_delta = 0;
  }
}

static int
get_flag(struct c6t_element* el, struct type_info* type)
{

  if (el->value[0] == zero)
  {
    if (type->flag_1 == 4) return in_keep_list(el);
    else return type->flag_1;
  }
  if (el->n_values < 7) return type->flag_2;
  else return (el->value[6] == 0 ? type->flag_2 : type->flag_3);
}

#if 0 // not used
static struct c6t_element*
get_from_ellist(char* name, char* type)
{
  int i, j;

#ifdef _call_tree_
  puts("+++++++ get_from_ellist");
#endif
  for (i = 0; i < types.curr; i++)
  {
    if (strcmp(types.member[i]->base_name, type) == 0)
    {
      for (j = 0; j < types.member[i]->curr; j++) /* loop over el. in type */
      {
        if (strcmp(types.member[i]->elem[j]->name, name) == 0)
          return types.member[i]->elem[j];
      }
    }
  }
  return NULL;
}
#endif

static void
get_multi_refs(void)
{
  int i;
  for (i = 0; i < types.curr; i++)  /* loop over base types */
  {
    if (strcmp(types.member[i]->base_name, "multipole") == 0)
    {
      multi_type = i;  break;
    }
  }
}

static int
get_next_name(char* name, char acro)
{
  int j, k = -1;

  for (j = 0; j < acro_occ; j++)
    if (acro_list[j] == acro)  k = j;

  if (k < 0) {
    k = acro_occ++; acro_list[k] = acro; acro_cnt[k] = 0;
  }

  sprintf(name, "%c_c6t_%d", acro, ++acro_cnt[k]);

  return 1;
}

#if 0 // not used
static void
gnu_file(struct c6t_element* el)
{
  double el_start, el_end;

  el_start = el->position - el->value[0] / two;
  el_end   = el->position + el->value[0] / two;
  printf("%s %e 0.\n", el->name, el_start);
  printf("%s %e 1.\n", el->name, el_start);
  printf("%s %e 1.\n", el->name, el_end);
  printf("%s %e 0.\n", el->name, el_end);
}
#endif

static void
grow_ellist( /* doubles object list size */
  struct c6t_el_list* p)
{
  struct c6t_element** p_loc = p->elem;
  int j, new = 2*p->max;
  const char *rout_name = "c6t:grow_ellist";

#ifdef _call_tree_
  puts("+++++++ grow_ellist");
#endif
  p->max = new;
  p->elem = mycalloc(rout_name, new, sizeof *p->elem);
  for (j = 0; j < p->curr; j++) p->elem[j] = p_loc[j];
  myfree(rout_name, p_loc);
}

static int
ident_el(struct c6t_element* el1, struct c6t_element* el2)
{
  double s, tolerance = eps_9;
  int j, m = el1->n_values < el2->n_values ? el1->n_values : el2->n_values;
  if (el1->mult_order != el2->mult_order)  return 1;
  if (el1->ref_radius != el2->ref_radius)  return 2;
  if (el1->ref_delta  != el2->ref_delta)   return 2;
  if (el1->keep_in  != el2->keep_in)       return 6;
  for (j = 0; j < m; j++)
  {
    s = fabs(el1->value[j]) + fabs(el2->value[j]);
    if (s > zero
        && fabs(el1->value[j] - el2->value[j])/s > tolerance) return 3;
  }
  if (m != el1->n_values)
  {
    for (j = m; j < el1->n_values; j++)
      if (el1->value[j] != zero) return 4;
  }
  else if (m != el2->n_values)
  {
    for (j = m; j < el2->n_values; j++)
      if (el2->value[j] != zero) return 5;
  }
  return 0;
}

static int
ident_zero(struct c6t_element* el)
{
  int j;
  for (j = 12; j < el->n_values; j++)
    if (el->value[j] != zero) return 1;
  return 0;
}

static int
in_keep_list(struct c6t_element* el)
{
  static char keep_these[MM_KEEP][24] = {"ip", "mt_"};
  char temp[24];
  int j;

  if (markall_flag) return 2;

  strcpy(temp, el->name); stolower(temp);
  for (j = 0; j < MM_KEEP; j++)
  {
    if (strncmp(temp, keep_these[j], strlen(keep_these[j])) == 0) return 2;
  }
  return 0;
}

static void
invert_normal(int count, double array[])
{
  int i;
  for (i = 0; i < (count+1)/2; i++)  array[2*i] = -array[2*i];
}

#if 0 // not used
static void
invert_skew(int count, double array[])
{
  int i;
  for (i = 0; i < count/2; i++)  array[2*i+1] = -array[2*i+1];
}
#endif

static void
link_c6t_in_front(struct c6t_element* new, struct c6t_element* el)
{
  if (el->previous == NULL) first_in_sequ = new;
  else el->previous->next = new;
  new->previous = el->previous; new->next = el;
  el->previous = new;
}

static void
link_behind(struct c6t_element* new, struct c6t_element* el)
{
  if (el->next == NULL)
  {
//    last_in_sequ = new; // not used
    last_in_sequ_org = new;
  }
  else el->next->previous = new;
  new->previous = el; new->next = el->next;
  el->next = new;
}

/* removed, replaced by stolower from mad_str.h
void lower(char* s)
{
  char* cp = s;
  while(*cp != '\0')
  {
    *cp = (char) tolower((int)*cp); cp++;
  }
}
*/

static struct c6t_element*
make_c6t_element(struct node* p, int ncombined)
{
  struct c6t_element *tmp_element;
  if ((tmp_element = convert_madx_to_c6t(p, ncombined)))
  {
    prev_element = current_element;
    current_element = tmp_element;
    if (elem_cnt++ == 0) first_in_sequ = current_element;
    else                 prev_element->next = current_element;
    current_element->previous = prev_element;
    current_element->next = NULL;
  }
  return tmp_element;
}

/* this is taken from doom but all it does is malloc the structure and fill
   it in */
static struct object*
make_obj(   /* creates a new object */
  const char* key,
  int vlint,       /* length of integer array */
  int vldble,      /* length of double array */
  int vlchar,      /* length of char array */
  int vlpobj)      /* length of object pointer array */
{
  struct object* p;
  const char *rout_name = "c6t:make_obj";

#ifdef _call_tree_
  put_info("+++++++ make_object","");
#endif
  p = mycalloc(rout_name, 1, sizeof *p);
  mycpy(p->key, key);
  if ((p->l_int  = vlint ) > 0) p->a_int  = mymalloc_atomic(rout_name, p->l_int  * sizeof *p->a_int );
  if ((p->l_dble = vldble) > 0) p->a_dble = mymalloc_atomic(rout_name, p->l_dble * sizeof *p->a_dble);
  if ((p->l_char = vlchar) > 0) p->a_char = mymalloc_atomic(rout_name, p->l_char * sizeof *p->a_char);
  if ((p->l_obj  = vlpobj) > 0) {
    p->p_obj = mycalloc(rout_name, p->l_obj, sizeof *p->p_obj);
    p->names = mycalloc(rout_name, p->l_obj, sizeof *p->names);
  }
  p->parent = NULL;
  /*      my_time(); */
  /*      p->ma_time = major_time; p->mi_time = minor_time; */
  return p;
}

static void
make_multipole(struct c6t_element* el)
{
  if (el->force > 0) /* multiply forces with length */
    app_factor(el->value[0], &el->value[12], el->n_values-12);
  el->value[0] = zero; /* set element length to zero */
  el->flag = 2;
  remove_from_ellist(el);
  strcpy(el->base_name, "multipole");
  add_to_ellist(el);
}

static void
mod_errors(void)
{
  current_element = first_in_sequ;
  while (current_element != NULL)
  {
    if (current_element->nf_err > 0)
      invert_normal(current_element->nf_err, current_element->p_fd_err->a_dble);
    current_element = current_element->next;
  }
}

static void
mod_lcavity(struct c6t_element* p)
{
  total_voltage += p->value[1];   /* accumulate voltage */
}

static void
mod_multipole(struct c6t_element* p)
{
  supp_small_comp(p);
}

static void
mod_rfmultipole(struct c6t_element* p)
{
  supp_small_comp(p);
}

static void
mod_octupole(struct c6t_element* p)
{
  supp_small_comp(p);
}

static void
mod_quadrupole(struct c6t_element* p)
{
  supp_small_comp(p);
}

static void
mod_rbend(struct c6t_element* p)
{
  supp_small_comp(p);
  /* commented out in the original c6t (MEH)
     set length to arc length
     if (fabs(p->a_dble[0]) > eps_9 && fabs(p->a_dble[10]) > eps_9)
     p->a_dble[0] *= p->a_dble[10] / (two * sin(p->a_dble[10]/two));
  */
  /* change 1/r to value for straight line */
  p->value[1] = p->value[10] / p->value[0];
}

static void
mod_rfcavity(struct c6t_element* p)
{
  total_voltage += p->value[1];  /* accumulate voltage */
}

static void
mod_crabcavity(struct c6t_element* p)
{
  total_voltage += p->value[1];  /* accumulate voltage */
}

static void
mod_sextupole(struct c6t_element* p)
{
  supp_small_comp(p);
  /* supress tilt angle (not component !) */
  /*    p->value[6] = zero; */
}

static void
multi_loop(void)
{
  int i, j, nup;
  struct c6t_element* el;
  for (i = 0; i < types.curr; i++)  /* loop over base types */
  {
    if (strcmp(types.member[i]->base_name, "multipole") == 0)
    {
      nup = types.member[i]->curr;
      for (j = 0; j < nup; j++) /* loop over mutipoles */
      {
        el = types.member[i]->elem[j];
        pre_multipole(el);
      }
    }
  }
}

static struct block*
new_block(void)
{
  struct block* p;
  const char *rout_name = "c6t:new_block";
  p = mycalloc(rout_name, 1, sizeof *p);
  sprintf(p->name, "BLOC%d", block_count++);
  return p;
}

static struct c6t_element*
new_c6t_element(int size, const char* name, const char* base)
{
  struct c6t_element* p;
  const char *rout_name = "c6t:new_c6t_element";
  p = mycalloc(rout_name, 1, sizeof *p);
  strcpy(p->name, name);
  p->equiv = p;
  strcpy(p->base_name, base);
  ++size;
  p->value = mycalloc_atomic(rout_name, size, sizeof *p->value);
  p->n_values = size;
  p->do_not_free = 0;
  return p;
}

static void
post_multipoles(void) /* post equiv. treatment of multipoles */
{
  int i, j;
  struct c6t_element *el, *eln;
  struct object *p;

  if (multi_type > -1) /* there are multipoles */
  {
    for (j = 0; j < types.member[multi_type]->curr; j++)
    {
      el = types.member[multi_type]->elem[j]; eln = el->equiv;
      if (el->nf_err > 0)
      {
        eln->mult_order = el->mult_order;
        eln->ref_radius = el->ref_radius;
        if (eln->p_fd_err == NULL)
        {
          eln->p_fd_err = p_err_zero;
          eln->nf_err = FIELD_MAX;
        }
        if (eln->nf_err < el->nf_err)
        {
          strcpy(tmp_name, eln->p_fd_err->key);
          p = eln->p_fd_err;
          eln->p_fd_err = make_obj(tmp_name, 0, el->nf_err, 0, 0);
          /* first initialise */
          for (i = 0; i < el->nf_err; i++)
            eln->p_fd_err->a_dble[i] = 0.0;
          for (i = 0; i < eln->nf_err; i++)
            eln->p_fd_err->a_dble[i] = p->a_dble[i];
          eln->nf_err = el->nf_err;
        }
      }
    }
  }
}

static double
power_of(double d, int i)
{
  int j;
  double tmp = 1;
  if (i > 0)
  {
    for (j = 0; j < i; j++)  tmp *= d;
    return tmp;
  }
  else if(i < 0)
  {
    for (j = 0; j < -i; j++)  tmp *= d;
    return 1./tmp;
  }
  else return 1.;
}

static void
pre_multipole(struct c6t_element* el) /* pre-process multipoles */
{
  /*
    1. first count multipole components < decapole (cnt)
    if   cnt == 1:
    if    dipole, set el->nc_pos
    else  make n_pole, insert in front, zero el->component
    2. add all comp. > dipole + errors
    if all zero:
    straight quad, yank element
    else  set error count to zero, keep n-pole
    else
    put sum into errors, zero all el->components
  */
  int i, last_nzero = -1, nz_cnt = 0, cnt = 0, s_pole = 0, ndmax;
  // int n_pole, // not used
  int low, new_el_t;
  struct c6t_element *new_el = NULL;
  char t_list[5][12] = {"dipole", "quadrupole", "sextupole", "octupole",
                        "decapole"};

  ndmax = el->n_values > 22 ? 22 : el->n_values;
  for (i = 12; i < ndmax; i++)
  {
    if (el->value[i] != zero)
    {
      s_pole = i; cnt++;
    }
  }
  if ((cnt == 1) || (el->value[12]!=zero) || (el->value[13]!=zero))
  {
    if (el->value[12]!=zero) { s_pole=12; cnt=1; }
    if (el->value[13]!=zero) { s_pole=13; cnt=1; }
    if ((new_el_t = (s_pole-12)/2) == 0)  el->nc_pos = s_pole;
    else
    {
      get_next_name(tmp_name, t_list[new_el_t][0]);
      new_el = new_c6t_element(s_pole+1, tmp_name, t_list[new_el_t]);
      new_el->do_not_free = 1;
      for (i = 0; i <= s_pole; i++) new_el->value[i] = el->value[i];
      for (i = 12; i <= s_pole; i++) el->value[i] = 0;
      new_el->flag = s_pole > 13 ? 2 : 1; new_el->npole_sign = 1;
      new_el->keep_in = el->keep_in;
      new_el->position = el->position;
      new_el->twtab_row = el->twtab_row;
      new_el->na_err = el->na_err; /* el->na_err = 0; */
      new_el->p_al_err = el->p_al_err; /*  el->p_al_err = NULL; */
      new_el->tilt_err = el->tilt_err; /*  keep tilt info */
      link_c6t_in_front(new_el, el);
      add_to_ellist(new_el);
      strcpy(tmp_name, el->name); strcpy(el->name, new_el->name);
      strcpy(new_el->name, tmp_name);
    }
  }
  for (i = 0; i < FIELD_MAX; i++) tmp_buff[i] = zero;
  low = cnt == 1 ? 2 : 0;
  for (i = low; i < el->n_values-14; i++)
    tmp_buff[i] = el->value[12+i];
  for (i = 0; i < el->nf_err; i++) tmp_buff[i] += el->p_fd_err->a_dble[i];
  for (i = 0; i < FIELD_MAX; i++) if (tmp_buff[i] != zero)
  {
    last_nzero = i; nz_cnt++;
  }
  // n_pole = last_nzero / 2; // not used
  el->rad_length = el->value[11];
  if (last_nzero < 0)  /* all quad+ components and all errors = zero */
  {
    el->nf_err = 0;
    if (el->keep_in == 0 && (s_pole == 0 || s_pole > 13)) yank(el);
    else                            el->nc_pos = s_pole;
  }
  else  /* element becomes multipole, all comp. above dipole -> errors */
  {
    el->nc_pos = s_pole > 13 ? 0 : s_pole;
    if (el->ref_delta == zero)
    {
      el->ref_delta = c1p3;
      el->ref_radius = ref_def;
      el->mult_order = 1;
    }
    if (++last_nzero > el->nf_err)
    {
      if (el->p_fd_err != NULL) strcpy(tmp_name, el->p_fd_err->key);
      else  snprintf(tmp_name, sizeof tmp_name, "%.42s_arfa", el->name);
      el->nf_err = last_nzero;
      el->p_fd_err = make_obj(tmp_name, 0, el->nf_err, 0, 0);
    }
    for (i = 0; i < el->nf_err; i++) el->p_fd_err->a_dble[i] = tmp_buff[i];
    for (i = 14; i < el->n_values; i++)  el->value[i] = zero;
  }
}

static void
pro_elem(struct node* cnode,  int ncombined)
/* processes one element, makes linked list */
/* converts MADX linked list to c6t internal linked list */
{
  int i;
  char t_key[KEY_LENGTH];
  struct c6t_element *tag_element, *tmp_element;
  double tmp_vk,tmp_hk;

  tag_aperture.apply=0;
  /* do the fiddly conversion but skip element if not needed */
  if (make_c6t_element(cnode, ncombined) == NULL) return;

  if      (strcmp(cnode->base_name, "rbend") == 0)       mod_rbend(current_element);
  else if (strcmp(cnode->base_name, "lcavity") == 0)     mod_lcavity(current_element);
  else if (strcmp(cnode->base_name, "multipole") == 0)   mod_multipole(current_element);
  else if (strcmp(cnode->base_name, "octupole") == 0)    mod_octupole(current_element);
  else if (strcmp(cnode->base_name, "quadrupole") == 0)  mod_quadrupole(current_element);
  else if (strcmp(cnode->base_name, "sextupole") == 0)   mod_sextupole(current_element);
  else if (strcmp(cnode->base_name, "rfcavity") == 0)    mod_rfcavity(current_element);
  else if (strcmp(cnode->base_name, "crabcavity") == 0)  mod_crabcavity(current_element);
  else if (strcmp(cnode->base_name, "rfmultipole") == 0) mod_rfmultipole(current_element);

  if (strstr(cnode->base_name, "kicker") || strstr(cnode->base_name, "tkicker"))
  {
    if (cnode->p_elem)
    {
      tmp_hk = el_par_value("hkick",cnode->p_elem); current_element->value[12] += tmp_hk;
      tmp_vk = el_par_value("vkick",cnode->p_elem); current_element->value[13] += tmp_vk;
    }
    current_element->value[12] += cnode->chkick;
    current_element->value[13] += cnode->cvkick;
  }

  strcpy(current_element->org_name, current_element->name);
  current_element->occ_cnt = cnode->occ_cnt;
  if (cnode->occ_cnt > 1)  /* add occurence count to name */
  {
    snprintf(t_key, sizeof t_key, "%.45s+%d", current_element->name,cnode->occ_cnt);
    strcpy(current_element->name, t_key);
  }
  current_element->position = cnode->position;

  /* errors in MADX are stored in the same way as c6t */

  if (cnode->p_fd_err) {
    field_cnt++;
    current_element->nf_err = cnode->p_fd_err->curr;
    current_element->p_fd_err = make_obj("FDDUM",0,FIELD_MAX,0,0);
    current_element->p_fd_err->c_dble = cnode->p_fd_err->curr;
    for (i=0;i<cnode->p_fd_err->curr;i++)
      current_element->p_fd_err->a_dble[i] = cnode->p_fd_err->a[i];
    for (i=12;((i<current_element->n_values) && (current_element->value[i]==0));i+=2);
    if (i>current_element->n_values) i-=2;
    current_element->mult_order = i-12;
    current_element->ref_radius = ref_def;
    get_error_refs(current_element);
  }
  else if(cnode->p_fd_err==NULL && (strcmp(cnode->base_name, "multipole") == 0) && el_par_value("angle",cnode->p_elem)!=0){
      double knl_tmp [FIELD_MAX];
      element_vector(cnode->p_elem, "knl", knl_tmp);
      //printf("angggleee", knl[0],  )
      if(fabs(knl_tmp[0] - el_par_value("angle",cnode->p_elem))>1e-7){
        field_cnt++;
        current_element->nf_err = 1;
        current_element->p_fd_err = make_obj("FDDUM",0,FIELD_MAX,0,0);
        current_element->p_fd_err->c_dble = 1;
        current_element->p_fd_err->a_dble[0]=-999; //maybe 0 works to be seen.. .
        current_element->mult_order = 0;
        current_element->ref_radius = ref_def;
        get_error_refs(current_element);
      } 
  }

  if (cnode->p_al_err) {
    align_cnt++;
    current_element->na_err = cnode->p_al_err->curr;
    current_element->p_al_err = make_obj("ALDUM",0,ALIGN_MAX,0,0);
    current_element->p_al_err->c_dble = cnode->p_al_err->curr;
    for (i=0;i<cnode->p_al_err->curr;i++)
      current_element->p_al_err->a_dble[i] = cnode->p_al_err->a[i];
  }

  /* if we have a tilt set the flag */
  ///// AL: WARNING: NOT ALL ELEMENTS STORE "TILT" in value[6]
  current_element->tilt_err = 0;
  if (current_element->n_values >= 7 && fabs(current_element->value[6]) > zero) {
    align_cnt++;
    current_element->tilt_err = 1;
  }

  // store the aperture type as keyword
  char keyword[3]="00";
  if (tag_aperture.apply == 1 ) {
    if      (0 == strcmp(tag_aperture.style,"circle"))      strcpy(keyword, "CR");
    else if (0 == strcmp(tag_aperture.style,"ellipse"))     strcpy(keyword, "EL");
    else if (0 == strcmp(tag_aperture.style,"rectangle"))   strcpy(keyword, "RE");
    else if (0 == strcmp(tag_aperture.style,"rectcircle") ||
       0 == strcmp(tag_aperture.style,"lhcscreen"))   strcpy(keyword, "RC");
    else if (0 == strcmp(tag_aperture.style,"rectellipse")) strcpy(keyword, "RL");
    else if (0 == strcmp(tag_aperture.style,"racetrack"))   strcpy(keyword, "RT");
    else if (0 == strcmp(tag_aperture.style,"octagon"))     strcpy(keyword, "OC");
    else warning("general aperture element not supported in sixtrack",tag_aperture.name);
 }

  // 2015-Oct-13  16:10:42  ghislain: before we add the current element to element list
  // check whether the element is long and has an aperture, and therefore
  // whether we should insert an aperture marker before the current element
  if (strcmp(keyword,"00") != 0 && current_element->value[0] > 0.) {


    tmp_element = current_element;
    tag_element = create_aperture(tag_aperture.name, keyword,
            tag_aperture.value[1], tag_aperture.value[2],
            tag_aperture.value[3], tag_aperture.value[4],
            tag_aperture.value[5], tag_aperture.value[6],
            tag_aperture.value[7],
            cnode->p_al_err);

    tag_element->previous = prev_element;
    prev_element->next = tag_element;
    tag_element->position = cnode->position - current_element->value[0] / 2.;

    current_element = tag_element;
    add_to_ellist(current_element);

    tmp_element->previous = current_element;
    current_element->next = tmp_element;
    tmp_element->position = cnode->position;

    prev_element = current_element;
    current_element = tmp_element;

  }


  add_to_ellist(current_element);


  /* add aperture element if necessary */
  if (strcmp(keyword,"00") != 0) 
    {
    if(tag_aperture.value[1] > 0 || tag_aperture.value[2] > 0 || tag_aperture.value[3] > 0)
    {
      tag_element = create_aperture(tag_aperture.name, keyword,
      tag_aperture.value[1], tag_aperture.value[2],
      tag_aperture.value[3], tag_aperture.value[4],
      tag_aperture.value[5], tag_aperture.value[6],
      tag_aperture.value[7],
      cnode->p_al_err);
      tag_element->previous = current_element;
      tag_element->next = current_element->next;
      current_element->next = tag_element;

      prev_element = current_element;
      current_element = tag_element;
      // 2015-Oct-13  15:45:58  ghislain: add half the prev_element length to account for thick elements
      current_element->position = cnode->position + prev_element->value[0] / 2.;
      add_to_ellist(current_element);
    }
    else
    {
      warning("An aperture type has been defined without any settings. It will NOT be converted:",
      tag_aperture.name);
    }
  }
}


static void
read_sequ(void)
{
  struct node* cnode;
  if ((current_sequ->n_nodes) > 0)  sequ_start = current_sequ->ex_start->position;
  cnode = current_sequ->ex_start;
  while (cnode && cnode != current_sequ->ex_end)
  {
    int ncombined = 0;
    if(strcmp(cnode->base_name, "wire") == 0){
      double len = el_par_value("l", cnode->p_elem);
      if(fabs(len) > 0)  mad_error("Wire elements length needs to be 0","Makethin will save you! ");    
      double inorm [20];
      ncombined = element_vector(cnode->p_elem, "current", inorm);
      for(int i=0; i<ncombined; i++){
        printf("nnwires %d \n", i);
        pro_elem(cnode, i);
      }
      cnode = cnode->next;
    }
    else{
      if (strstr(cnode->name,"$")==NULL) pro_elem(cnode, 0);
      cnode = cnode->next;
    }
  }

  sequ_end = current_sequ->ex_end->position;
  sequ_length = sequ_end - sequ_start;
//  last_in_sequ = current_element; // not used
  last_in_sequ_org = current_element;
  put_info("MADX sequence converted to c6t internal.","");
}

/* removes element from correct object list */
static void
remove_from_ellist(struct c6t_element* p_elem)
{
  int i, j;

#ifdef _call_tree_
  puts("+++++++ remove_from_ellist");
#endif
  for (i = 0; i < types.curr; i++)
  {
    if (strcmp(types.member[i]->base_name, p_elem->base_name) == 0)
    {
      for (j = 0; j < types.member[i]->curr; j++) /* loop over el. in type */
      {
        if (types.member[i]->elem[j] == p_elem)
        {
          types.member[i]->elem[j]
            = types.member[i]->elem[--types.member[i]->curr];
          return;
        }
      }
    }
  }
}

static void
replace_c6t(struct c6t_element* old, struct c6t_element* new)
{
  if (old->previous != NULL)  old->previous->next = new;
  new->previous = old->previous;
  if (old->next != NULL)      old->next->previous = new;
  new->next = old->next;
  old->flag = 0;
}

static void
split(void)
{
  if (split_list != NULL)
  {
    int i;
    for (i = 0; i < split_list->curr; i++)
    {
      struct c6t_element *el = split_list->elem[i];
      if (el->flag == 1
          && (split_flag != 0 || el->nf_err > 0)) split_special(el);
      else if (el->flag == 2 || el->flag == 3)  split_other(el);
      else if (el->split == 3)  split_kicker(el);
    }
  }
}


static void
split_kicker(struct c6t_element* el)
{
  struct c6t_element *k1, *k2;
  char c[24];

  if (el->value[0] > zero) split_other(el);
  if (el->flag == 6) /*split kicker into h + v */
  {
    get_next_name(c, 'h');
    k1 = new_c6t_element(13, c, "hkicker");
    get_next_name(c, 'v');
    k2 = new_c6t_element(14, c, "vkicker");
    k1->force = k2->force = el->force;
    k1->value[12] = el->value[12];
    k2->value[13] = el->value[13];
    k1->flag = k2->flag = 5;
    k1->position = el->position;
    k2->position = el->position;
    replace_c6t(el, k1);
    link_behind(k2, k1);
    add_to_ellist(k1); add_to_ellist(k2);
  }
}

static void
split_other(struct c6t_element* el)
{
  /* -> two drifts with non-lin. thin lens at centre */
  struct c6t_element *d1, *d2;
  double length = el->value[0] / two;
  char c[24];

  get_next_name(c, 'd');
  d1 = new_c6t_element(1, c, "drift");
  get_next_name(c, 'd');
  d2 = new_c6t_element(1, c, "drift");
  d1->value[0] = d2->value[0] = length;
  d1->flag = d2->flag = 1;
  d1->position = el->position - d1->value[0] / two;
  d2->position = el->position + d2->value[0] / two;
  treat_split(el);
  link_c6t_in_front(d1, el);
  link_behind(d2, el);
  add_to_ellist(d1); add_to_ellist(d2);
}

static void
split_special(struct c6t_element* el)
/* -> two lin. halves with multipole at centre */
{
  struct c6t_element *d1, *mt;
  double length = el->value[0] / two, mt_position = el->position;
  char c[24];
  int j;

  if (el->nf_err > 0)
    app_factor(el->value[0], el->p_fd_err->a_dble, el->nf_err);
  get_next_name(c, *el->base_name);
  d1 = new_c6t_element(el->n_values, c, el->base_name);
  d1->value[0] = el->value[0] = length;
  for (j = 1; j < el->n_values; j++)
    d1->value[j] = el->value[j];
  d1->flag = el->flag = 1;
  d1->force = el->force = 1;
  d1->position = el->position + d1->value[0] / two;
  el->position = el->position - el->value[0] / two;
  get_next_name(c, 'm');
  mt = new_c6t_element(MULTI_MAX, c, "multipole");
  mt->force = 1;
  mt->flag = 2; mt->position = mt_position;
  mt->n_values = el->n_values; /* multipole forces remain zero */
  mt->p_fd_err = el->p_fd_err; el->p_fd_err = NULL;
  mt->nf_err = el->nf_err; el->nf_err = 0;
  mt->p_al_err = el->p_al_err; mt->na_err = el->na_err;
  /* el->p_al_err = NULL; el->na_err = 0; */
  mt->keep_in = split_flag;
  add_to_ellist(mt);
  add_to_ellist(d1);
  link_behind(mt, el);
  link_behind(d1, mt);
}

static void
supp_elem(void)
{
  struct c6t_element *d1, *el;
  char c[24];
  int af;

  current_element = first_in_sequ;
  while (current_element != NULL)
  {
    el = current_element;
    if (el->value[0] == zero)  /* zero length */
    {
      if (el->flag == 0)  yank(el);
      else if (el->keep_in == 0 && el->npole_sign == 0
               && (el->flag == 1 || el->flag > 4)
               && ident_zero(el) == 0) yank(el);
      else if (el->flag == 3) /* cavity */
      {
        cavity_count++;
        if (cavall_flag == 0 && cavity_count > 1) yank(el);
      }
    }
    else if(el->c_drift > 0)
    {
      af = get_next_name(c, 'd');
      d1 = new_c6t_element(1, c, "drift");
      d1->value[0] = el->value[0];
      d1->flag = 1; d1->position = el->position;
      replace_c6t(el, d1);
      if (af != 0)  add_to_ellist(d1);
    }
    current_element = current_element->next;
  }
}

static void
supp_small_comp(struct c6t_element* p)
{
  int i;
  for (i = 12; i < p->n_values-1; i+=2)
  {
    if (fabs(p->value[i]) > fabs(p->value[i+1]))
    {
      if (fabs(p->value[i+1]) / fabs(p->value[i]) < eps_6)
        p->value[i+1] = zero;
    }
    else if (fabs(p->value[i+1]) > fabs(p->value[i]))
    {
      if (fabs(p->value[i]) / fabs(p->value[i+1]) < eps_6)
        p->value[i] = zero;
    }
  }
}

static void
treat_split(struct c6t_element* el)
{
  if (el->flag == 2)  el->keep_in = 1;
  if (el->nf_err > 0)
  {
    make_multipole(el);
  }
  else
  {
    if (el->force != 0)
      app_factor(el->value[0], &el->value[12], el->n_values-12);
    el->value[0] = zero; /* set element length to zero */
  }
}

static void
yank(struct c6t_element* el)
{
  if (el->previous != NULL)  el->previous->next = el->next;
  else                       first_in_sequ      = el->next;
  if (el->next != NULL)      el->next->previous = el->previous;
//  else                       last_in_sequ       = el->previous; // not used
  el->flag = 0;
}

static void
write_all_el(void)
{
  const char title[] =
    "SINGLE ELEMENTS---------------------------------------------------------";
  f2 = fopen("fc.2", "w");
  fprintf(f2, "%s\n", title);
  current_element = first_in_sequ;
  while (current_element != NULL)
  {
    if (current_element->flag > 0
        && current_element == current_element->equiv
        && current_element->w_flag == 0)
      write_c6t_element(current_element);
    current_element = current_element->next;
  }
  fprintf(f2, "NEXT\n");
}

static void write_rfmultipole(struct c6t_element* el)
{
  if(six_version < general_rf_req)
  {
    const double knl[] = {
      el->value[4],
      el->value[5],
      el->value[6],
      el->value[7]
    };
    const double ksl[] = {
      el->value[18],
      el->value[12],
      el->value[13],
      el->value[14]
    };
    double tilt = el->value[2];
    double freq = el->value[3];
    char name[48];

    if (fabs(knl[0])>eps_9) {
      double lag = 0.25-el->value[8];
      double pc0 = get_value("beam", "pc"); // GeV/c
      el->out_1 = fabs(tilt - M_PI/2)<eps_9 ? -23 : 23; // ID
      el->out_2 = knl[0] * pc0 * 1e3; // rad * GeV/c * 1e3 == rad * MeV/c => MV
      el->out_3 = freq; // freq
      el->out_4 = 2.0 * M_PI * lag; // rad
      strcpy(name, el->name);
      strcat(name, "d");
      fprintf(f2, name_format,
        name, el->out_1, el->out_2, el->out_3, el->out_4, el->out_5, el->out_6, el->out_7);
    }
    if (fabs(knl[1])>eps_9) {
      double lag = -el->value[9];
      el->out_1 = 26; // ID
      el->out_2 = -knl[1]; // 1/m
      el->out_3 = freq; // freq
      el->out_4 = 2.0 * M_PI * lag; // rad
      strcpy(name, el->name);
      strcat(name, "q");
      fprintf(f2, name_format,
        name, el->out_1, el->out_2, el->out_3, el->out_4, el->out_5, el->out_6, el->out_7);
    }
    if (fabs(knl[2])>eps_9) {
      double lag = -el->value[10];
      el->out_1 = 27; // ID
      el->out_2 = -knl[2] / 2.0; // 1/m^2
      el->out_3 = freq; // freq
      el->out_4 = 2.0 * M_PI * lag; // rad
      strcpy(name, el->name);
      strcat(name, "s");
      fprintf(f2, name_format,
        name, el->out_1, el->out_2, el->out_3, el->out_4, el->out_5, el->out_6, el->out_7);
    }
    if (fabs(knl[3])>eps_9) {
      double lag = -el->value[11];
      el->out_1 = 28; // ID
      el->out_2 = -knl[3] / 6.0; // 1/m^3
      el->out_3 = freq; // freq
      el->out_4 = 2.0 * M_PI * lag; // rad
      strcpy(name, el->name);
      strcat(name, "o");
      fprintf(f2, name_format,
        name, el->out_1, el->out_2, el->out_3, el->out_4, el->out_5, el->out_6, el->out_7);
    }
    if (fabs(ksl[0])>eps_9) {
      double lag = -0.25-el->value[19];
      double pc0 = get_value("beam", "pc"); // GeV/c
      el->out_1 = -23; // ID
      el->out_2 = ksl[0] * pc0 * 1e3; // rad * GeV/c * 1e3 == rad * MeV/c => MV
      el->out_3 = freq; // freq
      el->out_4 = 2.0 * M_PI * lag; // rad
      strcpy(name, el->name);
      strcat(name, "ds");
      fprintf(f2, name_format,
        name, el->out_1, el->out_2, el->out_3, el->out_4, el->out_5, el->out_6, el->out_7);
    }
    if (fabs(ksl[1])>eps_9) {
      double lag = -el->value[15];
      el->out_1 = -26; // ID
      el->out_2 = ksl[1]; // 1/m
      el->out_3 = freq; // freq
      el->out_4 = 2.0 * M_PI * lag; // rad
      strcpy(name, el->name);
      strcat(name, "qs");
      fprintf(f2, name_format,
        name, el->out_1, el->out_2, el->out_3, el->out_4, el->out_5, el->out_6, el->out_7);
    }
    if (fabs(ksl[2])>eps_9) {
      double lag = -el->value[16];
      el->out_1 = -27; // ID
      el->out_2 = ksl[2] / 2.0; // 1/m^2
      el->out_3 = freq; // freq
      el->out_4 = 2.0 * M_PI * lag; // rad
      strcpy(name, el->name);
      strcat(name, "ss");
      fprintf(f2, name_format,
        name, el->out_1, el->out_2, el->out_3, el->out_4, el->out_5, el->out_6, el->out_7);
    }
    if (fabs(ksl[3])>eps_9) {
      double lag = -el->value[17];
      el->out_1 = -28; // ID
      el->out_2 = ksl[3] / 6.0; // 1/m^3
      el->out_3 = freq; // freq
      el->out_4 = 2.0 * M_PI * lag; // rad
      strcpy(name, el->name);
      strcat(name, "os");
      fprintf(f2, name_format,
        name, el->out_1, el->out_2, el->out_3, el->out_4, el->out_5, el->out_6, el->out_7);
    }
  }
  else{
    el->out_1 = 41; // ID
    fprintf(f2, name_format,
    el->name, el->out_1, el->out_2, el->out_3, el->out_4, el->out_5, el->out_6, el->out_7);
    write_f3_rfmultipoles(el);
  }
}

static void
write_c6t_element(struct c6t_element* el)
{
  if (strcmp(el->name, "CAV") != 0) {
    if (strcmp(el->base_name, "rfmultipole")==0) {
      write_rfmultipole(el);
    } else {

      fprintf(f2, name_format,
        el->name, el->out_1, el->out_2, el->out_3, el->out_4, el->out_5, el->out_6, el->out_7);
    }
  }
  el->w_flag = 1;
}

static void
write_blocks(void)
{
  char title1[] =
    "BLOCK DEFINITIONS-------------------------------------------------------";
  char title2[] = "1  1";
  struct block* p = first_block;
  int i, lc = 0, nbct = 0;
  double sum = 0;

  fprintf(f2, "%s\n", title1); fprintf(f2, "%s\n", title2);
  while (p != NULL)
  {
    if (p->equiv == p)
    {
      if (p->flag != 0)
      {
        sprintf(p->name, "BLOC%d", ++nbct);
        fprintf(f2, name_format_short, p->name); lc++;
        for (i = 0; i < p->elements->curr; i++)
        {
          if (lc++ == LINES_MAX)
          {
            fprintf(f2,"\n"); fprintf(f2,"                  "); lc = 2;
          }
          fprintf(f2, name_format_short,p->elements->elem[i]->equiv->name);
        }
      }
      if (lc > 0)
      {
        fprintf(f2,"\n"); lc = 0;
      }
    }
    sum += p->length;
    p = p->next;
  }
  printf("\ntotal block length: %f\n", sum);
  fprintf(f2, "NEXT\n");
}

static void
write_f8_errors(void)
{
  double tiltval;
  if (align_cnt == 0)  return;
  current_element = first_in_sequ;
  while (current_element != NULL)
  {
    if (current_element->tilt_err > 0)
    {
      tiltval = current_element->value[6];
    }
    else {tiltval=0.0;}
    if (current_element->na_err > 0)
    {
      if (f8_cnt++ == 0)    f8 = fopen("fc.8", "w");
      fprintf(f8, name_format_4,current_element->equiv->name,
              1000*current_element->p_al_err->a_dble[0],
              1000*current_element->p_al_err->a_dble[1],
              1000*(current_element->p_al_err->a_dble[5]+tiltval));
    }
    else if (current_element->tilt_err > 0)
    {
      if (f8_cnt++ == 0)    f8 = fopen("fc.8", "w");
      fprintf(f8, name_format_4,current_element->equiv->name,
              0.0,
              0.0,
              1000*tiltval);
    }
    current_element = current_element->next;
  }
}

static void
write_f16_errors(void)
{
  int i;
  double factor;

  current_element = first_in_sequ;
  while (current_element != NULL)
  {
    if (current_element->nf_err > 0 && current_element->ref_delta != zero)
    {
      if (f16_cnt++ == 0)    f16 = fopen("fc.16", "w");
      if (current_element->equiv->f3_flag++ == 0)
        write_f3_entry("multipole", current_element->equiv);
      fprintf(f16,"%s\n", current_element->equiv->name);
      for (i = 0; i < current_element->nf_err; i++)
        tmp_buff[i] = current_element->p_fd_err->a_dble[i];

      if(fabs(current_element->value[10])>0){

      if(tmp_buff[0]==999)
        tmp_buff[0] = -(current_element->value[12] - current_element->value[10]); // SixTrack is opposite from MAD-X -> Change of sign compared to TWISS
      else
        tmp_buff[0] = tmp_buff[0] - (current_element->value[12] - current_element->value[10]);
      }

      for (i = current_element->nf_err; i < FIELD_MAX; i++)
        tmp_buff[i] = zero;
      factor = c1p3 / current_element->ref_delta;
      /* commented out because FIELD_MAX in MAD-X is larger than in original c6t */
      /*          for (i = 0; i < FIELD_MAX/2; i++) */
      for (i = 0; i < FIELD_MAX/2-1; i++)
      {
        fprintf(f16, "%23.15e", factor*tmp_buff[2*i]);
        factor *= current_element->ref_radius / (i+1);
        if ((i+1)%3 == 0) fprintf(f16,"\n");
      }
      if (i%3 != 0) fprintf(f16,"\n");
      factor = c1p3 / current_element->ref_delta;
      /*          for (i = 0; i < FIELD_MAX/2; i++) */
      for (i = 0; i < FIELD_MAX/2-1; i++)
      {
        fprintf(f16, "%23.15e", factor*tmp_buff[2*i+1]);
        factor *= current_element->ref_radius / (i+1);
        if ((i+1)%3 == 0) fprintf(f16,"\n");
      }
      if (i%3 != 0) fprintf(f16,"\n");
    }
    current_element = current_element->next;
  }
}

static void
write_f34_special(void)
{
  int i, j, n, flags[FIELD_MAX];
  double values[FIELD_MAX];
  char* t_list[NT34];
  char t_name[NAME_L];
  char* cp;
  double spos=zero,betx=zero,bety=zero,mux=zero,muy=zero;
  int err;

  t_list[0] = &mpole_names[1][0];
  t_list[1] = &mpole_names[2][0];
  t_list[2] = &mpole_names[3][0];
  t_list[3] = &mpole_names[4][0];
  t_list[4] = &mpole_names[5][0];

  if (special_flag == 0)  return;

  n = 0;
  if(f34_cnt++ == 0)    f34 = fopen("fc.34", "w");
  current_element = first_in_sequ;
  while (current_element != NULL)
  {
    for (i = 0; i < NT34; i++)
    {
      if (strcmp(current_element->base_name, t_list[i]) == 0)
      {
        n = f34_values(current_element, flags, values);
        for (j = 0; j < n; j++)
        {
          strcpy(t_name, current_element->name);
          if ((cp = strchr(t_name, '+')) != NULL) *cp = '\0';
          if ((err=double_from_table_row("twiss","s",&(current_element->twtab_row),&spos)))
            printf ("Not found double_from table = %i\n",err);
          if ((err=double_from_table_row("twiss","betx",&(current_element->twtab_row),&betx)))
            printf ("Not found double_from table = %i\n",err);
          if ((err=double_from_table_row("twiss","bety",&(current_element->twtab_row),&bety)))
            printf ("Not found double_from table = %i\n",err);
          if ((err=double_from_table_row("twiss","mux",&(current_element->twtab_row),&mux)))
            printf ("Not found double_from table = %i\n",err);
          if ((err=double_from_table_row("twiss","muy",&(current_element->twtab_row),&muy)))
            printf ("Not found double_from table = %i\n",err);
          fprintf(f34,
                  name_format_error,
                  spos,t_name,flags[j],values[j],betx,bety,mux,muy);
        }
      }
    }
    current_element = current_element->next;
  }
  if (last_in_sequ_org->twtab_row > 0)
  {
    if ((err=double_from_table_row("twiss","s",&(last_in_sequ_org->twtab_row),&spos)))
      printf ("Not found double_from table = %i\n",err);
    if ((err=double_from_table_row("twiss","betx",&(last_in_sequ_org->twtab_row),&betx)))
      printf ("Not found double_from table = %i\n",err);
    if ((err=double_from_table_row("twiss","bety",&(last_in_sequ_org->twtab_row),&bety)))
      printf ("Not found double_from table = %i\n",err);
    if ((err=double_from_table_row("twiss","mux",&(last_in_sequ_org->twtab_row),&mux)))
      printf ("Not found double_from table = %i\n",err);
    if ((err=double_from_table_row("twiss","muy",&(last_in_sequ_org->twtab_row),&muy)))
      printf ("Not found double_from table = %i\n",err);
  }
  fprintf(f34,
          name_format_error,
          spos,"end_marker",100,zero,betx,bety,mux,muy);
}

static void
write_f3_aper(void)
{
  char keyword[3];
  int f3aper_cnt = 0;
  current_element = first_in_sequ;

  while (current_element != NULL)
  {
    if (strstr(current_element->name,"_AP") != NULL
        && (current_element->equiv == current_element))
    {
      if (f3aper_cnt++ == 0)
      {
        f3aper  = fopen("fc.3.aper", "w");
        // fprintf(f3aper,"LIMI\n");
      }
      if      (current_element->value[1] == 1) strcpy(keyword, "CR") ;
      else if (current_element->value[1] == 2) strcpy(keyword, "RE") ;
      else if (current_element->value[1] == 3) strcpy(keyword, "EL") ;
      else if (current_element->value[1] == 4) strcpy(keyword, "RC") ;
      else if (current_element->value[1] == 5) strcpy(keyword, "RL") ;
      else if (current_element->value[1] == 6) strcpy(keyword, "RT") ;
      else if (current_element->value[1] == 7) strcpy(keyword, "OC") ;
      else strcpy(keyword, "UK") ; // unknown aperture type

      fprintf(f3aper,name_format_aper,
        current_element->name, keyword,
        current_element->value[2], current_element->value[3],
        current_element->value[4], current_element->value[5],
              current_element->value[6], current_element->value[7],
        current_element->value[8]);
    }
    current_element = current_element->next;
  }
  // if (f3aper_cnt > 0) fprintf(f3aper,"NEXT\n");
}

static void
write_f3_aux(void)
{
  double aux_val[4] = {-1.e20, -1.e20, -1.e20, -1.e20};
  double tw_alfa;
  int row=1;
  if ((double_from_table_row("summ","q1", &row, &(aux_val[0])) !=0) ||
      (double_from_table_row("summ","q2",  &row, &(aux_val[1])) !=0) ||
      (double_from_table_row("summ","dq1", &row, &(aux_val[2])) !=0) ||
      (double_from_table_row("summ","dq2", &row, &(aux_val[3])) !=0))
  {
    printf("c6t error: tunes or chromaticities not found!\n");
  }
  if (current_beam != NULL)
  {
    if (f3aux_cnt++ == 0) f3aux  = fopen("fc.3.aux", "w");
    if (double_from_table_row("summ","alfa", &row, &tw_alfa) !=0)
      printf("c6t warning: alfa not found in twiss\n");
    fprintf(f3aux, "SYNC\n");
    fprintf(f3aux,"%12.0f%10.6f%10.3f 0.  %12.6f%12.6f  1\n",
            harmon, tw_alfa, total_voltage, sequ_length,
            c1p3*command_par_value("mass", current_beam));
    fprintf(f3aux,"      1.        1.\n");
    fprintf(f3aux, "NEXT\n");
    fprintf(f3aux, "BEAM\n");
    fprintf(f3aux, "%12.4e%14.6g%14.6g%12.4e%12.4e  1  0\n",
            command_par_value("npart", current_beam),
      1e6*command_par_value("exn", current_beam), // um.rad in sixtrack
      1e6*command_par_value("eyn", current_beam),
            command_par_value("sigt", current_beam),
            command_par_value("sige", current_beam));
    fprintf(f3aux, "NEXT\n");
  }
  if (aux_val[0] > -1.e10 && aux_val[1] > -1.e10)
  {
    fprintf(f3aux, "TUNE\n");
    fprintf(f3aux, "QF%23.15f\n", aux_val[0]);
    fprintf(f3aux, "QD%23.15f\n", aux_val[1]);
    fprintf(f3aux, "NEXT\n");
  }
  if (aux_val[2] > -1.e10 && aux_val[3] > -1.e10)
  {
    fprintf(f3aux, "CHRO\n");
    fprintf(f3aux, "SXF%23.15f\n", aux_val[2]);
    fprintf(f3aux, "SXD%23.15f\n", aux_val[3]);
    fprintf(f3aux, "NEXT\n");
  }
}

static void
write_f3_sixmarker(struct element* el){
  int nl;
  char* tmp;
  char tmp2 [500], str_att_value[50]; 
  double attrtemp [20] ;
  char nstr [5], nstr2[5];
  tmp = mymalloc("sixmarker", 500*sizeof *tmp);


  if (!f3) f3 = fopen("fc.3", "w");
    
  tmp =command_par_string("f3string",el->def);
  nl = element_vector(el, "f3vector", attrtemp);
  
  for(int i=0; i<nl; i++){

    sprintf(str_att_value, "%12.8e", attrtemp[i]);
    
    strcpy(nstr2,"{");
    sprintf(nstr, "%d", i);
    strcat(nstr2,nstr);

    strcat(nstr2,"}");
    myrepl(nstr2, str_att_value, tmp, tmp2);
    strcpy(tmp, tmp2);
    
  }

  myrepl("{newline}", "\n", tmp2, tmp);
  fprintf(f3, "%s", tmp);
  fprintf(f3, "%s", "\n");

}

static void
write_f3_matrix(void)
{
  int i, i_max = 43, dim=6;
  current_element = first_in_sequ;
  
  
  double beta, value;

 
  beta= get_value("beam ","beta ");
 
  if (!f3) f3 = fopen("fc.3", "w");

  while (current_element != NULL)
  {
    if (strcmp(current_element->base_name, "matrix") == 0)
    {
      fprintf(f3,"TROM\n");
      fprintf(f3,"%-48s\n",current_element->name);
    
      for (i = 1; i < i_max; i++) {
        value=current_element->value[i];
        // The if statemenst are to go from pt to psigma and from t to sigma.     
        if((i+1)%dim==0){
          value=value*beta;
        }
        if(i%dim==0){
          value=value/beta;
        }
        if(i>(dim+24) && i <(31+dim)){
          value = value/beta;
        }
        if(i>(dim+30) && i < (37+dim)){
          value = value*beta;
        }
        if(i<(dim+1)){
      value = value * 1000;
        }
    

        fprintf(f3,"%23.15e", value);
        if (i%3 == 0) fprintf(f3,"\n");
      }

      fprintf(f3,"NEXT\n");
    }
    current_element = current_element->next;
  }
}

static void
write_f3_wire(void)
{

  current_element = first_in_sequ;
  if (!f3) f3 = fopen("fc.3", "w");
  int isfirst = 0;
  while (current_element != NULL)
  {
    if (strcmp(current_element->base_name, "wire") == 0)
    {

      if(isfirst==0) {
        fprintf(f3,"WIRE\n");
        isfirst =1;
      }

      fprintf(f3,name_format_short,current_element->name);
      fprintf(f3, "%d", (int)current_element->value[1]);
      for(int i=2; i < 9; i++) fprintf(f3,name_format_6, current_element->value[i]);
      fprintf(f3,"\n"); 
    }
     current_element = current_element->next;
  }
    if(isfirst >0) fprintf(f3,"NEXT\n");
}
   



static void
write_f3_entry(const char* option, struct c6t_element* el)
{
  if (f3_cnt++ == 0 && !f3) f3 = fopen("fc.3", "w");
  if (strcmp(option, "multipole") == 0) write_f3_mult(el);
}

static void
write_f3_mult(struct c6t_element* el)
{
  int i, j, i_max = -1;
  struct c6t_element* eln;
  if (multi_type < 0)  return;
  fprintf(f3,"MULT\n");
  fprintf(f3,name_format_3, el->name, c1p3*el->ref_radius,
          el->ref_delta);
  /* find non-zero errors in all elements equiv. to this, print error matrix */
  for (i = 0; i < FIELD_MAX; i++) error_matrix[i] = zero;
  for (j = 0; j < types.member[multi_type]->curr; j++)
  {
    eln = types.member[multi_type]->elem[j];
    if (eln->equiv == el)
    {
      for (i = 0; i < eln->nf_err; i++)
      {
        if (mult_auto_off)
        {
          error_matrix[i] = 1.;
        }
        else
        {
          if (eln->p_fd_err->a_dble[i] != zero)
          {
            i_max = i; error_matrix[i] = 1.;
          }
        }
      }
    }
  }
  if (mult_auto_off)
  {
    i_max = max_mult_ord * 2;
  }
  else
  {
    if (++i_max > 0)  i_max += i_max%2;
  }
  for (i = 0; i < i_max; i++)
  {
    fprintf(f3,"%4.0f.%4.0f.", 0., error_matrix[i]);
    if ((i+1)%2 == 0) fprintf(f3,"\n");
  }
  fprintf(f3,"NEXT\n");
}

static void
write_f3_rfmultipoles(struct c6t_element* current_element)
{

  if (!f3) f3 = fopen("fc.3", "w");

  if (strcmp(current_element->base_name, "rfmultipole") == 0)
  {
    fprintf(f3,"RFMULTIPOLE\n");
    fprintf(f3, "%s %f \n", current_element->name,current_element->value[2]);
    
    for (int i=0; i < current_element->value[3]; i++){
      fprintf(f3, name_format_5, current_element->value[i*4+7], current_element->value[i*4+8],
        current_element->value[i*4+9], current_element->value[i*4+10]);
    }
    fprintf(f3,"NEXT\n");
  }
}
static void
rfmultipole_name(char *name, struct c6t_element* el)
{
  const double knl[] = {
    el->value[4],
    el->value[5],
    el->value[6],
    el->value[7]
  };
  const double ksl[] = {
    el->value[18],
    el->value[12],
    el->value[13],
    el->value[14]
  };
  char tmp[256] = "";
  int n = 0;
  if (fabs(knl[0])>eps_9) {
    strcpy(tmp, el->name);
    strcat(tmp, "d");
    n += sprintf(name+n, name_format_short, tmp);
  }
  if (fabs(knl[1])>eps_9) {
    strcpy(tmp, el->name);
    strcat(tmp, "q");
    n += sprintf(name+n, name_format_short, tmp);
  }
  if (fabs(knl[2])>eps_9) {
    strcpy(tmp, el->name);
    strcat(tmp, "s");
    n += sprintf(name+n, name_format_short, tmp);
  }
  if (fabs(knl[3])>eps_9) {
    strcpy(tmp, el->name);
    strcat(tmp, "o");
    n += sprintf(name+n, name_format_short, tmp);
  }
  if (fabs(ksl[0])>eps_9) {
    strcpy(tmp, el->name);
    strcat(tmp, "ds");
    n += sprintf(name+n, name_format_short, tmp);
  }
  if (fabs(ksl[1])>eps_9) {
    strcpy(tmp, el->name);
    strcat(tmp, "qs");
    n += sprintf(name+n, name_format_short, tmp);
  }
  if (fabs(ksl[2])>eps_9) {
    strcpy(tmp, el->name);
    strcat(tmp, "ss");
    n += sprintf(name+n, name_format_short, tmp);
  }
  if (fabs(ksl[3])>eps_9) {
    strcpy(tmp, el->name);
    strcat(tmp, "os");
    n += sprintf(name+n, name_format_short, tmp);
  }
}

static void
write_struct(void)
{
  struct block* p = first_block;
  int lc = 0;
  char title[] =
    "STRUCTURE INPUT---------------------------------------------------------";

  fprintf(f2, "%s\n", title);
  if(multicol_flag == 1) {
    fprintf(f2, "%s\n", "MULTICOL");
  }
  while (p != NULL)
  {
    char name[256] = "";

    if (p->flag == 0) {
      if (strcmp(p->first->equiv->base_name,"rfmultipole") == 0 && six_version < general_rf_req)
        rfmultipole_name(name, p->first->equiv);
      else
        strcpy(name, p->first->equiv->name);
    }
    else
      strcpy(name,p->equiv->name);

    if(multicol_flag == 1) {
      if(long_names_flag==1) {
        fprintf(f2, "%-48s %-48s %17.9f\n", p->first->name, name, p->first->position);
      } else {
        fprintf(f2, "%-20s %-20s %17.9f\n", p->first->name, name, p->first->position);
      }
    } else {
      if (lc++ == LINES_MAX) {
        fprintf(f2,"\n"); lc = 1;
      }
      if(long_names_flag==1) {
        fprintf(f2, "%-48s ", name);
      } else {
        fprintf(f2, "%-17s ", name);
      }
    }
    
    p = p->next;
  }
  if (lc > 0)
  {
    fprintf(f2,"\n");
  }
  fprintf(f2, "NEXT\n");
}

static int
my_table_row(struct table* table, char* name)
{
  int i, j, ret = 0;
  char t_name[255];
  char* cp;
  for (i = 0; i < table->num_cols; i++)
    if(table->columns->inform[i] == 3) break;
  if (i < table->num_cols && last_row < table->curr)
  {
    for (j = last_row; j < table->curr; j++)
    {
      strcpy(t_name, table->s_cols[i][j]);
      if ((cp = strchr(t_name, ':')) != NULL) *cp = '\0';
      if (strcmp(name, t_name) == 0) break;
    }
    if (j < table->curr)
    {
      ret = j+1;
      last_row = j+1;
    }
  }
  return ret;
}

static void
c6t_finish(void)
{
  const char *rout_name = "c6t_finish";
  int i,j;
  struct block* p;
  /* remove elements and elements list */
  for(i=0; i<types.curr; i++)
  {
    for(j=0; j<types.member[i]->curr; j++)
    {
      if (types.member[i]->elem[j]->value)
        myfree(rout_name, types.member[i]->elem[j]->value);
      if (types.member[i]->elem[j]->p_al_err &&
          types.member[i]->elem[j]->do_not_free != 1)
      {
        if (types.member[i]->elem[j]->p_al_err->a_dble)
          myfree(rout_name, types.member[i]->elem[j]->p_al_err->a_dble);
        myfree(rout_name, types.member[i]->elem[j]->p_al_err);
        types.member[i]->elem[j]->p_al_err = NULL;
      }
      if (types.member[i]->elem[j]->p_fd_err &&
          types.member[i]->elem[j]->do_not_free != 1)
      {
        if (types.member[i]->elem[j]->p_fd_err->a_dble)
          myfree(rout_name, types.member[i]->elem[j]->p_fd_err->a_dble);
        myfree(rout_name, types.member[i]->elem[j]->p_fd_err);
        types.member[i]->elem[j]->p_fd_err = NULL;
      }
      myfree(rout_name, types.member[i]->elem[j]);
      types.member[i]->elem[j]=NULL;
    }
    myfree(rout_name, types.member[i]);
  }
  types.curr=0; first_in_sequ = NULL; last_in_sequ_org = NULL; // last_in_sequ = NULL; // not used
  current_element=NULL;
  /* remove blocks */
  p = first_block;
  while (p != NULL)
  {
    p = p->next;
    if (p) myfree(rout_name, p->previous);
  }
  first_block = NULL; prev_block=NULL; // last_block=NULL; not used
  current_block = NULL;
  /* remove split_list */
  if (split_list)
  {
    myfree(rout_name, split_list); split_list = NULL;
  }
  /* clear acro_cnt and acro_list */
  for(i=0; i<20; i++)
  {
    acro_list[i]='\0';
    acro_cnt[i]=0;
  }
  /* remember that this is not the first time we run */
  virgin_c6t=0;

  /* added by LD 2011-10-18 */
  if (f2) { fclose(f2); f2 = 0; }
  if (f3) { fclose(f3); f3 = 0; }
  if (f3aux) { fclose(f3aux); f3aux = 0; }
  if (f3aper) { fclose(f3aper); f3aper = 0; }
  if (f8) { fclose(f8); f8 = 0; }
  if (f16) { fclose(f16); f16 = 0; }
  if (f34) { fclose(f34); f34 = 0; }
}

static void
c6t_init(void)
{
  int j;
  const char *rout_name = "c6t_init";

  if (virgin_c6t)
  {
    p_err_zero = make_obj("zero_errors", 0, FIELD_MAX, 0, 0);
    for (j = 0; j < FIELD_MAX; j++)
      p_err_zero->a_dble[j]=0.0;

    for (j = 0; j < N_TYPES; j++) {
      t_info[j] = mymalloc(rout_name, sizeof *t_info[0]);
      sscanf(el_info[j],"%s%d%d%d%d%d%d",t_info[j]->name, &t_info[j]->flag_1,
             &t_info[j]->flag_2, &t_info[j]->flag_3, &t_info[j]->flag_4,
             &t_info[j]->flag_5, &t_info[j]->flag_6);
    }
  }
  if (current_sequ == NULL)           fatal_error("c6t - no current sequence.","");
  if (current_sequ->ex_start == NULL) fatal_error("c6t - sequence not expanded.","");
  if (current_sequ->tw_table == NULL) fatal_error("c6t - twiss table not found.","");
  if (attach_beam(current_sequ) == 0) fatal_error("c6t - sequence without beam command.","");

  /* initialise everything */
  block_count = 0;     /* current block count for naming */
  elem_cnt = 0;        /* element count */
  acro_occ = 0;        /* acro list occupation */
  align_cnt = 0;       /* element with align errors count */
  field_cnt = 0;       /* element with field errors count */
  f3_cnt = 0;          /* f3 write flag */
  f3aux_cnt = 0;       /* f3aux write flag */
  f8_cnt = 0;          /* f8 write count */
  f16_cnt = 0;         /* f16 write count */
  f34_cnt = 0;         /* f34 write count */
  special_flag = 1;    /* produce special output file from twiss */
  aperture_flag = 0;   /* if 1 insert apertures into structure */
  cavall_flag = 0;     /* if 0 dump all cavities into first */
  markall_flag = 0;    /* if 0 dump all marker into first */
  long_names_flag = 0; 
//  radius_flag = 0; // not used    /* change the default reference radius */
  split_flag = 0;      /* if 1 keep zero multipoles after split */
  mangle_flag = 0;     /* if 1 truncate to 14 chars and mangle names */
  multi_type = -1;     /* is set to multipole type if any found */
  cavity_count = 0;    /* count cavities in output */

  total_voltage = 0;
  harmon = 0;

  /* added by LD 2011-10-18 */
  f2 = 0;
  f3 = 0;
  f3aux = 0;
  f3aper = 0;
  f8 = 0;
  f16 = 0;
  f34 = 0;
}

static void
process_c6t(void)  /* steering routine */
{
  read_sequ();   /* from db read sequence, store all elements */
  add_c6t_drifts();
  conv_elem();   /* tag elements */
  split();       /* convert to thin */
  multi_loop();
  supp_elem();   /* suppress/replace zero force, most markers,
                    and possibly some cavities */
  concat_drifts();
  get_multi_refs();  /* get multipole flag */
  equiv_elem();  /* find first equivalent for all elements */
  post_multipoles();  /* give errors to all equiv. multipoles */
  block_it();    /* group linear elements into blocks */
  assign_att();  /* assign attributes + errors to all single elements */
  mod_errors();  /* flip normal components */

  setup_output_string(); /* sets the output format used for the out put files */


  write_all_el();
  write_blocks();
  write_struct();
  write_f16_errors();
  write_f34_special();
  write_f3_aux();
  write_f3_matrix();
  write_f3_wire();
  write_f3_aper();
  write_f8_errors();
}

// public interface
static void
setup_output_string(void)
{
  if(long_names_flag==1){
    strcpy(name_format,"%-48s %3d  %23.15e %23.15e  %23.15e  %23.15e  %23.15e  %23.15e\n");
    strcpy(name_format_short,"%-48s" );
    strcpy(name_format_error, " %23.15e  %-48s %3d %23.15e %23.15e %23.15e %23.15e %23.15e\n");
    strcpy(name_format_aper, "%-48s   %s %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n");
    strcpy(name_format_3,  "%-48s%20.10e%20.10e\n");
    strcpy(name_format_4, "%-48s  %14.6e%14.6e%17.9e\n");
    strcpy(name_format_5, "%23.15e %23.15e %23.15e %23.15e\n");
    strcpy(name_format_6, "%23.15e");
  }
    else{
    strcpy(name_format,"%-16s %3d  %16.9e %17.9e  %17.9e  %17.9e  %17.9e  %17.9e\n");
    strcpy(name_format_short, "%-18s");
    strcpy(name_format_error, " %20.13e  %-16s %3d %20.13e %20.13e %20.13e %20.13e %20.13e\n");
    strcpy(name_format_aper, "%-16s   %s %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n");
    strcpy(name_format_3,"%-16s%20.10e%20.10e\n");
    strcpy(name_format_4,"%-16s  %14.6e%14.6e%17.9e\n");
    strcpy(name_format_5, "%17.9e %17.9e %17.9e %17.9e\n");
    strcpy(name_format_6, "%17.9e");

  }
}
void
conv_sixtrack(struct in_cmd* mycmd) /* writes sixtrack input files from MAD-X */
{
  last_row = 0;

  puts("  ++++++++++++++++++++++++++++"  );
  puts("  +   c6t version 2.0        +"  );
  puts("  ++++++++++++++++++++++++++++\n");

  c6t_init();
  get_args(mycmd);
  process_c6t();
  printf("\nc6t terminated - total number of elements: %d\n", elem_cnt);
  printf("                    field errors    MAD-X: %d\n", field_cnt);
  printf("                    field errors SixTrack: %d\n", f16_cnt);
  printf("                 alignment errors   MAD-X: %d\n", align_cnt);
  printf("                alignment errors SixTrack: %d\n", f8_cnt);
  printf("                          sequence length: %f [m]\n", sequ_length);
  c6t_finish();
}