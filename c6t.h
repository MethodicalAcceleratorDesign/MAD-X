#define BASE_TYPES 100    /* maximum no. of element types allowed */
#define EL_COUNT   100    /* initial length of element type list */

/* MADX name     : circle, ellipse, rectangle, lhcscreen */
/* internal code :    1       2         3          4     */
struct aper_struct {
  int apply;
  char name[255];
  char style[255];
  double value[3];
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
  double* value;           /* element strength etc. values */
  struct object* p_al_err; /* pointer to alignment error object */
  struct object* p_fd_err; /* pointer to field error object */
  int tilt_err;            /* allow write_f8 to dump tilt as well */
  int do_not_free;         /* avoid free crash */
};

struct c6t_el_list /* contains list of element pointers */
{
  int max,                /* max. pointer array size */
      curr;               /* current occupation */
  char base_name[48];
  struct c6t_element** elem; /* element pointer list */
};

struct block
{
  char name[48];
  double length;
  int flag;              /* if 0 take element, else block */
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
  int      flag_1,  /* for length = 0 */
           flag_2,  /* for length > 0, normal */
           flag_3,  /* for length > 0, skew */
           flag_4,  /* if > 0: make drift, print warning when encountered */
           flag_5,  /* if > 0: magnet (for k0n * l) */
           flag_6;  /* if length > 0: 0 = no split
                       1 = split; if flag_2(_3) = 1: two identical + zero m.
                                  if flag_2(_3) = 2: two drift + full m. */
};

/* (some) constants and structure definitions for DOOM */

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

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

  int* a_int;             /* integer array */
  double* a_dble;         /* d.p. array */
  char* a_char;           /* string */
  struct object* parent;  /* pointer to parent object */
  struct object** p_obj;  /* object pointer array */
  char** names;           /* name pointers into a_char */
};

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif
