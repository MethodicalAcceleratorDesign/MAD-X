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

/* already defined as 42 in fulll.h */
/*#define FIELD_MAX 40*/        /* field error array length */
#define KEY_LENGTH 48       /* from DOOM */
#define MM_KEEP 2           /* no. of element name starts to keep */
#define N_TYPES 27          /* no. of valid element types */
#define MULTI_MAX 24        /* element array length for multipoles */
#define NT34 5              /* no. of element types in special fort.34 */
#define LINES_MAX 3         /* structure output line max. names */
#define SEQ_DUMP_LEVEL 0    /* chooses amount of dumped output */

void add_c6t_drifts();
void add_split_list(struct c6t_element*);
void add_to_ellist(struct c6t_element*);
void app_factor(double, double*, int);
void arr_print(double*, int);
void assign_att();
void att_aperture(struct c6t_element*);
void att_beambeam(struct c6t_element*);
void att_colli(struct c6t_element*);
void att_decapole(struct c6t_element*);
void att_drift(struct c6t_element*);
int f34_values(struct c6t_element*, int*, double*);
void att_hkicker(struct c6t_element*);
void att_kicker(struct c6t_element*);
void att_lcavity(struct c6t_element*);
void att_marker(struct c6t_element*);
void att_matrix(struct c6t_element*);
void att_multipole(struct c6t_element*);
void att_octupole(struct c6t_element*);
void att_quadrupole(struct c6t_element*);
void att_rbend(struct c6t_element*);
void att_rfcavity(struct c6t_element*);
void att_sbend(struct c6t_element*);
void att_sextupole(struct c6t_element*);
void att_vkicker(struct c6t_element*);
void att_undefined(struct c6t_element*);
void clean_c6t_element(struct c6t_element*);
struct c6t_element* create_aperture(char* ,char* ,int , int , struct double_array*);
void concat_drifts();
void conv_elem();
void c6t_finish();
void c6t_init();
struct c6t_element* convert_madx_to_c6t(struct node*);
void dump_c6t_element(struct c6t_element*);
void dump_c6t_sequ(int);
void dump_types(int);
void equiv_elem();
struct block* get_block_equiv(struct block*);
void get_args(struct in_cmd*);
void get_error_refs(struct c6t_element*);
int get_flag(struct c6t_element*, struct type_info*);
struct c6t_element* get_from_ellist(char*, char*);
void get_multi_refs();
int get_next_name(char*, char);
void gnu_file(struct c6t_element*);
void grow_ellist(struct c6t_el_list*);
int ident_el(struct c6t_element*, struct c6t_element*);
int ident_zero(struct c6t_element*);
int in_keep_list(struct c6t_element*);
void invert_normal(int, double*);
void invert_skew(int, double*);
void link_behind(struct c6t_element*, struct c6t_element*);
void link_c6t_in_front(struct c6t_element*, struct c6t_element*);
void lower(char*);
struct c6t_element* make_c6t_element(struct node*);
struct object* make_obj(char*, int, int, int, int);
void make_multipole(struct c6t_element*);
void mod_errors();
void mod_lcavity(struct c6t_element*);
void mod_multipole(struct c6t_element*);
void mod_octupole(struct c6t_element*);
void mod_quadrupole(struct c6t_element*);
void mod_rbend(struct c6t_element*);
void mod_rfcavity(struct c6t_element*);
void mod_sextupole(struct c6t_element*);
void multi_loop();
struct c6t_element* new_c6t_element(int, char*, char*);
struct block* new_block();
void post_multipoles();
double power_of(double, int);
void pre_multipole(struct c6t_element*);
void pro_elem(struct node*);
void process_c6t();
void read_sequ();
void remove_from_ellist(struct c6t_element*);
void replace_c6t(struct c6t_element*, struct c6t_element*);
void split();
void split_kicker(struct c6t_element*);
void split_other(struct c6t_element*);
void split_special(struct c6t_element*);
void supp_elem();
void supp_small_comp(struct c6t_element*);
void treat_split(struct c6t_element*);
void yank(struct c6t_element*);
void write_all_el();
void write_blocks();
void write_c6t_element(struct c6t_element*);
void write_f16_errors();
void write_f8_errors();
void write_f3_aper();
void write_f3aux();
void write_f3_entry(char*, struct c6t_element*);
void write_f3_mult(struct c6t_element*);
void write_f34_special();
void write_struct();

/* routines used from makethin.c */
double el_par_value_recurse(char*, struct element*);
struct command_parameter* return_param_recurse(char*, struct element* );
struct command_parameter* return_param(char* , struct element* );

struct li_list types;

struct type_info* t_info[N_TYPES];

struct block   *first_block, *last_block;
struct block*   prev_block;
struct block*   current_block = NULL;

int virgin_c6t = 1;

struct c6t_element *first_in_sequ, *last_in_sequ;
struct c6t_element* prev_element;
struct c6t_element* current_element = NULL;
struct c6t_element* debug_element = NULL;
struct c6t_el_list* split_list = NULL;
struct aper_struct tag_aperture;

struct object *p_err_zero;  /* pointer to error object with all zeroes */
              
char el_info[N_TYPES][60] = /* see type_info definition */
/*           l=0 l>0,normal l>0,skew ->drift make_k*l split */
{"aperture    2       2       2       0       0       0",
"beambeam     2       2       2       0       0       0",
"beamint      0       1       1       1       0       0",
"drift        0       1       1       0       0       0",
"decapole     2       2       2       0       1       2", 
"ecollimator  2       2       1       0       0       0",
"elseparator  0       1       1       1       0       0",
"gbend        1       1       1       2       1       1",
"hkicker      5       5       5       1       0       3",
"hmonitor     0       1       1       1       0       0",
"instrument   0       1       1       1       0       0",
"kicker       6       6       6       1       0       3",
"lcavity      3       3       3       0       0       2",
"marker       4       0       0       0       0       0",
"matrix       2       2       2       0       0       0",
"monitor      0       1       1       1       0       0",
"multipole    2       2       2       0       0       0",
"octupole     2       2       2       0       1       2",
"quadrupole   2       1       2       0       1       1",
"rbend        2       1       1       0       1       1",
"rcollimator  2       2       1       0       0       0",
"rfcavity     3       3       3       0       0       2",
"sbend        2       1       1       0       1       1",
"sextupole    2       2       2       0       1       2",
"solenoid     0       1       1       2       1       0",
"vkicker      5       5       5       1       0       3",
"vmonitor     0       1       1       1       0       0"};

char keep_these[MM_KEEP][24] = {"ip", "mt_"};
char mpole_names[][16] = {"dipole", "quadrupole", "sextupole",
                          "octupole", "decapole", "multipole"};
char acro_list[20];   /* list for name starts */
int acro_cnt[20];    /* counters for name starts */
char tmp_name[KEY_LENGTH];

int           block_count = 0,     /* current block count for naming */
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
              cavall_flag = 0,     /* if 0 lump all cavities into first */
              aperture_flag = 0,   /* if 1 insert apertures into structure */
              radius_flag = 0,     /* change the default reference radius */
              split_flag = 0,      /* if 1 keep zero multipoles after split */
              multi_type = -1,     /* is set to multipole type if any found */
              cavity_count = 0;    /* count cavities in output */

double        sequ_length,         /* length of  sequence */
              sequ_start, 
              sequ_end,
              total_voltage = 0,
              harmon = 0,
              error_matrix[FIELD_MAX],
              tmp_buff[FIELD_MAX];

const double ten   = 10;
const double c1p3 = 1.e3;
const double eps_6 = 1.e-6;
const double eps_9 = 1.e-9;
const double eps_12 = 1.e-12;
double ref_def = 0.017;

FILE *f2, *f3, *f3aux, *f3aper, *f8, *f16, *f34;

void conv_sixtrack(struct in_cmd* mycmd) /* writes sixtrack input files from MAD-X */
{
  puts("  ++++++++++++++++++++++++++++");
  puts("  +   c6t version 2.0        +");
  puts("  ++++++++++++++++++++++++++++\n");

  c6t_init();
  get_args(mycmd);
  process_c6t();
  printf("\nc6t terminated - total number of elements: %d\n", elem_cnt);
  printf("                    with alignment errors: %d\n", align_cnt);
  printf("                    with field     errors: %d\n", field_cnt);
  printf("                          sequence length: %f [m]\n", sequ_length);
  c6t_finish();
}

void add_c6t_drifts()
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
     else if (dl > eps_9)
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

void add_split_list(struct c6t_element* el)
{
  int i;
  char rout_name[] = "c6t:add_split_list";
  if (split_list == NULL) 
    {
     split_list = (struct c6t_el_list*) mycalloc(rout_name,1, sizeof(struct c6t_el_list));
     split_list->elem = 
       (struct c6t_element**) mycalloc(rout_name,EL_COUNT, sizeof(struct elem*));
     split_list->max = EL_COUNT;
    }
  else if (split_list->curr == split_list->max) grow_ellist(split_list);
  for (i = 0; i < split_list->curr; i++) if (split_list->elem[i] == el) return;
  split_list->elem[split_list->curr++] = el;
}

void add_to_ellist( /* adds element to correct object list */
		 struct c6t_element* p_elem)
{
  int j;
  char rout_name[] = "c6t:add_to_ellist";

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
  types.member[types.curr] 
       = (struct c6t_el_list*) mycalloc(rout_name,1,sizeof(struct c6t_el_list));
  types.member[types.curr]->elem
       = (struct c6t_element**) mycalloc(rout_name,EL_COUNT, sizeof(struct c6t_element*));
  types.member[types.curr]->elem[types.member[types.curr]->curr++] = p_elem;
  types.member[types.curr]->max = EL_COUNT;
  strcpy(types.member[types.curr]->base_name, p_elem->base_name);
  types.curr++;
}

void app_factor(double fact, double* array, int count)
{
  int i;
  for (i = 0; i < count; i++) array[i] *= fact;
}

void arr_print(double array[], int occ)
{
  int i;
  for (i = 0; i < occ; i++)
    {
     printf(" %12.4e", array[i]); if ((i+1)%5 == 0) printf("\n");
    }
  printf("\n");
}

void assign_att()
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
    	 if (strcmp(el->base_name, "aperture") == 0) att_aperture(el);
	 else if (strcmp(el->base_name, "beambeam") == 0) att_beambeam(el);
    	 else if (strcmp(el->base_name, "decapole") == 0) att_decapole(el);
    	 else if (strcmp(el->base_name, "drift") == 0) att_drift(el);
    	 else if (strcmp(el->base_name, "ecollimator") == 0) att_colli(el);
    	 else if (strcmp(el->base_name, "hkicker") == 0) att_hkicker(el);
    	 else if (strcmp(el->base_name, "kicker") == 0) att_kicker(el);
    	 else if (strcmp(el->base_name, "lcavity") == 0) att_lcavity(el);
    	 else if (strcmp(el->base_name, "marker") == 0) att_marker(el);
    	 else if (strcmp(el->base_name, "matrix") == 0) att_matrix(el);
    	 else if (strcmp(el->base_name, "multipole") == 0) att_multipole(el);
    	 else if (strcmp(el->base_name, "octupole") == 0) att_octupole(el);
    	 else if (strcmp(el->base_name, "quadrupole") == 0) att_quadrupole(el);
    	 else if (strcmp(el->base_name, "rbend") == 0) att_rbend(el);
    	 else if (strcmp(el->base_name, "rcollimator") == 0) att_colli(el);
    	 else if (strcmp(el->base_name, "rfcavity") == 0) att_rfcavity(el);
    	 else if (strcmp(el->base_name, "sbend") == 0) att_sbend(el);
    	 else if (strcmp(el->base_name, "sextupole") == 0) att_sextupole(el);
    	 else if (strcmp(el->base_name, "vkicker") == 0) att_vkicker(el);
    	 else att_undefined(el);
       }
     }
   }
}

void att_aperture(struct c6t_element* el)
{
  el->out_1 = 3;
  el->out_2 = 1e-8;
  el->out_3 = 0.0;
  el->out_4 = 0.0;
}

void att_beambeam(struct c6t_element* el)
{
  double beamx,beamy;
  if (double_from_table("twiss","x",&(el->twtab_row),&beamx) != 0 || 
      double_from_table("twiss","y",&(el->twtab_row),&beamy) != 0) {
    warning("c6t: beambeam element not found in twiss table","");
  }
  el->out_1 = 20; 
  el->out_2 = c1p3*(el->value[12] - beamx);
  el->out_3 = c1p3*(el->value[13] - beamy);
  el->out_4 = el->value[16];
}

void att_colli(struct c6t_element* el)
     /* ecollim. + rcollim. - make drift, do not concatenate */
{
  el->out_1 = 0; el->out_4 = el->value[0];
}

void att_decapole(struct c6t_element* el)
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

void att_drift(struct c6t_element* el)
{
  el->out_4 = el->value[0];
}

void att_hkicker(struct c6t_element* el)
{
 el->out_1 = 1; el->out_2 = el->value[12];
}

void att_kicker(struct c6t_element* el)
{
}

void att_lcavity(struct c6t_element* el)
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

void att_marker(struct c6t_element* el)
{
}

void att_matrix(struct c6t_element* el)
{
}

void att_multipole(struct c6t_element* el)
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

void att_octupole(struct c6t_element* el)
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

void att_quadrupole(struct c6t_element* el)
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

void att_rbend(struct c6t_element* el)
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

void att_rfcavity(struct c6t_element* el)
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

void att_sbend(struct c6t_element* el)
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

void att_sextupole(struct c6t_element* el)
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

void att_vkicker(struct c6t_element* el)
{
 el->out_1 = -1; el->out_2 = el->value[13];
}

void att_undefined(struct c6t_element* el)
{
  el->out_4 = el->value[0];
}

void block_it()
{
  struct c6t_element* el;
  char rout_name[] = "c6t:block_it";

  current_element = first_in_sequ;
  while ((el = current_element) != NULL)
    {
     current_block = new_block();
     current_block->previous = prev_block;
     current_block->next = NULL;
     if (prev_block == NULL) first_block = current_block;
     else                    prev_block->next = current_block;
     current_block->elements
       = (struct c6t_el_list*) mycalloc(rout_name,1,sizeof(struct c6t_el_list));
     current_block->elements->elem
       = (struct c6t_element**) mycalloc(rout_name,EL_COUNT, sizeof(struct c6t_element*));
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
   last_block = current_block;
}

void concat_drifts()
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

void conv_elem()
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

void c6t_finish()
{
  int i,j;
  struct block* p;
  /* remove elements and elements list */
  for(i=0; i<types.curr; i++) {
    for(j=0; j<types.member[i]->curr; j++) {
      if (types.member[i]->elem[j]->value) 
	free(types.member[i]->elem[j]->value);
      if (types.member[i]->elem[j]->p_al_err) {
	if (types.member[i]->elem[j]->p_al_err->a_dble)
	  free(types.member[i]->elem[j]->p_al_err->a_dble);
	free(types.member[i]->elem[j]->p_al_err);
      }
      if (types.member[i]->elem[j]->p_fd_err) {
	if (types.member[i]->elem[j]->p_fd_err->a_dble)
	  free(types.member[i]->elem[j]->p_fd_err->a_dble);
	free(types.member[i]->elem[j]->p_fd_err);
      }
      free(types.member[i]->elem[j]);
      types.member[i]->elem[j]=NULL;
    }
    free(types.member[i]);
  }
  types.curr=0; first_in_sequ = NULL; last_in_sequ = NULL;
  current_element=NULL;
  /* remove blocks */
  p = first_block;
  while (p != NULL)  
    {
     p = p->next;
     if (p) free(p->previous);
    }
  first_block = NULL; last_block=NULL; prev_block=NULL;
  current_block = NULL;
  /* remove split_list */
  if (split_list) {
    free(split_list); split_list = NULL;
  }
  /* clear acro_cnt and acro_list */
  for(i=0; i<20; i++) {
    acro_list[i]='\0';
    acro_cnt[i]=0;
  }
  /* remember that this is not the first time we run */
  virgin_c6t=0;
}

void c6t_init()
{
  int j;
  char rout_name[] = "c6t_init";

  if (virgin_c6t) {
    p_err_zero = make_obj("zero_errors", 0, FIELD_MAX, 0, 0);
    for (j = 0; j < FIELD_MAX; j++)
      {
	p_err_zero->a_dble[j]=0.0;
      }

    for (j = 0; j < N_TYPES; j++)
      {
	t_info[j] = (struct type_info*) mymalloc(rout_name,sizeof(struct type_info));
	sscanf(el_info[j],"%s%d%d%d%d%d%d",t_info[j]->name, &t_info[j]->flag_1,
	       &t_info[j]->flag_2, &t_info[j]->flag_3, &t_info[j]->flag_4,
	       &t_info[j]->flag_5, &t_info[j]->flag_6);
      }
  }
  if (current_sequ == NULL) 
    fatal_error("c6t - no current sequence.","");
  if (current_sequ->ex_start == NULL) 
    fatal_error("c6t - sequence not expanded.","");
  if (current_sequ->tw_table == NULL) 
    fatal_error("c6t - twiss table not found.","");
  if (attach_beam(current_sequ) == 0) 
    fatal_error("c6t - sequence without beam command.","");

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
  cavall_flag = 0;     /* if 0 lump all cavities into first */
  radius_flag = 0;     /* change the default reference radius */
  split_flag = 0;      /* if 1 keep zero multipoles after split */
  multi_type = -1;     /* is set to multipole type if any found */
  cavity_count = 0;    /* count cavities in output */
  
  total_voltage = 0;
  harmon = 0;
}

void clean_c6t_element(struct c6t_element* cleanme)
{
  int i;
  for(i=0; i<cleanme->n_values; i++) { cleanme->value[i]=0; }
}

struct c6t_element* create_aperture(char* name,char* type,int a, int b, struct double_array* p_al_err)
{
  struct c6t_element* aper_element;
  aper_element = new_c6t_element(4,name,"aperture");
  clean_c6t_element(aper_element);
  strcpy(aper_element->org_name,name);
  aper_element->value[0] = 0.0;
  aper_element->value[1] = a;
  aper_element->value[2] = b;
  if (strcmp(type,"RE")==0) {
    aper_element->value[3] = 1;
  } else {
    aper_element->value[3] = 2;
  }
  aper_element->keep_in=1;
  /* alignment errors of aperture are to be copied toalignment errors of element */
  if (p_al_err && p_al_err->curr>11) {
    align_cnt++;
    aper_element->na_err = p_al_err->curr;
    aper_element->p_al_err = make_obj("ALDUM",0,ALIGN_MAX,0,0);
    aper_element->p_al_err->c_dble = p_al_err->curr;
    aper_element->p_al_err->a_dble[0] = p_al_err->a[10];
    aper_element->p_al_err->a_dble[1] = p_al_err->a[11];
  }
  return aper_element;
}

struct c6t_element* convert_madx_to_c6t(struct node* p)
{
  struct command_parameter *kn_param = NULL,*ks_param = NULL,*aper_param = NULL;
  int i,j,maxkn,maxks;
  struct c6t_element* c6t_elem = NULL;
  char t_name[255];
  char* cp;
  int index=-1;

  strcpy(t_name, p->name);
  if ((cp = strchr(t_name, ':')) != NULL) *cp = '\0';
  if ((strcmp(p->base_name,"rbend") == 0)      ||
      (strcmp(p->base_name,"sbend") == 0)      ||
      (strcmp(p->base_name,"quadrupole") == 0) ||
      (strcmp(p->base_name,"sextupole") == 0)  ||
      (strcmp(p->base_name,"octupole") == 0)   ||
      (strcmp(p->base_name,"vkicker") == 0)    ||
      (strcmp(p->base_name,"hkicker") == 0)    ||
      (strcmp(p->base_name,"kicker") == 0)) {
    c6t_elem = new_c6t_element(19,t_name,p->base_name);
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
    if (c6t_elem->value[12] == zero && c6t_elem->value[0] > zero) {
      c6t_elem->value[12] = c6t_elem->value[10]/c6t_elem->value[0];
    }
    c6t_elem->value[13] = el_par_value_recurse("k0s",p->p_elem);
    c6t_elem->value[14] = el_par_value_recurse("k1",p->p_elem);
    c6t_elem->value[15] = el_par_value_recurse("k1s",p->p_elem);
    c6t_elem->value[16] = el_par_value_recurse("k2",p->p_elem);
    c6t_elem->value[17] = el_par_value_recurse("k2s",p->p_elem);
    c6t_elem->value[18] = el_par_value_recurse("k3",p->p_elem);
    c6t_elem->value[19] = el_par_value_recurse("k3s",p->p_elem);
   } else if ((strcmp(p->base_name,"multipole") == 0)) {
    maxkn=0;maxks=0;
/*      if ((kn_param = return_param_recurse("knl",p->p_elem))) maxkn=kn_param->double_array->curr; */
/*      if ((ks_param = return_param_recurse("ksl",p->p_elem))) maxks=ks_param->double_array->curr; */
    if ((index = name_list_pos("knl",p->p_elem->def->par_names))>-1) {
      kn_param = p->p_elem->def->par->parameters[index];
      maxkn=kn_param->double_array->curr;
    }
    if ((index = name_list_pos("ksl",p->p_elem->def->par_names))>-1) {
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
    c6t_elem->value[11] = el_par_value_recurse("lrad",p->p_elem);
    for (i=0; i<j; i++){
      if (i<maxkn) {c6t_elem->value[i*2+12] = kn_param->double_array->a[i]; }
      if (i<maxks) {c6t_elem->value[i*2+13] = ks_param->double_array->a[i]; }
    }
  } else if ((strcmp(p->base_name,"rfcavity") == 0)) {
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
  } else if ((strcmp(p->base_name,"marker") == 0)   ||
	     (strcmp(p->base_name,"instrument") == 0)    ||
	     (strcmp(p->base_name,"hmonitor") == 0) ||
	     (strcmp(p->base_name,"vmonitor") == 0) ||
	     (strcmp(p->base_name,"monitor") == 0)) {
    c6t_elem = new_c6t_element(0,t_name,p->base_name);
    clean_c6t_element(c6t_elem);
    strcpy(c6t_elem->org_name,t_name);
    c6t_elem->value[0] = el_par_value_recurse("l",p->p_elem);
  } else if ((strcmp(p->base_name,"rcollimator") == 0) ||
	     (strcmp(p->base_name,"ecollimator") == 0)){
    c6t_elem = new_c6t_element(13,t_name,p->base_name);
    clean_c6t_element(c6t_elem);
    strcpy(c6t_elem->org_name,t_name);
    c6t_elem->value[0] = el_par_value_recurse("l",p->p_elem);
    c6t_elem->value[12] = el_par_value_recurse("xsize",p->p_elem);
    c6t_elem->value[13] = el_par_value_recurse("ysize",p->p_elem);
  } else if((strcmp(p->base_name,"beambeam") == 0  )){
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
  } else if((strcmp(p->base_name,"elseparator") == 0  )){
    c6t_elem = new_c6t_element(3,t_name,p->base_name);
    clean_c6t_element(c6t_elem);
    strcpy(c6t_elem->org_name,t_name);
    c6t_elem->value[0] = el_par_value_recurse("l",p->p_elem);
    c6t_elem->value[2] = el_par_value_recurse("ex",p->p_elem);
    c6t_elem->value[3] = el_par_value_recurse("ey",p->p_elem);
  } else if (strcmp(p->base_name,"drift") == 0) {
    c6t_elem = new_c6t_element(0,t_name,p->base_name);
    clean_c6t_element(c6t_elem);
    strcpy(c6t_elem->org_name,t_name);
    c6t_elem->value[0] = el_par_value_recurse("l",p->p_elem);
  } else if (strcmp(p->base_name,"solenoid") == 0) {
    warning("Solenoid converted as drift :",t_name);
    c6t_elem = new_c6t_element(0,t_name,"drift");
    clean_c6t_element(c6t_elem);
    strcpy(c6t_elem->org_name,t_name);
    c6t_elem->value[0] = el_par_value_recurse("l",p->p_elem);
  } else {
    printf("Element not convertible! name= %s, basename = %s\n",p->name,p->base_name);
  }


  if (c6t_elem) {
    for (j = 0; j < c6t_elem->n_values; j++) 
      if (fabs(c6t_elem->value[j]) < eps_12) 
	c6t_elem->value[j] = 0.0;
    /* check to see if this has an aperture assigned, check for aperture flag */
    if ((aperture_flag) 
	&& (aper_param = return_param_recurse("apertype", p->p_elem))) {
      tag_aperture.apply=1;
      strcpy(tag_aperture.style,aper_param->string);
      strcpy(tag_aperture.name,t_name);
      strcat(tag_aperture.name,"_AP");
      if ((aper_param = return_param_recurse("aperture", p->p_elem))) {
        if (aper_param->expr_list != NULL) 
	  update_vector(aper_param->expr_list, aper_param->double_array);
	j=3; 
	if (aper_param->double_array->curr<3) j=aper_param->double_array->curr;
	for(i=0;i<j;i++) {
	  tag_aperture.value[i] = aper_param->double_array->a[i];
	}
      }
    }
    
    /* name used has to be without occ_cnt as this is added 
       (only 1) in tab_name_code */
    c6t_elem->twtab_row = table_row(current_sequ->tw_table,t_name);
  }

  return c6t_elem;
}

void dump_c6t_element(struct c6t_element* el)
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

void dump_c6t_sequ(int level)
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

void dump_types(int flag)
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

void equiv_elem()
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
            && strcmp(el->base_name,"marker") != 0) /* do not touch markers */
	      {
	       for (k = j+1; k < types.member[i]->curr; k++)
		 {
		  eln = types.member[i]->elem[k];
		  if (eln->flag > 0 
                      && eln->equiv == eln 
                      && ident_el(el, eln) == 0
                      && strcmp(eln->base_name,"marker") != 0
                      && strstr(eln->base_name,"colli") == NULL)
                    eln->equiv = el;
		 }
	      }
	   }
	}
    }
}

int f34_values(struct c6t_element* el, int* flags, double* values)
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

struct block* get_block_equiv(struct block* current)
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

void get_args(struct in_cmd* my_cmd)
{
  double tmp_ref_def;
  if ((aperture_flag = command_par_value("aperture", my_cmd->clone)))
    put_info("c6t - aperture flag selected","");
  if ((cavall_flag = command_par_value("cavall", my_cmd->clone)))
    put_info("c6t - cavall flag selected","");
  if ((split_flag = command_par_value("split", my_cmd->clone)))
    put_info("c6t - split flag selected","");
  if ((tmp_ref_def = command_par_value("radius", my_cmd->clone))>0.) {
    radius_flag = 1;
    ref_def = tmp_ref_def;
    printf("Reference radius set to : %f\n",ref_def);
  }
}

void get_error_refs(struct c6t_element* el)
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
  if (el->mult_order==0) { 
    el->ref_delta = c1p3 * tmp * power_of(el->ref_radius, el->mult_order);
  }else {
    el->ref_delta = 0;
  }
}

int get_flag(struct c6t_element* el, struct type_info* type)
{

  if (el->value[0] == zero)
    {
     if (type->flag_1 == 4) return in_keep_list(el);
     else return type->flag_1;
    }
  if (el->n_values < 7) return type->flag_2;
  else return (el->value[6] == 0 ? type->flag_2 : type->flag_3);
}

struct c6t_element* get_from_ellist(char* name, char* type)
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

void get_multi_refs()
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

int get_next_name(char* name, char acro)
{
  char s[2];
  int j, k = -1;

  for (j = 0; j < acro_occ; j++) if (acro_list[j] == acro)  k = j;
  if (k < 0)
    {
     k = acro_occ++; acro_list[k] = acro; acro_cnt[k] = 0;
    }
  s[0] = acro; s[1] = '\0';
  sprintf(name, "%s_c6t_%d", s, ++acro_cnt[k]);
  return 1;
}

void gnu_file(struct c6t_element* el)
{
  double el_start, el_end;

  el_start = el->position - el->value[0] / two;
  el_end   = el->position + el->value[0] / two;
  printf("%s %e 0.\n", el->name, el_start);
  printf("%s %e 1.\n", el->name, el_start);
  printf("%s %e 1.\n", el->name, el_end);
  printf("%s %e 0.\n", el->name, el_end);
}

void grow_ellist( /* doubles object list size */
		 struct c6t_el_list* p)
{
  struct c6t_element** p_loc = p->elem;
  int j, new = 2*p->max;
  char rout_name[] = "c6t:grow_ellist";

#ifdef _call_tree_
       puts("+++++++ grow_ellist");
#endif
  p->max = new;
  p->elem = (struct c6t_element**) mycalloc(rout_name,new, sizeof(struct c6t_element*));
  for (j = 0; j < p->curr; j++) p->elem[j] = p_loc[j];
  free(p_loc);
}

int ident_el(struct c6t_element* el1, struct c6t_element* el2)
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

int ident_zero(struct c6t_element* el)
{
  int j;
  for (j = 12; j < el->n_values; j++)
    if (el->value[j] != zero) return 1;
  return 0;
}

int in_keep_list(struct c6t_element* el)
{
  char temp[24];
  int j;

  strcpy(temp, el->name); lower(temp);
  for (j = 0; j < MM_KEEP; j++)
    {
     if (strncmp(temp, keep_these[j], strlen(keep_these[j])) == 0) return 2;
    }
  return 0;
}

void invert_normal(int count, double array[])
{
  int i;
  for (i = 0; i < (count+1)/2; i++)  array[2*i] = -array[2*i];
}

void invert_skew(int count, double array[])
{
  int i;
  for (i = 0; i < count/2; i++)  array[2*i+1] = -array[2*i+1];
}

void link_c6t_in_front(struct c6t_element* new, struct c6t_element* el)
{
  if (el->previous == NULL) first_in_sequ = new;
  else el->previous->next = new;
  new->previous = el->previous; new->next = el;
  el->previous = new;
}

void link_behind(struct c6t_element* new, struct c6t_element* el)
{
  if (el->next == NULL) last_in_sequ = new;
  else el->next->previous = new;
  new->previous = el; new->next = el->next;
  el->next = new;
}

void lower(char* s)
{
  char* cp = s;
  while(*cp != '\0') 
    {
     *cp = (char) tolower((int)*cp); cp++;
    }
}

struct c6t_element* make_c6t_element(struct node* p)
{
  struct c6t_element *tmp_element;
  if ((tmp_element = convert_madx_to_c6t(p))) {
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
struct object* make_obj(   /* creates a new object */
               char* key,
               int vlint,       /* length of integer array */
               int vldble,      /* length of double array */
               int vlchar,      /* length of char array */
               int vlpobj)      /* length of object pointer array */
{
  struct object* p;
  char rout_name[] = "c6t:make_obj";

#ifdef _call_tree_
       put_info("+++++++ make_object","");
#endif
  p = (struct object*)  mycalloc(rout_name, 1, sizeof(struct object));
  mycpy(p->key, key);
  if ((p->l_int = vlint) > 0) 
       p->a_int = (int*) mymalloc(rout_name, p->l_int * sizeof(int));
  if ((p->l_dble = vldble) > 0)
       p->a_dble = (double*) mymalloc(rout_name, p->l_dble * sizeof(double));
  if ((p->l_char = vlchar) > 0)
       p->a_char = (char*) mymalloc(rout_name, p->l_char);
  if ((p->l_obj  = vlpobj) > 0)
    {
     p->p_obj = (struct object**) mycalloc(rout_name, p->l_obj, sizeof(struct object*));
     p->names = (char**) mycalloc(rout_name, p->l_obj, sizeof(char*));
    }
  p->parent = NULL;
/*      my_time(); */
/*      p->ma_time = major_time; p->mi_time = minor_time; */
  return p;
}

void make_multipole(struct c6t_element* el)
{
  if (el->force > 0) /* multiply forces with length */
  app_factor(el->value[0], &el->value[12], el->n_values-12);
  el->value[0] = zero; /* set element length to zero */
  el->flag = 2;
  remove_from_ellist(el);
  strcpy(el->base_name, "multipole");
  add_to_ellist(el);
}


void mod_errors()
{
  current_element = first_in_sequ;
  while (current_element != NULL)
   {
    if (current_element->nf_err > 0) 
     invert_normal(current_element->nf_err, current_element->p_fd_err->a_dble);
    current_element = current_element->next;
   }
}

void mod_lcavity(struct c6t_element* p)
{
    total_voltage += p->value[1];   /* accumulate voltage */
}

void mod_multipole(struct c6t_element* p)
{
  supp_small_comp(p);
}

void mod_octupole(struct c6t_element* p)
{
  supp_small_comp(p);
}

void mod_quadrupole(struct c6t_element* p)
{
  supp_small_comp(p);
}

void mod_rbend(struct c6t_element* p)
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

void mod_rfcavity(struct c6t_element* p)
{
  total_voltage += p->value[1];  /* accumulate voltage */
}

void mod_sextupole(struct c6t_element* p)
{
  supp_small_comp(p); 
  /* supress tilt angle (not component !) */
/*    p->value[6] = zero; */
}

void multi_loop()
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

struct block* new_block()
{
  struct block* p;
  char rout_name[] = "c6t:new_block";
  p = (struct block*) mycalloc(rout_name, 1, sizeof(struct block));
  sprintf(p->name, "BLOC%d", block_count++);
  return p;
}

struct c6t_element* new_c6t_element(int size, char* name, char* base)
{
  struct c6t_element* p;
  char rout_name[] = "c6t:new_c6t_element";
  p = (struct c6t_element*) mycalloc(rout_name, 1, sizeof(struct c6t_element));
  strcpy(p->name, name);
  p->equiv = p;
  strcpy(p->base_name, base);
  p->value = (double*) mycalloc(rout_name,++size,sizeof(double));
  p->n_values = size;
  return p;
}

void post_multipoles() /* post equiv. treatment of multipoles */
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

double power_of(double d, int i)
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

void pre_multipole(struct c6t_element* el) /* pre-process multipoles */
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
  int i, last_nzero = -1, nz_cnt = 0, cnt = 0, n_pole, s_pole = 0, ndmax;
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
        for (i = 0; i <= s_pole; i++) new_el->value[i] = el->value[i];
        for (i = 12; i <= s_pole; i++) el->value[i] = 0;
        new_el->flag = s_pole > 13 ? 2 : 1; new_el->npole_sign = 1;
        new_el->keep_in = el->keep_in; 
        new_el->position = el->position;
        new_el->twtab_row = el->twtab_row;
        new_el->na_err = el->na_err; /* el->na_err = 0; */
        new_el->p_al_err = el->p_al_err; /*  el->p_al_err = NULL; */
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
  n_pole = last_nzero / 2;
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
        else  sprintf(tmp_name,"%s_arfa", el->name);
        el->nf_err = last_nzero;
        el->p_fd_err = make_obj(tmp_name, 0, el->nf_err, 0, 0);
       }
     for (i = 0; i < el->nf_err; i++) el->p_fd_err->a_dble[i] = tmp_buff[i];
     for (i = 14; i < el->n_values; i++)  el->value[i] = zero;
    }
}

void process_c6t()  /* steering routine */
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

  write_all_el();
  write_blocks();
  write_struct();
  write_f16_errors();
  write_f34_special();
  write_f3aux();
  write_f3_aper();
  write_f8_errors();
}

void pro_elem(struct node* cnode) 
     /* processes one element, makes linked list */
     /* converts MADX linked list to c6t internal linked list */
{
  int i;
  char t_key[KEY_LENGTH];
  char ap_name[255];
  struct c6t_element *tag_element;
  double tmp_vk,tmp_hk;

  tag_aperture.apply=0;
  /* do the fiddly conversion but skip element if not needed */
  if (make_c6t_element(cnode) == NULL) return;

  if (strcmp(cnode->base_name, "rbend") == 0) mod_rbend(current_element);
  else if (strcmp(cnode->base_name, "lcavity") == 0) mod_lcavity(current_element);
  else if (strcmp(cnode->base_name, "multipole") == 0) mod_multipole(current_element);
  else if (strcmp(cnode->base_name, "octupole") == 0) mod_octupole(current_element);
  else if (strcmp(cnode->base_name, "quadrupole") == 0) mod_quadrupole(current_element);
  else if (strcmp(cnode->base_name, "sextupole") == 0) mod_sextupole(current_element);
  else if (strcmp(cnode->base_name, "rfcavity") == 0) mod_rfcavity(current_element);
  if (strstr(cnode->base_name, "kicker")) {
    if (cnode->p_elem) {
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
      sprintf(t_key, "%s+%d", current_element->name,cnode->occ_cnt);
      strcpy(current_element->name, t_key);
    }
  current_element->position = cnode->position;
  add_to_ellist(current_element);

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
  if (cnode->p_al_err) {
    align_cnt++;
    current_element->na_err = cnode->p_al_err->curr;
    current_element->p_al_err = make_obj("ALDUM",0,ALIGN_MAX,0,0);
    current_element->p_al_err->c_dble = cnode->p_al_err->curr;
    for (i=0;i<cnode->p_al_err->curr;i++)
      current_element->p_al_err->a_dble[i] = cnode->p_al_err->a[i];
  }
  /* if we have a tilt set the flag */
  if (current_element->n_values >= 7 && current_element->value[6] > zero) {
    current_element->tilt_err = 1;
  } else {
    current_element->tilt_err = 0;
  }
  /* add aperture element if necessary */
  if (tag_aperture.apply==1) {
    if (strstr(tag_aperture.style,"circle")!=NULL) {
      tag_element = create_aperture(tag_aperture.name,"EL",
				    tag_aperture.value[0],tag_aperture.value[0],cnode->p_al_err);
    } else if (strstr(tag_aperture.style,"ellipse")!=NULL) {
      tag_element = create_aperture(tag_aperture.name,"EL",
				    tag_aperture.value[0],tag_aperture.value[1],cnode->p_al_err);
    } else if (strstr(tag_aperture.style,"rectangle")!=NULL) {
      tag_element = create_aperture(tag_aperture.name,"RE",
				    tag_aperture.value[0],tag_aperture.value[1],cnode->p_al_err);
    } else if (strstr(tag_aperture.style,"lhcscreen")!=NULL) {
      strcpy(ap_name,tag_aperture.name); strcat(ap_name,"1");
      tag_element = create_aperture(ap_name,"EL",
				    tag_aperture.value[0],tag_aperture.value[0],cnode->p_al_err);
      tag_element->previous = current_element;
      tag_element->next = current_element->next;
      current_element->next = tag_element;
      prev_element = current_element;
      current_element = tag_element;
      current_element->position = cnode->position;
      add_to_ellist(current_element);

      strcpy(ap_name,tag_aperture.name); strcat(ap_name,"2");
      tag_element = create_aperture(ap_name,"RE",
				    tag_aperture.value[1],tag_aperture.value[2],cnode->p_al_err);
    } else {
      warning("general aperture element not supported in sixtrack",tag_aperture.name);
    }

    tag_element->previous = current_element;
    tag_element->next = current_element->next;
    current_element->next = tag_element;
    prev_element = current_element;
    current_element = tag_element;
    current_element->position = cnode->position;
    add_to_ellist(current_element);
  }
}

void read_sequ()
{
  struct node* cnode;
  if ((current_sequ->n_nodes) > 0)  sequ_start = current_sequ->ex_start->position;
  cnode=current_sequ->ex_start;
  while(cnode && cnode!=current_sequ->ex_end)
    {
      if (strstr(cnode->name,"$")==NULL) pro_elem(cnode);
      cnode=cnode->next;
    }
  sequ_end = current_sequ->ex_end->position;
  sequ_length = sequ_end - sequ_start;
  last_in_sequ = current_element;
  put_info("MADX sequence converted to c6t internal.","");
}

/* removes element from correct object list */
void remove_from_ellist(struct c6t_element* p_elem)
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

void replace_c6t(struct c6t_element* old, struct c6t_element* new)
{
  if (old->previous != NULL)  old->previous->next = new;
  new->previous = old->previous;
  if (old->next != NULL)      old->next->previous = new;
  new->next = old->next;
  old->flag = 0;
}

void split()
{
  int i;
  struct c6t_element* el;
  if (split_list != NULL)
    {
     for (i = 0; i < split_list->curr; i++)
       {
	el = split_list->elem[i];
        if (el->flag == 1 
            && (split_flag != 0 || el->nf_err > 0)) split_special(el);
        else if (el->flag == 2 || el->flag == 3)  split_other(el);
        else if (el->split == 3)  split_kicker(el);
       }
    }
}

void split_kicker(struct c6t_element* el)
{
  struct c6t_element *k1, *k2;
  char c[24];
  int af;
  if (el->value[0] > zero) split_other(el);
  if (el->flag == 6) /*split kicker into h + v */
    {
     af = get_next_name(c, 'h');
     k1 = new_c6t_element(13, c, "hkicker");
     af = get_next_name(c, 'v');
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

void split_other(struct c6t_element* el)
{
/* -> two drifts with non-lin. thin lens at centre */
  struct c6t_element *d1, *d2;
  double length = el->value[0] / two; 
  char c[24];
  int af;

  af = get_next_name(c, 'd');
  d1 = new_c6t_element(1, c, "drift");
  af = get_next_name(c, 'd');
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

void split_special(struct c6t_element* el) 
/* -> two lin. halves with multipole at centre */
{
  struct c6t_element *d1, *mt;
  double length = el->value[0] / two, mt_position = el->position; 
  char c[24];
  int j, af;

  if (el->nf_err > 0) 
    app_factor(el->value[0], el->p_fd_err->a_dble, el->nf_err);
  af = get_next_name(c, *el->base_name);
  d1 = new_c6t_element(el->n_values, c, el->base_name);
  d1->value[0] = el->value[0] = length;
  for (j = 1; j < el->n_values; j++) 
      d1->value[j] = el->value[j];
  d1->flag = el->flag = 1;
  d1->force = el->force = 1;
  d1->position = el->position + d1->value[0] / two;
  el->position = el->position - el->value[0] / two;
  af = get_next_name(c, 'm');
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

void supp_elem()
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

void supp_small_comp(struct c6t_element* p)
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

void treat_split(struct c6t_element* el)
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

void yank(struct c6t_element* el)
{
  if (el->previous != NULL)  el->previous->next = el->next;
  else                       first_in_sequ      = el->next;
  if (el->next != NULL)      el->next->previous = el->previous;
  else                       last_in_sequ       = el->previous;
  el->flag = 0;
}

void write_all_el()
{
  char title[] =
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

void write_c6t_element(struct c6t_element* el)
{
  if (strcmp(el->name, "CAV") != 0)
      fprintf(f2, "%-16s %2d  %16.9e %17.9e  %17.9e\n", 
          el->name, el->out_1, el->out_2, el->out_3, el->out_4);
  el->w_flag = 1;
}

void write_blocks()
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
           fprintf(f2, "%-18s", p->name); lc++;
           for (i = 0; i < p->elements->curr; i++)
	     {
              if (lc++ == LINES_MAX)
	       {
                fprintf(f2,"\n"); fprintf(f2,"                  "); lc = 2;
	       }
	      fprintf(f2, "%-18s",p->elements->elem[i]->equiv->name);
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

void write_f8_errors()
{
  double tiltval;
  if (align_cnt == 0)  return;
  current_element = first_in_sequ;
  while (current_element != NULL)
    {
      if (current_element->tilt_err > 0) {
	tiltval = current_element->value[6]; 
      } else {tiltval=0.0;}
      if (current_element->na_err > 0)
	{
	  if (f8_cnt++ == 0)    f8 = fopen("fc.8", "w");
	  fprintf(f8, "%-16s  %14.6e%14.6e%17.9e\n",current_element->equiv->name,
		  1000*current_element->p_al_err->a_dble[0],
		  1000*current_element->p_al_err->a_dble[1],
		  1000*(current_element->p_al_err->a_dble[5]+tiltval));
	} else if (current_element->tilt_err > 0) {
	  if (f8_cnt++ == 0)    f8 = fopen("fc.8", "w");
	  fprintf(f8, "%-16s  %14.6e%14.6e%17.9e\n",current_element->equiv->name,
		  0.0,
		  0.0,
		  1000*tiltval);
	}
      current_element = current_element->next;
    }
}

void write_f16_errors()
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

void write_f34_special()
{
  int i, j, n, flags[FIELD_MAX];
  double values[FIELD_MAX];
  char* t_list[NT34];
  char t_name[NAME_L];
  char* cp;
  double spos,betx,bety,mux,muy;
  int err;
  
  t_list[0] = &mpole_names[1][0];
  t_list[1] = &mpole_names[2][0];
  t_list[2] = &mpole_names[3][0];
  t_list[3] = &mpole_names[4][0];
  t_list[4] = &mpole_names[5][0];
  
  if (special_flag == 0)  return;

  current_element = first_in_sequ;
  while (current_element != NULL)
    {
     for (i = 0; i < NT34; i++)
       {
        if (strcmp(current_element->base_name, t_list[i]) == 0)
	  {
	   n = f34_values(current_element, flags, values);
           if (f34_cnt++ == 0)    f34 = fopen("fc.34", "w");
           for (j = 0; j < n; j++)
	     {
	      strcpy(t_name, current_element->name);
	      if ((cp = strchr(t_name, '+')) != NULL) *cp = '\0';
	      if ((err=double_from_table("twiss","s",&(current_element->twtab_row),&spos)))
		printf ("Not found double_from table = %i\n",err);
	      if ((err=double_from_table("twiss","betx",&(current_element->twtab_row),&betx)))
		printf ("Not found double_from table = %i\n",err);
	      if ((err=double_from_table("twiss","bety",&(current_element->twtab_row),&bety)))
		printf ("Not found double_from table = %i\n",err);
	      if ((err=double_from_table("twiss","mux",&(current_element->twtab_row),&mux)))
		printf ("Not found double_from table = %i\n",err);
	      if ((err=double_from_table("twiss","muy",&(current_element->twtab_row),&muy)))
		printf ("Not found double_from table = %i\n",err);
              fprintf(f34, 
              " %20.13e  %-16s %3d %20.13e %20.13e %20.13e %20.13e %20.13e\n",
		      spos,t_name,flags[j],values[j],betx,bety,mux,muy);
	     }
	  }
       }
     current_element = current_element->next;
    }
  if (last_in_sequ->twtab_row > 0) {
    if ((err=double_from_table("twiss","s",&(last_in_sequ->twtab_row),&spos)))
      printf ("Not found double_from table = %i\n",err);
    if ((err=double_from_table("twiss","betx",&(last_in_sequ->twtab_row),&betx)))
      printf ("Not found double_from table = %i\n",err);
    if ((err=double_from_table("twiss","bety",&(last_in_sequ->twtab_row),&bety)))
      printf ("Not found double_from table = %i\n",err);
    if ((err=double_from_table("twiss","mux",&(last_in_sequ->twtab_row),&mux)))
      printf ("Not found double_from table = %i\n",err);
    if ((err=double_from_table("twiss","muy",&(last_in_sequ->twtab_row),&muy)))
      printf ("Not found double_from table = %i\n",err);
  }
  fprintf(f34, 
	  " %20.13e  %-16s %3d %20.13e %20.13e %20.13e %20.13e %20.13e\n",
	  spos,"end_marker",100,zero,betx,bety,mux,muy);
}

void write_f3_aper()
{
  int f3aper_cnt = 0;
  current_element = first_in_sequ;
  while (current_element != NULL)
    {
      if (strstr(current_element->name,"_AP")!=NULL
	  && (current_element->equiv == current_element)) {
	if (f3aper_cnt++ == 0) {
	  f3aper  = fopen("fc.3.aper", "w");
	  fprintf(f3aper,"LIMI\n");
	}
	if (current_element->value[3] == 1) {
	  fprintf(f3aper,"%s %s %f %f\n",current_element->name,"RE",
		 current_element->value[1], current_element->value[2]);
	} else {
	  fprintf(f3aper,"%s %s %f %f\n",current_element->name,"EL",
		 current_element->value[1], current_element->value[2]);	
	}
      }
     current_element = current_element->next;
    }
  if (f3aper_cnt > 0) fprintf(f3aper,"NEXT\n");
}

void write_f3_entry(char* option, struct c6t_element* el)
{
  if (f3_cnt++ == 0)     f3  = fopen("fc.3", "w");
  if (strcmp(option, "multipole") == 0) write_f3_mult(el);
}

void write_f3aux()
{
  double aux_val[4] = {-1.e20, -1.e20, -1.e20, -1.e20};
  double tw_alfa;
  int row=1;
  if ((double_from_table("summ","q1", &row, &(aux_val[0])) !=0) ||
      (double_from_table("summ","q2",  &row, &(aux_val[1])) !=0) ||
      (double_from_table("summ","dq1", &row, &(aux_val[2])) !=0) ||
      (double_from_table("summ","dq2", &row, &(aux_val[3])) !=0)) {
    printf("c6t error: tunes or chromaticities not found!\n");
  }      
  if (current_beam != NULL)
    {
     if (f3aux_cnt++ == 0)     f3aux  = fopen("fc.3.aux", "w");
     if (double_from_table("summ","alfa", &row, &tw_alfa) !=0) 
       printf("c6t warning: alfa not found in twiss\n");
     fprintf(f3aux, "SYNC\n");
     fprintf(f3aux,"%12.0f%10.6f%10.3f 0.%10.3f%12.6f  1\n",
     	     harmon, tw_alfa, total_voltage, sequ_length, 
	     c1p3*command_par_value("mass", current_beam));
     fprintf(f3aux,"      1.        1.\n");
     fprintf(f3aux, "NEXT\n");
     fprintf(f3aux, "BEAM\n");
     fprintf(f3aux, "%12.4e%14.6g%14.6g%12.4e%12.4e  1  0\n",
	     command_par_value("npart", current_beam),
	     0.25e6*command_par_value("exn", current_beam),
	     0.25e6*command_par_value("eyn", current_beam),
	     command_par_value("sigt", current_beam),
	     command_par_value("sige", current_beam));
     fprintf(f3aux, "NEXT\n");
    }
  if (aux_val[0] > -1.e10 && aux_val[1] > -1.e10)
    {
     fprintf(f3aux, "TUNE\n");
     fprintf(f3aux, "QF%12.5f\n", aux_val[0]);
     fprintf(f3aux, "QD%12.5f\n", aux_val[1]);
     fprintf(f3aux, "NEXT\n");
    }
  if (aux_val[2] > -1.e10 && aux_val[3] > -1.e10)
    {
     fprintf(f3aux, "CHRO\n");
     fprintf(f3aux, "SXF%12.5f\n", aux_val[2]);
     fprintf(f3aux, "SXD%12.5f\n", aux_val[3]);
     fprintf(f3aux, "NEXT\n");
    }
}

void write_f3_mult(struct c6t_element* el)
{
  int i, j, i_max = -1;
  struct c6t_element* eln;
  if (multi_type < 0)  return;
  fprintf(f3,"MULT\n");
  fprintf(f3,"%-16s%20.10e%20.10e\n", el->name, c1p3*el->ref_radius,
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
             if (eln->p_fd_err->a_dble[i] != zero)
	       {
                i_max = i; error_matrix[i] = 1.;
	       }
	    }
  	}
     }
  if (++i_max > 0)  i_max += i_max%2;
  for (i = 0; i < i_max; i++) 
    {
     fprintf(f3,"%4.0f.%4.0f.", 0., error_matrix[i]);
     if ((i+1)%2 == 0) fprintf(f3,"\n");
    }
  fprintf(f3,"NEXT\n");
}

void write_struct()
{
  struct block* p = first_block;
  int lc = 0;
  char* out;
  char title[] =
  "STRUCTURE INPUT---------------------------------------------------------";

  fprintf(f2, "%s\n", title);
  while (p != NULL)  
    {
     if (p->flag == 0) out = p->first->equiv->name;
     else              out = p->equiv->name;
     if (lc++ == LINES_MAX)
      {
       fprintf(f2,"\n"); lc = 1;
      }
     fprintf(f2, "%-18s", out);
     p = p->next;
    }
  if (lc > 0)
    {
     fprintf(f2,"\n");
    }
  fprintf(f2, "NEXT\n");
}
