/* makethin.c
   Thick to thin lens converter. Helmut Burkhardt
   Early versions in 2001, 2002 by Mark Hayes
*/

#ifdef _WRAP_FORTRAN_CALLS
#include "fortran_wrappers.h"
#endif
#ifdef _WRAP_C_CALLS
#include "c_wrappers.h"
#endif


/* define bool like in C++ */
#ifndef bool_for_c
#define bool_for_c
typedef unsigned char bool;
#define true 1
#define false 0
#endif

/* forward declarations of some routines which are only used within makethin */
struct element* create_thin_obj(struct element*,int slice_no);
struct sequence* seq_diet(struct sequence*);
double at_shift(int ,int );
double q_shift(int ,int );

/* this structure is used to store a lookup table of thick to thin
   element conversions already done */
struct thin_lookup
{
    struct element *thick_elem;
    struct element *thin_elem;
    int slice;
    struct thin_lookup *next;
};
struct thin_lookup *my_list = NULL;

/* this structure is used to store a lookup table of thick to thin
   sequence conversions already done */
struct thin_sequ_lookup
{
    struct sequence *thick_sequ;
    struct sequence *thin_sequ;
    struct thin_sequ_lookup *next;
};
struct thin_sequ_lookup *my_sequ_list = NULL;

/* this is used when choosing the style of slicing */
char* thin_style = NULL;
char collim_style[] = "collim";
/* this is used to return how to slice the selected elements */
struct el_list *thin_select_list = NULL;

/* code starts here **************************************************/

void force_consistent_slices(void)
/* hbu 10/2005
   loop over all elements and check that #slices of child and parent agree
   if not, use the maximum for both
*/
{
  struct element* el_i;
  struct command_parameter *child,*parent;
  int i,el_i_slice_pos,slices,slices_parent;
  for(i=0; i< element_list->curr; i++) /* loop over element_list */
  {
    el_i = element_list->elem[i];
    el_i_slice_pos = name_list_pos("slice",el_i->def->par_names);
    if(el_i_slice_pos>0 && el_i->parent!=NULL && el_i != el_i->parent )
    {
      child=el_i->def->par->parameters[el_i_slice_pos];
      parent=el_i->parent->def->par->parameters[el_i_slice_pos];
      slices=child->double_value;
      slices_parent=parent->double_value;
      if(slices != slices_parent)
      {
        if(slices>slices_parent) slices_parent=slices; else slices=slices_parent;
        child->double_value=slices;
        parent->double_value=slices_parent;
      }
    }
  }
}

void dump_slices(void)
/* Loops over all current elements and prints the number of slices. Used for debug and info */
{
  struct element* el_i;
  int i,el_i_slice_pos,slices,slices_parent,n_elem_with_slice=0,n_elem_with_slice_gt_1=0;
  char* parent_name;
  printf("++++++ dump_slices");
  printf("            name #slices  derived from #slices\n");
  for(i=0; i< element_list->curr; i++) /* loop over element_list */
  {
    el_i = element_list->elem[i];
    el_i_slice_pos = name_list_pos("slice",el_i->def->par_names);
    if(el_i_slice_pos>0)
    {
      n_elem_with_slice++;
      slices=el_i->def->par->parameters[el_i_slice_pos]->double_value;
      /* look also at parent if existing */
      slices_parent=0;
      parent_name="no parent";
      if(el_i->parent!=NULL)
      {
        slices_parent=el_i->parent->def->par->parameters[el_i_slice_pos]->double_value;
        parent_name=el_i->parent->name;
      }
      if(slices>1) n_elem_with_slice_gt_1++;
    }
  }
  printf("------ end of dump slices. There were %4d elements, %3d with slice numbers and %2d with slice numbers>1\n\n",element_list->curr,n_elem_with_slice,n_elem_with_slice_gt_1);
}

int get_slices_from_elem(struct element* elem)
{
  int elem_slice_pos=0,slices=1;
  elem_slice_pos = name_list_pos("slice",elem->def->par_names);
  if(elem_slice_pos > 0)
  {
    slices=elem->def->par->parameters[elem_slice_pos]->double_value;
  }
  if (slices==0) slices = 1; /* must always slice to thin */
  return slices;
}

/* Has this element already been dieted? returns NULL for NO.*/
struct element* get_thin(struct element* thick_elem, int slice)
{
  struct thin_lookup *cur;
  if (my_list)
  {
    cur = my_list;
    while (cur)
    {
      if (cur->thick_elem == thick_elem && cur->slice == slice)
      {
        return cur->thin_elem;
      }
      cur = cur->next;
    }
  }
  return NULL;
}

/* Enter a newly dieted element */
void put_thin(struct element* thick_elem, struct element* thin_elem, int slice)
{
  struct thin_lookup *p,*cur;
  char rout_name[] = "makethin:put_thin";
  p = (struct thin_lookup*) mycalloc(rout_name,1, sizeof(struct thin_lookup));
  p->thick_elem = thick_elem;
  p->thin_elem = thin_elem;
  p->slice = slice;
  p->next = NULL;
  if (my_list)
  {
    cur = my_list;
    while (cur->next) cur = cur->next;
    cur->next = p;
  }
  else
  {
    my_list = p;
  }
  return;
}

/* Has this sequence already been dieted? returns NULL for NO.*/
struct sequence* get_thin_sequ(struct sequence* thick_sequ)
{
  struct thin_sequ_lookup *cur;
  if (my_sequ_list)
  {
    cur = my_sequ_list;
    while (cur)
    {
      if (cur->thick_sequ == thick_sequ)
      {
        return cur->thin_sequ;
      }
      cur = cur->next;
    }
  }
  return NULL;
}

/* Enter a newly dieted element */
void put_thin_sequ(struct sequence* thick_sequ, struct sequence* thin_sequ)
{
  struct thin_sequ_lookup *p,*cur;
  char rout_name[] = "makethin:put_thin_sequ";
  p = (struct thin_sequ_lookup*) mycalloc(rout_name,1, sizeof(struct thin_sequ_lookup));
  p->thick_sequ = thick_sequ;
  p->thin_sequ = thin_sequ;
  p->next = NULL;
  if (my_sequ_list)
  {
    cur = my_sequ_list;
    while (cur->next) cur = cur->next;
    cur->next = p;
  }
  else
  {
    my_sequ_list = p;
  }
  return;
}

/* makes node name from element name and slice number*/
char* make_thin_name(char* e_name, int slice)
{
  char c_dummy[128];
  sprintf(c_dummy,"%s..%d", e_name, slice);
  return buffer(c_dummy);
}

struct expression* compound_expr(struct expression* e1, double v1,
                                 char* oper, struct expression* e2, double v2)
/* make one out of two expressions, using oper to connect them
   hbu 9/2005 moved from madxn.c to makethin.c as only used here
   and increased precision   sprintf(tmp, "%e"  ->   sprintf(tmp, "%.14g" */
{
  char** toks = tmp_l_array->p;
  struct expression* expr = NULL;
  char tmp[30];
  int n;
  char lb[] = "(", rb[] = ")";
  if (e1 != NULL || e2 != NULL)
  {
    if (e1 != NULL)
    {
      if (e2 != NULL)
      {
        toks[0] = lb; toks[1] = e1->string; toks[2] = rb;
        toks[3] = oper;
        toks[4] = lb; toks[5] = e2->string; toks[6] = rb;
      }
      else
      {
        sprintf(tmp, "%.14g", v2); /* hbu */
        toks[0] = lb; toks[1] = e1->string; toks[2] = rb;
        toks[3] = oper;
        toks[4] = lb; toks[5] = tmp; toks[6] = rb;
      }
    }
    else
    {
      sprintf(tmp, "%.14g", v1);  /* hbu */
      toks[0] = lb; toks[1] = tmp; toks[2] = rb;
      toks[3] = oper;
      toks[4] = lb; toks[5] = e2->string; toks[6] = rb;
    }
    join(toks, 7);
    pre_split(c_join->c, l_wrk, 0);
    n = mysplit(l_wrk->c, tmp_l_array);
    expr = make_expression(n, toks);
  }
  return expr;
}

/* scale an expression by a number - or leave it NULL */
struct expression* scale_expr(struct expression* expr,double scale)
{
  if (expr) return compound_expr(expr,0,"*",NULL,scale);
  return NULL;
}

/* combine two parameters using compound expression */
struct expression* comb_param(struct command_parameter* param1,
                              char* op, struct command_parameter* param2)
{
  return compound_expr(param1->expr,param1->double_value,op,param2->expr,param2->double_value);
}

/* returns parameter if it has been modified, otherwise NULL */
struct command_parameter* return_param(char* par, struct element* elem)
{
  int index;
  /* don't return base type definitions */
  if (elem==elem->parent) return NULL;

  if ((index = name_list_pos(par,elem->def->par_names))>-1
      && elem->def->par_names->inform[index] > 0)
    return elem->def->par->parameters[index];
  return NULL;
}

/* returns parameter if it has been modified, otherwise NULL  - recursively */
struct command_parameter* return_param_recurse(char* par, struct element* elem)
{
  struct command_parameter* param;
  param = return_param(par,elem);

  if (param) return param;
  if (elem!=elem->parent)
    return return_param_recurse(par,elem->parent);
  return NULL;
}

/* returns first parameter value found and recusively checks sub_elements */
double el_par_value_recurse(char* par, struct element* elem)
{
  if (return_param(par, elem)) return el_par_value(par,elem);
  if (elem != elem->parent)
    return el_par_value_recurse(par,elem->parent);
  return 0;
}

void add_cmd_parameter_clone(struct command* cmd,struct command_parameter *param,char* par_name,int inf) /*hbu add an identical copy (clone) of param to cmd */
{
  if(param)
  {
    cmd->par->parameters[cmd->par->curr] = clone_command_parameter(param); /* set current to identical copy (clone) of param */
    add_to_name_list(par_name,inf,cmd->par_names);
    cmd->par->curr++;
  }
}

void add_cmd_parameter_new(struct command* cmd,double par_value,char* par_name,int inf) /*hbu add a new param with one value to cmd */
{
  cmd->par->parameters[cmd->par->curr] = new_command_parameter(par_name, 2);
  cmd->par->parameters[cmd->par->curr]->double_value = par_value;
  add_to_name_list(par_name,inf,cmd->par_names);
  cmd->par->curr++;
}

/* multiply the k by length and divide by slice */
struct command_parameter* scale_and_slice(struct command_parameter *kn_param,
                                          struct command_parameter *length_param,
                                          int slices, int slice_no,
                                          int angle_conversion, int kl_flag)
{
  int last_non_zero=-1,i;
  struct expression *kn_i_expr;
  double kn_i_val;
  if (kn_param == NULL) return NULL;

  for (i=0; i<kn_param->expr_list->curr; i++)
  {
    kn_i_expr = kn_param->expr_list->list[i];
    kn_i_val  = kn_param->double_array->a[i];
    if ((kn_i_expr != NULL && zero_string(kn_i_expr->string)==0)  || kn_i_val!=0)
    {
      last_non_zero=i;
      if (kl_flag == 0 && (angle_conversion==0||i>0)) /*hbu apply the angle_conversion==0 check only to zero order multipole */
      {
        if ((length_param->expr) || (kn_i_expr))
        {
          kn_i_expr = compound_expr(kn_i_expr,kn_i_val,"*",length_param->expr,length_param->double_value); /* multiply expression with length */
        }
        else
        { /* multiply value with length */
          kn_i_val *= length_param->double_value;
        }
      }
      if (slices > 1)
      { /* give the correct weight by slice (multiply with the inverse of the number of slices) */
        if (kn_i_expr)
        {
          kn_i_expr = compound_expr(kn_i_expr,kn_i_val,"*",NULL,q_shift(slices,slice_no));
        }
        else
        {
          kn_i_val *= q_shift(slices,slice_no);
        }
      }
    }
    if(kn_i_expr) kn_param->expr_list->list[i] = kn_i_expr;
    kn_param->double_array->a[i] = kn_i_val;
  } /* for i ..*/
  if (last_non_zero==-1)
  {
    delete_command_parameter(kn_param); kn_param=NULL;
  }
  return kn_param;
}

/* translate k0,k1,k2,k3 & k0s,k1s,k2s,k3s to kn{} and ks{} */
/* 26/11/01 - removed tilt param */
int translate_k(struct command_parameter* *kparam,
                struct command_parameter* *ksparam,
                struct command_parameter *angle_param,
                struct command_parameter *kn_param,
                struct command_parameter *ks_param)
{
  int i,angle_conversion=0;
/*    char *zero[1]; */
/*    zero[0] = buffer("0"); */

  if ((kparam == NULL) && (ksparam == NULL))
    fatal_error("translate_k: no kparams to convert","");

  /* if we have a angle we ignore any given k0 */
  if (angle_param)
  {
    kparam[0] =  new_command_parameter("k0", 2);
    angle_conversion=1; /* note we do not divide by length, just to multiply again afterwards */
    if (angle_param->expr)
    {
      kparam[0]->expr =  clone_expression(angle_param->expr);
    }
    kparam[0]->double_value = angle_param->double_value;
  }

  for (i=0; i<4; i++)
  {
    /* zero all the parameters */
    kn_param->expr_list->list[i] = NULL; kn_param->double_array->a[i] = 0;
    ks_param->expr_list->list[i] = NULL; ks_param->double_array->a[i] = 0;
    /* copy across the k's */
    if (kparam[i])
    {
      if (kparam[i]->expr)
      {
        kn_param->expr_list->list[i] = clone_expression(kparam[i]->expr);
      }
      kn_param->double_array->a[i] = kparam[i]->double_value;
    }
    if (ksparam[i])
    {
      if (ksparam[i]->expr)
      {
        ks_param->expr_list->list[i] = clone_expression(ksparam[i]->expr);
      }
      ks_param->double_array->a[i] = ksparam[i]->double_value;
    }
    /* update the number of k's in our arrays */
    kn_param->expr_list->curr++; kn_param->double_array->curr++;
    ks_param->expr_list->curr++; ks_param->double_array->curr++;
  }

  return angle_conversion;
}

/* adds a node to the end of a sequence */
void seq_diet_add(struct node* node, struct sequence* sequ)
{
  if (sequ->start == NULL)
  { /* first node in new sequence? */
    sequ->start = node;
    sequ->end = node;
    node->next = NULL;
    node->previous = NULL;
  }
  else
  { /* no? then add to end */
    sequ->end->next = node;
    node->previous  = sequ->end;
    sequ->end = node;
  }
  add_to_node_list(node, 0, sequ->nodes);

  return;
}

/* adds a sequence to a sequence */
void seq_diet_add_sequ(struct node* thick_node, struct sequence* sub_sequ, struct sequence* sequ)
{
  struct node* node = new_sequ_node(sub_sequ, thick_node->occ_cnt); /* 1 is the occ_cnt*/
  node->length = 0;
  node->at_value = thick_node->at_value;
  if (node->at_expr)
  {
    node->at_expr = clone_expression(thick_node->at_expr);
  }
  seq_diet_add(node,sequ);
  return;
}

void add_lrad(struct command* cmd,struct command_parameter *length_param,int slices)
{
  struct command_parameter *l_par;
  if(length_param)
  {
    add_cmd_parameter_new(cmd,0.,"l",1); /* new parameter l with value of 0 */
    l_par = cmd->par->parameters[cmd->par->curr] = clone_command_parameter(length_param); /* keep what was l */
    strcpy(l_par->name,"lrad"); /* but rename to lrad and slice : */
    if (slices > 1) /* divide numbers or expressions by the number of slices */
    {
      if (l_par->expr) l_par->expr = compound_expr(l_par->expr,0.,"/",NULL,slices);
      else l_par->double_value /= slices;
    }
    add_to_name_list("lrad",1,cmd->par_names);
    cmd->par->curr++;
  }
}

/* creates the thin magnetic element - recursively for classes from which dericed (parent) */
struct element* create_thin_multipole(struct element* thick_elem, int slice_no)
{
  struct command_parameter *angle_param, *length_param, *kparam[4], *ksparam[4], *kn_param, *ks_param, *at_param, *fint_param;
  struct element *thin_elem_parent, *thin_elem;
  struct command* cmd;
  char *thin_name;
  int angle_conversion = 0;
  int slices, minimizefl;
  int knl_flag = 0,ksl_flag = 0;

  /* next is new to handle parent with possibly different slice number than child */
  slices = get_slices_from_elem(thick_elem);
  at_param = return_param("at",thick_elem);

  if (thick_elem == thick_elem->parent) return NULL; /* no further parent to consider */
  else
  {
    thin_elem_parent = create_thin_multipole(thick_elem->parent,slice_no); /* slice also the parent */
  }

  minimizefl=get_option("minimizeparents") && !at_param && thick_elem == thick_elem->parent;
  if(minimizefl)
  {
    slice_no=slices=1; /* do not slice this one */
  }
  if(slice_no > slices && thick_elem!=thick_elem->parent ) /* check, but not for base classes */
  {
    slice_no=1;
  }

  /* check to see if we've already done this one */
  thin_elem = get_thin(thick_elem,slice_no);
  if (thin_elem) return thin_elem;

  /* issue a warning in case of element parameter combinations not suitable for slicing */
  fint_param   = return_param_recurse("fint",thick_elem);
  if(fint_param)
  {
    printf("    *** warning %s is a thick %s with fringe fields. These will be lost in the translation to a multipole. Use dipedge.\n",
           thick_elem->name,thick_elem->parent->name);
  }

  length_param = return_param_recurse("l",thick_elem);
  angle_param  = return_param_recurse("angle",thick_elem);
  kparam[0]    = return_param_recurse("k0",thick_elem);
  kparam[1]    = return_param_recurse("k1",thick_elem);
  kparam[2]    = return_param_recurse("k2",thick_elem);
  kparam[3]    = return_param_recurse("k3",thick_elem);
  ksparam[0]   = return_param_recurse("k0s",thick_elem);
  ksparam[1]   = return_param_recurse("k1s",thick_elem);
  ksparam[2]   = return_param_recurse("k2s",thick_elem);
  ksparam[3]   = return_param_recurse("k3s",thick_elem);
  kn_param     = return_param_recurse("knl",thick_elem);
  ks_param     = return_param_recurse("ksl",thick_elem);
  if (kn_param) {kn_param = clone_command_parameter(kn_param); knl_flag++;}
  if (ks_param) {ks_param = clone_command_parameter(ks_param); ksl_flag++;}

  /* translate k0,k1,k2,k3,angle */
  if ((kparam[0] || kparam[1] || kparam[2] || kparam[3] || angle_param
       || ksparam[0] || ksparam[1] || ksparam[2] || ksparam[3])
      && (kn_param==NULL && ks_param==NULL))
  {
    kn_param = new_command_parameter("knl", 12);
    kn_param->expr_list = new_expr_list(10);
    kn_param->double_array = new_double_array(10);
    ks_param = new_command_parameter("ksl", 12);
    ks_param->expr_list = new_expr_list(10);
    ks_param->double_array = new_double_array(10);
    angle_conversion = translate_k(kparam,ksparam,angle_param,kn_param,ks_param);
  }

  kn_param = scale_and_slice(kn_param,length_param,slices,slice_no,
                             angle_conversion,knl_flag+ksl_flag);
  ks_param = scale_and_slice(ks_param,length_param,slices,slice_no,
                             angle_conversion,knl_flag+ksl_flag);
  /* set up new multipole command */
  cmd = new_command(buffer("thin_multipole"), 11, 11, /* max num names, max num param */
                    buffer("element"), buffer("none"), 0, 8); /* 0 is link, multipole is 8 */
  add_cmd_parameter_new(cmd,1.,"magnet",0); /* parameter magnet with value of 1 and inf=0 */
  if(!minimizefl)
  {
    add_cmd_parameter_clone(cmd,return_param("at"  ,thick_elem),"at"  ,1);
    add_cmd_parameter_clone(cmd,return_param("from",thick_elem),"from",1);
    add_lrad(cmd,length_param,slices);
    add_cmd_parameter_clone(cmd,kn_param,"knl",1);
    add_cmd_parameter_clone(cmd,ks_param,"ksl",1);
  }
  add_cmd_parameter_clone(cmd,return_param_recurse("apertype",thick_elem),"apertype",1);
  add_cmd_parameter_clone(cmd,return_param_recurse("aperture",thick_elem),"aperture",1);
  add_cmd_parameter_clone(cmd,return_param("bv",thick_elem),"bv",1);
  add_cmd_parameter_clone(cmd,return_param_recurse("tilt",thick_elem),"tilt",1);
  add_cmd_parameter_clone(cmd,return_param_recurse("kmax",    thick_elem),"kmax",    1);
  add_cmd_parameter_clone(cmd,return_param_recurse("kmin",    thick_elem),"kmin",    1);
  add_cmd_parameter_clone(cmd,return_param_recurse("calib",   thick_elem),"calib",   1);
  add_cmd_parameter_clone(cmd,return_param_recurse("polarity",thick_elem),"polarity",1);
  add_cmd_parameter_clone(cmd,return_param_recurse("mech_sep",thick_elem),"mech_sep",1);
  /* create element with this command */
  if (slices==1 && slice_no==1) thin_name=buffer(thick_elem->name);
  else
  {
    thin_name = make_thin_name(thick_elem->name,slice_no);
  }

  if (thin_elem_parent)
  {
    thin_elem = make_element(thin_name,thin_elem_parent->name,cmd,-1);
  }
  else
  {
    thin_elem = make_element(thin_name,"multipole",cmd,-1);
  }
  thin_elem->length = 0;
  thin_elem->bv = el_par_value("bv",thin_elem);
  if (thin_elem_parent && thin_elem_parent->bv)
  {
    thin_elem->bv = thin_elem_parent->bv;
  }
  put_thin(thick_elem,thin_elem,slice_no);
  return thin_elem;
}

struct element* create_thin_solenoid(struct element* thick_elem, int slice_no)
/*hbu create thin solenoid element, similar to create_thin_multipole */
{
  struct command_parameter *length_param, *ks_param, *at_param ,*ks_par;
  struct element *thin_elem_parent, *thin_elem;
  struct command* cmd;
  char *thin_name;
  int slices,minimizefl;

  if (thick_elem == thick_elem->parent) return NULL;
  else
  {
    thin_elem_parent = create_thin_solenoid(thick_elem->parent,slice_no); /*hbu move up to parent */
  }
  thin_elem = get_thin(thick_elem,slice_no);
  if (thin_elem) return thin_elem; /* is already thin */
  slices = get_slices_from_elem(thick_elem);
  /* get parameters from the thick solenoid element */
  length_param  = return_param_recurse("l",thick_elem);
  ks_param      = return_param_recurse("ks",thick_elem);
  at_param      = return_param("at",thick_elem);

  minimizefl=get_option("minimizeparents") && !at_param && thick_elem == thick_elem->parent;
  if(minimizefl)
  {
    slice_no=slices=1; /* do not slice this one */
  }

  /* set up new solenoid command */
  cmd = new_command(buffer("thin_solenoid"), 11, 11, /* max num names, max num param */
                    buffer("element"), buffer("none"), 0, 9); /* 0 is link, solenoid is 9 */  /*hbu trial */
  add_cmd_parameter_new(cmd,1.,"magnet",0); /* parameter magnet with value of 1 and inf=0 */


  if(!minimizefl)
  {
    add_cmd_parameter_clone(cmd,return_param("at"  ,thick_elem),"at"  ,1);
    add_cmd_parameter_clone(cmd,return_param("from",thick_elem),"from",1);
    add_lrad(cmd,length_param,slices);
  }
  add_cmd_parameter_clone(cmd,ks_param,"ks",1); /* keep ks */
  if(!minimizefl)
  {
    if (length_param && ks_param) /* in addition provide   ksi = ks * l /slices */
    {
      ks_par = cmd->par->parameters[cmd->par->curr] = clone_command_parameter(ks_param); /* start from clone of ks */
      strcpy(ks_par->name,"ksi"); /* change name to ksi */
      if (length_param->expr && ks_par->expr) /* first step is ks * l calculation, expression or value */
      {
        ks_par->expr = compound_expr(ks_par->expr,ks_par->double_value,"*",length_param->expr,length_param->double_value); /* multiply expression with length */
      }
      else ks_par->double_value *= length_param->double_value; /* multiply value with length */
      if (slices > 1) /* 2nd step, divide by slices, expression or number */
      {
        if (ks_par->expr) ks_par->expr = compound_expr(ks_par->expr,0.,"/",NULL,slices);
        else ks_par->double_value /= slices;
      }
      add_to_name_list("ksi",1,cmd->par_names);
      cmd->par->curr++;
    }
  }
  add_cmd_parameter_clone(cmd,return_param_recurse("apertype",thick_elem),"apertype",1);
  add_cmd_parameter_clone(cmd,return_param_recurse("aperture",thick_elem),"aperture",1);
  add_cmd_parameter_clone(cmd,return_param("bv",thick_elem),"bv",1);
  add_cmd_parameter_clone(cmd,return_param("tilt",thick_elem),"tilt",1);
  /* create element with this command */
  if (slices==1 && slice_no==1) thin_name=buffer(thick_elem->name);
  else thin_name = make_thin_name(thick_elem->name,slice_no);
  if (thin_elem_parent)
  {
    thin_elem = make_element(thin_name,thin_elem_parent->name,cmd,-1);
  }
  else
  {
    thin_elem = make_element(thin_name,"solenoid",cmd,-1);
  }
  thin_elem->length = 0;
  thin_elem->bv = el_par_value("bv",thin_elem);
  if (thin_elem_parent && thin_elem_parent->bv)
  {
    thin_elem->bv = thin_elem_parent->bv;
  }
  put_thin(thick_elem,thin_elem,slice_no);
  return thin_elem;
}

struct element* create_thin_elseparator(struct element* thick_elem, int slice_no)
/*hbu create thin elseparator element, similar to and made from create_thin_solenoid */
{
  struct command_parameter *length_param, *ex_param, *ey_param, *tilt_param, *at_param ,*ex_par, *ey_par;
  struct element *thin_elem_parent, *thin_elem;
  struct command* cmd;
  char *thin_name;
  int slices,minimizefl;

  if (thick_elem == thick_elem->parent) return NULL;
  else
  {
    thin_elem_parent = create_thin_elseparator(thick_elem->parent,slice_no); /*hbu move up to parent */
  }
  thin_elem = get_thin(thick_elem,slice_no);
  if (thin_elem) return thin_elem; /* is already thin */
  slices = get_slices_from_elem(thick_elem);
  /* get parameters from the thick elseparator element */
  length_param  = return_param_recurse("l",thick_elem);
  ex_param      = return_param_recurse("ex",thick_elem);
  ey_param      = return_param_recurse("ey",thick_elem);
  tilt_param    = return_param_recurse("tilt",thick_elem);
  at_param      = return_param("at",thick_elem);

  minimizefl=get_option("minimizeparents") && !at_param && thick_elem == thick_elem->parent;
  if(minimizefl)
  {
    slice_no=slices=1; /* do not slice this one */
  }

  /* set up new solenoid command */
  cmd = new_command(buffer("thin_elseparator"), 11, 11, /* max num names, max num param */
                    buffer("element"), buffer("none"), 0, 11); /* 0 is link, elseparator is 11 */  /*hbu trial */
  add_cmd_parameter_new(cmd,1.,"magnet",0); /* parameter magnet with value of 1 and inf=0 */


  if(!minimizefl)
  {
    add_cmd_parameter_clone(cmd,return_param("at"  ,thick_elem),"at"  ,1);
    add_cmd_parameter_clone(cmd,return_param("from",thick_elem),"from",1);
    add_lrad(cmd,length_param,slices);
  }
  add_cmd_parameter_clone(cmd,ex_param,"ex",1); /* keep ex */
  add_cmd_parameter_clone(cmd,ey_param,"ey",1); /* keep ey */
  add_cmd_parameter_clone(cmd,tilt_param,"tilt",1); /* keep tilt */
  if(!minimizefl)
  {
    /* create ex_l from ex */
    if (length_param && ex_param) /* in addition provide   ex_l = ex * l /slices */
    {
      ex_par = cmd->par->parameters[cmd->par->curr] = clone_command_parameter(ex_param); /* start from clone of ex */
      strcpy(ex_par->name,"ex_l"); /* change name to ex_l */
      if (length_param->expr && ex_par->expr) /* first step is ex * l calculation, expression or value */
      {
        ex_par->expr = compound_expr(ex_par->expr,ex_par->double_value,"*",length_param->expr,length_param->double_value); /* multiply expression with length */
      }
      else ex_par->double_value *= length_param->double_value; /* multiply value with length */
      if (slices > 1) /* 2nd step, divide by slices, expression or number */
      {
        if (ex_par->expr) ex_par->expr = compound_expr(ex_par->expr,0.,"/",NULL,slices);
        else ex_par->double_value /= slices;
      }
      add_to_name_list("ex_l",1,cmd->par_names);
      cmd->par->curr++;
    }
    /* create ey_l from ey */
    if (length_param && ey_param) /* in addition provide   ey_l = ey * l /slices */
    {
      ey_par = cmd->par->parameters[cmd->par->curr] = clone_command_parameter(ey_param); /* start from clone of ey */
      strcpy(ey_par->name,"ey_l"); /* change name to ey_l */
      if (length_param->expr && ey_par->expr) /* first step is ey * l calculation, expression or value */
      {
        ey_par->expr = compound_expr(ey_par->expr,ey_par->double_value,"*",length_param->expr,length_param->double_value); /* multiply expression with length */
      }
      else ey_par->double_value *= length_param->double_value; /* multiply value with length */
      if (slices > 1) /* 2nd step, divide by slices, expression or number */
      {
        if (ey_par->expr) ey_par->expr = compound_expr(ey_par->expr,0.,"/",NULL,slices);
        else ey_par->double_value /= slices;
      }
      add_to_name_list("ey_l",1,cmd->par_names);
      cmd->par->curr++;
    }
  }
  add_cmd_parameter_clone(cmd,return_param_recurse("apertype",thick_elem),"apertype",1);
  add_cmd_parameter_clone(cmd,return_param_recurse("aperture",thick_elem),"aperture",1);
  add_cmd_parameter_clone(cmd,return_param("bv",thick_elem),"bv",1);
  add_cmd_parameter_clone(cmd,return_param("tilt",thick_elem),"tilt",1);
  /* create element with this command */
  if (slices==1 && slice_no==1) thin_name=buffer(thick_elem->name);
  else thin_name = make_thin_name(thick_elem->name,slice_no);
  if (thin_elem_parent)
  {
    thin_elem = make_element(thin_name,thin_elem_parent->name,cmd,-1);
  }
  else
  {
    thin_elem = make_element(thin_name,"elseparator",cmd,-1);
  }
  thin_elem->length = 0;
  thin_elem->bv = el_par_value("bv",thin_elem);
  if (thin_elem_parent && thin_elem_parent->bv)
  {
    thin_elem->bv = thin_elem_parent->bv;
  }
  put_thin(thick_elem,thin_elem,slice_no);
  return thin_elem;
}

/* put in one of those nice marker kind of things */
struct node* new_marker(struct node *thick_node, double at, struct expression *at_expr)
{
  struct node* node=NULL;
  struct element* elem=NULL;

  int pos;
  struct command* p;
  struct command* clone;

  if (thick_node->p_elem)
  {
    pos = name_list_pos("marker", defined_commands->list);
    /* clone = clone_command(defined_commands->commands[pos]); */
    p = defined_commands->commands[pos];
    clone = new_command(p->name, 10, 10, p->module, p->group, p->link_type,p->mad8_type);
    /* 10 is the maximum number of par names and parvalues */
    add_cmd_parameter_clone(clone,return_param_recurse("at",      thick_node->p_elem),"at",      1);
    add_cmd_parameter_clone(clone,return_param_recurse("from",    thick_node->p_elem),"from",    1);
    add_cmd_parameter_clone(clone,return_param_recurse("apertype",thick_node->p_elem),"apertype",1);
    add_cmd_parameter_clone(clone,return_param_recurse("aperture",thick_node->p_elem),"aperture",1);
    add_cmd_parameter_clone(clone,return_param_recurse("aper_tol",thick_node->p_elem),"aper_tol",1);
    add_cmd_parameter_clone(clone,return_param_recurse("kmax",    thick_node->p_elem),"kmax",    1);
    add_cmd_parameter_clone(clone,return_param_recurse("kmin",    thick_node->p_elem),"kmin",    1);
    add_cmd_parameter_clone(clone,return_param_recurse("calib",   thick_node->p_elem),"calib",   1);
    add_cmd_parameter_clone(clone,return_param_recurse("polarity",thick_node->p_elem),"polarity",1);
    add_cmd_parameter_clone(clone,return_param_recurse("mech_sep",thick_node->p_elem),"mech_sep",1);
    elem = make_element(thick_node->p_elem->name, "marker", clone,-1);
    node = new_elem_node(elem, thick_node->occ_cnt);
    strcpy(node->name, thick_node->name);
    node->occ_cnt = thick_node->occ_cnt;
    node->at_value = at;
    if (at_expr) node->at_expr = clone_expression(at_expr);
    node->from_name = thick_node->from_name;
  }
  else
  {
    fatal_error("Oh dear, this is not an element!",thick_node->name);
  }

  return node;
}

/* adds a thin elem in sliced nodes to the end of a sequence */
void seq_diet_add_elem(struct node* node, struct sequence* to_sequ)
{
  struct command_parameter *at_param, *length_param;
  struct expression *l_expr = NULL, *at_expr = NULL;
  struct node* thin_node;
  struct element* elem;
  double length = 0, at = 0;
  int i,middle=-1,slices = 1;
  char* old_thin_style;

  old_thin_style = NULL;
  if (strstr(node->base_name,"collimator"))
  {
    elem = create_thin_obj(node->p_elem,1);
    old_thin_style = thin_style;
    thin_style = collim_style;
  }
  else if (strstr(node->base_name,"solenoid"))
  {
    elem = create_thin_solenoid(node->p_elem,1);  /* create the first thin solenoid slice */
  }
  else if (strstr(node->base_name,"elseparator"))
  {
    elem = create_thin_elseparator(node->p_elem,1);  /* create the first thin elseparator slice */
  }
  else
  {
    elem = create_thin_multipole(node->p_elem,1); /* get info from first slice */
  }
  slices = get_slices_from_elem(node->p_elem); /*hbu June 2005 */

  at_param = return_param_recurse("at",elem);
  length_param = return_param_recurse("l",node->p_elem); /*get original length*/
  if (length_param) l_expr  = length_param->expr;
  if (at_param)     at_expr = at_param->expr;

  at     = el_par_value_recurse("at", elem);
  length = el_par_value_recurse("l",node->p_elem);

  if (node->at_expr) at_expr = node->at_expr;
  if (node->at_value != zero) at = node->at_value;
  if (node->length   != zero) length = node->length;
  /* note that a properly created clone node will contain the length of the element */
  /* this will override all other definitions and hence the already sliced element length
     is irrelevant */

  if (slices>1)
  { /* sets after which element I should put the marker */
    middle = abs(slices/2);
  }

  for (i=0; i<slices; i++)
  {
    if (strstr(node->base_name,"collimator"))
    {
      elem = create_thin_obj(node->p_elem,i+1);
    }
    else if (strstr(node->base_name,"solenoid"))
    {
      elem = create_thin_solenoid(node->p_elem,i+1);
    }
    else if (strstr(node->base_name,"elseparator"))
    {
      elem = create_thin_elseparator(node->p_elem,i+1);
    }
    else
    {
      elem = create_thin_multipole(node->p_elem,i+1);
    }
    thin_node = new_elem_node(elem, node->occ_cnt);
    thin_node->length   = 0.0;
    thin_node->from_name = buffer(node->from_name);
    if (fabs(at_shift(slices,i+1))>0.0)
    {
      if (at_expr || l_expr)
      {
        thin_node->at_expr =
          compound_expr(at_expr,at,"+",scale_expr(l_expr,at_shift(slices,i+1)),
                        length*at_shift(slices,i+1));
      }
    }
    else
    {
      if (at_expr) thin_node->at_expr = clone_expression(at_expr);
    }
    thin_node->at_value = at + length*at_shift(slices,i+1);
    if (i==middle) seq_diet_add(new_marker(node,at,at_expr),to_sequ);
    seq_diet_add(thin_node,to_sequ);
  }
  if (strstr(node->base_name,"collimator")) thin_style=old_thin_style;
  return;
}

/* creates the thin non-magnetic element - recursively */
struct element* create_thin_obj(struct element* thick_elem, int slice_no)
{
  struct element *thin_elem_parent = NULL , *thin_elem = NULL;
  struct command* cmd = NULL;
  struct command_parameter*  length_param= NULL;
  int length_i = -1,lrad_i = -1,slices=1;
  char* thin_name = NULL;

  if (thick_elem == thick_elem->parent)
  {
    return NULL;
  }
  else
  {
    thin_elem_parent = create_thin_obj(thick_elem->parent,slice_no);
  }

  /* check to see if we've already done this one */
  thin_elem = get_thin(thick_elem,slice_no);
  if (thin_elem) return thin_elem;

  /* set up new multipole command */
  cmd = clone_command(thick_elem->def);
  length_param = return_param_recurse("l",thick_elem);
  length_i = name_list_pos("l",thick_elem->def->par_names);
  lrad_i   = name_list_pos("lrad",thick_elem->def->par_names);
  if (length_param)
  {
    if (lrad_i > -1 && thick_elem->def->par_names->inform[lrad_i]>0)
    {
      /* already exists so replace lrad */
      cmd->par->parameters[lrad_i]->double_value = cmd->par->parameters[length_i]->double_value;
      if (cmd->par->parameters[length_i]->expr)
      {
        if (cmd->par->parameters[lrad_i]->expr)
          delete_expression(cmd->par->parameters[lrad_i]->expr);
        cmd->par->parameters[lrad_i]->expr =
          clone_expression(cmd->par->parameters[length_i]->expr);
      }
    }
    else
    { /* doesn't exist */
      if (name_list_pos("lrad",thick_elem->base_type->def->par_names)>-1)
      {
        /* add lrad only if allowed by element */
        if (cmd->par->curr == cmd->par->max) grow_command_parameter_list(cmd->par);
        if (cmd->par_names->curr == cmd->par_names->max)
          grow_name_list(cmd->par_names);
        cmd->par->parameters[cmd->par->curr] = clone_command_parameter(length_param);
        add_to_name_list("lrad",1,cmd->par_names);
        cmd->par->parameters[name_list_pos("lrad",cmd->par_names)]->expr =
          clone_expression(cmd->par->parameters[length_i]->expr);
        cmd->par->curr++;
      }
    }
  }

  if (length_i > -1)
  {
    cmd->par->parameters[length_i]->double_value = 0;
    cmd->par->parameters[length_i]->expr = NULL;
  }
  if (strstr(thick_elem->base_type->name,"collimator"))
  {
    slices = get_slices_from_elem(thick_elem);
  }
  if (slices==1 && slice_no==1)
  {
    thin_name=buffer(thick_elem->name);
  }
  else
  {
    thin_name=make_thin_name(thick_elem->name,slice_no);
  }

  if (thin_elem_parent)
  {
    thin_elem = make_element(thin_name,thin_elem_parent->name,cmd,-1);
  }
  else
  {
    thin_elem = make_element(thin_name,thick_elem->base_type->name,cmd,-1);
  }
  thin_elem->length = 0;
  thin_elem->bv = el_par_value("bv",thin_elem);

  put_thin(thick_elem,thin_elem,slice_no);
  return thin_elem;
}

/* this copies an element node and sets the length to zero
   and radiation length to the length
   to be used for "copying" optically neutral elements */
struct node* copy_thin(struct node* thick_node)
{
  struct node* thin_node = NULL;

  thin_node = clone_node(thick_node, 0);
  thin_node->length=0;
  thin_node->p_elem->length=0;
  /* if we have a non zero length then an lrad has to be created */
  if (el_par_value("l",thick_node->p_elem)>zero)
    thin_node->p_elem = create_thin_obj(thick_node->p_elem,1);

  return thin_node;
}

/* this decides how to split an individual node and
   sends it onto the thin_sequ builder */
void seq_diet_node(struct node* thick_node, struct sequence* thin_sequ)
{
  struct node* thin_node;
  if (thick_node->p_elem)
  { /* this is an element to split and add */
    if (el_par_value("l",thick_node->p_elem)==zero) /* if it's already thin copy it directly*/
    {
      seq_diet_add(thin_node = copy_thin(thick_node),thin_sequ);
    }
    else if(strcmp(thick_node->base_name,"matrix") == 0)
    { /*hbu. Take matrix as it is, including any length */
      seq_diet_add(thick_node,thin_sequ);
    }
    else
    { /* we have to slim it down a bit...*/
      if (strcmp(thick_node->base_name,"marker") == 0    ||
          strcmp(thick_node->base_name,"instrument") == 0  ||
          strcmp(thick_node->base_name,"hmonitor") == 0    ||
          strcmp(thick_node->base_name,"vmonitor") == 0    ||
          strcmp(thick_node->base_name,"monitor") == 0     ||
          strcmp(thick_node->base_name,"vkicker") == 0     ||
          strcmp(thick_node->base_name,"hkicker") == 0     ||
          strcmp(thick_node->base_name,"kicker") == 0      ||
          strcmp(thick_node->base_name,"rfcavity") == 0    ||
	  strcmp(thick_node->base_name,"crabcavity") == 0
        )
      {
        seq_diet_add(thin_node = copy_thin(thick_node),thin_sequ);
        /*   delete_node(thick_node); */
        /* special cavity list stuff */
        if (strcmp(thin_node->p_elem->base_type->name, "rfcavity") == 0 &&
            find_element(thin_node->p_elem->name, thin_sequ->cavities) == NULL)
          add_to_el_list(&thin_node->p_elem, 0, thin_sequ->cavities, 0);
	/* special crab cavity list stuff */
        if (strcmp(thin_node->p_elem->base_type->name, "crabcavity") == 0 &&
            find_element(thin_node->p_elem->name, thin_sequ->crabcavities) == NULL)
          add_to_el_list(&thin_node->p_elem, 0, thin_sequ->crabcavities, 0);
      }
      else if (strcmp(thick_node->base_name,"rbend") == 0 ||
               strcmp(thick_node->base_name,"sbend") == 0       ||
               strcmp(thick_node->base_name,"quadrupole") == 0  ||
               strcmp(thick_node->base_name,"sextupole") == 0   ||
               strcmp(thick_node->base_name,"octupole") == 0    ||
               strcmp(thick_node->base_name,"solenoid") == 0    || /*hbu */
               strcmp(thick_node->base_name,"multipole") == 0
               || /* special spliting required. */
               strcmp(thick_node->base_name,"rcollimator") == 0 ||
               strcmp(thick_node->base_name,"ecollimator") == 0 ||
               strcmp(thick_node->base_name,"elseparator") == 0
        )
      {
        seq_diet_add_elem(thick_node,thin_sequ);
        /*   delete_node(thick_node); */
      }
      else if (strcmp(thick_node->base_name,"drift") == 0)
      {
        /* ignore this as it makes no sense to slice */
      }
      else
      {
        fprintf(prt_file, "Found unknown basename %s, doing copy with length set to zero.\n",thick_node->base_name);
        seq_diet_add(copy_thin(thick_node),thin_sequ);
        /*        delete_node(thick_node); */
      }
    }
  }
  else if (thick_node->p_sequ)
  { /* this is a sequence to split and add */
    seq_diet_add_sequ(thick_node,seq_diet(thick_node->p_sequ),thin_sequ);
  }
  else
  { /* we have no idea what this is - serious error */
    fatal_error("node is not element or sequence",thick_node->base_name);
  }
}

/* slim down this sequence - this is the bit to be called recursively */
/* this actually creates the thin sequence */
struct sequence* seq_diet(struct sequence* thick_sequ)
{
  struct node *thick_node = NULL;
  struct sequence* thin_sequ;
  char name[128];
  int pos;

  /* first check to see if it had been already sliced */
  if ((thin_sequ=get_thin_sequ(thick_sequ))) return thin_sequ;

  strcpy(name,thick_sequ->name);
  fprintf(prt_file, "makethin: slicing sequence : %s\n",name);
  thin_sequ = new_sequence(name, thick_sequ->ref_flag);
  thin_sequ->start = NULL;
  thin_sequ->share = thick_sequ->share;
  thin_sequ->nested = thick_sequ->nested;
  thin_sequ->length = sequence_length(thick_sequ);
  thin_sequ->refpos = buffer(thick_sequ->refpos);
  thin_sequ->ref_flag = thick_sequ->ref_flag;
  thin_sequ->beam = thick_sequ->beam;
  if (thin_sequ->cavities != NULL)  thin_sequ->cavities->curr = 0;
  else thin_sequ->cavities = new_el_list(100);
  if (thin_sequ->crabcavities != NULL)  thin_sequ->crabcavities->curr = 0;
  else thin_sequ->crabcavities = new_el_list(100);
  thick_node = thick_sequ->start;
  while(thick_node != NULL)
  { /* loop over current sequence */
    /* the nodes are added to the sequence in seq_diet_add() */
    seq_diet_node(thick_node,thin_sequ);
    if (thick_node == thick_sequ->end)  break;
    thick_node = thick_node->next;
  }
  thin_sequ->end->next = thin_sequ->start;
  /* now we have to move the pointer in the sequences list
     to point to our thin sequence */
  if ((pos = name_list_pos(name, sequences->list)) < 0)
  {
    fatal_error("unknown sequence sliced:", name);
  }
  else
  {
    sequences->sequs[pos]= thin_sequ;
    /* delete_sequence(thick_sequ) */
  }

  /* add to list of sequences sliced */
  put_thin_sequ(thick_sequ,thin_sequ);

  return thin_sequ;
}

/* This converts the MAD-X command to something I can use
   if a file has been specified we send the command to exec_save
   which writes the file for us */
void makethin(struct in_cmd* cmd)
{
  struct sequence *thick_sequ = NULL ,*thin_sequ = NULL;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  char *name = NULL;
  int pos,pos2;
  int k=0;
/*    time_t start; */

/*    start = time(NULL); */
  pos = name_list_pos("style", nl);
  if (nl->inform[pos] && (name = pl->parameters[pos]->string))
  {
    thin_style = buffer(pl->parameters[pos]->string);
    fprintf(prt_file, "makethin: style chosen : %s\n",thin_style);
  }

  /* first check makethin parameters which influence the selection */

  pos = name_list_pos("minimizeparents", nl);
  /* k = true; */  /* Use this to set minimizeparents to true by default. */
  if( pos > -1 && nl->inform[pos])
  {
    k=pl->parameters[pos]->double_value;
  }
  set_option("minimizeparents", &k);

  pos = name_list_pos("makeconsistent", nl);
  if( pos > -1 && nl->inform[pos])
  {
    k=pl->parameters[pos]->double_value;
    set_option("makeconsistent", &k);
  }

  if (slice_select->curr > 0)
  {
    set_selected_elements(); /* makethin selection */
    thin_select_list = selected_elements;
  }
  if (thin_select_list == NULL)
  {
    warning("makethin: no selection list,","slicing all to one thin lens.");
  }
  else if (thin_select_list->curr == 0)
  {
    warning("makethin selection list empty,","slicing all to one thin lens.");
  }
  if(get_option("makeconsistent"))
  {
    force_consistent_slices();
  }
  pos = name_list_pos("sequence", nl);
  if (nl->inform[pos] && (name = pl->parameters[pos]->string))
  {
    if ((pos2 = name_list_pos(name, sequences->list)) >= 0)
    {
      thick_sequ = sequences->sequs[pos2];
      thin_sequ = seq_diet(thick_sequ);
      disable_line(thin_sequ->name, line_list);
    }
    else warning("unknown sequence ignored:", name);
  }
  else warning("makethin without sequence:", "ignored");

  /* fprintf(prt_file, "makethin: finished in %f seconds.\n",difftime(time(NULL),start)); */
  thin_select_list = NULL;
}

/*************************************************************************/
/* these are the routines to determine the method of splitting */
/* note slice number is counted from 1 NOT 0 */

/* return at relative shifts from center of unsliced magnet */
double simple_at_shift(int slices,int slice_no)
{
  double at = 0;
  at = ((double) 2*slice_no-1)/((double) 2*slices)-0.5;
  return at;
}

double teapot_at_shift(int slices,int slice_no)
{
  double at = 0;
  switch (slices)
  {
    case 1:
      at = 0.;
      break;
    case 2:
      if (slice_no == 1) at = -1./3.;
      if (slice_no == 2) at = +1./3.;
      break;
    case 3:
      if (slice_no == 1) at = -3./8.;
      if (slice_no == 2) at = 0.;
      if (slice_no == 3) at = +3./8.;
      break;
    case 4:
      if (slice_no == 1) at = -2./5.;
      if (slice_no == 2) at = -2./15.;
      if (slice_no == 3) at = +2./15.;
      if (slice_no == 4) at = +2./5.;
      break;
  }
  /* return the simple style if slices > 4 */
  if (slices > 4) at = simple_at_shift(slices,slice_no);
  return at;
}

double collim_at_shift(int slices,int slice_no)
{
  double at = 0;
  if (slices==1)
  {
    at = 0.0;
  }
  else
  {
    at = (slice_no-1.0)/(slices-1.0)-0.5;
  }
  return at;
}

/* return at relative strength shifts from unsliced magnet */
double teapot_q_shift(int slices,int slice_no)
{
  return 1./slices;
}

double simple_q_shift(int slices,int slice_no)
{
  return 1./slices;
}

double collim_q_shift(int slices,int slice_no)
{ /* pointless actually, but it pleases symmetrically */
  return 1./slices;
}


/* return at relative shifts from center of unsliced magnet */
double at_shift(int slices,int slice_no)
{
  if (thin_style == NULL || strcmp(thin_style,"teapot")==0)
  {
    return teapot_at_shift(slices,slice_no);
  }
  else if (strcmp(thin_style,"simple")==0)
  {
    return simple_at_shift(slices,slice_no);
  }
  else if (strcmp(thin_style,"collim")==0)
  {
    return collim_at_shift(slices,slice_no);
  }
  else
  {
    fatal_error("makethin: Style chosen not known:",thin_style);
  }
  return 0;
}

/* return at relative strength shifts from unsliced magnet */
double q_shift(int slices,int slice_no)
{
  if (thin_style == NULL || strcmp(thin_style,"teapot")==0)
  {
    return teapot_q_shift(slices,slice_no);
  }
  else if (strcmp(thin_style,"simple")==0)
  {
    return simple_q_shift(slices,slice_no);
  }
  else if (strcmp(thin_style,"collim")==0)
  {
    return collim_q_shift(slices,slice_no);
  }
  else
  {
    fatal_error("makethin: Style chosen not known:",thin_style);
  }
  return 0;
}

void set_selected_elements()
{ /*hbu June 2005.  New set_selected_elements */
  struct element* el_j;
  struct name_list* nl;
  struct command_parameter_list* pl;
  int i, j, pos_slice, pos_full, pos_range, slice, el_j_slice_pos;
  bool full_fl,range_fl,slice_fl;
  struct node* c_node;    /* for range check.  current node */
  struct node* nodes[2];  /* for range check.  first and last in range */
  /* Init curr and list->curr in global el_list structure.  selected_elements is passed to add_to_el_list and used at the end as thin_select_list
     selected_elements  is only used in makethin (set here and read in and could be named thin_select_list
  */
  selected_elements->curr = 0;
  selected_elements->list->curr = 0;  /* Reset list->curr in global el_list structure.   selected_elements is passed to add_to_el_list */
  if (current_sequ == NULL || current_sequ->ex_start == NULL) /* check that there is an active sequence, otherwise crash in get_ex_range */
  {
    warning("makethin selection without active sequence,", "ignored");
    return;
  }
  /* default is full sequence from start to end */
  nodes[0] = current_sequ->ex_start;
  nodes[1] = current_sequ->ex_end;
  for (i = 0; i < slice_select->curr; i++) /* loop over "select,flag=makethin" commands */
  {
    nl = slice_select->commands[i]->par_names;
    pl = slice_select->commands[i]->par;
    pos_full  = name_list_pos("full", nl);
    full_fl   = pos_full  > -1 && nl->inform[pos_full];  /* selection with full */
    pos_range = name_list_pos("range", nl);
    range_fl  = pos_range > -1 && nl->inform[pos_range]; /* selection with range */
    pos_slice = name_list_pos("slice", nl);              /* position of slice parameter in select command list */
    slice_fl  = pos_slice > -1 && nl->inform[pos_slice]; /* selection with slice */
    if (slice_fl) slice = pl->parameters[pos_slice]->double_value; /* Parameter has been read. Slice number from select command */
    else slice = 1;
    if(full_fl) /* use full sequence from start to end, the default */
    {
      nodes[0] = current_sequ->ex_start;
      nodes[1] = current_sequ->ex_end;
    }
    if(range_fl)
    {
      if (current_sequ == NULL || current_sequ->ex_start == NULL) /* check that there is an active sequence, otherwise crash in get_ex_range */
      {
        warning("makethin range selection without active sequence,", "ignored");
        return;
      }
      if( get_ex_range(pl->parameters[pos_range]->string, current_sequ, nodes) == 0) /* set start nodes[0] and end notes[1] depending on the range string */
      {
        printf("    +++ warning, empty range");
        continue;
      }
    }
    if(slice_fl) /* Set slice number in elements. Add to list of selected_elements */
    {
      if(range_fl) /* now elements in the sequence in the range */
      {
        c_node = nodes[0];
        while (c_node != NULL) /* loop over nodes in range,  set slice number in elements */
        {
          el_j = c_node->p_elem;
          el_j_slice_pos = name_list_pos("slice",el_j->def->par_names); /* position of slice parameter in element list */
          if (pass_select(el_j->name, slice_select->commands[i]) != 0) /* selection on class and pattern done in pass_select. element el_j selected */
          { /* the element el_j passes the selection */
            if(el_j_slice_pos > 0) el_j->def->par->parameters[el_j_slice_pos]->double_value=slice; /* Set the element slice number to the number of slices given in the select statement. */
            if( name_list_pos(el_j->name, selected_elements->list) < 0) /* el_j not yet in selected_elements */
            {
              add_to_el_list(&el_j, slice, selected_elements, 0);
            } /* new selection */
          } /* selection */
          if (c_node == nodes[1]) break; /* done with last node */
          c_node = c_node->next;
        } /* end of while loop over nodes in range */
      } /* range_fl */
      else /* no range_fl */
      {
        for(j=0; j< element_list->curr; j++) /* loop over element_list */
        {
          el_j = element_list->elem[j];
          el_j_slice_pos = name_list_pos("slice",el_j->def->par_names);
          if (pass_select(el_j->name, slice_select->commands[i]) != 0) /* selection on class and pattern done in pass_select. element el_j selected */
          { /* the element el_j passes the selection */
            if(el_j_slice_pos > 0) el_j->def->par->parameters[el_j_slice_pos]->double_value=slice; /* Set the element slice number to the number of slices given in the select statement. */
            if( name_list_pos(el_j->name, selected_elements->list) < 0) /* el_j not yet in selected_elements */
            {
              add_to_el_list(&el_j, slice, selected_elements, 0);
            } /* new selection */
          } /* selection */
        } /* loop over element_list */
      } /* range_fl */
    } /* slice_fl */
  } /* end of loop over select slice commands */
}

/*************************************************************************/
