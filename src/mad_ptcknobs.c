#include "madx.h"

/*______________________________________________________________
  matchptcknobs.c
  Piotr Skowronski (CERN) 2006
  MAD-X module responsible for matching with help of PTC knobs (parameters)

  Feb.2007: added proper support for variable limits

  It implements MAD-X the command:
  match, useptcknobs=true;
  (...)
  endmatch;

  It is basically a macro that performs automatically a set of commands
  depicted above as (...) within a rather complicated loop (points 3-17 below).
  It enables fast matching of any order parameters with PTC.

  The algorithm:

  1. Buffer the key commands (ptc_varyknob, constraint, ptc_setswitch, ptc_twiss or ptc_normal, etc)
  appearing between
  match, useptcknobs=true;
  and any of matching actions calls (migrad,lmdif,jacobian, etc)

  2.  When matching action appears,
  a) set "The Current Variables Values" (TCVV) to zero
  b) perform THE LOOP, i.e. points 3-17

  3.  Prepare PTC environment (ptc_createuniverse, ptc_createlayout)

  4.  Set the user defined knobs (with ptc_knob).

  5.  Set TCVV using ptc_setfieldcomp command

  6.  Run a PTC command (twiss or normal)

  7.  Run a runtime created script that performs a standard matching
  All the user defined knobs are variables of this matching.

  8.  Evaluate constraints expressions to get the matching function vector (I)

  9.  Add the matched values to TCVV

  10. End PTC session (run ptc_end)

  11. If the matched values are are not close enough to zeroes then goto 3

  12. Prepare PTC environment (ptc_createuniverse, ptc_createlayout)

  13. Set TCVV using ptc_setfieldcomp command
  ( --- please note that knobs are not set in this case  -- )
  14. Run a PTC command (twiss or normal)

  15. Evaluate constraints expressions to get the matching function vector (II)

  16. Evaluate a penalty function that compares matching function vectors (I) and (II)
  See points 7 and 14

  17  If the matching function vectors are not similar to each other
  within requested precision then goto 3

  18. Print TCVV, which are the matched values.

  ______________________________________________________________*/

#define MAX_CONTRAINS  100
#define MAX_KNOBS  100
#define COMM_LENGTH  500

struct madx_mpk_variable
{
    char*  name; /*variable name, if initial parameter just copy with mpk_ prefix, if field comp then {elname}_{"kn"/"ks"}{number}*/
    char*  namecv;
    double upper;
    double lower;
    double trustrange;
    double step;
    int    knobidx; /*index of the knob that makes paramter of this variable */
    double currentvalue;
    double oldvalue;
    int    kn; /*if a field component then different from -1*/
    int    ks;
    char   IsIniCond; /**/
};

struct madx_mpk_knob
{
    char* elname;  /*if a field property element name */
    char* initial; /*if initial condition, name of a variable */
    int*  KN; /*array with filed components*/
    int   nKN;
    int*  KS;
    int   nKS;
    int   exactnamematch;
};

// private variables

// static int                        madx_mpk_debug;

// constraints
static int                        madx_mpk_Nconstraints;
static char*                      madx_mpk_constraints[MAX_CONTRAINS];

// knobs
static int                        madx_mpk_Nknobs;
static struct madx_mpk_knob       madx_mpk_knobs[MAX_KNOBS];
static char*                      madx_mpk_setknobs[MAX_KNOBS]; /*this is ptc_setknobs*/
// static struct in_cmd*             madx_mpk_varyknob_cmds[MAX_KNOBS]; /*this is ptc_setknobs*/

// other
static int                        madx_mpk_Nvariables;  /*total number if matching variables >= number of knob commands (one knob command may define more then one knob/variable)*/
static struct madx_mpk_variable   madx_mpk_variables[MAX_KNOBS];

static struct in_cmd* madx_mpk_comm_createuniverse;
static struct in_cmd* madx_mpk_comm_createlayout;
static struct in_cmd* madx_mpk_comm_setswitch;
static struct in_cmd* madx_mpk_comm_calculate;/*ptc_twiss or ptc_normal*/

// static char twisscommand[]="ptc_twiss, table=ptc_twiss, icase=6, no=2, betx=10, alfx=.0,  bety=10., alfy=0, betz=10., alfz=0;";

// private functions

static int
factorial(int v)
{
  if (v <= 1 )
  {
    return  1;
  }
  else
  {
    return v*factorial(v-1);
  }
}

static int
findsetknob(char* ename, int exactnamematch, char* initialpar)
{
  /*
    finds a knob that corresponds to this name, and detects eventual errors
    if not found returns fresh empty knob
  */
  int i;
  int cmpres;
  int result = madx_mpk_Nknobs; /*a new empty one*/
  char* bconta;
  char* acontb;

  if (ename)
  {
    for (i = 0; i < madx_mpk_Nknobs; i++)
    {
      if (madx_mpk_knobs[i].elname == 0x0) continue;
      cmpres = strcmp(ename,madx_mpk_knobs[i].elname);
      if (cmpres == 0)
      {
        if ( exactnamematch == madx_mpk_knobs[i].exactnamematch )
        {
          result = i; /*the same element or family*/
          continue; /*to find out if there is any setting problem*/
        }
        else
        {
          mad_error("findsetknob","A knob for such named element(s) found, but name matching flag does not agree.");
          return -i;
        }
      }

      /*
        if a not-exact and b exact,  b can not contain a
        if a not-exact and b non-exact,  b can not contain a and a can not contain b
        if a exact, and b not-exact, a can not contain b
      */
      acontb = strstr(ename, madx_mpk_knobs[i].elname);
      bconta = strstr(madx_mpk_knobs[i].elname, ename);
      if ( (acontb && (exactnamematch == 0)) || (bconta && (madx_mpk_knobs[i].exactnamematch  == 0)) )
      {
        mad_error("findsetknob",
              "This variable (name %s, exactmatch %d) can cause ambiguity with another already defined variable (name %s, exactmatch %d)",
              ename, exactnamematch, madx_mpk_knobs[i].elname, madx_mpk_knobs[i].exactnamematch);
        return -i;
      }
    }
  }

  if (initialpar)
  {
    for (i = 0; i < madx_mpk_Nknobs; i++)
    {
      if (madx_mpk_knobs[i].initial)
      {
        if (( strcmp(initialpar,madx_mpk_knobs[i].initial) == 0 ))
        {

          mad_error("findsetknob","Such initial parameter is already defined");
          return -i;
        }
      }
    }
  }

  if (result == madx_mpk_Nknobs)
  {
    madx_mpk_Nknobs++;  /*we have not found a setknob to add this vary knob*/
  }

  return result;
}

static int
mapptctomad(char* ptcname, char* madxname)
{

  if( strcmp(ptcname,"beta11") == 0 )
  {
    strcpy(madxname,"betx");
    return 0;
  }

  if( strcmp(ptcname,"beta22") == 0 )
  {
    strcpy(madxname,"bety");
    return 0;
  }

  if( strcmp(ptcname,"beta33") == 0 )
  {
    strcpy(madxname,"betz");
    return 0;
  }

  if( strcmp(ptcname,"alfa11") == 0 )
  {
    strcpy(madxname,"alfx");
    return 0;
  }

  if( strcmp(ptcname,"alfa22") == 0 )
  {
    strcpy(madxname,"alfy");
    return 0;
  }

  if( strcmp(ptcname,"alfa33") == 0 )
  {
    strcpy(madxname,"alfz");
    return 0;
  }


  if( strcmp(ptcname,"disp1") == 0 )
  {
    strcpy(madxname,"dx");
    return 0;
  }

  if( strcmp(ptcname,"disp2") == 0 )
  {
    strcpy(madxname,"dpz");
    return 0;
  }

  if( strcmp(ptcname,"disp3") == 0 )
  {
    strcpy(madxname,"dy");
    return 0;
  }

  if( strcmp(ptcname,"disp4") == 0 )
  {
    strcpy(madxname,"dpy");
    return 0;
  }


  return 1;
}

static int
readstartvalues(void)
{
/*
  reads the starting values of variable to initializ currentvalue
*/

  int i;
  int n;
  int ncomp;
  double nvalue;
  struct node*                node = 0x0;
  struct madx_mpk_variable*      v = 0x0;
  struct madx_mpk_knob*         kn = 0x0;
  char  buff[50];
  char* p;

  if (debuglevel)  printf("\n\n\n  READING INITIAL VALUES \n\n\n");

  for (i = 0; i < madx_mpk_Nvariables; i++)
  {
    v = &(madx_mpk_variables[i]);
    kn = &(madx_mpk_knobs[v->knobidx]);

    if (v->IsIniCond)
    {
      mapptctomad(kn->initial,buff);
      v->currentvalue = command_par_value(buff, madx_mpk_comm_calculate->clone);
      if (debuglevel)  printf("Initialized current value for %s to %f\n", kn->initial,v->currentvalue);
    }
    else if(kn->exactnamematch == 0)
    {
      if (debuglevel)  printf("Family here\n");
      n =  0;
      node = current_sequ->range_start;

      while (node != 0x0)
      {
        strcpy(buff,node->name);

        p = strstr(buff,kn->elname);
        if ( p == buff )
        {
          break;
        }

        n++;
        node = node->next;

        if ( node == current_sequ->range_start )
        {
          fatal_error("readstartvalues: Can not find element: ",kn->elname);
          return 1;
        }
      }

      if (v->kn >=0)
      {
        ncomp = v->kn;
        w_ptc_getnfieldcomp_(&n,&ncomp, &nvalue);
        v->currentvalue = nvalue;
      }
      else
      {
        ncomp = v->ks;
        w_ptc_getsfieldcomp_(&n,&ncomp, &nvalue);
        v->currentvalue = nvalue;
      }
      if (debuglevel)  printf("Got first element %s of family %s, field is %f\n",kn->elname,buff,v->currentvalue);

      n++;
      node = node->next;

      while (node != 0x0)
      {
        strcpy(buff,node->name);

        p = strstr(buff,kn->elname);
        if ( p == buff )
        {
          if (debuglevel)  printf("Got another element %s of the family %s\n",node->name,kn->elname);

          if (v->kn >=0)
          {
            ncomp = v->kn;
            w_ptc_getnfieldcomp_(&n,&ncomp, &nvalue);
          }
          else
          {
            ncomp = v->ks;
            w_ptc_getsfieldcomp_(&n,&ncomp, &nvalue);
          }

          if (v->currentvalue != nvalue)
          {
            warningnew("matchptcknobs",
                       "Element %s has incoherent field %f strngth with its family %f.\n",
                       node->name, nvalue,v->currentvalue);
          }
        }

        n++;
        node = node->next;

        if ( node == current_sequ->range_start )
        {
          break;
        }
      }
    }
    else
    {
      n =  0;
      node = current_sequ->range_start;

      while (node != 0x0)
      {
        strcpy(buff,node->name);
        p = strstr(buff,":");
        if (p) *p = 0;

        if ( strcmp(buff,kn->elname) == 0 )
        {
          break;
        }

        n++;

        node = node->next;

        if ( node == current_sequ->range_start )
        {
          fatal_error("readstartvalues: Can not find element: ",kn->elname);
          return 1;
        }
      }

      if (v->kn >=0)
      {
        ncomp = v->kn;
        w_ptc_getnfieldcomp_(&n,&ncomp, &nvalue);
        v->currentvalue = nvalue;
      }
      else
      {
        ncomp = v->ks;
        w_ptc_getsfieldcomp_(&n,&ncomp, &nvalue);
        v->currentvalue = nvalue;
      }

      if (debuglevel)  printf("Got %f for %s\n",v->currentvalue,kn->elname);

    }
  }

  return 0;
}

static int
madx_mpk_scalelimits(int nv)
{
  /*User specifies uppper and lower limits of the matched k-values
    here they are converted to the PTC units */

/*
  double l;
  struct element* el = 0x0;
  char*  p;
  struct node*                node = 0x0;
*/
  int    fieldorder;/*PTC nomenclature, 1 dipol, 2 quad ...*/
  float  fact = 0.0;

  struct madx_mpk_variable*     v = 0x0;
//  struct madx_mpk_knob*         kn = 0x0; // not used

  if ( (nv < 0) || (nv >= MAX_KNOBS) )
  {
    mad_error("madx_mpk_scalelimits","Passed variable out of range");
    return 1;
  }


  v = &madx_mpk_variables[nv];
  // kn = &(madx_mpk_knobs[v->knobidx]); // not used

  fieldorder = 1 + (v->kn >= 0)?v->kn:v->ks;/*PTC nomenclature, 1 dipol, 2 quad ...*/

  fact = factorial(fieldorder);

  if ( ( v->kn == 0) || (v->ks == 0) )
  {
    printf("madx_mpk_scalelimits: Dipol limits don't need to be scaled\n");
    return 1;
  }

/* If lenght of the magnet will be needed, thin elements ???
 *   if(kn->exactnamematch)
 *    {
 *      el = find_element(kn->elname, element_list);
 *    }
 *   else
 *    {
 *      node = current_sequ->range_start;
 *      while (node != 0x0)
 *       {
 *         p = strstr(node->name,kn->elname);
 *         if ( p == node->name )
 *          {
 *            el = node->p_elem;
 *            break;
 *          }
 *
 *         node = node->next;
 *         if ( node == current_sequ->range_start )
 *          {
 *            fatal_error("madx_mpk_scalelimits: Can not find element starting with: ",kn->elname);
 *            return 1;
 *          }
 *       }
 *
 *
 *    }
 *   if (el == 0x0)
 *    {
 *      mad_error("madx_mpk_scalelimits","Can not find element named %s in the current sequence",madx_mpk_knobs[v->knobidx].elname);
 *      return 1;
 *    }
 *
 *   if (debuglevel)  printf("Element %s\n",el->name);
 *
 *   l = el_par_value("l",el);
 *   if (l <= 0.0)
 *    {
 *      printf("Element %s has zero lenght\n",el->name);
 *      l = 1.0; zero lentgh elements has field defined such, compatible with madx_ptc_module.f90:SUMM_MULTIPOLES_AND_ERRORS
 *
 */

  v->upper = v->upper/fact;
  v->lower = v->lower/fact;

  if (debuglevel)  printf("Set limits to field order (PTC) %d, fact=%f: %f %f\n",
                          fieldorder, fact, v->lower, v->upper );

  return 0;

}

static void
madx_mpk_addfieldcomp(struct madx_mpk_knob* knob, int kn, int ks)
{
/*
  adds a new field component to a knob
*/
  if (knob == 0x0) {
    warning("madx_mpk_addfieldcomp","knob parameter is null");
    return;
  }

  if (kn >= 0) {
    knob->nKN++;
    knob->KN = myrealloc("madx_mpk_addfieldcomp", knob->KN, knob->nKN * sizeof *knob->KN);
    knob->KN[knob->nKN - 1] = kn;
  }

  if (ks >= 0) {
    knob->nKS++;
    knob->KS = myrealloc("madx_mpk_addfieldcomp", knob->KS, knob->nKS * sizeof *knob->KS);
    knob->KS[knob->nKS - 1] = ks;
  }
}

static void
makestdmatchfile(char* fname, char* matchactioncommand)
{
  FILE* f = 0x0;
  struct madx_mpk_variable*      v = 0x0;
  int  i;
  unsigned int anumber = abs(time(0)*rand());
  double lower, upper;

  if (debuglevel)  printf("I am in makestdmatchfile, file name is >%s<\n", fname);

  while(f == 0x0)
  {/*repeat until generate unique file name*/
    if (fname[0] == 0)
    { /*if string is empty*/
      sprintf(fname,"/tmp/match_ptcknobs_%u.madx",anumber);
    }
    f = fopen(fname,"w");
    if (f == 0x0)
    {
      warningnew("makestdmatchfile","Could not open file %s",fname);
    }
  }

  /*if (debuglevel < 2) fprintf(f,"assign, echo=/tmp/mpk_stdmatch.out;\n");*/

  fprintf(f,"match, use_macro;\n");

  for (i = 0; i<madx_mpk_Nvariables; i++)
  {

    v = &(madx_mpk_variables[i]);

    /*create vary command as string*/
    if (debuglevel)  printf("\ncurrentvalue=%f trustrange=%f lower=%f upper=%f\n",
                            v->currentvalue, v->trustrange , v->lower ,v->upper);

    if ( (v->currentvalue - v->trustrange) < v->lower )
    {
      lower = v->currentvalue - v->lower;
    }
    else
    {
      lower = - v->trustrange;
    }

    if ( (v->currentvalue + v->trustrange) > v->upper )
    {
      upper = v->upper - v->currentvalue;
    }
    else
    {
      upper =  v->trustrange;
    }

    fprintf(f,"   vary, name=%s, step= %g, lower= %g, upper = %g;\n",
            v->name,  v->step, lower,     upper);

  }
  fprintf(f,"   \n");

  fprintf(f,"   m1: macro =  \n");
  fprintf(f,"     {\n");

  for (i = 0; i<madx_mpk_Nvariables; i++)
  {
    if (madx_mpk_variables[i].IsIniCond )
    {
      fprintf(f,"      ptc_setknobvalue ,element=%s, value=%s, refreshtables=false;\n",
              madx_mpk_knobs[ madx_mpk_variables[i].knobidx ].initial,
              madx_mpk_variables[i].name);
    }
    else
    {
      fprintf(f,"      ptc_setknobvalue ,element=%s, kn=%d ,ks=%d, value=%s, refreshtables=false;\n",
              madx_mpk_knobs[ madx_mpk_variables[i].knobidx ].elname,
              madx_mpk_variables[i].kn,
              madx_mpk_variables[i].ks,
              madx_mpk_variables[i].name);
    }
    if (debuglevel)  fprintf(f,"      value , %s, %s;\n", madx_mpk_variables[i].name, madx_mpk_variables[i].namecv );
  }

  fprintf(f,"      ptc_refreshpartables;\n");

  fprintf(f,"     };\n");


  for (i = 0; i<madx_mpk_Nconstraints; i++)
  {
    fprintf(f,"   %s;\n",madx_mpk_constraints[i]);
  }


  fprintf(f,"  %s\n",matchactioncommand);

  fprintf(f,"endmatch;\n");

  /*if (debuglevel < 2) fprintf(f,"assign, echo=terminal;\n");*/

  fclose(f);

}

static int
run_ptccalculation(int setknobs, char* readstartval)
{
  int i,n;
  char buff[500];
  char comd[600];
  char* iniparname;

  char **toks=madx_mpk_comm_calculate->tok_list->p;
  int ntoks = madx_mpk_comm_calculate->tok_list->curr;
  int start;

  struct madx_mpk_variable*      v = 0x0;
  struct madx_mpk_knob*         kn = 0x0;
  struct node*                node = 0x0;
  char* p;


  this_cmd = madx_mpk_comm_createuniverse;
  current_command =  madx_mpk_comm_createuniverse->clone;
  process();

  this_cmd = madx_mpk_comm_createlayout;
  current_command =  madx_mpk_comm_createlayout->clone;
  process();

  if (madx_mpk_comm_setswitch)
  {
    this_cmd = madx_mpk_comm_createlayout;
    current_command =  madx_mpk_comm_createlayout->clone;
    process();
  }

  if (*readstartval == 0)
  {
    for(i=0;i<madx_mpk_Nvariables;i++)
    {

      v = &(madx_mpk_variables[i]);
      kn = &(madx_mpk_knobs[v->knobidx]);
      set_variable_(v->namecv, &(v->currentvalue));

      if (v->IsIniCond)
      { /*Set the initial parameter in ptc_twiss*/

        iniparname = kn->initial;
        mapptctomad(iniparname,buff);

        for(start=0; start<ntoks; start++)
        {
          if (strcmp(toks[start],buff)==0)
          {
            break;
          }
        }

        set_command_par_value( buff, madx_mpk_comm_calculate->clone ,v->currentvalue);

        if (debuglevel)   printf("Setting Initial %s to CV %f, now it is %f\n",
                                 buff,v->currentvalue,
                                 command_par_value(buff, madx_mpk_comm_calculate->clone));

      }
      else
      {

        if(kn->exactnamematch == 0)
        {
          n =  0;
          node = current_sequ->range_start;

          while (node != 0x0)
          {
            strcpy(buff,node->name);

            p = strstr(buff,kn->elname);
            if ( p == buff )
            {
              p = strstr(buff,":");
              if (p) *p = 0;


              sprintf(comd,"ptc_setfieldcomp, element=%s, kn=%d, ks=%d, value=%s;",
                      buff, v->kn,v->ks,v->namecv);
              if (debuglevel)  printf("%s\n",comd);
              pro_input_(comd);

            }

            n++;
            node = node->next;

            if ( node == current_sequ->range_start )
            {
              break;
            }
          }

        }
        else
        {
          sprintf(comd,"ptc_setfieldcomp, element=%s, kn=%d, ks=%d, value=%s;",
                  kn->elname, v->kn,v->ks,v->namecv);
          if (debuglevel)  printf("%s\n",comd);
          pro_input_(comd);
        }
      }
    }
  }

  if(setknobs)
  {

    for(i=0;i<madx_mpk_Nknobs;i++)
    {
/*        printf("Setting knob %d: \n%s\n",i,madx_mpk_setknobs[i]);*/
      pro_input_(madx_mpk_setknobs[i]);
    }
  }
  else
  {
    if (debuglevel)  printf("Knob Setting Is not requested this time.\n");
  }

/*  pro_input_(twisscommand);*/

  if (debuglevel)  printf("Running ptc_twiss or ptc_normal.\n");

  this_cmd = madx_mpk_comm_calculate;
  current_command =  madx_mpk_comm_calculate->clone;
  current_twiss = current_command;

  pro_ptc_twiss();
/*  process();*/

  if ( *readstartval )
  {
    readstartvalues();
    *readstartval = 0;
  }

  return 0;
}

static double
match2_summary(void)
{
  int i,j;

  double penalty;

  printf("\n");
  printf("MATCH SUMMARY\n\n");
/*  fprintf(prt_file, "Macro Constraint            Value                     Penalty\n");*/
  printf("--------------------------------------------------------------------\n");
  penalty=0;
  for(i=0;i<MAX_MATCH_MACRO;i++)
  {
    if(match2_macro_name[i]==NULL) break;
    printf("macro: %-20s\n",match2_macro_name[i]);
    for(j=0;j<MAX_MATCH_CONS;j++)
    {
      if (match2_cons_name[i][j]==NULL) break;
      printf("  constraint: %-40s\n",match2_cons_name[i][j]);
      printf("  values:     %+12.5e%c%+12.5e\n",
             match2_cons_value_lhs[i][j],
             match2_cons_sign[i][j],
             match2_cons_value_rhs[i][j]);
      printf("  weight:     %+12.5e\n", match2_cons_weight[i][j]);
      printf("  penalty:    %+12.5e\n\n",match2_cons_value[i][j]);
      penalty+=pow(match2_cons_value[i][j],2);
    }
  }

  return penalty;
}

static int
madx_mpk_init(void)
{
  /*performs initialization,
    prepares commands for exacution out of gethered user information*/

  /*Create ptc_command: just copy paste the valid parameters - no error checking*/

  char                           vary[COMM_LENGTH];
  int i,k;

  for(i=0; i< madx_mpk_Nknobs; i++)
  {
    sprintf(vary,"ptc_knob, ");

    if (madx_mpk_knobs[i].elname)
    {
      sprintf(&(vary[strlen(vary)]),"element=%s, ",madx_mpk_knobs[i].elname);

      if (madx_mpk_knobs[i].nKN > 0)
      {

        sprintf(&(vary[strlen(vary)]),"kn=");

        for (k=0; k < madx_mpk_knobs[i].nKN; k++)
        {
          sprintf(&(vary[strlen(vary)]),"%d,",madx_mpk_knobs[i].KN[k]);
        }
      }

      if (madx_mpk_knobs[i].nKS > 0)
      {

        sprintf(&(vary[strlen(vary)]),"ks=");

        for (k=0; k < madx_mpk_knobs[i].nKS; k++)
        {
          sprintf(&(vary[strlen(vary)]),"%d,",madx_mpk_knobs[i].KS[k]);
        }
      }

      /*this is the last one anyway if ename is specified ; at the end*/
      if (madx_mpk_knobs[i].exactnamematch)
      {
        sprintf(&(vary[strlen(vary)]),"exactmatch=%s; ","true");
      }
      else
      {
        sprintf(&(vary[strlen(vary)]),"exactmatch=%s; ","false");
      }

    }

    if (madx_mpk_knobs[i].initial)
    {
      sprintf(&(vary[strlen(vary)]),"initial=%s; ",madx_mpk_knobs[i].initial);
    }

    int len = strlen(vary)+1;
    madx_mpk_setknobs[i] = mymalloc_atomic("madx_mpk_addvariable", len * sizeof *madx_mpk_setknobs[0]);
    strcpy(madx_mpk_setknobs[i],vary);

    if (debuglevel)  printf("madx_mpk_setknobs[%d]= %s\n",i,madx_mpk_setknobs[i]);
  }

  return 0;
}

// public interface
 
void
madx_mpk_run(struct in_cmd* cmd)
{
/*the main routine of the module, called after at matching action (migrad, lmdif. etc...) */
  const char *rout_name = "madx_mpk_run";
  int i;
  double  tolerance;
  int     calls = 0, maxcalls;
  double  ktar, penalty, penalty2, oldtar = 0.0, tar;
  double  var;
/*  double* matchedvalues;*/
  double* function_vector1 = 0x0;
  double* function_vector2 = 0x0;
  char    ptcend[20];
  char    callmatchfile[400];
  char    matchfilename[300];
  char    matchactioncommand[400];
  char    maxNCallsExceeded = 0;
  char    knobsatmatchedpoint = 0; /*flag indicating that matched values of knobs are with precision close to zero (expansion if valid only close to )*/
  char    matched2 = 0;
  char    matched = 0;
  char    readstartval = 1;
//  int     retval; // not used
  int     fieldorder;
  FILE*   fdbg = fopen("currentvalues.txt","a+");

  matchfilename[0] = 0; /*sign for makestdmatchfile to generate a new file name*/

  tolerance = command_par_value("tolerance",cmd->clone);
  if (tolerance < 0)
  {
    warningnew("matchknobs.c: madx_mpk_run","Tolerance is less then 0. Command ignored.");
    return;
  }
  maxcalls = command_par_value("calls",cmd->clone);

  if(madx_mpk_comm_createuniverse == 0x0)
  {
    warningnew("matchknobs.c: madx_mpk_run","ptc_createuniverse command is missing.");
    return;
  }

  if(madx_mpk_comm_createlayout == 0x0)
  {
    warningnew("matchknobs.c: madx_mpk_run","ptc_createlayout is missing.");
    return;
  }

  if(madx_mpk_comm_calculate == 0x0)
  {
    warningnew("matchknobs.c: madx_mpk_run","Neither ptc_twiss nor ptc_normal seen since \"match, use_ptcknob;\"");
    return;
  }

  if (madx_mpk_Nknobs == 0)
  {
    warningnew("matchknobs.c: madx_mpk_run","no knobs seen yet.");
    return;
  }

  if (madx_mpk_Nvariables == 0)
  {
    warningnew("matchknobs.c: madx_mpk_run","no variables seen yet.");
    return;
  }


  sprintf(matchactioncommand,"%s,",cmd->tok_list->p[0]); /*there is no comma after the command name*/
  i = 1;
  while(cmd->tok_list->p[i] != 0x0)
  {
    sprintf(&(matchactioncommand[strlen(matchactioncommand)])," %s",cmd->tok_list->p[i]);
    if (debuglevel)  printf("%s",cmd->tok_list->p[i]);
    i++;
  }

  sprintf(&(matchactioncommand[strlen(matchactioncommand)]),";");

  madx_mpk_init(); // retval = not used



  strcpy(ptcend,"ptc_end;");
  /*check if ptc starting commands are fine */
  match_is_on = kMatch_NoMatch;

  match2_keepexpressions = 1;

  fprintf(fdbg,"\n###############################################\n");
  for (i=0;i<madx_mpk_Nvariables;i++)
  {
    fprintf(fdbg,"%16s   ",madx_mpk_variables[i].name);
  }
  fprintf(fdbg,"\n"); fflush(0);

  do
  {
    calls++;

    maxNCallsExceeded =  ( (maxcalls > 0)&&(calls > maxcalls) );

    printf("###########################################################################################\n");
    printf("\n\n\n Call %d\n",calls);
    printf("###########################################################################################\n");

    printf("\n CURRENT VALUES \n");
    for (i=0;i<madx_mpk_Nvariables;i++)
    {
      printf("%16s   ",madx_mpk_variables[i].name);
    }
    printf("\n");
    for (i=0;i<madx_mpk_Nvariables;i++)
    {
      printf("%16f   ",madx_mpk_variables[i].currentvalue);
    }
    printf("\n");

    fflush(0);

    if (debuglevel)  printf("\n\nCalling TWISS with KNOBS\n");
    run_ptccalculation(1,&readstartval);

    if (geterrorflag())
    {
      mad_error("Matching With Knobs","PTC calculation ended with an error. Check your setting and matching limits.");
      pro_input_(ptcend);
      goto cleaning;
    }


    /*set all variables to 0*/
    var = 0.0;
    for (i=0;i<madx_mpk_Nvariables;i++)
    {
      set_variable_(madx_mpk_variables[i].name, &var);
    }

    if (debuglevel)  printf("\n\nDoing Match on KNOBS\n");
    makestdmatchfile(matchfilename,matchactioncommand);
    if (debuglevel)  printf("matchfilename is %s\n",matchfilename);

    sprintf(callmatchfile,"call, file=\"%s\" ;",matchfilename);
    pro_input_(callmatchfile);

    if (function_vector1 == 0x0) {
      function_vector1 = mymalloc_atomic("madx_mpk_run", (total_const+1) * sizeof *function_vector1);
      function_vector2 = mymalloc_atomic("madx_mpk_run", (total_const+1) * sizeof *function_vector2);
    }

    i = match2_evaluate_exressions(0,0,function_vector1);
    if (debuglevel) printf("match2_evaluate_exressions returned fun_vector 1 of length %d\n",i);


    pro_input_(ptcend);

    /* END OF STD MATCHING*/

    /* CALCULATE KTAR*/

    ktar = 0.0;
    fprintf(fdbg,"\n%d V ",calls);
    printf("MPK STEP %d",calls);
    for (i=0;i<madx_mpk_Nvariables;i++)
    {
      var = get_variable_(madx_mpk_variables[i].name);
      fprintf(fdbg," %20.16f ", var);
      ktar += var*var;
    }
    fprintf(fdbg,"\n");


    /* BIG  KTAR*/
    /*BIG TAR*/
    /*SCALE AND BRING BACK OLD VALUES*/

    tar = get_variable_("tar");

    if ( (oldtar <= 0) || (tar < 10.0*oldtar)) /*correct within order of magnitude - we can not truest the oldtar */
    {                                         /*because it also bears error of approximation*/
      for (i=0;i<madx_mpk_Nvariables;i++)
      {
        madx_mpk_variables[i].oldvalue = madx_mpk_variables[i].currentvalue;

        var = get_variable_(madx_mpk_variables[i].name);
        madx_mpk_variables[i].currentvalue += var;
      }
    }
    else
    {
      fprintf(fdbg,"Narrowing trust ranges Bringing back the previous step values\n");
      for (i=0;i<madx_mpk_Nvariables;i++)
      {
        madx_mpk_variables[i].trustrange /= 2.0;
        madx_mpk_variables[i].step       /= 2.0;

        madx_mpk_variables[i].currentvalue = madx_mpk_variables[i].oldvalue;

      }

    }

    /*SMALL TAR*/
    /*ADD DELTAS*/

    /*CONTINUE*/

    /* SMALL KTAR*/
    /*VERIFY*/


    fprintf(fdbg,"%d CV",calls);
    for (i=0;i<madx_mpk_Nvariables;i++)
    {
      fprintf(fdbg," %20.16f ", madx_mpk_variables[i].currentvalue);
    }
    fprintf(fdbg,"\n");

    fprintf(fdbg," KTAR=%20.16E TAR=%20.16E OLDTAR=%20.16E \n",ktar, tar, oldtar);
    printf(" KTAR=%20.16E TAR=%20.16E OLDTAR=%20.16E \n",ktar, tar, oldtar);

    if (debuglevel)  printf("KTAR=%E \n", ktar );

    fflush(0x0);


    oldtar = tar;

    if (ktar < tolerance)
    {
      knobsatmatchedpoint = 1;
      if (debuglevel)  printf("Matched values of knobs are close to zeroes within the tolerance\n");
    }
    else
    {

      knobsatmatchedpoint = 0;
      if (debuglevel)  printf("Matched values of knobs are NOT close to zeroes within the tolerance\n");
      if (!maxNCallsExceeded) continue;
    }



    if (debuglevel)  printf("\n\nCalling TWISS with NO knobs\n");
    run_ptccalculation(0,&readstartval);

    i = match2_evaluate_exressions(0,0,function_vector2);
    if (debuglevel)  printf("match2_evaluate_exressions returned fun_vector 2 of length %d\n",i);

    pro_input_(ptcend);

    penalty  = 0.0;
    penalty2 = 0.0;

    for (i=0; i< total_const; i++)
    {
      penalty += function_vector2[i]*function_vector2[i];
      var = function_vector2[i] - function_vector1[i];
      /*if (debuglevel>1) */
      if (debuglevel)  printf("Constraint%d: v1=%f v2=%f diff=%f \n",i,function_vector1[i],function_vector2[i],var);
      penalty2 += var*var;
    }

    if (debuglevel > 1) match2_summary();

    fprintf(fdbg,"Verification: penalty=%E  penalty2=%E \n",penalty, penalty2);

    if (debuglevel)  printf("Closenes of the function vectors with and without knobs %e\n",penalty2);

    if (penalty2 < tolerance)
    {
      matched2 = 1;
      if (debuglevel)  printf("Constraints values with and without knobs within the tolerance\n");
    }
    else
    {
      matched2 = 0;
      if (debuglevel)  printf("Constraints values with and without knobs NOT within the tolerance\n");
    }

    matched = matched2 && knobsatmatchedpoint;


    if (debuglevel)  printf("matched=%d, maxNCallsExceeded=%d\n",matched,maxNCallsExceeded);

  }
  while( ( !(matched || maxNCallsExceeded ))  &&  (ktar != 0.0) );
  /*if ktar is 0.0 we can not make better, because next iteration will be the same*/
  match2_keepexpressions = 0;
  match_is_on = kMatch_PTCknobs;

  fflush(0);
  tar = match2_summary();
  printf("\n\n");
  printf("================================================================================== \n");
  printf("== Matching finished successfully after %3.3d calls                               == \n",calls);
  printf("== Final penalty function TAR = %E                                    == \n", tar);

  printf("== Precision of Taylor expansion at the result point = %E             ==\n",ktar);
  printf("==------------------------------------------------------------------------------== \n");
  printf("==           Variable      Final Value                                          == \n");
  printf("\n");


  for (i=0;i<madx_mpk_Nvariables;i++)
  {
    set_variable_(madx_mpk_variables[i].name,&(madx_mpk_variables[i].currentvalue));
    printf("== %2d %16s    %+12.8E",
           i,madx_mpk_variables[i].name, madx_mpk_variables[i].currentvalue);

    if (madx_mpk_variables[i].IsIniCond)
    {
      printf("                                      == \n");
    }
    else
    {
      fieldorder = 0;
      if (madx_mpk_variables[i].kn >= 0)
      {
        fieldorder = factorial(madx_mpk_variables[i].kn);
      }
      if (madx_mpk_variables[i].ks >= 0)
      {
        fieldorder = factorial(madx_mpk_variables[i].ks);
      }

      if (fieldorder >= 0)
      {

        printf("  PTC units -> MADX: %+12.8E  == \n",
               fieldorder*madx_mpk_variables[i].currentvalue);

      }
      else
      {
        printf("                                    == \n");
      }
    }
    printf("\n");
  }

  printf("==========================================================================\n");

  cleaning:
  if (function_vector1) myfree(rout_name,function_vector1);
  if (function_vector2) myfree(rout_name,function_vector2);

  function_vector1 = 0x0;
  function_vector2 = 0x0;

  fclose(fdbg);
/*  remove(matchfilename);    */
}

void
madx_mpk_prepare(void)
{
/* prepares the internal variables*/
  int i; /**/

  madx_mpk_Nconstraints = 0;
  madx_mpk_Nknobs = 0;
  madx_mpk_Nvariables = 0;

  madx_mpk_comm_createuniverse = 0x0;
  madx_mpk_comm_createlayout = 0x0;
  madx_mpk_comm_setswitch = 0x0;
  madx_mpk_comm_calculate = 0x0;

  for (i = 0; i< MAX_KNOBS; i++)
  {

    madx_mpk_knobs[i].elname  = 0x0;
    madx_mpk_knobs[i].initial = 0x0;
    madx_mpk_knobs[i].KN      = 0x0;
    madx_mpk_knobs[i].nKN     = 0;
    madx_mpk_knobs[i].KS      = 0x0;
    madx_mpk_knobs[i].nKS     = 0;
    madx_mpk_knobs[i].exactnamematch = 1;

    madx_mpk_variables[i].name          = 0x0;
    madx_mpk_variables[i].namecv        = 0x0;
    madx_mpk_variables[i].upper         = 0.0;
    madx_mpk_variables[i].lower         = 0.0;
    madx_mpk_variables[i].trustrange    = 0.0;
    madx_mpk_variables[i].step          = 0.0;
    madx_mpk_variables[i].knobidx       = -1;
    madx_mpk_variables[i].currentvalue  = 0.0;

    madx_mpk_variables[i].kn            = -1;
    madx_mpk_variables[i].ks            = -1;
    madx_mpk_variables[i].IsIniCond     = -1;
  }
}

void
madx_mpk_addconstraint(const char* constr)
{
  char* buff;
  int l;
  if (constr == 0x0) return;
  l = strlen(constr);
  if (l<1) return;
  l++;
  buff = mymalloc_atomic("madx_mpk_addconstraint", l * sizeof *buff);
  strcpy(buff,constr);
  madx_mpk_constraints[madx_mpk_Nconstraints++] = buff;
}

void
madx_mpk_addvariable(struct in_cmd* cmd)
{
  char                           vary[COMM_LENGTH];
  int                            slen;
  char*                          string;
  char*                          ename;
  char*                          initialpar;
  int                            exactnamematch;/*0 for families with name starting with ename, 1 for for element having name ename*/
  int                            kn,ks;
/*  int                            int_arr[100];*/
/*  int                            i;*/
  int                            knobidx; /*index of setknob command for this varyknob*/

  struct madx_mpk_variable*      v = 0x0;

  if (cmd == 0x0) return;

  /*get a name for the variable*/
  ename  = command_par_string("element",cmd->clone);
  initialpar = command_par_string("initial",cmd->clone);
  if (( ename == 0x0 ) && ( initialpar == 0x0 ))
  {
    mad_error("matchknobs.c: madx_mpk_addvariable",
          "Neither element nor initial parameter specified. Command ignored!");
    return;
  }

  if ( ename && initialpar)
  {
    mad_error("matchknobs.c: madx_mpk_addvariable",
          "Single command may define only one of two, field component or initial parameter. Command ignored!");
    return;
  }

  kn = command_par_value("kn",cmd->clone);
  ks = command_par_value("ks",cmd->clone);

  if ( ename && (kn >= 0) && (ks >=0) )
  {
    mad_error("matchknobs.c: madx_mpk_addvariable",
          "Single command may define only one field component, not ks and kn together. Command ignored.");
    return;
  }

  exactnamematch = command_par_value("exactmatch",cmd->clone);

  knobidx = findsetknob(ename,exactnamematch,initialpar);

  if (knobidx < 0)
  {
    mad_error("madx_mpk_addvariable","Error occured while adding this variable.");
    return;
  }

  /*so we make a new variable...*/

  v = &(madx_mpk_variables[madx_mpk_Nvariables]);

  v->upper = command_par_value("upper",cmd->clone);
  v->lower = command_par_value("lower",cmd->clone);
  v->trustrange = command_par_value("trustrange",cmd->clone);
  v->step = command_par_value("step",cmd->clone);
  v->knobidx = knobidx;

  if (initialpar)
  {
    int len = strlen(initialpar)+1;
    madx_mpk_knobs[knobidx].initial = mymalloc_atomic("madx_mpk_addvariable", len * sizeof *madx_mpk_knobs[0].initial);
    strcpy(madx_mpk_knobs[knobidx].initial, initialpar);

    len = sprintf(vary,"mpk_%s", initialpar)+1;
    v->name = mymalloc_atomic("madx_mpk_addvariable", len * sizeof *v->name);
    strcpy(v->name,vary);

    len = sprintf(vary,"mpk_%s_0",initialpar)+1;
    v->namecv = mymalloc_atomic("madx_mpk_addvariable", len * sizeof *v->namecv);
    strcpy(v->namecv,vary);
    v->IsIniCond = 1;
    v->kn = -1;
    v->ks = -1;
  }

  if (ename)
  {
    int len = strlen(ename)+1;
    madx_mpk_knobs[knobidx].elname = mymalloc_atomic("madx_mpk_addvariable", len * sizeof *madx_mpk_knobs[0].elname);
    strcpy(madx_mpk_knobs[knobidx].elname, ename);

    if (exactnamematch != 0)
    {
      madx_mpk_knobs[knobidx].exactnamematch = 1;
    }
    else
    {
      madx_mpk_knobs[knobidx].exactnamematch = 0;
    }

    if (kn >=0)
    {
      /*create variable name, use vary as a temp buffer */
      sprintf(vary,"%s_kn%d",ename,kn);
      madx_mpk_addfieldcomp(&(madx_mpk_knobs[knobidx]), kn, -1);
    }

    if (ks >=0)
    {
      /*create variable name, use vary as a temp buffer */
      sprintf(vary,"%s_ks%d",ename,kn);
      madx_mpk_addfieldcomp(&(madx_mpk_knobs[knobidx]), -1, ks);
    }


    slen = strlen(vary);
    string = mymalloc_atomic("madx_mpk_addvariable", (slen+1) * sizeof *string);
    strcpy(string,vary);
    v->name = string;

    vary[slen  ] = '_';
    vary[slen+1] = '0';
    vary[slen+2] =  0;/*end of string*/
    slen  = slen+2;

    string = mymalloc_atomic("madx_mpk_addvariable", (slen+1) * sizeof *string);
    strcpy(string,vary);
    v->namecv = string;

    v->kn = kn;
    v->ks = ks;
    v->IsIniCond = 0;

    madx_mpk_scalelimits(madx_mpk_Nvariables);

  }

  madx_mpk_Nvariables++;

  if (debuglevel)  printf("Added new knobs: now there is\n  knobs: %d\n  variables: %d\n",
                          madx_mpk_Nknobs,madx_mpk_Nvariables);
}

void
madx_mpk_end(void)
{
/*
  remove_from_command_list(this_cmd->clone->name,stored_commands);
*/
  match_is_on = kMatch_NoMatch;
}

void
madx_mpk_setcreateuniverse(struct in_cmd* cmd)
{
  cmd->clone_flag = 1; /* do not delete for the moment*/
  cmd->label = mymalloc_atomic("madx_mpk_setcreateuniverse", 24 * sizeof *cmd->label);
  strcpy(cmd->label,"matchptcknob_ptc_CU");

  madx_mpk_comm_createuniverse = cmd;
}

void
madx_mpk_setcreatelayout(struct in_cmd* cmd)
{
  cmd->clone_flag = 1; /* do not delete for the moment*/
  cmd->label = mymalloc_atomic("madx_mpk_setcreatelayout", 24 * sizeof *cmd->label);
  strcpy(cmd->label,"matchptcknob_ptc_CL");
  madx_mpk_comm_createlayout = cmd;
}

void
madx_mpk_setsetswitch(struct in_cmd* cmd)
{
  cmd->clone_flag = 1; /* do not delete for the moment*/
  cmd->label = mymalloc_atomic("madx_mpk_setsetswitch", 24 * sizeof *cmd->label);
  strcpy(cmd->label,"matchptcknob_ptc_SSW");
  madx_mpk_comm_setswitch = cmd;
}

void
madx_mpk_setcalc(struct in_cmd* cmd)
{
  cmd->clone_flag = 1; /* do not delete for the moment*/
  cmd->label = mymalloc_atomic("madx_mpk_setcalc", 24 * sizeof *cmd->label);
  strcpy(cmd->label,"matchptcknob_ptc_CMD");
  madx_mpk_comm_calculate = cmd;
}

