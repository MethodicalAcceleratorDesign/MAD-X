#include "matchptcknobs.h"
/*______________________________________________________________
matchptcknobs.c
Piotr Skowronski (CERN) 2006
MAD-X module responsible for matching with help of PTC knobs (parameters)

It implements MAD-X the command:
match, useptcknobs=true;
(...)
endmatch;

It is basically a macro that performs automatically a set of commands
depicted above as (...) within a rather complicated loop (points 3-17 below). 
It enables fast matching of any order parameters with PTC.

The algorithm:

1. Buffer the key commands (ptc_knob, constraint, ptc_setswitch, ptc_twiss or ptc_normal, etc) 
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
#include "stdlib.h"
#include "time.h"
#include "string.h"

#define MAX_CONTRAINS  100
#define MAX_KNOBS  20
#define COMM_LENGTH  500


/************************************/
/*          DATA structure          */
/************************************/

int            madx_mpk_Nconstraints;
char*          madx_mpk_constraints[MAX_CONTRAINS];
int            madx_mpk_Nknobs;
struct in_cmd* madx_mpk_knobs[MAX_KNOBS]; /*this is ptc_setknobs*/
int            madx_mpk_Nvariables;  /*total number if matching variables >= number of knob commands (one knob command may define more then one knob/variable)*/
char*          madx_mpk_variables[MAX_KNOBS]; /*vary command*/
char*          madx_mpk_variable_names[MAX_KNOBS]; /*names of each variable*/
char*          madx_mpk_variable_elnames[MAX_KNOBS]; /*element names for each variable*/
int            madx_mpk_variable_KNs[MAX_KNOBS]; /*kn componment */
int            madx_mpk_variable_KSs[MAX_KNOBS]; /*ks component */

struct in_cmd* madx_mpk_comm_createuniverse;
struct in_cmd* madx_mpk_comm_createlayout;
struct in_cmd* madx_mpk_comm_setswitch;
struct in_cmd* madx_mpk_comm_calculate;/*ptc_twiss or ptc_normal*/

char twisscommand[]="ptc_twiss, table=ptc_twiss, icase=6, no=2, betx=10, alfx=.0,  bety=10., alfy=0, betz=10., alfz=0;";

void makestdmatchfile(char* fname, struct in_cmd* cmd);
int run_ptccalculation(double* currentvalues,int setknobs);

/*******************************************/
/*MAD-X Global Variables used in this file */
/*******************************************/

extern int              debuglevel;
extern int              match_is_on;
extern struct in_cmd*   this_cmd;
extern struct command*  current_command;
extern char             match2_keepexpressions;
extern double           match2_cons_value[10][30];
extern double           match2_cons_value_rhs[10][30];
extern double           match2_cons_value_lhs[10][30];

/************************************/
/*MAD-X Functions used in this file */
/************************************/

extern void             pro_input(char* statement);
extern void             process();
extern struct command*  clone_command(struct command*);
extern double           command_par_value(char*, struct command*);
extern char*            command_par_string(char*, struct command*);
extern void             comm_para_(char*, int*, int*, int*, int*, double*, char*, int*);
extern double           get_variable_(char*);
extern void             set_variable_(char*, double*);
extern void*            mymalloc(char* caller, size_t size);
extern void             myfree(char* rout_name, void* p);
extern void             warningnew(char* t1, char* fmt, ...); 
extern int              match2_evaluate_exressions(int i, double* fun_vec);
extern int              name_list_pos(char*, struct name_list*);
extern void             pro_ptc_knob(struct in_cmd* cmd);
/*_________________________________________________________________________*/
/*_________________________________________________________________________*/
/*_________________________________________________________________________*/
/*******************************/
/* Implementations starts here */
/*******************************/

void madx_mpk_run(struct in_cmd* cmd)
{
  char rout_name[] = "madx_mpk_run";
  int i;
  double  tolerance;
  int     calls = 0, maxcalls;
  double  ktar, penalty, penalty2;
  double  var;
  double* matchedvalues;
  double* function_vector1;
  double* function_vector2;
  char    ptcend[20];
  char    callmatchfile[200];
  char    matchfilename[100];
  char    maxNCallsExceeded = 0;
  char    knobsatmatchedpoint = 0; /*flag indicating that matched values of knobs are with precision close to zero (expansion if valid only close to )*/
  char    matched2 = 0;
  char    matched = 0;
  
  
  tolerance = command_par_value("tolerance",cmd->clone);
  if (tolerance < 0)
   {
     warningnew("matchknobs.c: madx_mpk_run","Tolerance is less then 0. Command ignored.");
     return;
   } 
  maxcalls = (int)command_par_value("calls",cmd->clone);
  
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
  
  
  
  makestdmatchfile(matchfilename,cmd);
  sprintf(callmatchfile,"call, file=\"%s\" ;",matchfilename);

  matchedvalues = (double*)mymalloc("madx_mpk_run",madx_mpk_Nvariables*sizeof(double));
  memset(matchedvalues,0,madx_mpk_Nvariables*sizeof(double));
  
  function_vector1 = (double*)mymalloc("madx_mpk_run",madx_mpk_Nconstraints*sizeof(double));
  function_vector2 = (double*)mymalloc("madx_mpk_run",madx_mpk_Nconstraints*sizeof(double));
  
  
  strcpy(ptcend,"ptc_end;");
  /*check if ptc starting commands are fine */
  match_is_on = kMatch_NoMatch;   

  match2_keepexpressions = 1;
  
  do
   {
     calls++;

     maxNCallsExceeded =  ( (maxcalls > 0)&&(calls > maxcalls) );
      
     printf("\n\n\n\n\n Call %d\n",calls);
     
     printf("\n CURRENT VALUES \n");
     for (i=0;i<madx_mpk_Nvariables;i++)
      {
        printf("%16s   ",madx_mpk_variable_names[i]);
      }
     printf("\n"); 
     for (i=0;i<madx_mpk_Nvariables;i++)
      {
        printf("%16f   ",matchedvalues[i]);
      }
     
    
     printf("\n\nCalling TWISS with KNOBS\n");
     run_ptccalculation(matchedvalues,1);
     
     /*set all variables to 0*/
     var = 0.0;
     for (i=0;i<madx_mpk_Nvariables;i++)
      {
        set_variable_(madx_mpk_variable_names[i],&var);
      }
     
     printf("\n\nDoing Match on KNOBS\n");
     pro_input(callmatchfile);
     
     printf("\n\nPeeking Function_Vector 1\n");
     i = match2_evaluate_exressions(0,function_vector1);
     printf("match2_evaluate_exressions returned fun_vector of length %d\n",i);

     
     pro_input(ptcend);

     printf("\n\n\n");
     printf("\nFunction_Vector 1:\n");

     ktar = 0.0;
     for (i=0;i<madx_mpk_Nvariables;i++)
      {
    /*      set_variable_(madx_mpk_variable_names[i],&penalty);*/
        var = get_variable_(madx_mpk_variable_names[i]);
        printf("%s=%f  ",madx_mpk_variable_names[i],var);
        matchedvalues[i] += var;
        ktar += var*var;
      }
     printf("\n");
     printf("KTAR=%E \n", ktar ); 
     
     if (ktar < tolerance)
      { 
        knobsatmatchedpoint = 1;
        printf("Matched values of knobs are close to zeroes within the tolerance\n");
      }  
     else
      {
        knobsatmatchedpoint = 0;
        printf("Matched values of knobs are NOT close to zeroes within the tolerance\n");
        if (!maxNCallsExceeded) continue;
      }

     
     printf("\n\nCalling TWISS with NO knobs\n");
     run_ptccalculation(matchedvalues,0);
     
     printf("\n\nPeeking Function_Vector 2\n");
     i = match2_evaluate_exressions(0,function_vector2);

     pro_input(ptcend);
    
     
     
     penalty  = 0.0;
     penalty2 = 0.0;

     printf("\n Function_Vector 2:\n");
     for (i=0; i< madx_mpk_Nconstraints; i++)
      {
        penalty += function_vector2[i]*function_vector2[i];
        var = function_vector2[i] - function_vector1[i];
        if (debuglevel>1) printf("Constraint%d: %f diff=%f \n",i,function_vector2[i],var);
        penalty2 += var*var;
      }
     
     printf("Closenes of the function vectors with and without knobs %e\n",penalty2);
     
     if (penalty2 < tolerance)
      { 
        matched2 = 1;
        printf("Constraints values with and without knobs within the tolerance\n");
      }  
     else
      {
        matched2 = 0;
        printf("Constraints values with and without knobs NOT within the tolerance\n");
      }
    
     matched = matched2 && knobsatmatchedpoint;
   
     
     printf("matched=%d, maxNCallsExceeded=%d\n",matched,maxNCallsExceeded);
   }
  while( !(matched || maxNCallsExceeded) );

  match2_keepexpressions = 0;
  match_is_on = kMatch_PTCknobs;   
  
  fflush(0);
  printf("\n\n");
  printf("========================================================================== \n");
  printf("== Matching finished successfully after %3.3d calls                       == \n",calls);
  
  printf("== Precision of Taylor expansion the result point = %E        ==\n",ktar);
  printf("==----------------------------------------------------------------------== \n");
  printf("==      Variable                   Final Value                          == \n");
  printf("\n");
  for (i=0;i<madx_mpk_Nvariables;i++)
   {
     set_variable_(madx_mpk_variable_names[i],&matchedvalues[i]);
     printf("== %d %16s          %+f                                  == \n",
             i,madx_mpk_variable_names[i],matchedvalues[i]);
     printf("\n");           
   }
  
  printf("==========================================================================\n");

  myfree(rout_name,matchedvalues);
  myfree(rout_name,function_vector1);
  myfree(rout_name,function_vector2);
/*  remove(matchfilename);    */
}
/*_________________________________________________________________________*/

void madx_mpk_prepare()
{
  printf("I am in MATCH PTC WITH KNOBS\n");

  madx_mpk_Nconstraints = 0;
  madx_mpk_Nknobs = 0;
  madx_mpk_Nvariables = 0;
  
  madx_mpk_comm_createuniverse = 0x0;
  madx_mpk_comm_createlayout = 0x0;
  madx_mpk_comm_setswitch = 0x0;
  madx_mpk_comm_calculate = 0x0;
  
  return;
  
}


/*_________________________________________________________________________*/

void madx_mpk_addconstraint(const char* constr)
{
  char* buff;
  int l;
  if (constr == 0x0) return;
  l = strlen(constr);
  if (l<1) return;
  l++;
  buff = (char*)mymalloc("madx_mpk_addconstraint",l*sizeof(char));
  strcpy(buff,constr);
  printf("madx_mpk_addconstraint: got %s\n",buff);
  madx_mpk_constraints[madx_mpk_Nconstraints++] = buff;
  
}
/*_________________________________________________________________________*/

void madx_mpk_addvariable(struct in_cmd* cmd)
{
  struct command_parameter_list* c_parameters;
  struct name_list*              c_parnames;
  char                           vary[COMM_LENGTH];
  int                            pos;
  char*                          string;
  char*                          ename;
  static char*                   verypars[3]={"step","upper","lower"};
  int                            i;
  int                            int_arr[100];
  
  
  if (cmd == 0x0) return;
  
  sprintf(vary,"matchptcknob_variable%d",madx_mpk_Nknobs);
  cmd->clone_flag = 1; /* do not delete for the moment*/
  cmd->label =  (char*)mymalloc("madx_mpk_addvariable",strlen(vary)*sizeof(char));
  strcpy(cmd->label,vary);
  
  madx_mpk_knobs[madx_mpk_Nknobs] = cmd;
  

  c_parameters = cmd->clone->par;
  c_parnames   = cmd->clone->par_names;
  
  pos   = name_list_pos("element", c_parnames);
  if (pos < 0)
  {
    warningnew("matchknobs.c: madx_mpk_addvariable","element parameter does not exist");
    return;
  }

  ename  = c_parameters->parameters[pos]->string;
  if ( ename == 0x0 )
  {
    warningnew("matchknobs.c: madx_mpk_addvariable"," no element name: ignored");
    return;
  }
  
  

  
  /*i is dummy... we know there is no strings or doubles */
  comm_para_("kn", &pos, &i, &i, int_arr, 0x0, 0x0, 0x0);

  printf("pos is %d\n",pos);
  for (i=0;i<pos;i++)
   {
     if (int_arr[i] < 0) break;
     printf("vary, name=%s_kn%d variable no. %d \n",ename,int_arr[i],madx_mpk_Nvariables);
     
     /*create variable name, use vary as a temp buffer */
     sprintf(vary,"%s_kn%d",ename,int_arr[i]);
     string = (char*)mymalloc("madx_mpk_addvariable",strlen(vary)*sizeof(char));  
     strcpy(string,vary);
     madx_mpk_variable_names[madx_mpk_Nvariables] = string;
     
     /*create vary command as string*/
     sprintf(vary,"vary, name=%s, %s= %g , %s= %g , %s= %g ",madx_mpk_variable_names[madx_mpk_Nvariables],
                   verypars[0],command_par_value(verypars[0],cmd->clone),
                   verypars[1],command_par_value(verypars[1],cmd->clone),
                   verypars[2],command_par_value(verypars[2],cmd->clone));

     
     string = (char*)mymalloc("madx_mpk_addvariable",strlen(vary)*sizeof(char));  
     strcpy(string,vary);

     madx_mpk_variables[madx_mpk_Nvariables]      = string;
     
     /*remember element name*/ 
     string = (char*)mymalloc("madx_mpk_addvariable",strlen(ename)*sizeof(char));  
     strcpy(string,ename);
     madx_mpk_variable_elnames[madx_mpk_Nvariables] = string;
     
     /*and field coefficients of this knob*/
     madx_mpk_variable_KNs[madx_mpk_Nvariables]   = int_arr[i];
     madx_mpk_variable_KSs[madx_mpk_Nvariables]   = -1;

     madx_mpk_Nvariables++;
   }

  comm_para_("ks", &pos, &i, &i, int_arr, 0x0, 0x0, 0x0);
  printf("pos is %d\n",pos);
  for (i=0;i<pos;i++)
   {
     if (int_arr[i] < 0) break;
     printf("vary, name=%s_ks%d variable no. %d \n",ename,int_arr[i],madx_mpk_Nvariables);
     
     /*create variable name, use vary as a temp buffer */
     sprintf(vary,"%s_ks%d",ename,int_arr[i]);
     string = (char*)mymalloc("madx_mpk_addvariable",strlen(vary)*sizeof(char));  
     strcpy(string,vary);
     madx_mpk_variable_names[madx_mpk_Nvariables] = string;
     
     /*create vary command as string*/
     sprintf(vary,"vary, name=%s, %s= %g , %s= %g , %s= %g ",madx_mpk_variable_names[madx_mpk_Nvariables],
                   verypars[0],command_par_value(verypars[0],cmd->clone),
                   verypars[1],command_par_value(verypars[1],cmd->clone),
                   verypars[2],command_par_value(verypars[2],cmd->clone));


     string = (char*)mymalloc("madx_mpk_addvariable",strlen(vary)*sizeof(char));  
     strcpy(string,vary);

     madx_mpk_variables[madx_mpk_Nvariables]      = string;
     
     /*remember element name*/ 
     string = (char*)mymalloc("madx_mpk_addvariable",strlen(ename)*sizeof(char));  
     strcpy(string,ename);
     madx_mpk_variable_elnames[madx_mpk_Nvariables] = string;
     
     /*and field coefficients of this knob*/
     madx_mpk_variable_KSs[madx_mpk_Nvariables]   = int_arr[i];
     madx_mpk_variable_KNs[madx_mpk_Nvariables]   = -1;
     
     madx_mpk_Nvariables++;

   }
  
  madx_mpk_Nknobs++;
  
  printf("Added new knobs: now there is\n  knobs: %d\n  variables: %d\n",
          madx_mpk_Nknobs,madx_mpk_Nvariables);
}
/*_________________________________________________________________________*/

void madx_mpk_setcreateuniverse(struct in_cmd* cmd)
{
  cmd->clone_flag = 1; /* do not delete for the moment*/
  cmd->label =  (char*)mymalloc("madx_mpk_setcreateuniverse",24*sizeof(char));
  strcpy(cmd->label,"matchptcknob_ptc_CU");

  madx_mpk_comm_createuniverse = cmd;
}
/*_________________________________________________________________________*/
void madx_mpk_setcreatelayout(struct in_cmd* cmd)
{
  cmd->clone_flag = 1; /* do not delete for the moment*/
  cmd->label =  (char*)mymalloc("madx_mpk_setcreatelayout",24*sizeof(char));
  strcpy(cmd->label,"matchptcknob_ptc_CL");
  madx_mpk_comm_createlayout = cmd;
}
/*_________________________________________________________________________*/
void madx_mpk_setsetswitch(struct in_cmd* cmd)
{
  cmd->clone_flag = 1; /* do not delete for the moment*/
  cmd->label =  (char*)mymalloc("madx_mpk_setsetswitch",24*sizeof(char));
  strcpy(cmd->label,"matchptcknob_ptc_SSW");
  madx_mpk_comm_setswitch = cmd;
}
/*_________________________________________________________________________*/

void madx_mpk_setcalc(struct in_cmd* cmd)
{
  cmd->clone_flag = 1; /* do not delete for the moment*/
  cmd->label =  (char*)mymalloc("madx_mpk_setcalc",24*sizeof(char));
  strcpy(cmd->label,"matchptcknob_ptc_CMD");
  madx_mpk_comm_calculate = cmd;
}
/*_________________________________________________________________________*/

void madx_mpk_end()
{

/*
  remove_from_command_list(this_cmd->clone->name,stored_commands);


*/
  match_is_on = kMatch_NoMatch;

}


/*_________________________________________________________________________*/
/*_________________________________________________________________________*/
/*_________________________________________________________________________*/

void makestdmatchfile(char* fname, struct in_cmd* cmd)
{
  FILE* f = 0x0;
  int  i;
  unsigned int anumber = time(0x0)*rand();
  
  printf("I am in makestdmatchfile\n");
  
  while(f == 0x0)
   {/*repeat until generate unique file name*/
     sprintf(fname,"/tmp/match_ptcknobs_%d.madx",anumber);
     f = fopen(fname,"w");
     if (f == 0x0)
      {
        warningnew("madx_mpk_run","Could not open file %s",fname);
      }
     else
      {
        printf("madx_mpk_run: Knob Matching file in %s\n",fname); 
      } 
   }

  printf("Std Match file name is %s\n",fname);
  
  fprintf(f,"match, use_macro;\n");
  
  for (i = 0; i<madx_mpk_Nvariables; i++)
   {
     fprintf(f,"   %s ; \n",madx_mpk_variables[i]);
     
   }
  fprintf(f,"   \n");

  fprintf(f,"   m1: macro =  \n");
  fprintf(f,"     {\n");

  printf("Std Match file name is %s\n",fname);
  printf("There is %d variables.\n",madx_mpk_Nvariables);
  
  for (i = 0; i<madx_mpk_Nvariables; i++)
   {
     
     fprintf(f,"      ptc_setknobvalue ,element=%s, kn=%d ,ks=%d, value=%s;\n",
                    	madx_mpk_variable_elnames[i], 
		madx_mpk_variable_KNs[i],
		madx_mpk_variable_KSs[i],
                    	madx_mpk_variable_names[i]);
     
     if (debuglevel)
      {
        fprintf(f,"      value , %s;\n", madx_mpk_variable_names[i]);
      }
   }

  
  fprintf(f,"     };\n");


  for (i = 0; i<madx_mpk_Nconstraints; i++)
   {
     fprintf(f,"   %s;\n",madx_mpk_constraints[i]);
   }
   

  i = 0;
  while(cmd->tok_list->p[i] != 0x0)
   {
     fprintf(f," %s",cmd->tok_list->p[i]);
     i++;
   }
  fprintf(f,";\n");

  fprintf(f,"endmatch;\n");

  fclose(f);

}

int run_ptccalculation(double* currentvalues,int setknobs)
{
  int i;
  char buff[500];
  
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
  
  for(i=0;i<madx_mpk_Nvariables;i++)
   {
     sprintf(buff,"ptc_setfieldcomp, element=%s, kn=%d, ks=%d, value=%f;",
                                     madx_mpk_variable_elnames[i], 
	                 madx_mpk_variable_KNs[i],
	                 madx_mpk_variable_KSs[i],
                                     currentvalues[i]);
     printf("%s\n",buff);
     pro_input(buff);
   }

  if(setknobs)
   {
    for(i=0;i<madx_mpk_Nknobs;i++)
     {
        this_cmd = madx_mpk_knobs[i];
        current_command =  madx_mpk_knobs[i]->clone;
        pro_ptc_knob(madx_mpk_knobs[i]);
     }
   }
   
/*  pro_input(twisscommand);*/
  
  this_cmd = madx_mpk_comm_calculate;
  current_command =  madx_mpk_comm_calculate->clone;
  process();
  
  return 0;
} 
