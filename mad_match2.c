#include "madx.h"

// public variables

int MAX_MATCH_CONS =  0; /*these are set to proper values at the initialization of the match2 module*/
int MAX_MATCH_MACRO = 0; /*zero values means that it is not initialized yet*/

char match2_keepexpressions = 0; /*do not delete expressions at the end matching used by match with PTC knobs*/

char   **match2_macro_name;
char*  **match2_cons_name;
double **match2_cons_value;
double **match2_cons_value_rhs;
double **match2_cons_value_lhs;
double **match2_cons_weight;
char   **match2_cons_sign;
struct expression* **match2_cons_rhs;
struct expression* **match2_cons_lhs;

int match2_cons_curr[3];

/************************************************************************/
/* The functions below belongs to matchc2, but were moved here to make
   mpars compile */

/************************************************************************/

int
match2_augmentnmacros(void)
{
  /*makes place in the working arrays for a new macro
    Piotr Skowronski Mar 2007
  */
  int i,j;
  char fn[]={"match2_augmentnmacros"};
  char**   new_match2_macro_name;
  char* ** new_match2_cons_name;
  double** new_match2_cons_value;
  double** new_match2_cons_value_rhs;
  double** new_match2_cons_value_lhs;
  double** new_match2_cons_weight;
  char**   new_match2_cons_sign;
  struct expression* ** new_match2_cons_rhs;
  struct expression* ** new_match2_cons_lhs;

  if(MAX_MATCH_MACRO == 0)
  {
    error("match2_augmentnconstraints","match with use_maco was not initialized");
    return 1;
  }

  new_match2_macro_name     = (char**)  mycalloc(fn,MAX_MATCH_MACRO+1,sizeof(char*));
  new_match2_cons_name      = (char* **)mycalloc(fn,MAX_MATCH_MACRO+1,sizeof(char**));
  new_match2_cons_value     = (double**)mycalloc(fn,MAX_MATCH_MACRO+1,sizeof(double*));
  new_match2_cons_value_rhs = (double**)mycalloc(fn,MAX_MATCH_MACRO+1,sizeof(double*));
  new_match2_cons_value_lhs = (double**)mycalloc(fn,MAX_MATCH_MACRO+1,sizeof(double*));
  new_match2_cons_weight    = (double**)mycalloc(fn,MAX_MATCH_MACRO+1,sizeof(double*));
  new_match2_cons_sign      = (char**)  mycalloc(fn,MAX_MATCH_MACRO+1,sizeof(char*));
  new_match2_cons_rhs       = (struct expression* **) mycalloc(fn,MAX_MATCH_MACRO+1,sizeof(struct expression**));
  new_match2_cons_lhs       = (struct expression* **) mycalloc(fn,MAX_MATCH_MACRO+1,sizeof(struct expression**));

  /*copy old pointers to arrays*/
  for(i=0;i<MAX_MATCH_MACRO;i++)
  {
    new_match2_cons_name[i]      = match2_cons_name[i];

    new_match2_cons_value[i]     = match2_cons_value[i];

    new_match2_cons_value_rhs[i] = match2_cons_value_rhs[i];
    new_match2_cons_value_lhs[i] = match2_cons_value_lhs[i];
    new_match2_cons_weight[i]    = match2_cons_weight[i];
    new_match2_cons_sign[i]      = match2_cons_sign[i];

    new_match2_cons_rhs[i]       = match2_cons_rhs[i];
    new_match2_cons_lhs[i]       = match2_cons_lhs[i];

    new_match2_macro_name[i]     = match2_macro_name[i];
  }

  /*free the old arrays*/
  myfree(fn,match2_cons_name);
  myfree(fn,match2_cons_value);
  myfree(fn,match2_cons_value_rhs);
  myfree(fn,match2_cons_value_lhs);
  myfree(fn,match2_cons_weight);
  myfree(fn,match2_cons_sign);
  myfree(fn,match2_cons_rhs);
  myfree(fn,match2_cons_lhs);
  myfree(fn,match2_macro_name);


  /*assign freed pointers to the new arrays*/

  match2_cons_name = new_match2_cons_name  ;
  match2_cons_value = new_match2_cons_value  ;
  match2_cons_value_rhs = new_match2_cons_value_rhs  ;
  match2_cons_value_lhs = new_match2_cons_value_lhs  ;
  match2_cons_weight = new_match2_cons_weight  ;
  match2_cons_sign = new_match2_cons_sign  ;
  match2_cons_rhs = new_match2_cons_rhs  ;
  match2_cons_lhs = new_match2_cons_lhs  ;
  match2_macro_name = new_match2_macro_name  ;



  /*make arrays in the new row*/

  match2_cons_name[MAX_MATCH_MACRO]      = (char**)mycalloc(fn,MAX_MATCH_CONS,sizeof(char*));
  match2_cons_value[MAX_MATCH_MACRO]     = (double*)mycalloc(fn,MAX_MATCH_CONS,sizeof(double));

  match2_cons_value_rhs[MAX_MATCH_MACRO] = (double*)mycalloc(fn,MAX_MATCH_CONS,sizeof(double));
  match2_cons_value_lhs[MAX_MATCH_MACRO] = (double*)mycalloc(fn,MAX_MATCH_CONS,sizeof(double));
  match2_cons_weight[MAX_MATCH_MACRO]    = (double*)mycalloc(fn,MAX_MATCH_CONS,sizeof(double));
  match2_cons_sign[MAX_MATCH_MACRO]      = (char*)mycalloc(fn,MAX_MATCH_CONS,sizeof(char));

  match2_cons_rhs[MAX_MATCH_MACRO]       = (struct expression**) mycalloc(fn,MAX_MATCH_CONS,sizeof(struct expression*));
  match2_cons_lhs[MAX_MATCH_MACRO]       = (struct expression**) mycalloc(fn,MAX_MATCH_CONS,sizeof(struct expression*));


  /*initializes arrays in the last row*/

  match2_macro_name[MAX_MATCH_MACRO]=NULL;

  for(j=0;j<MAX_MATCH_CONS;j++)
  {
    match2_cons_name     [MAX_MATCH_MACRO][j]=0x0;

    match2_cons_value    [MAX_MATCH_MACRO][j]=0.0;
    match2_cons_value_lhs[MAX_MATCH_MACRO][j]=0.0;
    match2_cons_value_rhs[MAX_MATCH_MACRO][j]=0.0;
    match2_cons_weight   [MAX_MATCH_MACRO][j]=0.0;
    match2_cons_sign     [MAX_MATCH_MACRO][j]='n';

    match2_cons_rhs      [MAX_MATCH_MACRO][j]=0x0;
    match2_cons_lhs      [MAX_MATCH_MACRO][j]=0x0;
  }


  return ++MAX_MATCH_MACRO;

}
/************************************************************************/

void
match2_delete_arrays(void)
{
  /*clean the stuff;*/
  int i;
  char fn[]={"match2_delete_arrays"};

  if(MAX_MATCH_MACRO <= 0) return;

  for(i=0;i<MAX_MATCH_MACRO;i++)
  {
    if(match2_cons_name[i] == 0x0) break;
    myfree(fn,match2_cons_name     [i]);
    myfree(fn,match2_cons_value    [i]);
    myfree(fn,match2_cons_value_lhs[i]);
    myfree(fn,match2_cons_value_rhs[i]);
    myfree(fn,match2_cons_weight   [i]);
    myfree(fn,match2_cons_sign     [i]);
    myfree(fn,match2_cons_rhs      [i]);
    myfree(fn,match2_cons_lhs      [i]);
  }

  myfree(fn,match2_cons_name);
  myfree(fn,match2_cons_value);
  myfree(fn,match2_cons_value_rhs);
  myfree(fn,match2_cons_value_lhs);
  myfree(fn,match2_cons_weight);
  myfree(fn,match2_cons_sign);
  myfree(fn,match2_cons_rhs);
  myfree(fn,match2_cons_lhs);
  myfree(fn,match2_macro_name);


  match2_cons_name = 0x0;
  match2_cons_value = 0x0;
  match2_cons_value_rhs = 0x0;
  match2_cons_value_lhs = 0x0;
  match2_cons_weight = 0x0;
  match2_cons_sign = 0x0;
  match2_cons_rhs = 0x0;
  match2_cons_lhs = 0x0;
  match2_macro_name = 0x0;

  /*for security so we cannot add more constraints if the module is not initialized*/
  MAX_MATCH_CONS =  0;
  MAX_MATCH_MACRO = 0;

}
/************************************************************************/


void
match2_delete_expressions(void)
{
  char rout_name[] = "match2_delete_expressions";

  int i,j;

  for(i=0;i<MAX_MATCH_MACRO;i++)
  {
    if ( match2_cons_name[i][0] == 0x0) break;
    for(j=0;j<MAX_MATCH_CONS;j++)
    {
      if ( match2_cons_name[i][j] == 0x0) break;
      myfree(rout_name,match2_cons_name[i][j]);
      delete_expression(match2_cons_rhs[i][j]);
      delete_expression(match2_cons_lhs[i][j]);
      match2_cons_rhs[i][j] = 0x0;
      match2_cons_lhs[i][j] = 0x0;
    }
  }

}
