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

// public interface

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

  match2_cons_name = new_match2_cons_name;
  match2_cons_value = new_match2_cons_value;
  match2_cons_value_rhs = new_match2_cons_value_rhs;
  match2_cons_value_lhs = new_match2_cons_value_lhs;
  match2_cons_weight = new_match2_cons_weight;
  match2_cons_sign = new_match2_cons_sign;
  match2_cons_rhs = new_match2_cons_rhs;
  match2_cons_lhs = new_match2_cons_lhs;
  match2_macro_name = new_match2_macro_name;

  /*make arrays in the new row*/

  match2_cons_name[MAX_MATCH_MACRO]      = mycalloc(fn,MAX_MATCH_CONS,sizeof(char*));
  match2_cons_value[MAX_MATCH_MACRO]     = mycalloc(fn,MAX_MATCH_CONS,sizeof(double));

  match2_cons_value_rhs[MAX_MATCH_MACRO] = mycalloc(fn,MAX_MATCH_CONS,sizeof(double));
  match2_cons_value_lhs[MAX_MATCH_MACRO] = mycalloc(fn,MAX_MATCH_CONS,sizeof(double));
  match2_cons_weight[MAX_MATCH_MACRO]    = mycalloc(fn,MAX_MATCH_CONS,sizeof(double));
  match2_cons_sign[MAX_MATCH_MACRO]      = mycalloc(fn,MAX_MATCH_CONS,sizeof(char));

  match2_cons_rhs[MAX_MATCH_MACRO]       = mycalloc(fn,MAX_MATCH_CONS,sizeof(struct expression*));
  match2_cons_lhs[MAX_MATCH_MACRO]       = mycalloc(fn,MAX_MATCH_CONS,sizeof(struct expression*));

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

void
match2_match(struct in_cmd* cmd)
{
  match_is_on = 2;
  total_const=0;

  if (MAX_MATCH_MACRO == 0)
  {
    MAX_MATCH_CONS =  100;
    MAX_MATCH_MACRO = 1;
    match2_alloc_arrays();
  }
  else
  {
    match2_delete_expressions();
  }

  match2_init_arrays();
}

int
match2_augmentnconstraints(void)
{
  /*makes place in the working arrays for a new macro
    Piotr Skowronski Mar 2007
  */
  int i,j;
  char fn[]={"match2_augmentnconstraints"};
  char*   new_match2_macro_name     = 0x0;
  char* * new_match2_cons_name      = 0x0;
  double* new_match2_cons_value     = 0x0;
  double* new_match2_cons_value_rhs = 0x0;
  double* new_match2_cons_value_lhs = 0x0;
  double* new_match2_cons_weight    = 0x0;
  char*   new_match2_cons_sign      = 0x0;
  struct expression* * new_match2_cons_rhs       = 0x0;
  struct expression* * new_match2_cons_lhs       = 0x0;

  if(MAX_MATCH_MACRO == 0)
  {
    error("match2_augmentnconstraints","match with use_maco was not initialized");
    return 1;
  }

  for(i=0;i<MAX_MATCH_MACRO;i++)
  {

    new_match2_macro_name     = mycalloc(fn,MAX_MATCH_CONS*2,sizeof(char));
    new_match2_cons_name      = mycalloc(fn,MAX_MATCH_CONS*2,sizeof(char*));
    new_match2_cons_value     = mycalloc(fn,MAX_MATCH_CONS*2,sizeof(double));
    new_match2_cons_value_rhs = mycalloc(fn,MAX_MATCH_CONS*2,sizeof(double));
    new_match2_cons_value_lhs = mycalloc(fn,MAX_MATCH_CONS*2,sizeof(double));
    new_match2_cons_weight    = mycalloc(fn,MAX_MATCH_CONS*2,sizeof(double));
    new_match2_cons_sign      = mycalloc(fn,MAX_MATCH_CONS*2,sizeof(char));
    new_match2_cons_rhs       = mycalloc(fn,MAX_MATCH_CONS*2,sizeof(struct expression*));
    new_match2_cons_lhs       = mycalloc(fn,MAX_MATCH_CONS*2,sizeof(struct expression*));

    /*copy old content to the new arrays*/
    for(j=0;j<MAX_MATCH_CONS;j++)
    {
      new_match2_cons_name     [j] = match2_cons_name     [i][j];

      new_match2_cons_value    [j] = match2_cons_value    [i][j];
      new_match2_cons_value_lhs[j] = match2_cons_value_lhs[i][j];
      new_match2_cons_value_rhs[j] = match2_cons_value_rhs[i][j];
      new_match2_cons_weight   [j] = match2_cons_weight   [i][j];
      new_match2_cons_sign     [j] = match2_cons_sign     [i][j];

      new_match2_cons_rhs      [j] = match2_cons_rhs      [i][j];
      new_match2_cons_lhs      [j] = match2_cons_lhs      [i][j];
    }

    /*initializes the new parts*/
    for(j=MAX_MATCH_CONS;j<MAX_MATCH_CONS*2;j++)
    {
      new_match2_cons_name     [j] = 0x0;

      new_match2_cons_value    [j] = 0.0;
      new_match2_cons_value_lhs[j] = 0.0;
      new_match2_cons_value_rhs[j] = 0.0;
      new_match2_cons_weight   [j] = 0.0;
      new_match2_cons_sign     [j] = 'n';

      new_match2_cons_rhs      [j] = 0x0;
      new_match2_cons_lhs      [j] = 0x0;
    }

    /*free the old arrays*/
    myfree(fn,match2_cons_name     [i]);
    myfree(fn,match2_cons_value    [i]);
    myfree(fn,match2_cons_value_lhs[i]);
    myfree(fn,match2_cons_value_rhs[i]);
    myfree(fn,match2_cons_weight   [i]);
    myfree(fn,match2_cons_sign     [i]);
    myfree(fn,match2_cons_rhs      [i]);
    myfree(fn,match2_cons_lhs      [i]);

    /*assign freed pointers to the new arrays*/
    match2_cons_name     [i] = new_match2_cons_name ;
    match2_cons_value    [i] = new_match2_cons_value ;
    match2_cons_value_lhs[i] = new_match2_cons_value_lhs ;
    match2_cons_value_rhs[i] = new_match2_cons_value_rhs ;
    match2_cons_weight   [i] = new_match2_cons_weight ;
    match2_cons_sign     [i] = new_match2_cons_sign ;
    match2_cons_rhs      [i] = new_match2_cons_rhs ;
    match2_cons_lhs      [i] = new_match2_cons_lhs ;
  }

  MAX_MATCH_CONS =  MAX_MATCH_CONS*2;

  return MAX_MATCH_CONS;
}

void
match2_end(struct in_cmd* cmd)
{
  int i,j;

  fprintf(prt_file, "\n");
  fprintf(prt_file, "MATCH SUMMARY\n\n");
/*  fprintf(prt_file, "Macro Constraint            Value                     Penalty\n");*/
  fprintf(prt_file, "--------------------------------------------------------------------\n");
  penalty=0;
  for(i=0;i<MAX_MATCH_MACRO;i++)
  {
    if(match2_macro_name[i]==NULL) break;
    fprintf(prt_file,"macro: %-20s\n",match2_macro_name[i]);
    for(j=0;j<MAX_MATCH_CONS;j++)
    {
      if (match2_cons_name[i][j]==NULL) break;
      fprintf(prt_file,"  constraint: %-40s\n",match2_cons_name[i][j]);
      fprintf(prt_file,"  values:     %+12.5e%c%+12.5e\n",
              match2_cons_value_lhs[i][j],
              match2_cons_sign[i][j],
              match2_cons_value_rhs[i][j]);
      fprintf(prt_file,"  weight:     %+12.5e\n", match2_cons_weight[i][j]);
      fprintf(prt_file,"  penalty:    %+12.5e\n\n",match2_cons_value[i][j]);
      penalty+=pow(match2_cons_value[i][j],2);
    }
  }

  fprintf(prt_file, "\n\n");
  fprintf(prt_file, "Final Penalty Function = %16.8e\n\n",penalty);
  match2_print_var(cmd);
  fprintf(prt_file, "END MATCH SUMMARY\n\n");
  current_gweight = delete_command(current_gweight);
  current_weight = delete_command(current_weight);
  vary_vect = delete_double_array(vary_vect);
  vary_dvect = delete_double_array(vary_dvect);
  fun_vect = delete_double_array(fun_vect);
  comm_constraints->curr = 0;
  stored_match_var = delete_command_list(stored_match_var);
  vary_cnt = 0;
  match_is_on = 0;
  current_call_lim = 0;
  current_calls = 0;
  set_option("twiss_print", &keep_tw_print);
  print_match_summary = 0;
  set_option("match_summary", &print_match_summary);


  fprintf(prt_file, "VARIABLE \"TAR\" SET TO %16.8e\n",penalty);
/*  sprintf(assign_cmd,"tar= %16.8e ;",penalty);*/
/*  pro_input(assign_cmd);*/
  set_variable("tar",&penalty);

  if (!match2_keepexpressions)
  {
    match2_delete_expressions();
    match2_delete_arrays();
    total_const = 0;
  }
}

void
match2_alloc_arrays(void)
{
  int i;
  char fn[]={"match2_alloc_arrays"};

  match2_macro_name     = (char**)  mycalloc(fn,MAX_MATCH_MACRO,sizeof(char*));
  match2_cons_name      = (char* **)mycalloc(fn,MAX_MATCH_MACRO,sizeof(char**));
  match2_cons_value     = (double**)mycalloc(fn,MAX_MATCH_MACRO,sizeof(double*));
  match2_cons_value_rhs = (double**)mycalloc(fn,MAX_MATCH_MACRO,sizeof(double*));
  match2_cons_value_lhs = (double**)mycalloc(fn,MAX_MATCH_MACRO,sizeof(double*));
  match2_cons_weight    = (double**)mycalloc(fn,MAX_MATCH_MACRO,sizeof(double*));
  match2_cons_sign      = (char**)  mycalloc(fn,MAX_MATCH_MACRO,sizeof(char*));
  match2_cons_rhs       = (struct expression* **) mycalloc(fn,MAX_MATCH_MACRO,sizeof(struct expression**));
  match2_cons_lhs       = (struct expression* **) mycalloc(fn,MAX_MATCH_MACRO,sizeof(struct expression**));

  for(i=0;i<MAX_MATCH_MACRO;i++)
  {
    match2_cons_name[i]      = (char**)mycalloc(fn,MAX_MATCH_CONS,sizeof(char*));

    match2_cons_value[i]     = (double*)mycalloc(fn,MAX_MATCH_CONS,sizeof(double));

    match2_cons_value_rhs[i] = (double*)mycalloc(fn,MAX_MATCH_CONS,sizeof(double));
    match2_cons_value_lhs[i] = (double*)mycalloc(fn,MAX_MATCH_CONS,sizeof(double));
    match2_cons_weight[i]    = (double*)mycalloc(fn,MAX_MATCH_CONS,sizeof(double));
    match2_cons_sign[i]      = (char*)mycalloc(fn,MAX_MATCH_CONS,sizeof(char));

    match2_cons_rhs[i]       = (struct expression**) mycalloc(fn,MAX_MATCH_CONS,sizeof(struct expression*));
    match2_cons_lhs[i]       = (struct expression**) mycalloc(fn,MAX_MATCH_CONS,sizeof(struct expression*));
  }

}

void
match2_init_arrays(void)
{
  /*clean the stuff;*/
  int i,j;

  for(i=0;i<MAX_MATCH_MACRO;i++)
  {
    match2_macro_name[i]=NULL;

    for(j=0;j<MAX_MATCH_CONS;j++)
    {
      match2_cons_name     [i][j]=0x0;

      match2_cons_value    [i][j]=0.0;
      match2_cons_value_lhs[i][j]=0.0;
      match2_cons_value_rhs[i][j]=0.0;
      match2_cons_weight   [i][j]=0.0;
      match2_cons_sign     [i][j]='n';

      match2_cons_rhs      [i][j]=0x0;
      match2_cons_lhs      [i][j]=0x0;
    }
  }

  for(i=0;i<3;i++) match2_cons_curr[i]=0;
}

void
match2_macro(struct in_cmd* cmd)
{
  int pos;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int i, idx = -1;


  pos = name_list_pos("name", nl);
  if (nl->inform[pos])
  {

    for(i=0; i < MAX_MATCH_MACRO;i++)
    {
      if (match2_macro_name[i]==NULL)
      {
        idx = i;
        break;
      }
    }

    if (idx < 0  )
    {
      printf("Max number of match macros reached. Augmenting.\n");
      match2_augmentnmacros();
      idx = MAX_MATCH_MACRO -1;
    }
/*    printf("%d\n",i);*/
    match2_macro_name[idx]=pl->parameters[pos]->string;
/*
  printf("%d: exec, %s;\n",idx,pl->parameters[pos]->string);
  printf("%s\n", execute);*/
    /*      pro_input(execute);*/
  }
}

void
match2_disasambleconstraint(struct in_cmd* cmd)
{
  /*
    Disassambles regular constraint with range into set of constraints with expr=...
    Piotr Skowronski Mar 2007
  */
  struct node* nodes[2];
  struct sequence* sequ;
  char* name;
  int k,  jj;
  struct command_parameter_list* pl = cmd->clone->par;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter* par = 0x0;
  char s;
  char buff[50];
  double w;
  sequ = current_sequ;
  name = command_par_string("range",cmd->clone);

  if ( strlen(name) > 0 )  /* parameter has been read */
  {

    k = get_ex_range(name, sequ, nodes);
    if (k == 0)
    {
      error("match2_disasambleconstraint","Bad range! Ignoring\n");
      return;
    }
  }
  else
  {
    printf("Range not specified explicitely, using FULL range.\n");
    nodes[0] = sequ->ex_start;
    nodes[1] = sequ->ex_end;
  }

  w = command_par_value("weight",cmd->clone);
  for (jj = 0; jj < pl->curr; jj++)
  {
    if (nl->inform[jj] && pl->parameters[jj]->type == 4)
    {
      par = pl->parameters[jj];
/*
  printf("Got constraint type %d name %s\n",par->c_type, par->name);
  printf("   min_expr: %#x  c_min=%f \n", par->min_expr, par->c_min);
  printf("   max_expr: %#x  c_max=%f \n", par->max_expr, par->c_max);
  printf("       expr: %#x  double_value %f \n\n", par->expr, par->double_value);
*/

      switch(par->c_type)
      {
        case 1: /* minimum */
        case 3: /* both */
          s = '>';
          if (par->min_expr == NULL)
          {
            sprintf(buff,"%f",par->c_min);
            match2_setconstrinrange(nodes,w, par->name,s,buff);
          }
          else
          {
            match2_setconstrinrange(nodes,w, par->name,s,par->min_expr->string);
          }
          if (par->c_type == 1) break;
        case 2: /* maximum */
          s = '<';
          if (par->max_expr == NULL)
          {
            sprintf(buff,"%f",par->c_max);
            match2_setconstrinrange(nodes,w, par->name,s,buff);
          }
          else
          {
            match2_setconstrinrange(nodes,w, par->name,s,par->max_expr->string);
          }
          break;
        case 4: /* value */
          s = '=';
          if (par->expr == NULL)
          {
            sprintf(buff,"%f",par->double_value);
            match2_setconstrinrange(nodes,w, par->name,s,buff);
          }
          else
          {
            match2_setconstrinrange(nodes,w, par->name,s,par->expr->string);
          }
      }
    }
  }
}

void
match2_setconstrinrange(struct node** nodes, double w, char* parname, char s, char* rexpr)
{
/*
  sets the same consraint for a range of elements defined bu nodes
  Piotr Skowronski Mar 2007
*/
  struct node* c_node;
  char tablecmd[500];
  char buff[500];
  char* p;

  c_node = nodes[0];
  do
  {
    strcpy(buff,c_node->name);

    p = strstr(buff,":");
    if ( p )
    {
      if (p[1] == '0') /*it means that this is a drift automatically added to a sequence*/
      {               /*this guy does not work with table command so we do not care about him*/
        c_node = c_node->next;
        continue;
      }
      p[0] = '[';
      p[2] = ']';
      p[3] =  0;
    }

    sprintf(tablecmd,"constraint, weight=%f, expr=table(twiss,%s,%s)%c%s ;",
            w,  buff, parname, s, rexpr);

    pro_input(tablecmd);

    if (nodes[1] == c_node) break; /*only one element in the range*/

    c_node = c_node->next;
  } while ( c_node && (c_node != nodes[1]) );
}

void
match2_constraint(struct in_cmd* cmd)
{
  int i,j,k,nitem,type;
  int start,end;
  char **toks=cmd->tok_list->p;
  int n = cmd->tok_list->curr;
  struct expression* expr = NULL;
  char* cname;
  char s;


  int exprfound = 0;

  i=0;j=0;s='n';


  for(i=0;i<MAX_MATCH_MACRO;i++)
  {
    if (match2_macro_name[i]==NULL)
    {
      break;
    }
  }

  i--;

  for(j=0;j<MAX_MATCH_CONS;j++)
  {
    if (match2_cons_lhs[i][j]==NULL)
    {
      break;
    }
  }


  if (j >= MAX_MATCH_CONS)
  {
    j=MAX_MATCH_CONS;
    printf("Max number of constraints %d reached. Increasing tables. Macro %d \n",MAX_MATCH_CONS, i);
    match2_augmentnconstraints();
  }

  exprfound = 0;
  for(start=0; start<n; start++)
  {
    if (strcmp(toks[start],"expr")==0)
    {
      exprfound = 1;
      break;
    }
  }

  /*the ckeck if "expr" is present should be here? /skowron/*/
  if (exprfound == 0)
  {
    match2_disasambleconstraint(cmd);
    return;
  }

  start=start+2;
/*  start=3; |+constraint [0], expr [1] = [2] start [3]+|*/
  /*  printf("%s\n",toks[start]);*/
  for (k = start; k < n; k++) {
    s=*toks[k];
    if (s == '<' || s == '>' || s == '=')
    {
      break;
    }
  }

  if (k>=n)
  {
    warning("match2_constraint: expr not present in this constraint","ignoring");
    return;
  }

  /*  printf("%d\n",k);*/
  if ((type = loc_expr(toks, k, start, &end)) > 0)
  {
    nitem = end + 1 - start;

    expr=make_expression(nitem,&toks[start]);


    /*      printf("%d %d\n",i,j);*/
    match2_cons_lhs[i][j]=expr;
    match2_cons_sign[i][j]=s;
    /*      comm->par->parameters[pos]->type=4;*/
    /*      comm->par->parameters[pos]->expr=expr;*/
  }
  else
  {
    warning("match2_constraint: no valid expression in constraint","ignoring");
  }


  start=end+2;
  for (k = start; k < n; k++)
  {
    if (*toks[k] == ';')
    {
      break;
    }
  }

  if ((type = loc_expr(toks, k, start, &end)) > 0)
  {
    nitem = end + 1 - start;
    expr=make_expression(nitem,&toks[start]);
    match2_cons_rhs[i][j]=expr;
  }
  else
  {
    match2_cons_lhs[i][j]=NULL;
    match2_cons_rhs[i][j]=NULL;
    warning("no valid expression in constraint","ignoring");
    return;
  }


  match2_cons_weight[i][j]=command_par_value("weight",cmd->clone);
  for(start=0; start<n; start++)
  {
    if (strcmp(toks[start],"expr")==0)
    {
      break;
    }
  }

  start=start+2;
  nitem = end-start+1;
  if ( (cname=command_par_string("name",cmd->clone) ) == NULL) {
/*    printf("not-given-name\n");*/
    cname=spec_join(&toks[start], nitem);
  } else {
/*    printf("given-name %s\n",cname);*/
  };
  match2_cons_name[i][j]=(char*) mymalloc("match2_constraint",strlen(cname)+1);
/*  strcpy(match2_cons_name[i][j],cname);*/
  n=0;
  for(k=0;k<strlen(cname);k++){
    if(cname[k]!=' ') {
      match2_cons_name[i][j][n]=cname[k];
      n++;
    }
  }
  match2_cons_name[i][j][n]='\0';
/*  strcpy(match2_cons_name[i][j],cname);*/
  total_const++;
/*
  printf("%d %d, ntot=%d: %s\n",i,j,total_const,match2_cons_name[i][j]);
  printf("%e %c %e\n",expression_value(match2_cons_lhs[i][j],type) , match2_cons_sign[i][j],
  expression_value(match2_cons_rhs[i][j],type) );
*/
}

int
match2_evaluate_exressions(int i, int k, double* fun_vec)
{
  int j;
  double rhs,lhs,r;/* RDM fork */
  char s;
  for(j=0; j < MAX_MATCH_CONS ;j++) {

    if (match2_cons_name[i][j]==NULL) break;

    rhs=expression_value(match2_cons_rhs[i][j],2);
    lhs=expression_value(match2_cons_lhs[i][j],2);
    s =match2_cons_sign[i][j];
    r=lhs - rhs;
    fun_vec[k]=match2_cons_weight[i][j]*r;
    if (s == '>' && r > 0) fun_vec[k]=0;
    else if (s == '<'  && r < 0) fun_vec[k]=0;
    match2_cons_value[i][j]=fun_vec[k];
    match2_cons_value_rhs[i][j]=rhs;
    match2_cons_value_lhs[i][j]=lhs;
    k++;
  }

  return k;
}

int
match2_print_var(struct in_cmd* cmd)
{
  int n,l;
  char *varname,*knobfilename,*knobname;
  double ivalue,fvalue;
  FILE *knobfile=NULL;
  knobfilename=command_par_string("knobfile",cmd->clone);
  if (knobfilename){
      knobfile=fopen(knobfilename,"w");
  };
  n=stored_match_var->curr;
  fprintf(prt_file, "\n\n");
  fprintf(prt_file, "%-24s %-12s %-12s %-12s %-12s\n",
                    "Variable", "Final Value",
                    "Initial Value","Lower Limit",
                    "Upper Limit");
  for(l=0;l<80;l++) fprintf(prt_file, "-");
  fprintf(prt_file,"\n");
  for(l=0;l<n;l++){
    varname=command_par_string("name",stored_match_var->commands[l]);
    ivalue=command_par_value("init",stored_match_var->commands[l]);
    fvalue=get_variable(varname);
    knobname=command_par_string("knob",stored_match_var->commands[l]);
    if (knobfilename){
      fprintf(knobfile, "%-12s :=%+15.8e%+15.8e*%s;\n",
                         varname,ivalue,fvalue-ivalue,knobname);
    }
    fprintf(prt_file,"%-24s",varname);
    fprintf(prt_file," %12.5e",fvalue);
    fprintf(prt_file," %12.5e",ivalue);
    fprintf(prt_file," %12.5e",command_par_value("lower",stored_match_var->commands[l]));
    fprintf(prt_file," %12.5e",command_par_value("upper",stored_match_var->commands[l]));
    fprintf(prt_file,"\n");
  }
  fprintf(prt_file,"\n");
  if ( knobfilename ) fclose(knobfile);
  return 0;
}
