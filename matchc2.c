
void match2_match(struct in_cmd* cmd)
{
  int i,j;
  match_is_on = 2;
  total_const=0;
  printf("Inside USE_MACRO mode\n");
  for(i=0;i<10;i++) match2_macro_name[i]=NULL;
  for(i=0;i<10;i++)  for(j=0;j<30;j++) match2_cons_name[i][j]=NULL;
  for(i=0;i<10;i++)  for(j=0;j<30;j++) match2_cons_rhs[i][j]=NULL;
  for(i=0;i<10;i++)  for(j=0;j<30;j++) match2_cons_lhs[i][j]=NULL;
  for(i=0;i<10;i++)  for(j=0;j<30;j++) match2_cons_value[i][j]=0;
  for(i=0;i<10;i++)  for(j=0;j<30;j++) match2_cons_value_lhs[i][j]=0;
  for(i=0;i<10;i++)  for(j=0;j<30;j++) match2_cons_value_rhs[i][j]=0;
  for(i=0;i<10;i++)  for(j=0;j<30;j++) match2_cons_sign[i][j]='n';
  for(i=0;i<10;i++)  for(j=0;j<30;j++) match2_cons_weight[i][j]=0;
  for(i=0;i<3;i++) match2_cons_curr[i]=0;
  return;
}

void match2_end(struct in_cmd* cmd)
{
  int i,j;


  fprintf(prt_file, "\n");
  fprintf(prt_file, "MATCH SUMMARY\n\n");
/*  fprintf(prt_file, "Macro Constraint            Value                     Penalty\n");*/
  fprintf(prt_file, "--------------------------------------------------------------------\n");
  penalty=0;
  for(i=0;match2_macro_name[i]!=NULL;i++) {
    printf("macro: %-20s\n",match2_macro_name[i]);
    for(j=0;match2_cons_name[i][j]!=NULL;j++) {
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
  
  if (!match2_keepexpressions) 
   {
     match2_delete_expressions();
   }  
   
  fprintf(prt_file, "\n\n");
  fprintf(prt_file, "Final Penalty Function = %16.8e\n\n",penalty);

  fprintf(prt_file, "\n\n");
  fprintf(prt_file, "Variable                   Final Value        Lower Limit        Upper Limit\n");
  fprintf(prt_file, "-------------------------------------------------------------------------------\n");
  print_match_summary=1;
  set_option("match_summary",&print_match_summary);
  mtgetc_(vary_vect->a, vary_dvect->a);
  fprintf(prt_file, "\n");
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
  total_const = 0;
  set_option("twiss_print", &keep_tw_print);
  print_match_summary = 0;
  set_option("match_summary", &print_match_summary);


  fprintf(prt_file, "EVALUATING \"tar= %16.8e;\"\n",penalty);
/*  sprintf(assign_cmd,"tar= %16.8e ;",penalty);*/
/*  pro_input(assign_cmd);*/
  set_variable("tar",&penalty);
  return;
}

void match2_macro(struct in_cmd* cmd)
{
  int pos;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int i;

  pos = name_list_pos("name", nl);
  if (nl->inform[pos]) {
    for(i=0;match2_macro_name[i]!=NULL;i++);
/*    printf("%d\n",i);*/
    match2_macro_name[i]=pl->parameters[pos]->string;
    printf("%d: exec, %s;\n",i,pl->parameters[pos]->string);
/*    printf("%s\n", execute);*/
    /*      pro_input(execute);*/
  }
  return;
}

void match2_constraint(struct in_cmd* cmd)
{
  int i,j,k,nitem,type;
  int start,end;
  char **toks=cmd->tok_list->p;
  int n = cmd->tok_list->curr;
  struct expression* expr = NULL;
  char* cname;
  char s;

  i=0;j=0;s='n';
  for(start=0; start<n; start++) {
    if (strcmp(toks[start],"expr")==0) break;
  }
  start=start+2;
/*  start=3; |+constraint [0], expr [1] = [2] start [3]+|*/
  /*  printf("%s\n",toks[start]);*/
  for (k = start; k < n; k++) {
    s=*toks[k];
    if (s == '<' || s == '>' || s == '=') break;
  }
  if (k==n) {
    warning("no valid expression in constraint","ignoring");
    return;
  }
  /*  printf("%d\n",k);*/
  if ((type = loc_expr(toks, k, start, &end)) > 0) {
    nitem = end + 1 - start;
    expr=make_expression(nitem,&toks[start]);
    for(i=0;match2_macro_name[i]!=NULL;i++);i--;
    for(j=0;match2_cons_lhs[i][j]!=NULL;j++);
    /*      printf("%d %d\n",i,j);*/
    match2_cons_lhs[i][j]=expr;
    match2_cons_sign[i][j]=s;
    /*      comm->par->parameters[pos]->type=4;*/
    /*      comm->par->parameters[pos]->expr=expr;*/
  } else {
    warning("no valid expression in constraint","ignoring");
  }
  start=end+2;
  for (k = start; k < n; k++) if (*toks[k] == ';') break;
  if ((type = loc_expr(toks, k, start, &end)) > 0) {
    nitem = end + 1 - start;
    expr=make_expression(nitem,&toks[start]);
    match2_cons_rhs[i][j]=expr;
  } else {
    match2_cons_lhs[i][j]=NULL;
    match2_cons_rhs[i][j]=NULL;
    warning("no valid expression in constraint","ignoring");
    return;
  }
  match2_cons_weight[i][j]=command_par_value("weight",cmd->clone);
  for(start=0; start<n; start++) {
    if (strcmp(toks[start],"expr")==0) break;
  }
  start=start+2;
  nitem = end-start+1;
  cname=spec_join(&toks[start], nitem);
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
  printf("%d %d: %s\n",i,j,match2_cons_name[i][j]);
/*  printf("%e\n",expression_value(match2_cons_expr[i][j],type));*/
  return;
}

int match2_evaluate_exressions(int i, int k, double* fun_vec)
{
  int j;
  
  double rhs,lhs,r;/* RDM fork */
  char s;
  for(j=0;match2_cons_name[i][j]!=NULL;j++) {
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

void match2_delete_expressions()
{
  char rout_name[] = "match2_delete_expressions";

  int i,j;
  
  for(i=0;match2_macro_name[i]!=NULL;i++) {
    for(j=0;match2_cons_name[i][j]!=NULL;j++) {
      myfree(rout_name,match2_cons_name[i][j]);
      delete_expression(match2_cons_rhs[i][j]);
      delete_expression(match2_cons_lhs[i][j]);
    }
  }

}
