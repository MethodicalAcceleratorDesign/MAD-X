#include "madx.h"

// private interface

static void
match_prepare_varypos(void)
  /* keeps constraints from nodes, reexpands, adds constraints to nodes */
{
  struct node* node = current_sequ->ex_start;
  struct constraint_list** tmplist = mymalloc("match_prepare_varypos", current_sequ->n_nodes * sizeof *tmplist);
  int i = 0;
  while (node != NULL) {
    tmplist[i] = current_sequ->all_nodes[i]->cl; i++;
    if (node == current_sequ->ex_end) break;
    node = node->next;
  }
  expand_curr_sequ(0);
  i = 0;
  node = current_sequ->ex_start;
  while (node != NULL) {
    current_sequ->all_nodes[i]->cl = tmplist[i]; i++;
    if (node == current_sequ->ex_end) break;
    node = node->next;
  }
  myfree("match_prepare_varypos", tmplist);
}

static void
mtjacprint(int m, int n,double* jac,struct in_cmd* cmd)
{
  int i,j,k,l,t;
  double *SV,*U,*VT;
  double tmp;
  char *knobfilename, *jacfilename;
  FILE *knobfile=NULL, *jacfile=NULL;
  k=0;
  l=0;
  jacfilename=command_par_string("jacfile",cmd->clone);
  knobfilename=command_par_string("knobfile",cmd->clone);
  if (jacfilename) jacfile=fopen(jacfilename,"w");
  fprintf(prt_file, "\n\nJACOBIAN:\n");
  fprintf(prt_file, "%4s %16s %10s %20s\n","Node","Constraint","Variable","Derivative");
  fprintf(prt_file, "---------------------------------------------------------------\n");
  for(i=0;i<MAX_MATCH_MACRO;i++)
  {
    if (match2_macro_name[i]==NULL) break;
    for(j=0;j<MAX_MATCH_CONS;j++)
    {
      if (match2_cons_name[i][j]==NULL) break;
      if ( strcmp(match2_cons_name[i][j],"0")!=0 )
        if (jacfilename)
          fprintf(jacfile, "%16s:=%+20.10e",match2_cons_name[i][j],match2_cons_value[i][j]);
      for(l=0;l<n;l++)
      {
        fprintf(prt_file, "%4s ",match2_macro_name[i]);
        fprintf(prt_file, "%16s ",match2_cons_name[i][j]);
        fprintf(prt_file, "%10s ",command_par_string("name",stored_match_var->commands[l]));
        fprintf(prt_file, "%20.10e",jac[l*m+k]);
        fprintf(prt_file, "\n");
      if (jacfilename) if ( strcmp(match2_cons_name[i][j],"0")!=0 )
        fprintf(jacfile, "%+20.10e*%16s",jac[l*m+k],
                command_par_string("name",stored_match_var->commands[l]));
      }
      if (jacfilename) if ( strcmp(match2_cons_name[i][j],"0")!=0 )
      fprintf(jacfile, ";\n");
      k++;
    }
  }
  if (jacfilename) fclose(jacfile);
  fprintf(prt_file, "\n\nSINGULAR VALUE DECOMPOSITION:\n");
  SV = mymalloc_atomic("match_match", (total_const+total_vars)  * sizeof *SV);
  VT = mymalloc_atomic("match_match", (total_vars *total_vars)  * sizeof *VT);
  U  = mymalloc_atomic("match_match", (total_const*total_const) * sizeof *U );
  mtsvd_(&total_const,&total_vars,jac,SV,U,VT);
  k=0;
  l=0;
  fprintf(prt_file, "%-25s %12s %-34s\n",
          "Variable vector","---> Sing. val.","* Node constraint vector");
  fprintf(prt_file, "--------------------------------------------------------------------\n");
  for(i=0;i<imax(n,m);i++){
    for(j=0;j<imax(n,m);j++){
      if ( (i<n)&&(j<n)) {
        fprintf(prt_file, "%-12s",command_par_string("name",stored_match_var->commands[j]));
        fprintf(prt_file, "%12.5g ",VT[j*n+i]);
      }
      else { fprintf(prt_file, "%24s",""); }

      if ( (i<imin(n,m)) &&(j==0))   {
        fprintf(prt_file, "%12.5g ",SV[i]);
      }
      else { fprintf(prt_file, "%12s ",""); }

      if ( (i<m)&&(j<m) ) {
        if (match2_cons_name[k][l]==NULL) {
          k++;
          l=0;
        }
        /*        fprintf(prt_file, "%i %i",k,l);*/
        fprintf(prt_file, "%12.5g ",U[i*m+j]);
        /*        fprintf(prt_file, "%10s ",match2_macro_name[k]);*/
        fprintf(prt_file, "%10s ",match2_cons_name[k][l]);
        l++;
      }
      else { fprintf(prt_file, "%22s",""); }
      fprintf(prt_file, "\n");
    }
    l=0;
    k=0;
    fprintf(prt_file, "\n");
    fprintf(prt_file, "\n");
  }
  /* Print knob file */
  if (knobfilename) {
    knobfile=fopen(knobfilename,"w");
    k=0;
    l=0;
    for(i=0;i<n;i++){
      k=0; l=0;
      fprintf(knobfile, "%-12s := %+15.8e ",
              command_par_string("name",stored_match_var->commands[i]),
              get_variable(command_par_string("name",stored_match_var->commands[i])));
      for(j=0;j<m;j++){
        if (match2_cons_name[k][l]==NULL) { k++; l=0; }
        /* Compute M-1=V SV-1 UT*/
        /* M-1[i,j] = sum_t^min(m,n)   V[i,t]  / SV[t] * U[j,t] ) */
        /*          = sum_t^min(m,n)  VT[t,i] / SV[t] * U[j,t]  ) */
        /* in fortran M[i,j]=M[i+j*n] */
        tmp=0;
        for(t=0;t<imin(m,n);t++) {
          if (SV[t]>jac_cond) { tmp+=VT[t+n*i]/SV[t]*U[j+t*m]; }
          /*        fprintf(prt_file,"VT %d,%d,%12.5e ",t,i,VT[t+n*i]);*/
          /*        fprintf(prt_file,"SV %d,%12.5e ",t,SV[t]);*/
          /*        fprintf(prt_file,"U  %d,%d,%12.5e\n",j,t,U[j+t*m]);*/
        }
        /*      fprintf(prt_file,"M-1 %d,%d,%e \n",i,j,tmp);*/
        if ( strcmp(match2_cons_name[k][l],"0")!=0 ) {
          fprintf(knobfile, "%+15.8e*%s",tmp,match2_cons_name[k][l]);
        };
        l++;
      }
      fprintf(knobfile, ";\n");
    }
    fclose(knobfile);
  }
  myfree("match_match", SV);
  myfree("match_match", VT);
  myfree("match_match", U);
}

static void
match_action(struct in_cmd* cmd)
{
  int i;
  int iseed, iprint;
  int izero = 0;
  int local_calls;
  int local_call_lim;

  if (match_is_on == kMatch_PTCknobs)
  {
    madx_mpk_run(cmd);
    return;
  }


  if (stored_match_var->curr == 0)
  {
    warning("no vary command seen,","match command terminated");
    match_is_on = 0;
    set_option("match_is_on", &izero);
    return;
  }
  total_vars = stored_match_var->curr;
  penalty = zero;

  match_calls = command_par_value("calls", cmd->clone);
  match_tol = command_par_value("tolerance", cmd->clone);
  fprintf(prt_file, "number of variables:    %d\n", total_vars);
  fprintf(prt_file, "user given constraints: %d\n", comm_constraints->curr);
  fprintf(prt_file, "total constraints:      %d\n\n", total_const);
  vary_vect = new_double_array(total_vars);
  vary_dvect = new_double_array(total_vars);
  fun_vect = new_double_array(total_const);

  current_call_lim += match_calls;

  local_call_lim = command_par_value("calls", cmd->clone);
  local_calls = 0;

  if (strcmp(cmd->tok_list->p[0], "lmdif") == 0 && total_vars > total_const)
  {
    print_match_summary = 0;
    warning("number of variables larger than number of constraints:", "match command ignored");
    return;
  }
  else if (strcmp(cmd->tok_list->p[0], "lmdif") == 0 && total_vars <= total_const)
  {
    print_match_summary = 1;
    match_work[0] = new_double_array(total_vars);
    match_work[1] = new_double_array(total_vars*total_const);
    match_work[2] = new_double_array(total_vars);
    match_work[3] = new_double_array(total_vars);
    match_work[4] = new_double_array(total_vars);
    match_work[5] = new_double_array(total_vars);
    match_work[6] = new_double_array(total_vars);
    match_work[7] = new_double_array(total_const);
    match_work[8] = new_double_array(total_const);
    match_work[9] = new_double_array(total_vars);
    fprintf(prt_file, "START LMDIF:\n\n");
    mtlmdf_(&total_const, &total_vars,
            &match_tol, &local_calls, &local_call_lim,
            vary_vect->a, vary_dvect->a, fun_vect->a, match_work[0]->a,
            match_work[1]->a, match_work[2]->a, match_work[3]->a,
            match_work[4]->a, match_work[5]->a, match_work[6]->a,
            match_work[7]->a,match_work[9]->a);
  }
  else if (strcmp(cmd->tok_list->p[0], "jacobian") == 0)
  {
    print_match_summary = 0;

    jac_strategy = command_par_value("strategy", cmd->clone);
    jac_cool = command_par_value("cool", cmd->clone);
    jac_repeat = command_par_value("repeat", cmd->clone);
    jac_balance = command_par_value("balance", cmd->clone);
    jac_random = command_par_value("random", cmd->clone);
    jac_bisec = command_par_value("bisec", cmd->clone);
    jac_cond = command_par_value("cond", cmd->clone);
    match_work[0] = new_double_array(total_vars*total_const);
    match_work[1] = new_double_array(total_const);
    match_work[2] = new_double_array(total_const);
    match_work[3] = new_double_array(total_vars);
    match_work[4] = new_double_array(total_vars);
    fprintf(prt_file, "START JACOBIAN:\n\n");
    // size of arrays allocated in mtjac_
    unsigned long mem_size = (2*total_const*total_vars + 1020 * (total_const+total_vars)) * sizeof(double);
    if (mem_size > 8000000) {
      fprintf(prt_file, "Problem size is too large (may cause stability problems): %d*%d (%lu bytes)\n\n", total_const, total_vars, mem_size);
//      fatal_error("stack overflow will occur in fortran routines,","good bye");
    }
    mtjac_(&total_const, &total_vars,
           &jac_strategy, &jac_cool,&jac_balance, &jac_random,
           &jac_repeat,&jac_bisec,&jac_cond,&match_is_on,
           &match_tol, &local_calls, &local_call_lim,
           vary_vect->a, vary_dvect->a, fun_vect->a,
           match_work[0]->a,match_work[1]->a,match_work[2]->a,
           match_work[3]->a,match_work[4]->a);
/*    if (jac_strategy==2 && match_is_on==2) {*/
    if (jac_strategy==2) {
      mtjacprint(total_const,total_vars,match_work[0]->a,cmd);
    }
  }
  else if (strcmp(cmd->tok_list->p[0], "migrad") == 0 && total_vars <= total_const)
  {
    print_match_summary = 1;
    mig_strategy = command_par_value("strategy", cmd->clone);
    match_work[0] = new_double_array(total_vars*total_vars);
    match_work[1] = new_double_array(7*total_vars);
    match_work[2] = new_double_array(total_vars);
    match_work[3] = new_double_array(total_vars);
    match_work[4] = new_double_array(total_vars);
    match_work[5] = new_double_array(total_vars*total_vars);
    match_work[6] = new_double_array(total_vars);
    match_work[7] = new_double_array(total_vars);
    fprintf(prt_file, "START MIGRAD:\n\n");
    mtmigr_(&total_const, &total_vars, &mig_strategy,
            &match_tol, &local_calls, &local_call_lim,
            vary_vect->a, vary_dvect->a, fun_vect->a,
            match_work[0]->a, match_work[1]->a, match_work[2]->a,
            match_work[3]->a, match_work[4]->a, match_work[5]->a,
            match_work[6]->a, match_work[7]->a);
  }
  else if (strcmp(cmd->tok_list->p[0], "simplex") == 0)
  {
    print_match_summary = 1;
    match_work[0] = new_double_array(total_vars*(total_vars+1));
    match_work[1] = new_double_array(total_vars+1);
    match_work[2] = new_double_array(4*total_vars);
    fprintf(prt_file, "START SIMPLEX:\n\n");
    mtsimp_(&total_const, &total_vars,
            &match_tol, &local_calls, &local_call_lim,
            vary_vect->a, vary_dvect->a, fun_vect->a,
            match_work[0]->a, match_work[1]->a, match_work[2]->a);
  }
  else if (strcmp(cmd->tok_list->p[0], "siman") == 0)
  {
    print_match_summary = 1;
    iseed = 1;
    iprint = 0;
    match_work[0] = new_double_array(total_vars+1);
    match_i_work[1] = new_int_array(total_vars+1);
    match_work[2] = new_double_array(total_vars+1);
    match_work[3] = new_double_array(total_vars+1);
    match_work[4] = new_double_array(total_vars+1);
    match_work[5] = new_double_array(total_vars+1);
    match_work[6] = new_double_array(total_vars+1);
    fprintf(prt_file, "START SIMAN:\n\n");
    mtsa_(&total_const, &total_vars,
          &match_tol, &local_calls, &local_call_lim,
          vary_vect->a, fun_vect->a,&iseed,&iprint,
          match_work[0]->a, match_i_work[1]->i, match_work[2]->a,
          match_work[3]->a, match_work[4]->a, match_work[5]->a,
          match_work[6]->a);
  }
  for (i = 0; i < MATCH_WORK; i++)
    match_work[i] = delete_double_array(match_work[i]);
}

#if 0 // not used...
static void
match_cell(struct in_cmd* cmd)
{
  (void)cmd;
}
#endif

static void
match_constraint(struct in_cmd* cmd)
{
  struct node* c_node;

  if(match_is_on==2)
  {
    match2_constraint(cmd);
  }
  else if (match_is_on == kMatch_PTCknobs)
  {
    madx_mpk_addconstraint(in->buffers[in->curr]->c_a->c);
  }
  else
  { /* old match */
    // TODO: use all sequences!
    comm_constraints->curr=0;
    fill_constraint_list(1, cmd->clone, comm_constraints);
    struct select_iter* it = start_iter_select(cmd->clone, match_sequs, NULL);
    while (fetch_node_select(it, &c_node, NULL)) {
      update_node_constraints(c_node, comm_constraints);
    }
  }
}

static void
match_couple(struct in_cmd* cmd)
{
  (void)cmd;
}

static void
match_end(struct in_cmd* cmd)
{
  int i;
  struct node* c_node;
  int izero = 0;
  int ione  = 1;
  if (fun_vect == NULL )
  {
    fprintf(prt_file, "WARNING: No matching method selected.\n");
    return;
  }


  if (match_is_on==kMatch_UseMacro) {
    match2_end(cmd);
    return;
  }

  if (match_is_on == kMatch_PTCknobs) {
    madx_mpk_end();
    return;
  }

  /* write out all final constraint values and vary parameters */
  penalty = zero;
  if (get_option("varylength") != zero) match_prepare_varypos();
  pro_twiss();
  current_const = 0;
  fprintf(prt_file, "\n");
  fprintf(prt_file, "MATCH SUMMARY\n\n");
  fprintf(prt_file, "Node_Name                  Constraint   Type  Target Value       Final Value        Penalty\n");
  fprintf(prt_file, "--------------------------------------------------------------------------------------------------\n");
  print_match_summary = 1;
  set_option("match_summary", &print_match_summary);
  set_option("chrom_match", &ione);
  for (i = 0; i < match_num_seqs; i++)
  {
    current_twiss = local_twiss[i]->clone;
    pro_twiss();
    if (twiss_success && print_match_summary == 1)
    {
      collect_(&current_const, &penalty, fun_vect->a);
    }
  }
  fprintf(prt_file, "\n\n");

  fprintf(prt_file, "Final Penalty Function = %16.8e\n\n",penalty);



  fprintf(prt_file, "\n\n");
/*  fprintf(prt_file, "Variable                   Final Value        Lower Limit        Upper Limit\n");*/
/*  fprintf(prt_file, "-------------------------------------------------------------------------------\n");*/
/*  if ((print_match_summary == 1) && vary_vect )*/
/*  {*/
/*    |+ drop first int* passed variable mtgetc_(&stored_match_var->curr, vary_vect->a, vary_dvect->a); +|*/
/*    mtgetc_(vary_vect->a, vary_dvect->a); |+ mtgeti->mtgetc +|*/
/*  }*/
/*  fprintf(prt_file, "\n");*/
  match2_print_var(cmd);
  fprintf(prt_file, "END MATCH SUMMARY\n\n");
  print_match_summary = 0;
  set_option("match_summary", &print_match_summary);
  set_option("chrom_match", &izero);
  i = 0; set_option("match_local", &i); /* flag */
  /* End summary of matching job */

  for (i = 0; i < match_num_seqs; i++)
  {
    c_node = match_sequs->sequs[i]->ex_start;
    while(c_node != NULL) /* drop all constraint lists at nodes */
    {
      if (c_node->cl) c_node->cl = delete_constraint_list(c_node->cl);
      if (c_node == match_sequs->sequs[i]->ex_end) break;
      c_node = c_node->next;
    }
  }
  current_gweight = delete_command(current_gweight);
  current_weight = delete_command(current_weight);
  for (i = 0; i < comm_constraints->curr; i++)
    delete_constraint(comm_constraints->constraints[i]);
  for (i = 0; i < sequences->curr; i++)
  {
    if (sequences->sequs[i]->cl != NULL)  sequences->sequs[i]->cl
                                            = delete_constraint_list(sequences->sequs[i]->cl);
  }
  vary_vect = delete_double_array(vary_vect);
  vary_dvect = delete_double_array(vary_dvect);
  fun_vect = delete_double_array(fun_vect);
  comm_constraints->curr = 0;
  stored_match_var = delete_command_list(stored_match_var);
  for (i = 0; i < match_num_seqs; i++) if (local_twiss[i])
  {
    if (local_twiss[i]->clone) delete_command(local_twiss[i]->clone);
    local_twiss[i] = delete_in_cmd(local_twiss[i]);
  }
  vary_cnt = 0;
  match_is_on = 0;
  set_option("match_is_on", &izero);
  match_sequs->curr = 0;
  match_sequs->list->curr = 0;
  current_call_lim = 0;
  current_calls = 0;
  total_const = 0;
  set_option("twiss_print", &keep_tw_print);


  fprintf(prt_file, "VARIABLE \"TAR\" SET TO %16.8e\n",penalty);
/*  sprintf(assign_cmd,"tar= %16.8e ;",penalty);*/
/*  pro_input(assign_cmd);*/
  set_variable("tar",&penalty);
}

static void
match_fix(struct in_cmd* cmd)
{
  (void)cmd;
}

static void
match_global(struct in_cmd* cmd)
{
  struct command_parameter* cp;
  struct name_list* nl = cmd->clone->par_names;
  struct sequence* sequ;
  int pos, n, low, up;
  pos = name_list_pos("sequence", nl);
  if((pos>=0)?nl->inform[pos]:0) /* sequence specified */
  {
    cp = cmd->clone->par->parameters[pos];
    for (n = 0; n < match_sequs->curr; n++)
    {
      if (strcmp(cp->string, match_sequs->sequs[n]->name) == 0) break;
    }
    if (n == match_sequs->curr)
    {
      warning(cp->string," :sequence not selected by MATCH, skipped");
      return;
    }
    low = up = n;
  }
  else
  {
    low = 0; up = match_sequs->curr - 1;
  }
  for (n = low; n <= up; n++)
  {
    sequ = match_sequs->sequs[n];
    if (sequ->cl == NULL) sequ->cl = new_constraint_list(10);
    comm_constraints->curr=0;
    fill_constraint_list(2, cmd->clone, comm_constraints);
    update_sequ_constraints(sequ, comm_constraints);
  }
  print_match_summary = 1;
}

static void
match_gweight(struct in_cmd* cmd)
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* plc = cmd->clone->par;
  struct command_parameter_list* pl = current_gweight->par;
  int j;

  (void)cmd;
  for (j = 0; j < pl->curr; j++)
  {
    if (nl->inform[j])
    {
      pl->parameters[j]->double_value = plc->parameters[j]->double_value;
      pl->parameters[j]->expr = plc->parameters[j]->expr;
    }
  }
}

static void
match_level(struct in_cmd* cmd)
{
  (void)cmd;
}

static void
match_match(struct in_cmd* cmd)
{
  struct name_list* tnl; /* local name list for TWISS input definition */
  double match_beta_value;
  struct command* comm;
  struct command_parameter* cp;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  struct sequence* sequ;
  int i, j, pos, n, spos, tpos, chrom_flg;
  int izero = 0, ione = 1, slow;
  char *dpp;

  if (match_is_on)
  {
    warning("already inside match:", "ignored");
    return;
  }
  set_option("match_is_on", &ione);
  keep_tw_print = get_option("twiss_print");
  set_option("twiss_print", &izero);

  pos=name_list_pos("use_macro", nl);
  if((pos>=0)?nl->inform[pos]:0) {
    match2_match(cmd);
    return;
  }
  pos=name_list_pos("use_ptcknob", nl);
  if((pos>=0)?nl->inform[pos]:0) {
    match_is_on = kMatch_PTCknobs;
    madx_mpk_prepare();
    return;
  }
  match_is_on = 1;
  pos = name_list_pos("sequence", nl);
  fprintf(prt_file, "START MATCHING\n\n");

  if((pos>=0)?nl->inform[pos]:0) /* sequence specified */
  {
    cp = cmd->clone->par->parameters[pos];
    if ((n = cp->m_string->curr) > match_sequs->max)
    {
      warning("excess sequence names", "skipped");
      n =  match_sequs->max;
    }
    for (i = 0; i < n; i++)
    {
      if ((pos = name_list_pos(cp->m_string->p[i], sequences->list)) < 0)
        warning("unknown sequence", "skipped");
      else
      {
        if ((sequ = sequences->sequs[pos]) == NULL
            || sequ->ex_start == NULL)
        {
          warning(sequ->name," :sequence not active, match killed");
          match_is_on = 0;
          return;
        }
        match_sequs->sequs[match_sequs->curr++] = sequ;
        add_to_name_list(sequ->name, 0, match_sequs->list);
      }
    }
    if (match_sequs->curr == 0)
    {
      warning("no active sequence,","match killed");
      match_is_on = 0;
      return;
    }
  }
  else {
    match_sequs->sequs[match_sequs->curr++] = current_sequ;
    add_to_name_list(current_sequ->name, 0, match_sequs->list);
  }
  pos = name_list_pos("vlength", nl);
  if((pos>=0)?nl->inform[pos]:0) i = pl->parameters[pos]->double_value;
  else i = 0;
  set_option("varylength", &i);

  /* START SLOW-CHECK */
  pos = name_list_pos("slow", nl);
  if((pos>=0)?nl->inform[pos]:0) slow = pl->parameters[pos]->double_value;
  else slow = 0;
  set_option("slow_match", &slow);
  /* END SLOW-CHECK */

  /* START CHK-SEQ */
  current_match = cmd->clone;
  pos = name_list_pos("sequence", nl);
  if((pos>=0)?nl->inform[pos]:0) /* sequence specified */
  {
    cp = cmd->clone->par->parameters[pos];
    match_num_seqs = cp->m_string->curr;
    fprintf(prt_file, "number of sequences: %d\n", match_num_seqs);
    for (i = 0; i < match_num_seqs; i++)
    {
      match_seqs[i] = buffer(cp->m_string->p[i]);
      /* check if the specified sequence is defined */
      if ((pos = name_list_pos(match_seqs[i], sequences->list)) > -1)
      {
        current_sequ = sequences->sequs[pos];
        fprintf(prt_file, "sequence name: %s\n", current_sequ->name);

        /* START defining a TWISS input command for each sequence */
        local_twiss[i] = new_in_cmd(10);
        local_twiss[i]->type = 0;
        local_twiss[i]->clone = local_twiss[i]->cmd_def
          = clone_command(find_command("twiss", defined_commands));
        tnl = local_twiss[i]->cmd_def->par_names;
        tpos = name_list_pos("sequence", tnl);

        local_twiss[i]->cmd_def->par->parameters[tpos]->string = match_seqs[i];
        local_twiss[i]->cmd_def->par_names->inform[tpos] = 1;
        if (slow)
	  {
           tpos = name_list_pos("save", tnl);
           local_twiss[i]->cmd_def->par->parameters[tpos]->string = NULL;
           local_twiss[i]->cmd_def->par_names->inform[tpos] = 1;
	  }
        /* END defining a TWISS input command for each sequence */
      }
      else
      {
        warning("unknown sequence ignored:", match_seqs[i]);
        return;
      }
    }
  }
  else
  {
    /* START defining a TWISS input command for default sequence */
    match_num_seqs = 1;
    local_twiss[0] = new_in_cmd(10);
    local_twiss[0]->type = 0;
    local_twiss[0]->clone = local_twiss[0]->cmd_def
      = clone_command(find_command("twiss", defined_commands));
    tnl = local_twiss[0]->cmd_def->par_names;
    tpos = name_list_pos("sequence", tnl);
    local_twiss[0]->cmd_def->par->parameters[tpos]->string = current_sequ->name;
    local_twiss[0]->cmd_def->par_names->inform[tpos] = 1;
    if (slow)
      {
       tpos = name_list_pos("save", tnl);
       local_twiss[0]->cmd_def->par->parameters[tpos]->string = NULL;
       local_twiss[0]->cmd_def->par_names->inform[tpos] = 1;
      }
     /* END defining a TWISS input command for default sequence */
  }
  if (current_sequ == NULL || current_sequ->ex_start == NULL)
  {
    warning("MATCH called without active sequence,", "ignored");
    return;
  }
  /* END CHK-SEQ */

  for (i = 0; i < match_num_seqs; i++)
  {
    for (j = 0; j < local_twiss[i]->cmd_def->par->curr; j++)
    {
      tnl = local_twiss[i]->cmd_def->par_names;
      tpos = name_list_pos("sequence", tnl);
      if (slow) spos = name_list_pos("save", tnl);
      else      spos = -1;
      if (j != tpos  && j != spos)
      local_twiss[i]->cmd_def->par_names->inform[j] = 0;
    }
  }

  /* START CHK-BETA-INPUT */
  /* START CHK-BETA0 */
  pos = name_list_pos("beta0", nl);
  if((pos>=0)?nl->inform[pos]:0) /* beta0 specified */
  {
    cp = cmd->clone->par->parameters[pos];
    match_num_beta = cp->m_string->curr;
    fprintf(prt_file, "number of beta0s: %d\n", match_num_beta);
    for (i = 0; i < match_num_beta; i++)
    {
      match_beta[i] = buffer(cp->m_string->p[i]);
      fprintf(prt_file, "BETA0 name: %s\n", match_beta[i]);
      /* START defining a TWISS input command for each sequence */
      tnl = local_twiss[i]->cmd_def->par_names;
      tpos = name_list_pos("beta0", tnl);
      local_twiss[i]->cmd_def->par_names->inform[tpos] = 1;
      local_twiss[i]->cmd_def->par->parameters[tpos]->string = match_beta[i];
      /* END defining a TWISS input command for each sequence */
    }
  }
  /* END CHK-BETA0 */

  /* START CHK-RANGE */
  pos = name_list_pos("range", nl);
  if((pos>=0)?nl->inform[pos]:0) /* range specified */
  {
    cp = cmd->clone->par->parameters[pos];
    match_num_range = cp->m_string->curr;
    for (i = 0; i < match_num_range; i++)
    {
      match_range[i] = buffer(cp->m_string->p[i]);
      /* START adding range to TWISS input command for each sequence */
      tnl = local_twiss[i]->cmd_def->par_names;
      tpos = name_list_pos("range", tnl);
      local_twiss[i]->cmd_def->par_names->inform[tpos] = 1;
      local_twiss[i]->cmd_def->par->parameters[tpos]->string
        = match_range[i];
      /* END adding range to TWISS input command for each sequence */
    }
  }
  /* END CHK-RANGE */

  /* START CHK-USEORBIT */
  pos = name_list_pos("useorbit", nl);
  if((pos>=0)?nl->inform[pos]:0) /* useorbit specified */
  {
    cp = cmd->clone->par->parameters[pos];
    for (i = 0; i < cp->m_string->curr; i++)
    {
      /* START adding useorbit to TWISS input command for each sequence */
      tnl = local_twiss[i]->cmd_def->par_names;
      tpos = name_list_pos("useorbit", tnl);
      local_twiss[i]->cmd_def->par_names->inform[tpos] = 1;
      local_twiss[i]->cmd_def->par->parameters[tpos]->string
        = buffer(cp->m_string->p[i]);
      /* END adding range to TWISS input command for each sequence */
    }
  }
  /* END CHK-USEORBIT */

  /* START CHK-KEEPORBIT */
  pos = name_list_pos("keeporbit", nl);
  if((pos>=0)?nl->inform[pos]:0) /* keeporbit specified */
  {
    cp = cmd->clone->par->parameters[pos];
    for (i = 0; i < cp->m_string->curr; i++)
    {
      /* START adding keeporbit to TWISS input command for each sequence */
      tnl = local_twiss[i]->cmd_def->par_names;
      tpos = name_list_pos("keeporbit", tnl);
      local_twiss[i]->cmd_def->par_names->inform[tpos] = 1;
      local_twiss[i]->cmd_def->par->parameters[tpos]->string
        = buffer(cp->m_string->p[i]);
      /* END adding range to TWISS input command for each sequence */
    }
  }
  /* END CHK-KEEPORBIT */

  /* START CHK-R-MATRIX */
  pos = name_list_pos("rmatrix", nl);
  if((pos>=0)?nl->inform[pos]:0) /* rmatrix specified */
  {
    cp = cmd->clone->par->parameters[pos];
    for (i = 0; i < match_num_seqs; i++)
    {
      /* START adding rmatrix to TWISS input command for each sequence */
      tnl = local_twiss[i]->cmd_def->par_names;
      tpos = name_list_pos("rmatrix", tnl);
      local_twiss[i]->cmd_def->par_names->inform[tpos] = 1;
      local_twiss[i]->cmd_def->par->parameters[tpos]->double_value
        = 1;
      /* END adding rmatrix to TWISS input command for each sequence */
    }
  }
  /* END CHK-R-MATRIX */

  /* START CHK-CHROM */
  pos = name_list_pos("chrom", nl);
  if (pos > 0)
   {
     chrom_flg = cmd->clone->par->parameters[pos]->double_value;
   }
  else
   {
     chrom_flg = 0;
   }

  if(chrom_flg) /* chrom specified */
  {
    for (i = 0; i < match_num_seqs; i++) {
      /* START adding chrom to TWISS input command for each sequence */
      tnl = local_twiss[i]->cmd_def->par_names;
      tpos = name_list_pos("chrom", tnl);
      local_twiss[i]->cmd_def->par_names->inform[tpos] = 1;
      local_twiss[i]->cmd_def->par->parameters[tpos]->double_value = 1;
      set_option("twiss_chrom",&ione);
      /* END adding chrom to TWISS input command for each sequence */
    }
  }
  else set_option("twiss_chrom",&izero);
  /* END CHK-CHROM */

  /* START CHK-SECTORMAP */
  pos = name_list_pos("sectormap", nl);
  if((pos>=0)?nl->inform[pos]:0) /* sectormap specified */
  {
    cp = cmd->clone->par->parameters[pos];
    for (i = 0; i < match_num_seqs; i++)
    {
      /* START adding sectormap to TWISS input command for each sequence */
      tnl = local_twiss[i]->cmd_def->par_names;
      tpos = name_list_pos("sectormap", tnl);
      local_twiss[i]->cmd_def->par_names->inform[tpos] = 1;
      local_twiss[i]->cmd_def->par->parameters[tpos]->double_value
        = 1;
      /* END adding sectormap to TWISS input command for each sequence */
    }
  }
  /* END CHK-CHROM */

  /* START CHK-DELTAP */
  pos = name_list_pos("deltap", nl);
  if((pos>=0)?nl->inform[pos]:0) /* deltap specified */
  {
    cp = cmd->clone->par->parameters[pos];
    for (i = 0; i < match_num_seqs; i++)
    {
      /* START adding deltap to TWISS input command for each sequence */
      tnl = local_twiss[i]->cmd_def->par_names;
      tpos = name_list_pos("deltap", tnl);
      local_twiss[i]->cmd_def->par_names->inform[tpos] = 1;
      dpp = mymalloc_atomic("match_match", 20 * sizeof *dpp);
      sprintf(dpp,"%e",cp->double_array->a[0]);
      local_twiss[i]->cmd_def->par->parameters[tpos]->string = dpp;
      /*       END adding deltap to TWISS input command for each sequence */
      /*      fprintf(prt_file, "entry value: %f\n", cp->double_array->a[0]);*/
      /*      twiss_deltas->curr=1;*/
      /*      twiss_deltas->a[0]=cp->double_array->a[0];*/
    }
  }
  /* END CHK-DELTAP */

  /* START CHK-ENTRIES of TYPE DOUBLE-REAL */
  for (j = 0; j < nl->curr; j++)
  {
    if(nl->inform[j]) /* initial conditions are specified
                         at the match command level */
    {
      cp = cmd->clone->par->parameters[j];
      match_num_beta = 0;
      if(cp->type == 2)
      {
        match_num_beta = 1;
      }
      if(cp->type == 12)
      {
        match_num_beta = cp->double_array->curr;
      }
      if(match_num_beta > 0)
      {
        fprintf(prt_file, "entry name: %s\n", cmd->clone->par_names->names[j]);
        fprintf(prt_file, "number of entries: %d\n", match_num_beta);
      }
      for (i = 0; i < match_num_beta; i++)
      {
        /* START defining a TWISS input command for each sequence */
        match_beta_value = cp->double_array->a[i];
        fprintf(prt_file, "entry value: %f\n", match_beta_value);
        tnl = local_twiss[i]->cmd_def->par_names;
        tpos = name_list_pos(cmd->clone->par_names->names[j], tnl);
        local_twiss[i]->cmd_def->par_names->inform[tpos] = 1;
        local_twiss[i]->cmd_def->par->parameters[tpos]->double_value
          = cp->double_array->a[i];
        /*        fprintf(prt_file, "entry value: %f %f\n", match_beta_value, cp->double_array->a[i]);*/
        /* END defining a TWISS input command for each sequence */
      }
    }
  }
  /* END CHK-ENTRIES of TYPE DOUBLE-REAL */
  /* END CHK-BETA-INPUT */

  /* START generating a TWISS table via 'pro_twiss' */
  for (i = 0; i < match_num_seqs; i++)
  {
    /* fprintf(prt_file, "%s %s\n", "call TWISS from MATCH: sequence =",
       match_sequs->sequs[i]->name); */
    current_twiss = local_twiss[i]->clone;
    pro_twiss();
  }
  /* END generating a TWISS table via 'pro_twiss' */

  /* for (i = 0; i < match_sequs->curr; i++)
     puts(match_sequs->sequs[i]->name); */
  if ((comm = find_command("weight", defined_commands)) != NULL)
    current_weight = clone_command(comm);
  else fatal_error("no weight command in dictionary,","good bye");
  if ((comm = find_command("gweight", defined_commands)) != NULL)
    current_gweight = clone_command(comm);
  else fatal_error("no gweight command in dictionary,","good bye");

}

static void
match_rmatrix(struct in_cmd* cmd)
{
  (void)cmd;
}

static void
match_tmatrix(struct in_cmd* cmd)
{
  (void)cmd;
}

static void
match_vary(struct in_cmd* cmd)
{
  int pos;
  double value;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;

  if (stored_match_var == NULL) stored_match_var = new_command_list("vary", 100);
  pos = name_list_pos("name", nl);
  if ((pos>=0)?nl->inform[pos]:0)
  {
    value=get_variable(pl->parameters[pos]->string);
    set_value("vary","init",&value);
    cmd->clone_flag = 1;
    add_to_command_list(pl->parameters[pos]->string,
                        cmd->clone, stored_match_var, 1);
  }
}

static void
match_weight(struct in_cmd* cmd)
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* plc = cmd->clone->par;
  struct command_parameter_list* pl = current_weight->par;
  int j;
  for (j = 0; j < pl->curr; j++)
  {
    if (nl->inform[j])
    {
      pl->parameters[j]->double_value = plc->parameters[j]->double_value;
      pl->parameters[j]->expr = plc->parameters[j]->expr;
    }
  }
}

// public interface

void
pro_match(struct in_cmd* cmd)
  /* controls the matching module */
{
  /* changed the sequence of if statements so that MAD
     can go through the whole matching sequence */

  if (strcmp(cmd->tok_list->p[0], "match") == 0)
  {
    match_match(cmd);
  }
  else if (match_is_on == 0)
  {
    warning("no MATCH command seen,","ignored");
    return;
  }
  else if (strcmp(cmd->tok_list->p[0], "endmatch") == 0)
  {
    match_end(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "migrad") == 0 ||
           strcmp(cmd->tok_list->p[0], "lmdif") == 0 ||
           strcmp(cmd->tok_list->p[0], "jacobian") == 0 ||
           strcmp(cmd->tok_list->p[0], "simplex") == 0)
  {
    match_action(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "constraint") == 0)
  {
    match_constraint(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "couple") == 0)
  {
    match_couple(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "fix") == 0)
  {
    match_fix(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "global") == 0)
  {
    match_global(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "level") == 0)
  {
    match_level(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "vary") == 0)
  {
    match_vary(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "weight") == 0)
  {
    match_weight(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "gweight") == 0)
  {
    match_gweight(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "rmatrix") == 0)
  {
    match_rmatrix(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "tmatrix") == 0)
  {
    match_tmatrix(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "global") == 0)
  {
    match_global(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "use_macro") == 0)
  {
    match2_macro(cmd);
  }
}

int
next_vary(char* name, int* name_l, double* lower, double* upper, double* step, int* slope, double* opt)
  /* returns the next variable to be varied during match;
     0 = none, else count */
{
  int pos;
  double l_step;
  const char* v_name;
  struct name_list* nl;
  struct command* comm;
  struct command_parameter_list* pl;
  if (vary_cnt == stored_match_var->curr)
  {
    vary_cnt = 0; return 0;
  }
  comm = stored_match_var->commands[vary_cnt];
  nl = comm->par_names;
  pl = comm->par;
  pos = name_list_pos("name", nl);
  v_name = (pos>=0)?pl->parameters[pos]->string:"";
  strfcpy(name, v_name, *name_l);
  *lower = command_par_value("lower", comm);
  *upper = command_par_value("upper", comm);
  if ((l_step = command_par_value("step", comm)) < ten_m_12) l_step = ten_m_12;
  *step = l_step;
  *slope = command_par_value("slope", comm);
  *opt = command_par_value("opt", comm);
  return ++vary_cnt;
}

// public interface used by fortran code

void
mtcond(int* print_flag, int* nf, double* fun_vec, int* stab_flag)
{
  int i,j,k=0;
  char execute[40];
  static int nconserrs = 0; /*number of call finihed with error*/

  if (match_is_on==2)
  {
    for(i=0; i<MAX_MATCH_MACRO; i++)
    {
      if (match2_macro_name[i]==NULL) break;

      sprintf(execute,"exec, %s;",match2_macro_name[i]);
      match_is_on=0;
      pro_input(execute);
      match_is_on=2;
      if (!geterrorflag())
      {
        *stab_flag=0;
        k=match2_evaluate_exressions(i,k,fun_vec);
        nconserrs = 0;
      }
      else
      {
        nconserrs++;
        if (nconserrs > 5)
        { /*return the error code only after 5 consecutive fails*/
          *stab_flag=1; return;
        }
        else
        { /*otherwise just put all the constraints to max double value*/
          *stab_flag=0;
          for(j=0; j < *nf; j++)
          {
            fun_vec[j] = DBL_MAX;
          }
        }
      }
    }
  }
  else
  { /* old match */
    current_const = 0;
    penalty = zero;
    set_option("match_print", print_flag);
    /* mtgeti_(&stored_match_var->curr, vary_vect->a, vary_dvect->a); */
    for (i = 0; i < match_num_seqs; i++)
    {
      /* fprintf(prt_file, "%s %s\n", "call TWISS from matching: sequence=",
         match_sequs->sequs[i]->name); */
      current_twiss = local_twiss[i]->clone;
      if (get_option("varylength") != zero) match_prepare_varypos();

      if (get_option("rmatrix") != zero) fprintf(prt_file, "%s\n", "call TWISS with RMATRIX");
      /* match with chrom */
      if (get_option("chrom") != zero) fprintf(prt_file, "%s\n", "call TWISS with CHROM");
      /* match with deltap */
      if (get_option("deltap") != zero) fprintf(prt_file, "%s\n", "call TWISS with DELTAP");
      /* match with deltap */
      if (get_option("sectormap") != zero) fprintf(prt_file, "%s\n", "call TWISS with SECTORMAP");

      pro_twiss();
      if (twiss_success && !geterrorflag())
      {
        *stab_flag = 0;
        collect_(&current_const, &penalty, fun_vec);
        /* Do not write the penalty function here.
          'lmdif' still needs to check if the iteration was succesfull!
          the penalty function should only be printed from the
          'lmdif' routine. */
        /* fprintf(prt_file, "penalty function: %e\n", penalty); */
      }
      else
      {
        *stab_flag = 1; break;  /* Twiss failed - give up */
      }
    }
  } /* old match */
}

int
mtputconsname(char* noden, int* nodei , char* consn, int* consi)
{
  int i,j;
  i = *nodei -1;
  j = *consi -1;
  match2_macro_name[i] = mymalloc_atomic("match_match", 20 * sizeof *match2_macro_name[0]);
  strncpy(match2_macro_name[i],noden,20);
  match2_macro_name[i][19]='\0';
  match2_cons_name[i][j] = mymalloc_atomic("match_match", 20 * sizeof *match2_cons_name[0][0]);
  strncpy(match2_cons_name[i][j],consn,20);
  match2_cons_name[i][j][19]='\0';
  return 0;
}
