void mtsa_(int*, int*, double*, int*, int*, double*, double*, int*, int*,
           double*, int*, double*, double*, double*, double*, double*);
void match_action(struct in_cmd* cmd)
{
  int i;
  int iseed, iprint;

  if (stored_match_var->curr == 0)
  {
    warning("no vary command seen,","match command terminated");
    match_is_on = 0;
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
            &match_tol, &current_calls, &current_call_lim,
            vary_vect->a, vary_dvect->a, fun_vect->a, match_work[0]->a,
            match_work[1]->a, match_work[2]->a, match_work[3]->a,
            match_work[4]->a, match_work[5]->a, match_work[6]->a,
            match_work[7]->a,match_work[8]->a,match_work[9]->a);
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
    match_work[0] = new_double_array(total_vars*total_const);
    match_work[1] = new_double_array(total_const);
    match_work[2] = new_double_array(total_const);
    match_work[3] = new_double_array(total_vars);
    match_work[4] = new_double_array(total_vars);
    fprintf(prt_file, "START JACOBIAN:\n\n");
    mtjac_(&total_const, &total_vars,
           &jac_strategy, &jac_cool,&jac_balance, &jac_random,
           &jac_repeat,&jac_bisec,&match_is_on,
           &match_tol, &current_calls, &current_call_lim,
           vary_vect->a, vary_dvect->a, fun_vect->a,
           match_work[0]->a,match_work[1]->a,match_work[2]->a,
           match_work[3]->a,match_work[4]->a);
   if (jac_strategy==2 && match_is_on==2) {
     printf("Print Jacobian\n");
     mtjacprint(total_const,total_vars,match_work[0]->a);
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
            &match_tol, &current_calls, &current_call_lim,
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
            &match_tol, &current_calls, &current_call_lim,
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
          &match_tol, &current_calls, &current_call_lim,
          vary_vect->a, fun_vect->a,&iseed,&iprint,
          match_work[0]->a, match_i_work[1]->i, match_work[2]->a,
          match_work[3]->a, match_work[4]->a, match_work[5]->a,
          match_work[6]->a);
  }
  for (i = 0; i < MATCH_WORK; i++)
    match_work[i] = delete_double_array(match_work[i]);
}

void mtcond(int* print_flag, int* nf, double* fun_vec, int* stab_flag)
{
  int i;
  int j,k=0;/* RDM fork */
  double rhs,lhs,r;/* RDM fork */
  char execute[40],s;/* RDM fork */

  if (match_is_on==2) { /* RDM fork */
    for(i=0;match2_macro_name[i]!=NULL;i++) {
      sprintf(execute,"exec, %s;",match2_macro_name[i]);
      pro_input(execute);
/*      if (twiss_success) {*/
        *stab_flag=0;
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
/*      } else {*/
/*        *stab_flag=1; return;*/
/*      }*/
    }
  } else { /* RDM old match */
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
    /* RDM 22.9.2005: match with chrom */
    if (get_option("chrom") != zero) fprintf(prt_file, "%s\n", "call TWISS with CHROM");
    /* RDM 12.12.2005: match with deltap */
    if (get_option("deltap") != zero) fprintf(prt_file, "%s\n", "call TWISS with DELTAP");
    /* RDM 16.12.2005: match with deltap */
    if (get_option("sectormap") != zero) fprintf(prt_file, "%s\n", "call TWISS with SECTORMAP");

    pro_twiss();
    if (twiss_success)
    {
      *stab_flag = 0;
      collect_(&current_const, &penalty, fun_vec);
      /*OB 5.3.2002: Do not write the penalty function here.
        'lmdif' still needs to check if the itteration was succesfull!
        the penalty function should only be printed from the
        'lmdif' routine. */
      /* fprintf(prt_file, "penalty function: %e\n", penalty); */
    }
    else
    {
      *stab_flag = 1; break;  /* Twiss failed - give up */
    }
  }
  } /*RDM old match */
}

void match_cell(struct in_cmd* cmd)
{
}

void match_constraint(struct in_cmd* cmd)
{
  char* name;
  struct command_parameter* cp;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  struct sequence* sequ;
  struct node* nodes[2];
  struct node* c_node;
  int pos, k, n, low, up;

  if(match_is_on==2) { /* RDM fork */
    match2_constraint(cmd);
  }
  else { /* RDM old match */
    pos = name_list_pos("sequence", nl);
    if(nl->inform[pos]) /* sequence specified */
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
      pos = name_list_pos("range", nl);
      if (pos > -1 && nl->inform[pos])  /* parameter has been read */
      {
        name = pl->parameters[pos]->string;
        if ((k = get_ex_range(name, sequ, nodes)) == 0)  return;
      }
      else
      {
        nodes[0] = sequ->ex_start; nodes[1] = sequ->ex_end;
      }
      c_node = nodes[0];
      comm_constraints->curr=0;
      fill_constraint_list(1, cmd->clone, comm_constraints);
      while (c_node)
      {
        if (pass_select(c_node->p_elem->name, cmd->clone) != 0)
          update_node_constraints(c_node, comm_constraints);
        if (c_node == nodes[1]) break;
        c_node = c_node->next;
      }
    }
    /* OB 12.2.2002 */
  }
}

void match_couple(struct in_cmd* cmd)
{
}

void match_end(struct in_cmd* cmd)
{
  int i;
  struct node* c_node;

  if (match_is_on==2) {
    match2_end(cmd);
    return;
  }
  
  /* OB 5.3.2002: write out all final constraint values and vary parameters */
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
  fprintf(prt_file, "Variable                   Final Value        Lower Limit        Upper Limit\n");
  fprintf(prt_file, "-------------------------------------------------------------------------------\n");
  if (print_match_summary == 1)
  {
    /* drop first int* passed variable mtgetc_(&stored_match_var->curr, vary_vect->a, vary_dvect->a); 05.02.2005 */
    mtgetc_(vary_vect->a, vary_dvect->a); /* mtgeti->mtgetc JMJ, 8/4/2003 */
  }
  fprintf(prt_file, "\n");
  fprintf(prt_file, "END MATCH SUMMARY\n\n");
  print_match_summary = 0;
  set_option("match_summary", &print_match_summary);
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
  match_sequs->curr = 0;
  current_call_lim = 0;
  current_calls = 0;
  total_const = 0;
  set_option("twiss_print", &keep_tw_print);


  fprintf(prt_file, "EVALUATING \"tar= %16.8e;\"\n",penalty);
/*  sprintf(assign_cmd,"tar= %16.8e ;",penalty);*/
/*  pro_input(assign_cmd);*/
  set_variable("tar",&penalty);
}

void match_fix(struct in_cmd* cmd)
{
}

void match_global(struct in_cmd* cmd)
{
  struct command_parameter* cp;
  struct name_list* nl = cmd->clone->par_names;
  struct sequence* sequ;
  int pos, n, low, up;
  pos = name_list_pos("sequence", nl);
  if(nl->inform[pos]) /* sequence specified */
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
  /* OB 15.3.2002 */
}

void match_gweight(struct in_cmd* cmd)
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* plc = cmd->clone->par;
  struct command_parameter_list* pl = current_gweight->par;
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

void match_level(struct in_cmd* cmd)
{
}

void match_match(struct in_cmd* cmd)
{
  struct name_list* tnl; /* OB 31.1.2002: local name list for TWISS input definition */
  double match_beta_value; /* OB 31.1.2002 */
  struct command* comm;
  struct command_parameter* cp;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  struct sequence* sequ;
  int i, j, pos, n, tpos;
  int izero = 0;
  /* RDM 19/12/2005*/
  int ione = 1;
  char *dpp;
  /* RDM 19/12/2005*/

  if (match_is_on)
  {
    warning("already inside match:", "ignored");
    return;
  }
  keep_tw_print = get_option("twiss_print");
  set_option("twiss_print", &izero);
  pos=name_list_pos("use_macro", nl);
  if(nl->inform[pos]) {/* RDM FORK */
    match2_match(cmd);
    return;
  }/* RDM FORK */
  else { /* RDM  old match */
    match_is_on = 1;
    pos = name_list_pos("sequence", nl);
    fprintf(prt_file, "START MATCHING\n\n");
    if(nl->inform[pos]) /* sequence specified */
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
        }
      }
      if (match_sequs->curr == 0)
      {
        warning("no active sequence,","match killed");
        match_is_on = 0;
        return;
      }
    }
    else match_sequs->sequs[match_sequs->curr++] = current_sequ;
    pos = name_list_pos("vlength", nl);
    if(nl->inform[pos]) i = pl->parameters[pos]->double_value;
    else i = 0;
    set_option("varylength", &i);
    /* START CHK-SEQ; OB 1.2.2002 */
    current_match = cmd->clone;
    pos = name_list_pos("sequence", nl);
    if(nl->inform[pos]) /* sequence specified */
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
          /* OB 1.2.2002 */
          /* START defining a TWISS input command for each sequence */
          local_twiss[i] = new_in_cmd(10);
          local_twiss[i]->type = 0;
          local_twiss[i]->clone = local_twiss[i]->cmd_def
            = clone_command(find_command("twiss", defined_commands));
          tnl = local_twiss[i]->cmd_def->par_names;
          tpos = name_list_pos("sequence", tnl);
          local_twiss[i]->cmd_def->par->parameters[tpos]->string = match_seqs[i];
          local_twiss[i]->cmd_def->par_names->inform[tpos] = 1;
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
      /* END defining a TWISS input command for default sequence */
    }
    if (current_sequ == NULL || current_sequ->ex_start == NULL)
    {
      warning("MATCH called without active sequence,", "ignored");
      return;
    }
    /* END CHK-SEQ; OB 1.2.2002 */

    for (i = 0; i < match_num_seqs; i++)
    {
      for (j = 0; j < local_twiss[i]->cmd_def->par->curr; j++)
      {
        tnl = local_twiss[i]->cmd_def->par_names;
        tpos = name_list_pos("sequence", tnl);
        if (j != tpos) local_twiss[i]->cmd_def->par_names->inform[j] = 0;
      }
    }

    /* START CHK-BETA-INPUT; OB 1.2.2002 */
    /* START CHK-BETA0; OB 23.1.2002 */
    pos = name_list_pos("beta0", nl);
    if(nl->inform[pos]) /* beta0 specified */
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
    /* END CHK-BETA0; OB 23.1.2002 */

    /* START CHK-RANGE; HG 12.11.2002 */
    pos = name_list_pos("range", nl);
    if(nl->inform[pos]) /* range specified */
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
    /* END CHK-RANGE; OB 12.11.2002 */

    /* START CHK-USEORBIT; HG 28.1.2003 */
    pos = name_list_pos("useorbit", nl);
    if(nl->inform[pos]) /* useorbit specified */
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
    /* END CHK-USEORBIT; HG 28.1.2003 */

    /* START CHK-KEEPORBIT; HG 28.1.2003 */
    pos = name_list_pos("keeporbit", nl);
    if(nl->inform[pos]) /* keeporbit specified */
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
    /* END CHK-KEEPORBIT; HG 28.1.2003 */

    /* START CHK-R-MATRIX; OB 6.10.2003 */
    pos = name_list_pos("rmatrix", nl);
    if(nl->inform[pos]) /* rmatrix specified */
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
    /* END CHK-R-MATRIX; OB 6.10.2003 */

    /* START CHK-CHROM; RDM 22.9.2005 */
    pos = name_list_pos("chrom", nl);
    if(nl->inform[pos]) /* chrom specified */
    {
      cp = cmd->clone->par->parameters[pos];
      for (i = 0; i < match_num_seqs; i++)
      {
        /* START adding chrom to TWISS input command for each sequence */
        tnl = local_twiss[i]->cmd_def->par_names;
        tpos = name_list_pos("chrom", tnl);
        local_twiss[i]->cmd_def->par_names->inform[tpos] = 1;
        local_twiss[i]->cmd_def->par->parameters[tpos]->double_value
          = 1;
        set_option("twiss_chrom",&ione);
        /* END adding chrom to TWISS input command for each sequence */
      }
    }
    /* END CHK-CHROM; RDM 22.9.2005 */

    /* START CHK-SECTORMAP; RDM 16.12.2005 */
    pos = name_list_pos("sectormap", nl);
    if(nl->inform[pos]) /* sectormap specified */
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
    /* END CHK-CHROM; RDM 16.12.2005 */

    /* START CHK-DELTAP; RDM 12.12.2005 */
    pos = name_list_pos("deltap", nl);
    if(nl->inform[pos]) /* deltap specified */
    {
      cp = cmd->clone->par->parameters[pos];
      for (i = 0; i < match_num_seqs; i++)
      {
        /* START adding deltap to TWISS input command for each sequence */
        tnl = local_twiss[i]->cmd_def->par_names;
        tpos = name_list_pos("deltap", tnl);
        local_twiss[i]->cmd_def->par_names->inform[tpos] = 1;
        dpp=malloc(20);
        sprintf(dpp,"%e",cp->double_array->a[0]);
        local_twiss[i]->cmd_def->par->parameters[tpos]->string
          = dpp;
        /*       END adding deltap to TWISS input command for each sequence */
        /*      fprintf(prt_file, "entry value: %f\n", cp->double_array->a[0]);*/
        /*      twiss_deltas->curr=1;*/
        /*      twiss_deltas->a[0]=cp->double_array->a[0];*/
      }
    }
    /* END CHK-DELTAP; RDM 12.12.2005 */

    /* START CHK-ENTRIES of TYPE DOUBLE-REAL; OB 23.1.2002 */
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
    /* END CHK-ENTRIES of TYPE DOUBLE-REAL; OB 23.1.2002 */
    /* END CHK-BETA-INPUT; OB 1.2.2002 */

    /* START generating a TWISS table via 'pro_twiss'; OB 1.2.2002 */
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
}

void match_prepare_varypos()
  /* keeps constraints from nodes, reexpands, adds constraints to nodes */
{
  struct node* node = current_sequ->ex_start;
  struct constraint_list** tmplist = (struct constraint_list**)
    mymalloc("match_prepare_varypos",current_sequ->n_nodes * sizeof(struct constraint_list*));
  int i = 0;
  while (node != NULL)
  {
    tmplist[i] = current_sequ->all_nodes[i]->cl; i++;
    if (node == current_sequ->ex_end) break;
    node = node->next;
  }
  expand_curr_sequ(0);
  i = 0;
  node = current_sequ->ex_start;
  while (node != NULL)
  {
    current_sequ->all_nodes[i]->cl = tmplist[i]; i++;
    if (node == current_sequ->ex_end) break;
    node = node->next;
  }
  myfree("match_prepare_varypos", tmplist);
}

void match_rmatrix(struct in_cmd* cmd)
{
}

void match_tmatrix(struct in_cmd* cmd)
{
}

void match_vary(struct in_cmd* cmd)
{
  int pos;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  if (stored_match_var == NULL) stored_match_var = new_command_list("vary", 100);
  pos = name_list_pos("name", nl);
  if (nl->inform[pos])
  {
    cmd->clone_flag = 1;
    add_to_command_list(pl->parameters[pos]->string,
                        cmd->clone, stored_match_var, 1);
  }
}

void match_weight(struct in_cmd* cmd)
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

void mtjacprint(int m, int n,double* jac){
  int i,j,k,l;
  k=0;
  printf("Macro Constraint Variable Derivative\n");
  printf("------------------------------------\n");
  for(i=0;match2_macro_name[i]!=NULL;i++) {
    for(j=0;match2_cons_name[i][j]!=NULL;j++) {
      for(l=0;l<n;l++){
        printf("%10s ",match2_macro_name[i]);
        printf("%10s ",match2_cons_name[i][j]);
        printf("%10s ",command_par_string("name",stored_match_var->commands[l]));
        printf("%e",jac[l*m+k]);
        printf("\n");
      }
      k++;
    }
  }
}
