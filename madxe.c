void pro_error(struct in_cmd* cmd)
{

 if (strcmp(cmd->tok_list->p[0], "eoption") == 0)
     {
      error_eoption(cmd);
      cmd->clone_flag = 1; /* do not drop */
      current_eopt = cmd->clone;
      return;
     }
  if (get_option("debug")) fprintf(prt_file, "enter ERROR module\n");
  if (current_sequ == NULL || current_sequ->ex_start == NULL)
    {
     warning("ERROR, but no active sequence:", "ignored");
     return;
    }
        setbuf(stdout,(char *)0);

 if (error_select->curr > 0) set_selected_errors();

 if (strcmp(cmd->tok_list->p[0], "ealign") == 0)
          {
           error_ealign(cmd);
          }
 else if (strcmp(cmd->tok_list->p[0], "efield") == 0)
          {
           error_efield(cmd);
          }
 else if (strcmp(cmd->tok_list->p[0], "efcomp") == 0)
          {
           error_efcomp(cmd);
          }
 else if (strcmp(cmd->tok_list->p[0], "eprint") == 0)
          {
           error_eprint(cmd);
          }
 else if (strcmp(cmd->tok_list->p[0], "esave") == 0)
          {
           error_esave(cmd);
          }
}

void error_esave(struct in_cmd* cmd)
{
    char *ef_table_file;
//  if(efield_table == NULL) {
       efield_table = make_table("efield", "efield", efield_table_cols,
                               efield_table_types, 10000);
       add_to_table_list(efield_table, table_register);
       pro_error_make_efield_table();
//  }
    ef_table_file = command_par_string("file",cmd->clone);
    out_table("efield",efield_table,ef_table_file);
}

void error_ealign(struct in_cmd* cmd)
{
  struct node *ndexe;
  struct node *nextnode;
  int i;
  int chcount[3] = {0,0,0};
  double val[ALIGN_MAX] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  static char att[ALIGN_MAX][7] = {"dx","dy","ds","dphi","dtheta","dpsi","mrex","mrey","mredx","mredy","arex","arey","mscalx","mscaly"};

  struct command_parameter_list* pl = current_error->par;
  struct sequence* mysequ = current_sequ;

  nextnode = mysequ->ex_start;
  ndexe = mysequ->ex_end;

  while (nextnode != ndexe) {

      if(nextnode->sel_err == 1) {
        if(nextnode->p_al_err == NULL) {
          chcount[0]++;
          nextnode->p_al_err = new_double_array(ALIGN_MAX);
          nextnode->p_al_err->curr = ALIGN_MAX;
        } else {
          if(add_error_opt == 1) {
              chcount[2]++;
          } else {
              chcount[1]++;
          }
        }
          for(i=0;i<ALIGN_MAX;i++){
               val[i] = pl->parameters[i]->double_value;
               if(pl->parameters[i]->expr != NULL) {
                  if(pl->parameters[i]->expr->status != 1) {
                      val[i] = command_par_value(att[i], cmd->clone);
                      pl->parameters[i]->expr->status = 0;
                  }
               }
               if(add_error_opt == 1) {
                  nextnode->p_al_err->a[i] += val[i];
               } else {
                  nextnode->p_al_err->a[i] = val[i];
               }
          }
      }  /* end of treatment of selected node */
      nextnode = nextnode->next;
  }  /* end of loop over all nodes */
  if(chcount[0] != 0)
  fprintf(prt_file, "Assigned alignment errors to %d elements\n",chcount[0]);
  if(chcount[1] != 0)
  fprintf(prt_file, "Replaced alignment errors for %d elements\n",chcount[1]);
  if(chcount[2] != 0)
  fprintf(prt_file, "Added alignment errors to %d elements\n",chcount[2]);

}

void error_eprint(struct in_cmd* cmd)
{
  struct node *ndexe;
  struct node *nextnode;
  static char pln_alig[ALIGN_MAX][7] = {"dx","dy","ds","dphi","dtheta","dpsi","mrex","mrey","mredx","mredy","arex","arey","mscalx","mscaly"};
  static float alig_fact[ALIGN_MAX] = {1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1.0,1.0};

  int i;
  struct sequence* mysequ = current_sequ;
  int mycount;

  ndexe = mysequ->ex_end;
  nextnode = mysequ->ex_start;

  mycount = 0;

  while (nextnode != ndexe) {

 //   if(nextnode->sel_err == 1) {
       if(nextnode->p_al_err != NULL) {
         fprintf(prt_file, "\n\nAlignment errors for element %s \n",nextnode->name);
         fprintf(prt_file,"\nDisplacements in [mm], rotations in [mrad] \n");
         fprintf(prt_file," %6s %10s %12s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n",
                  pln_alig[0], pln_alig[1],pln_alig[2],pln_alig[3],
                  pln_alig[4], pln_alig[5],pln_alig[6],pln_alig[7],
                  pln_alig[8], pln_alig[9],pln_alig[10],pln_alig[11],
                  pln_alig[12],pln_alig[13]);
         for(i=0;i<nextnode->p_al_err->curr;i++) {
            fprintf(prt_file, "%10.6f ",alig_fact[i]*nextnode->p_al_err->a[i]);
         }
         fprintf(prt_file, "\n");
       }

       if(nextnode->p_fd_err != NULL) {
         mycount++;
         if(mycount <= 50000) {
         fprintf(prt_file,"%s %d\n",nextnode->name,(int)nextnode->p_fd_err);
         fprintf(prt_file, "\n\nField errors for element %s \n",nextnode->name);
         fprintf(prt_file, "Multipole order:     Normal:           Skew: \n");
  /*     for(i=0;i<nextnode->p_fd_err->curr;i++) { */
         for(i=0;i<22;i++) {
            fprintf(prt_file, "%8d          %8e      %8e\n",i/2,
                   1000*nextnode->p_fd_err->a[i],
                   1000*nextnode->p_fd_err->a[i+1]);
            i++;
         }
         fprintf(prt_file, "\n");
         }
       }
//    }

      nextnode = nextnode->next;
  }

}


void error_efcomp(struct in_cmd* cmd)
{
  struct name_list* nl;
  int pos;

  struct node *ndexe;
  struct node *nextnode;
  int i,j,k;
  int    lvec;
  int    flgmgt = 0;
  int chcount[3] = {0,0,0};
  char rout_name[] = "error_efcomp";
  double norfac; /* factor for normalization at reference radius */
  int    n;     /* order of reference multipole */
  double rr, rrr;    /* reference radius for multipole error */
  struct double_array *ptr;
  double *nvec;
  double ref_str;
  double ref_len;
  double nlength;
  double nvec0, nvec1, nvec2, nvec3;
  double val[2] = {0,0};
  static char atts[2][7] = {"order","radius"};
  static char attv[4][7] = {"dkn","dks","dknr","dksr"};
  int         iattv[4] = {0,0,0,0};
  struct sequence* mysequ = current_sequ;

  nvec = (double *)mycalloc("error_efcomp",1000, sizeof(double));

  ndexe = mysequ->ex_end;
  nextnode = mysequ->ex_start;

  nl = cmd->clone->par_names;

// here comes a kludge, check which of the assignment vectors is there
/*
    i = 0;
    while((cmd->tok_list->p[i]) != NULL) {
      for(k=0;k<4;k++) {
          if(strcmp(cmd->tok_list->p[i],attv[k]) == 0) {
             iattv[k] = 1;
          }
      }
    i++;
    }
*/
    for(k=0; k<4; k++) {
         pos = name_list_pos(attv[k],nl);
         if(nl->inform[pos] > 0) {
             if (get_option("debug"))
             fprintf(prt_file, "set iattv %d for %s to 1\n",iattv[k],attv[k]);
             iattv[k] = 1;
         }
    }

  while (nextnode != ndexe) { /*loop over elements and get strengths in vector*/
           current_node = nextnode;
           flgmgt = (int)node_value("magnet");
        if((nextnode->sel_err == 1) && (flgmgt == 1))  {
            if(nextnode->p_fd_err == NULL) {
            chcount[0]++;
            nextnode->p_fd_err = new_double_array(FIELD_MAX);
            nextnode->p_fd_err->curr = FIELD_MAX;
          } else {
            if(add_error_opt == 1) {
              chcount[2]++;
            } else {
              chcount[1]++;
            }
          }
          if (get_option("debug"))
          fprintf(prt_file, "field for %s %s %d\n",
                 nextnode->name,nextnode->base_name,nextnode->sel_err);

/* now get order (n) and radius (rr) from command, if any */
          for(i=0;i<2;i++){
             val[i] = command_par_value(atts[i],cmd->clone);
             if(i==0) {
               n = val[i];
               /* debug printout */
               if (get_option("debug"))
               fprintf(prt_file, "order  is %d\n",n);
             } else if (i==1) {
               rrr = val[i];
               rr  = fabs(rrr);
               /* debug printout */
               if (get_option("debug"))
               fprintf(prt_file, "radius is %f\n",val[i]);
             }
          }
/* get length of node and check if magnet */
           ref_str = 0.0;
           nlength = node_value("l");
           ref_len = nlength;
           if (get_option("debug"))
           fprintf(prt_file, "original length is %f\n",nlength);

         if(strcmp(nextnode->base_name,"multipole") == 0) {
           if(rrr > 0 ) {
              get_node_vector("knl",&lvec,nvec);
           } else {
              get_node_vector("ksl",&lvec,nvec);
           }
           if (get_option("debug")) {
             for(i=0;i<4;i++) {
                fprintf(prt_file, "original field = %d is %f\n",i,nvec[i]);
             }
           }
           if (get_option("debug"))
           fprintf(prt_file, "====n====>>> %d %f %f \n\n",n,nvec[n],nlength);
           ref_str = nvec[n];
         } else if (strcmp(nextnode->base_name,"sbend") == 0) {
           nvec0 = node_value("k0");
           if (get_option("debug")) {
              fprintf(prt_file, "original field0 is %f\n",nvec0);
              fprintf(prt_file, "====0====>>> %d %f %f \n\n",n,nvec0,nlength);
           }
           ref_str = nvec0*nlength;
         } else if (strcmp(nextnode->base_name,"rbend") == 0) {
           nvec0 = node_value("k0");
           if (get_option("debug")) {
              fprintf(prt_file, "original field0 is %f\n",nvec0);
              fprintf(prt_file, "====0====>>> %d %f %f \n\n",n,nvec0,nlength);
           }
           ref_str = nvec0*nlength;
         } else if (strcmp(nextnode->base_name,"quadrupole") == 0) {
           nvec1 = node_value("k1");
           if (get_option("debug")) {
              fprintf(prt_file, "original field1 is %f\n",nvec1);
              fprintf(prt_file, "====1====>>> %d %f %f \n\n",n,nvec1,nlength);
           }
           ref_str = nvec1*nlength;
         } else if (strcmp(nextnode->base_name,"sextupole") == 0) {
           nvec2 = node_value("k2");
           if (get_option("debug")) {
              fprintf(prt_file, "original field2 is %f\n",nvec2);
              fprintf(prt_file, "====2====>>> %d %f %f \n\n",n,nvec2,nlength);
           }
           ref_str = nvec2*nlength;
         } else if (strcmp(nextnode->base_name,"octupole") == 0) {
           nvec3 = node_value("k3");
           if (get_option("debug")) {
              fprintf(prt_file, "original field3 is %f\n",nvec3);
              fprintf(prt_file, "====3====>>> %d %f %f \n\n",n,nvec3,nlength);
           }
           ref_str = nvec3*nlength;
         }

         /*  edbug print out field components , not done for production version
         */

         /* normal components -> 2j, skew components 2j+1 */
         if(flgmgt == 1) {
           for(i=0;i<4;i++)  {
             if (get_option("debug")) fprintf(prt_file, "%s %d\n",attv[i],iattv[i]);
             ptr = command_par_array(attv[i],cmd->clone);
             if((ptr != NULL) && (iattv[i] == 1)) {
               for(j=0;j<ptr->curr;j++){
                 if(i==0)  {
                   if(add_error_opt == 1) {
                      nextnode->p_fd_err->a[2*j]   += ptr->a[j];
                   } else {
                      nextnode->p_fd_err->a[2*j]    = ptr->a[j];
                   }
                 } else if(i==1) {
                   if(add_error_opt == 1) {
                   nextnode->p_fd_err->a[2*j+1] += ptr->a[j];
                   } else {
                   nextnode->p_fd_err->a[2*j+1]  = ptr->a[j];
                   }
                 } else if(i==2) {
                   if(fabs(rr) < 1.0E-6) {
                      printf("++++++ error: trying to assign relative field errors \n");
                      printf("       with no or zero reference radius specified\n");
                      exit(-1);
                   }
                   norfac = pow(rr,(n-j)) * (fact(j)/fact(n));
/*
                   if (get_option("debug"))
                   fprintf(prt_file, "norm(n): %d %d %f %f\n",n,j,rr,norfac);
*/
                   if(add_error_opt == 1) {
                   nextnode->p_fd_err->a[2*j]   += ptr->a[j]*ref_str*norfac;
                   } else {
                   nextnode->p_fd_err->a[2*j]    = ptr->a[j]*ref_str*norfac;
                   }
                 } else if(i==3) {
                   if(fabs(rr) < 1.0E-6) {
                      printf("++++++ error: trying to assign relative field errors \n");
                      printf("       with no or zero reference radius specified\n");
                      exit(-1);
                   }
                   norfac = pow(rr,(n-j)) * (fact(j)/fact(n));
/*
                   if (get_option("debug"))
                   fprintf(prt_file, "norm(s): %d %d %f %f\n",n,j,rr,norfac);
*/
                   if(add_error_opt == 1) {
                   nextnode->p_fd_err->a[2*j+1] += ptr->a[j]*ref_str*norfac;
                   } else {
                   nextnode->p_fd_err->a[2*j+1]  = ptr->a[j]*ref_str*norfac;
                   }
                 }

               }
            }
           }
         }
      }  /* end of treatment of selected node */
      nextnode = nextnode->next;
  }      /* end of loop over all nodes */
  if(chcount[0] != 0)
  fprintf(prt_file, "Assigned field errors to %d elements\n",chcount[0]);
  if(chcount[1] != 0)
  fprintf(prt_file, "Replaced field errors for %d elements\n",chcount[1]);
  if(chcount[2] != 0)
  fprintf(prt_file, "Added field errors to %d elements\n",chcount[2]);
  myfree(rout_name,nvec);

}

void error_efield(struct in_cmd* cmd)
{
  int i;
  /*
  struct node *ndexs, *ndexe;
  struct node *nextnode;
  struct name_list* nl = current_error->par_names;
  struct command_parameter_list* pl = current_error->par;
  struct sequence* mysequ = current_sequ;
  */

  fprintf(prt_file, "in efield routine\n");
  fprintf(prt_file, "efield command not yet implemented\n");

  if (get_option("debug"))
  {
   for(i=0;i<cmd->tok_list->curr;i++)
     fprintf(prt_file, "command(s): %s\n",cmd->tok_list->p[i]);
  }
}

void error_eoption(struct in_cmd* cmd)
{
  struct name_list* nl = cmd->clone->par_names;
  int i, debug;
  int val, pos, seed;

  if ((debug=get_option("debug")))  {
     fprintf(prt_file, "in eoption routine\n");
     for(i=0;i<cmd->tok_list->curr;i++) {
        fprintf(prt_file, "command(s): %s\n",cmd->tok_list->p[i]);
     }
  }
  if ((pos = name_list_pos("seed", nl)) > -1)
    {
     if (nl->inform[pos])
       {
        seed = command_par_value("seed", cmd->clone);
        init55(seed);
       }
    }
 val = command_par_value("add", cmd->clone);
 if(val == 0) {
      if (debug)  fprintf(prt_file, "add option not set\n");
      add_error_opt = 0;
 }else {
      if (debug)  fprintf(prt_file, "add option set\n");
      add_error_opt = 1;
 }

}

double fact(int i)
{
  int k;
  double nfact;

  if(i == 0) {
     return(1.0);
  } else if(i == 1) {
     return(1.0);
  } else if(i > 1) {
    nfact = 1.;
    for (k=1;k<=i;k++) {
        nfact = nfact*k;
    }
     return(nfact);
  } else if(i < 0){
     return(-1.0);
  } else {
     return(-1.0);
  }
}

void pro_error_make_efield_table()
{
  struct table *ttb = efield_table;
  struct node *nanf;
  struct node *nend;
  int         j;
  struct sequence* mysequ = current_sequ;

  setbuf(stdout,(char *)0);
  nanf = mysequ->ex_start;
  nend = mysequ->ex_end;

      while (nanf != nend) {
        if(nanf->sel_err == 1) {
           string_to_table("efield","name",nanf->name);
/* */
                  /*
           printf("=> %s %e %e %e\n",nanf->name,nanf->p_fd_err,nanf->p_al_err);
                  */
           if(nanf->p_fd_err != NULL) {
              for (j=1; j <= EFIELD_TAB; j++) {
                 ttb->d_cols[j][ttb->curr] = nanf->p_fd_err->a[j-1];
                  /*
                 printf("Field: %d %e\n",j,ttb->d_cols[j][ttb->curr]);
                  */
              }
           }
           if(nanf->p_al_err != NULL) {
              for (j=1; j < ttb->num_cols-EFIELD_TAB; j++) {
                 ttb->d_cols[j+EFIELD_TAB][ttb->curr] = nanf->p_al_err->a[j-1];
                  /*
                 printf("Align: %d %e\n",j,ttb->d_cols[j][ttb->curr]);
                  */
              }
           }
/* */
           augment_count("efield");
        }
        nanf = nanf->next;
      }
}

