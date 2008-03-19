#include "fortran_wrappers.h"

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
 else if (strcmp(cmd->tok_list->p[0], "seterr") == 0)
          {
           error_seterr(cmd);
          }
 else if (strcmp(cmd->tok_list->p[0], "esave") == 0)
          {
           error_esave(cmd);
          }
}


void error_seterr(struct in_cmd* cmd)
{

/* read the errors from a named table  and stores
   them in the nodes of the sequence.       
   Subsequent Twiss will use them correctly.
   ===> Must be preceded by a call to "read_table"
   ===> (unless table exists in memory !)
*/

  int i, ix;
  int j;

  struct node *ndexe;
  struct node *nextnode;

  char     name[NAME_L];
  char   slname[NAME_L];

  char     nname[NAME_L];
  char   slnname[NAME_L];

  char*    namtab;
  int      t1;

  struct   table  *err;

/* set up pointers to current sequence for later use */
  struct sequence* mysequ = current_sequ;
  nextnode = mysequ->ex_start;
  ndexe = mysequ->ex_end;

/* printf("Pointers: %d %d %d\n",mysequ,nextnode,ndexe); */

  if ((namtab = command_par_string("table",cmd->clone)) != NULL) {
       printf("Want to use named table: %s\n",namtab);
       if ((t1 = name_list_pos(namtab, table_register->names)) > -1) {
          printf("The table ==> %s <=== was found \n",namtab);
       } else {
          /* fatal_error("Error table requested, but not existing:",namtab); */
          /* exit(-77); */ 
          printf("No such error table in memory: %s\n",namtab);
          /* try now, better a clean exit afterwards ... */
          exit(-77);
       }

  } else {
       if (get_option("debug")) {
         printf("No table name requested\n");
         printf("Use default name\n");
       }
       strcpy(namtab,"error");
       if ((t1 = name_list_pos(namtab, table_register->names)) > -1) {
          printf("The default table ==> %s <=== was found \n",namtab);
       } else {
          /* fatal_error("Error table requested, but not existing:",namtab); */
          /* exit(-77); */ 
          printf("No default error table in memory: %s\n",namtab);
          /* try now, better a clean exit afterwards ... */
          exit(-77);
       }
  }
 
  err = table_register->tables[t1];


  i = 1; /* watch out ! i is the ROW number, not the C index !!! */
  ix=0;
  while(ix == 0) {
      ix =   str_from_table(namtab, "name", &i, name);
      if(ix == 0) {
             stolower(name);
             strcpy(slname,strip(name));
             supp_tb(slname);
          nextnode = mysequ->ex_start;
          while (nextnode != ndexe) {
            
             strcpy(nname,nextnode->name);
             stolower(nname);
             strcpy(slnname,strip(nname));
             supp_tb(slnname);
          
/*
             printf("seq and input (0): %s %d %s %d\n", nname,strlen(nname),  name,strlen(name));
             printf("seq d in (2): %s %d %s %d\n",slnname,strlen(slnname),slname,strlen(slname)); 
*/
       
             if(strcmp(slname,slnname) == 0) {

/*              printf("O.K.:  %s in sequence and input table\n",slname); */
                /*
                ===>  now we have the match of the elements ..
                */

                /* We have now the input and the node, generate array and selection flag */
                nextnode->sel_err = 1;
                nextnode->p_fd_err = new_double_array(FIELD_MAX);
                nextnode->p_fd_err->curr = FIELD_MAX;
                nextnode->p_al_err = new_double_array(ALIGN_MAX);
                nextnode->p_al_err->curr = ALIGN_MAX;

                for (j = 1; j < EFIELD_TAB; j++) {
             /*   printf("efield errors: %d %e\n",j,err->d_cols[j][i-1]); */
                  nextnode->p_fd_err->a[j-1] = err->d_cols[j][i-1];
                }
                for (j = 1; j < err->num_cols-EFIELD_TAB; j++) {
             /*   printf("ealign errors: %d %e\n",j,err->d_cols[j+EFIELD_TAB][i-1]); */
                  nextnode->p_al_err->a[j-1] = err->d_cols[j+EFIELD_TAB][i-1];
                }
                

                nextnode = ndexe;
             } else {
                nextnode = nextnode->next;
             }
          }
    
      } 
        i++;
  }


return;
}


void error_esave(struct in_cmd* cmd)
{
    char *ef_table_file;
/*  if(efield_table == NULL) { */
       efield_table = make_table("efield", "efield", efield_table_cols,
                               efield_table_types, 10000);
       add_to_table_list(efield_table, table_register);
       pro_error_make_efield_table();
/*  }                          */
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
  int fll;

  ndexe = mysequ->ex_end;
  nextnode = mysequ->ex_start;

  mycount = 0;
 
  fll = command_par_value("full", cmd->clone);

  while (nextnode != ndexe) {

    if((nextnode->sel_err == 1) || (fll > 0) ){
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
       } else {
         fprintf(prt_file, "\n\nAlignment errors for element %s \n",nextnode->name);
         fprintf(prt_file,"\nDisplacements in [mm], rotations in [mrad] \n");
         fprintf(prt_file," %6s %10s %12s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n",
                  pln_alig[0], pln_alig[1],pln_alig[2],pln_alig[3],
                  pln_alig[4], pln_alig[5],pln_alig[6],pln_alig[7],
                  pln_alig[8], pln_alig[9],pln_alig[10],pln_alig[11],
                  pln_alig[12],pln_alig[13]);
         fprintf(prt_file, "\n");
       }

       if(nextnode->p_fd_err != NULL) {
         mycount++;
         if(mycount <= 50000) {
           /* fprintf(prt_file,"%s %d\n",nextnode->name,(int)nextnode->p_fd_err); */
           fprintf(prt_file, "\n\nField errors for element %s \n",nextnode->name);
           fprintf(prt_file, "Multipole order:     Normal:           Skew: \n");
           for(i=0;i<22;i++) {
              fprintf(prt_file, "%8d          %8e      %8e\n",i/2,
                     nextnode->p_fd_err->a[i],
                     nextnode->p_fd_err->a[i+1]);
              i++;
           }
           fprintf(prt_file, "\n");
         }
       } else {
         mycount++;
         if(mycount <= 50000) {
           fprintf(prt_file, "\n\nField errors for element %s \n",nextnode->name);
           fprintf(prt_file, "Multipole order:     Normal:           Skew: \n");
         }
       }
    }

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
  int    hyst = 0;
  int    flgmgt = 0;
  int chcount[3] = {0,0,0};
  char rout_name[] = "error_efcomp";
  double norfac; /* factor for normalization at reference radius */
  int    n = -1;     /* order of reference multipole */
  double rr = 0.0;   /* reference radius for multipole error */
  double rrr = 0.0;  /* reference radius for multipole error */
  struct double_array *ptr;
  struct double_array *pcoef;
  double h_co_n[FIELD_MAX/2][4];
  double h_co_s[FIELD_MAX/2][4];
  double *hco_n;
  double *hco_s;
  double *nvec;
  double deer;  
  double ref_str;
  double ref_strn;
  double ref_len;
  double nlength;
  double nvec0, nvec1, nvec2, nvec3;
  double val[3] = {0,0,0};
  static char atts[3][7] = {"order","radius","hyster"};
  static char attv[4][7] = {"dkn","dks","dknr","dksr"};
  int         iattv[4] = {0,0,0,0};
  struct sequence* mysequ = current_sequ;

  nvec = (double *)mycalloc("error_efcomp",1000, sizeof(double));

  ndexe = mysequ->ex_end;
  nextnode = mysequ->ex_start;

  nl = cmd->clone->par_names;

/* here comes a kludge, check which of the assignment vectors is there */
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
    hco_n = &h_co_n[0][0];
    for(j=0;j<FIELD_MAX*2;j++) {
       *hco_n = 0.0;           
       hco_n++;
    }
    hco_s = &h_co_s[0][0];
    for(j=0;j<FIELD_MAX*2;j++) {
       *hco_s = 0.0;           
       hco_s++;
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

/* now get order (n), radius (rr) and hyster flag (hyst) from command, if any */
          for(i=0;i<3;i++){
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
             } else if (i==2) {
               hyst = val[i];
               /* debug printout */
               if (get_option("debug"))
               fprintf(prt_file, "hyster flag is %d\n",(int)val[i]);
             }
          }

/* now get coefficients for time memory effects in magnets */
          pcoef = command_par_array("hcoeffn",cmd->clone);
          hco_n = &h_co_n[0][0];
          if(pcoef != NULL) {
             for(j=0;j<pcoef->curr;j++) {
                *hco_n = pcoef->a[j];
                hco_n++;
             }
             if (get_option("debug")) {
               for(j=0;j<FIELD_MAX/2;j++) {
                printf("COEFF: %d %e %e %e %e\n",j,h_co_n[j][0],h_co_n[j][1],h_co_n[j][2],h_co_n[j][3]);
               }
             }
          }

          pcoef = command_par_array("hcoeffs",cmd->clone);
          hco_s = &h_co_s[0][0];
          if(pcoef != NULL) {
             for(j=0;j<pcoef->curr;j++) {
                *hco_s = pcoef->a[j];
                hco_s++;
             }
             if (get_option("debug")) {
               for(j=0;j<FIELD_MAX/2;j++) {
                printf("COEFF: %d %e %e %e %e\n",j,h_co_s[j][0],h_co_s[j][1],h_co_s[j][2],h_co_s[j][3]);
               }
             }
          }


/* get length of node and check if magnet */
           ref_str  = 0.0;
           ref_strn = 0.0;
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
           ref_strn = fabs(ref_str);
         } else if (strcmp(nextnode->base_name,"sbend") == 0) {
           nvec0 = node_value("k0");
           if (get_option("debug")) {
              fprintf(prt_file, "original field0 is %f\n",nvec0);
              fprintf(prt_file, "====0====>>> %d %f %f \n\n",n,nvec0,nlength);
           }
           ref_str = nvec0*nlength;
           ref_strn = fabs(nvec0);
         } else if (strcmp(nextnode->base_name,"rbend") == 0) {
           nvec0 = node_value("k0");
           if (get_option("debug")) {
              fprintf(prt_file, "original field0 is %f\n",nvec0);
              fprintf(prt_file, "====0====>>> %d %f %f \n\n",n,nvec0,nlength);
           }
           ref_str = nvec0*nlength;
           ref_strn = fabs(nvec0);
         } else if (strcmp(nextnode->base_name,"quadrupole") == 0) {
           nvec1 = node_value("k1");
           if (get_option("debug")) {
              fprintf(prt_file, "original field1 is %f\n",nvec1);
              fprintf(prt_file, "====1====>>> %d %f %f \n\n",n,nvec1,nlength);
           }
           ref_str = nvec1*nlength;
           ref_strn = fabs(nvec1);
         } else if (strcmp(nextnode->base_name,"sextupole") == 0) {
           nvec2 = node_value("k2");
           if (get_option("debug")) {
              fprintf(prt_file, "original field2 is %f\n",nvec2);
              fprintf(prt_file, "====2====>>> %d %f %f \n\n",n,nvec2,nlength);
           }
           ref_str = nvec2*nlength;
           ref_strn = fabs(nvec2);
         } else if (strcmp(nextnode->base_name,"octupole") == 0) {
           nvec3 = node_value("k3");
           if (get_option("debug")) {
              fprintf(prt_file, "original field3 is %f\n",nvec3);
              fprintf(prt_file, "====3====>>> %d %f %f \n\n",n,nvec3,nlength);
           }
           ref_str = nvec3*nlength;
           ref_strn = fabs(nvec3);
         }

         /*  edbug print out field components , not done for production version
         */

         /* normal components -> 2j, skew components 2j+1 */
         if(flgmgt == 1) {

           for(i=0;i<4;i++)  {   /* loop over possible commands */

             if (get_option("debug")) fprintf(prt_file, "%s %d\n",attv[i],iattv[i]);

             ptr = command_par_array(attv[i],cmd->clone);
             if((ptr != NULL) && (iattv[i] == 1))  { /* command [i] found ? */


               for(j=0;j<ptr->curr;j++)  { /* loop over all parameters */

                 /* start field error assignment */
                 /* NORMAL COMPONENTS, ABSOLUTE ERRORS */
                 if(i==0)  {
                   if(add_error_opt == 1) {
                      nextnode->p_fd_err->a[2*j]   += ptr->a[j];
                   } else {
                      nextnode->p_fd_err->a[2*j]    = ptr->a[j];
                   }

                 /* SKEW COMPONENTS, ABSOLUTE ERRORS */
                 } else if(i==1) {
                   if(add_error_opt == 1) {
                   nextnode->p_fd_err->a[2*j+1] += ptr->a[j];
                   } else {
                   nextnode->p_fd_err->a[2*j+1]  = ptr->a[j];
                   }

                 /* NORMAL COMPONENTS, RELATIVE ERRORS, MAY BE CORRECTED FOR MEMORY EFFECTS */
                 } else if(i==2) {
                   if(fabs(rr) < 1.0E-6) {
                      printf("++++++ error: trying to assign relative field errors \n");
                      printf("       with no or zero reference radius specified\n");
                      exit(-1);
                   }
                   norfac = pow(rr,(n-j)) * (fact(j)/fact(n));

                   /* if flag for hysteresis correction is set, use coefficients for correction */
                   deer = 0.0;
                   if(hyst == 1) {
                      deer = h_co_n[j][3]*pow(ref_strn,3) + h_co_n[j][2]*pow(ref_strn,2) + 
                             h_co_n[j][1]*pow(ref_strn,1) + h_co_n[j][0];
                      if (get_option("debug")) 
                      printf("after correction (n): %d %e %e %e %e\n",
                              j,ref_strn,ptr->a[j],deer,(ptr->a[j] + deer));
                   }
/*
                   if (get_option("debug"))
                   fprintf(prt_file, "norm(n): %d %d %f %f\n",n,j,rr,norfac);
*/
                   if(add_error_opt == 1) {
                   nextnode->p_fd_err->a[2*j]   += (ptr->a[j]+deer)*ref_str*norfac;
                   } else {
                   nextnode->p_fd_err->a[2*j]    = (ptr->a[j]+deer)*ref_str*norfac;
                   }

                 /* SKEW COMPONENTS, RELATIVE ERRORS, MAY BE CORRECTED FOR MEMORY EFFECTS */
                 } else if(i==3) {
                   if(fabs(rr) < 1.0E-6) {
                      printf("++++++ error: trying to assign relative field errors \n");
                      printf("       with no or zero reference radius specified\n");
                      exit(-1);
                   }
                   norfac = pow(rr,(n-j)) * (fact(j)/fact(n));

                   /* if flag for hysteresis correction is set, use coefficients for correction */
                   deer = 0.0;
                   if(hyst == 1) {
                      deer = h_co_s[j][3]*pow(ref_strn,3) + h_co_s[j][2]*pow(ref_strn,2) + 
                             h_co_s[j][1]*pow(ref_strn,1) + h_co_s[j][0];
                      if (get_option("debug")) 
                      printf("after correction (s): %d %e %e %e %e\n",
                              j,ref_strn,ptr->a[j],deer,(ptr->a[j] + deer));
                   }
/*
                   if (get_option("debug"))
                   fprintf(prt_file, "norm(s): %d %d %f %f\n",n,j,rr,norfac);
*/
                   if(add_error_opt == 1) {
                   nextnode->p_fd_err->a[2*j+1] += (ptr->a[j]+deer)*ref_str*norfac;
                   } else {
                   nextnode->p_fd_err->a[2*j+1]  = (ptr->a[j]+deer)*ref_str*norfac;
                   }
                 } /*  end  of field error assignment */

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
  int is, ia;
  static  int  ia_seen = 0;

  is = 0; ia = 0;
  
  i = 0;
  while(cmd->tok_list->p[i] != NULL) {
     if(strcmp("add",cmd->tok_list->p[i]) == 0) {
          ia = 1;
     }
     if(strcmp("seed",cmd->tok_list->p[i]) == 0) {
          is = 1;
     }
     i++;
  }
  if ((debug=get_option("debug"))) printf("FOUND: %d %d \n",ia,is);

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

  /* change only if present in command or not yet set */
  if ((ia == 1) || (ia_seen != 1))
  {
       val = command_par_value("add", cmd->clone);
       if(val == 0) {
         if (debug)  fprintf(prt_file, "add option not set\n");
         add_error_opt = 0;
       }else {
         if (debug)  fprintf(prt_file, "add option set\n");
         add_error_opt = 1;
       }
  }


  if(ia == 1) ia_seen = 1;

  if ((debug=get_option("debug"))) printf("err_add eoption: %d seen: %d\n",add_error_opt,ia_seen);

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

