void pro_correct(struct in_cmd* cmd)
{
  if (strcmp(cmd->tok_list->p[0], "correct") == 0)
    {
     correct_correct(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "usekick") == 0)
    {
     correct_usekick(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "usemonitor") == 0)
    {
     correct_usemonitor(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "getorbit") == 0)
    {
     correct_getorbit(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "putorbit") == 0)
    {
     correct_putorbit(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "coption") == 0)
    {
     correct_option(cmd);
    }
}
void correct_correct(struct in_cmd* cmd)
{
  int ix, im, ip, it, idrop;
  int i,j,err,nnnseq;
  int imon, icor;
  int ncorr, nmon;
  int niter;
  int resout;
  int twism;
  struct timeb tp;
  int ifail, sflag, svdflg;
  float  rms;
  double tmp1, tmp2, tmp3, tmp4;
  double  sigcut;            /* number of sigmas (normalized) for filter cut */
  char    *clist, *mlist;    /* file names for monitor and corrector output */
  double  *dmat = {NULL};    /* response matrix, double precision */
  double  *corvec, *monvec;  /* vectors to hold measured orbit and correctors */
  double  *resvec;           /* vector to hold corrected orbit */
  char    *conm;             /* vector to hold corrector names (for MICADO) */
  int     *sing;             /* array to store pointer to singular correctors */
  static int     *nm, *nx, *nc;
  struct id_mic  *corl;
  ip = pro_correct_getcommands(cmd);
  im = pro_correct_gettables(ip);
  ncorr = im%10000; nmon  = im/10000;
  printf("%d monitors and %d correctors found in input\n",nmon,ncorr);
  if(nmon == 0) {
    printf("No monitor found in input, no correction done\n");
    return;
  }
  if(ncorr == 0) {
    printf("No corrector found in input, no correction done\n");
    return;
  }
  /* For debugging set output buffer to zero */
  if (get_option("debug"))  setbuf(stdout,NULL);
  /* Prepare file descriptors for the output */
  if((resout = command_par_value("resout",cmd->clone)) > 0) {
     if(fddata == NULL) {
        if((fddata = fopen("corr.out","w")) == NULL)
           exit(99);
     }
     if(fcdata == NULL) {
        if((fcdata = fopen("stren.out","w")) == NULL)
           exit(99);
     }
  }
  /* If only Twiss summary is required prepare and write it */
  if((twism = command_par_value("twissum",cmd->clone)) > 0) {
     if(ftdata == NULL) {
        if((ftdata = fopen("twiss.summ","w")) == NULL)
           exit(99);
     }
     j = 1;
     if((nnnseq = get_variable("n")) == 0) {
         nnnseq = twism;
     }
     err = double_from_table("summ", "xcomax",&j,&tmp1);
     err = double_from_table("summ", "xcorms",&j,&tmp2);
     err = double_from_table("summ", "ycomax",&j,&tmp3);
     err = double_from_table("summ", "ycorms",&j,&tmp4);
     fprintf(ftdata," T: %d %e %e %e %e\n",nnnseq,tmp1,tmp2,tmp3,tmp4);
     return;
  }
  /* allocate vectors used by correction algorithms */
  nx  = (int *)mycalloc("correct_correct_nx",ncorr,sizeof(int));
  nc  = (int *)mycalloc("correct_correct_nc",ncorr,sizeof(int));
  nm  = (int *)mycalloc("correct_correct_nm",nmon,sizeof(int));
  sing = (int *)mycalloc("correct_correct_sing",ncorr*2,sizeof(int));
  corvec = (double *)mycalloc("correct_correct_corvec",ncorr,sizeof(double));
  monvec = (double *)mycalloc("correct_correct_monvec",nmon,sizeof(double));
  resvec = (double *)mycalloc("correct_correct_resvec",nmon,sizeof(double));
  conm = (char *)mycalloc("correct_correct_conm",ncorr*16,sizeof(char));
  /* get original settings of correctors from input Twiss-table */
  it = pro_correct_getcorrs(cmd);
  /* get input orbit, default is from input Twiss-table */
  it = pro_correct_getorbit(cmd);
  /* find and prepare enabled correctors and monitors, may be repeated */
  ix = pro_correct_getactive(ip, nm, nx, nc, corvec, monvec, conm);
  icor = ix%10000; imon  = ix/10000;
  printf("%d monitors and %d correctors enabled\n",imon,icor);
  /* normalized cut on beam position, if requested */
  if((sigcut = command_par_value("moncut",cmd->clone)) > 0) {
     idrop = pro_correct_filter(ip,sigcut);
     printf("Disabled %d monitors with %-2.2f sigma cut\n",idrop,sigcut);
     ix = pro_correct_getactive(ip, nm, nx, nc, corvec, monvec, conm);
     icor = ix%10000; imon  = ix/10000;
     printf("After filter of %-2.2f sigma:\n",sigcut);
     printf("%d monitors and %d correctors enabled\n",imon,icor);
  }
  /* set up response matrix for ring or line */
  corl = correct_orbit->cor_table;
  if(strcmp("ring",command_par_string("flag",cmd->clone)) == 0) {
    if(dmat != NULL) free(dmat);
    /* icor and imon used to set up correct matrix size !! */
    dmat = (double *)pro_correct_response_ring(ip,icor,imon);
    if((svdflg = command_par_value("cond",cmd->clone)) == 1) {
       printf("SVD conditioning requested ...\n");
       /* printf("Time before svd-comd:  %-6.3f\n",fextim());    */
       sflag=c_svddec(dmat,imon,icor,sing);
       /* printf("Time after svd-cond:  %-6.3f\n",fextim());     */
       /* printf("sflag: %d\n",sflag); */
       for(ix=0;ix<sflag;ix++) {
         corl[nx[sing[2*ix+0]]].enable = 0;
         printf("Removed:   %d %s\n",nx[sing[2*ix+0]],
                 corl[nx[sing[2*ix+0]]].p_node->name);
       }
       ix = pro_correct_getactive(ip, nm, nx, nc, corvec, monvec, conm);
       icor = ix%10000; imon  = ix/10000;
       printf("After SVD conditioning:             \n");
       printf("%d monitors and %d correctors enabled\n\n",imon,icor);
       if(dmat != NULL) free(dmat);
       /* icor and imon used to set up correct matrix size !! */
       dmat = (double *)pro_correct_response_ring(ip,icor,imon);
       sflag=c_svddec(dmat,imon,icor,sing);
       printf("Found %d singular values\n",sflag);
    }
  }
  else if(strcmp("line",command_par_string("flag",cmd->clone)) == 0) {
    if(dmat != NULL) free(dmat);
    dmat = (double *)pro_correct_response_line(ip,ncorr,nmon); }
  else { printf("INVALID MACHINE TYPE\n"); exit(-1);
  }
  if (get_option("debug")) {
    pro_correct_prtwiss();
    pro_correct_write_cocu_table();
  }
  /* LSQ correction, use all available correctors */
  if(strcmp("lsq",command_par_string("mode",cmd->clone)) == 0) {
    /*frs haveit_(dmat,monvec,corvec,resvec,nx,&imon,&icor); */
    /* printf("Time before lsq:  %-6.3f\n",fextim());   */
    c_haveit(dmat,monvec,corvec,resvec,nx,imon,icor);
    /* printf("Time after lsq:  %-6.3f\n",fextim());    */
    pro_correct_write_results(monvec, resvec, corvec, nx, nc, nm, imon, icor, ip);
  }
  /* SVD correction, use all available correctors */
  if(strcmp("svd",command_par_string("mode",cmd->clone)) == 0) {
    /*frs haveit_(dmat,monvec,corvec,resvec,nx,&imon,&icor); */
    /* printf("Time before svd-corr:  %-6.3f\n",fextim());   */
       sflag=c_svdcorr(dmat,monvec,corvec,resvec,nx,imon,icor);
    /* printf("Time after svd-corr:  %-6.3f\n",fextim());    */
    pro_correct_write_results(monvec, resvec, corvec, nx, nc, nm, imon, icor, ip);
  }
  /* MICADO correction, get desired number of correctors from command */
  corrl = command_par_value("corrlim",cmd->clone);
  set_variable("corrlim",&corrl);
  if(strcmp("micado",command_par_string("mode",cmd->clone)) == 0) {
    printf("enter MICADO correction ...\n");
    if((niter = command_par_value("ncorr",cmd->clone)) == 0) {
          printf("Requested %d correctors (???) set to %d\n",niter,icor);
          niter = icor;
    }
    else if((niter = command_par_value("ncorr",cmd->clone)) < 0) {
          printf("Requested %d correctors (???) set to 0\n",niter);
          niter = 0;
    }
    else if((niter = command_par_value("ncorr",cmd->clone)) > icor) {
          printf("Fewer correctors available than requested by ncorr\n");
          printf("you want %d,  you get %d\n",niter,icor);
          printf("ncorr reset to %d\n",icor);
          niter = icor;
    }
    rms  = 1000.0*command_par_value("error",cmd->clone);
    /*frs       micit_(dmat,monvec,corvec,resvec,nx,&rms,&imon,&icor,&niter); */
    /* printf("Time before micado:  %-6.3f\n",fextim());  */
    ifail = c_micit(dmat,conm,monvec,corvec,resvec,nx,rms,imon,icor,niter);
    /* printf("Time after micado:  %-6.3f\n",fextim());   */
    if(ifail != 0) {
       printf("MICADO correction completed with error code %d\n\n",ifail);
    }
    pro_correct_write_results(monvec, resvec, corvec, nx, nc, nm, imon, icor, ip);
  }
  /* write corrector output to tfs table */
  if ((clist = command_par_string("clist",cmd->clone)) != NULL) {
    out_table("corr",corr_table,clist);
  }
  /* write monitor output to tfs table */
  if ((mlist = command_par_string("mlist",cmd->clone)) != NULL) {
    out_table("mon",mon_table,mlist);
  }
  /* Clean up at the end of the module */
  free(nm);free(dmat);free(nx);free(nc);free(corvec);
  free(monvec);free(resvec); free(conm);
  return;
}
void pro_correct_option(struct in_cmd* cmd)
{
  struct name_list* nl = cmd->clone->par_names;
  int i, debug;
  int val, pos, seed;
  if ((debug=get_option("debug")))  {
     fprintf(prt_file, "in coption routine\n");
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
 val = command_par_value("print", cmd->clone);
 if(val == 0) {
      if (debug)  fprintf(prt_file, "print option not set\n");
      print_correct_opt = 0;
 }else {
      if (debug)  fprintf(prt_file, "print option set\n");
      print_correct_opt = val;
 }
}
int pro_correct_getcommands(struct in_cmd* cmd)
{
  static char att[10][8] = {"iterate","plane","ncorr","error","clist",
                            "mlist", "flag","mode","",""};
  static int iplane = 1;
  char plane[20];

  if (get_option("debug")) printf("enter CORRECT module\n");
  if (current_sequ == NULL || current_sequ->ex_start == NULL)
    {
      warning("CORRECT, but no active sequence:", "ignored");
      return(-1);
    }
  strcpy(plane,command_par_string(att[1],cmd->clone));
  if(strcmp("x",plane) == 0) {
    iplane = 1;
  } else if (strcmp("y",plane) == 0) {
    iplane = 2;
  } else if (strcmp("h",plane) == 0) {
    iplane = 1;
  } else if (strcmp("v",plane) == 0) {
    iplane = 2;
  } else {
    printf("No valid plane specified, x plane used \n");
    iplane = 1;
  }
  return(iplane);
}
int  pro_correct_gettables(int iplane)
{
  struct id_mic *cor_l;
  struct id_mic *mon_l;
  struct table *ttb;
  int j;
  int cntm = {0};
  int cntc = {0};
  static char atm[6][4] = {"hmon","vmon","moni","hkic","vkic","kick"};
  if(correct_orbit == NULL) {
    correct_orbit = (struct orb_cor*)mycalloc("pro_correct_gettables",1, sizeof(struct orb_cor));
  }
  if(corr_table == NULL) {
    corr_table = make_table("corr", "corr", corr_table_cols,
            corr_table_types, 5000);
    add_to_table_list(corr_table, table_register);
    pro_correct_make_corr_table();
  }
  if(mon_table == NULL) {
    mon_table = make_table("mon", "mon", mon_table_cols,
           mon_table_types, 5000);
    add_to_table_list(mon_table, table_register);
    pro_correct_make_mon_table();
  }

  if(correct_orbit->cor_table != NULL) free(correct_orbit->cor_table);
  if(correct_orbit->mon_table != NULL) free(correct_orbit->mon_table);
  correct_orbit->cor_table = (struct id_mic *)mycalloc("pro_correct_gettables_cor",5200, sizeof(struct id_mic));
  correct_orbit->mon_table = (struct id_mic *)mycalloc("pro_correct_gettables_mon",5200, sizeof(struct id_mic));
  ttb = twiss_table;
  correct_orbit->mon_table->previous = NULL;
  correct_orbit->mon_table->next = NULL;
  correct_orbit->cor_table->previous = NULL;
  correct_orbit->cor_table->next = NULL;
  mon_l = correct_orbit->mon_table;
  cor_l = correct_orbit->cor_table;
  for (j=0; j < ttb->curr; j++) {
   if((strncmp(atm[iplane-1],ttb->p_nodes[j]->base_name,4) == 0) ||
      (strncmp(atm[2],      ttb->p_nodes[j]->base_name,4) == 0))  {
      mon_l->id_ttb = j;
      mon_l->enable = ttb->p_nodes[j]->enable;
      mon_l->p_node = ttb->p_nodes[j];
      mon_l->next = mon_l;
      mon_l->next++; mon_l++;
      cntm++;
    }
    if((strncmp(atm[iplane+2],ttb->p_nodes[j]->base_name,4) == 0) ||
       (strncmp(atm[5],       ttb->p_nodes[j]->base_name,4) == 0))  {
      cor_l->id_ttb = j;
      cor_l->enable = ttb->p_nodes[j]->enable;
      cor_l->p_node = ttb->p_nodes[j];
      cor_l->next = cor_l;
      cor_l->next++; cor_l++;
      cntc++;
    }
  }
  mon_l--; mon_l->next = NULL;
  cor_l--; cor_l->next = NULL;
  return(10000*cntm + cntc);
}
int pro_correct_getorbit(struct in_cmd* cmd)
{
  struct name_list* nl;
  int i;
  int pos;
  struct id_mic *m;  /* access to tables for monitors and correctors */
  struct table *ttb;
  double **da1;
  double xlimit;
  ttb = twiss_table;
  da1 = ttb->d_cols;
  nl = cmd->clone->par_names;
  m = correct_orbit->mon_table;
  while(m) {
    /*    printf("m-list: %d %s %s\n",m->id_ttb,m->p_node->name,m->p_node->base_name); */
    m->val.before[0] = da1[ 9][m->id_ttb];
    m->val.before[1] = da1[11][m->id_ttb];
    m->val.before[0] = da1[ 9][m->id_ttb]*1000.;
    m->val.before[1] = da1[11][m->id_ttb]*1000.;
    pos = name_list_pos("monon", nl);
    if(nl->inform[pos] > 0) {
          xlimit = command_par_value("monon",cmd->clone);
          if(frndm() > xlimit) {
             m->enable = 0;
             printf("Monitor %s disabled\n",m->p_node->name);
          }
    }

    /* scaling error should come first, monitor alignment not scaled ... */
    pos = name_list_pos("monscale", nl);
    if(nl->inform[pos] > 0) {
      if((command_par_value("monscale",cmd->clone)) == 1) {
        if(m->p_node->p_al_err != NULL) {
          if (get_option("debug")) {
            printf("m-list: %d %s %s\n",m->id_ttb,m->p_node->name,m->p_node->base_name);
            printf("scales: %e %e \n",m->p_node->p_al_err->a[12], m->p_node->p_al_err->a[13]);
          }
          m->val.before[0] = m->val.before[0]*(1.0 + m->p_node->p_al_err->a[12]);
          m->val.before[1] = m->val.before[1]*(1.0 + m->p_node->p_al_err->a[13]);
        }
      }
    }

    /* monitor misalignment after all other reading manipulations ! */
    pos = name_list_pos("monerror", nl);
    if(nl->inform[pos] > 0) {
      if((command_par_value("monerror",cmd->clone)) == 1) {
        if(m->p_node->p_al_err != NULL) {
          if (get_option("debug")) {
            printf("m-list: %d %s %s\n",m->id_ttb,m->p_node->name,m->p_node->base_name);
            printf("errors: %e %e \n",m->p_node->p_al_err->a[6], m->p_node->p_al_err->a[7]);
          }
          m->val.before[0] += m->p_node->p_al_err->a[6]*1000.;
          m->val.before[1] += m->p_node->p_al_err->a[7]*1000.;
        }
      }
    }
    m = m->next;
  };
  i = 0;
  return(i);
}
int pro_correct_getcorrs(struct in_cmd* cmd)
{
  int i;
  struct id_mic *c;  /* access to tables for monitors and correctors */
  struct table *ttb;
  double **da1;
  ttb = twiss_table;
  da1 = ttb->d_cols;
  c = correct_orbit->cor_table;
  while(c) {
    c->val.before[0] = da1[59][c->id_ttb]*1000.;
    c->val.before[1] = da1[60][c->id_ttb]*1000.;
    /*
    printf("c-list: %d %s %s\n",c->id_ttb,c->p_node->name,c->p_node->base_name);
    printf("initial strengths: %e %e\n",c->val.before[0],c->val.before[1]);
    */
    c = c->next;
  };
  i = 0;
  return(i);
}
void pro_correct_prtwiss()
{
  int i,j;
  int pr_cols;
  struct table *ttb;
  double **da1;
  ttb = twiss_table;
  printf(" %d %d\n",ttb->curr,ttb->num_cols);
  for (i=0; i<ttb->curr; i++) {

    printf(" %s %s\n",ttb->s_cols[0][i],ttb->s_cols[1][i]);
  }
  da1 = ttb->d_cols;
  for (j=0; j < ttb->curr; j++) {
    printf("\n\n");
    printf("from table: %s \n",ttb->node_nm->p[j]);
    printf("from node:  %s \n",ttb->p_nodes[j]->name);
    printf(" %s %s\n",ttb->s_cols[0][j],ttb->s_cols[1][j]);
    pr_cols = ttb->num_cols;
    pr_cols = 19; /* print only for 20 columns */
    for (i=0; i<pr_cols; i++) {
      if(&da1[i][0] != NULL) {
      printf("%-8s %f\n",twiss_table_cols[i],da1[i][j]);
      }
    }
  }
  return;
}
void pro_correct_write_cocu_table()
{
  int i,j;
  int pr_cols;
  int cp[13] = {1,0,2,9,11,3,6,4,7,5,8,15,17};
  struct table *ttb;
  double **da1;
  FILE *fp1;
  fp1 = fopen ("cocu_in.opt","w");
  ttb = twiss_table;
  pr_cols = ttb->num_cols;
  pr_cols = 13; /* print only for 19 columns */
  fprintf(fp1,"*");
  for (i=0; i<pr_cols; i++) {
    fprintf(fp1,"%-8s ",twiss_table_cols[cp[i]]);
  }
  da1 = ttb->d_cols;
  for (j=0; j < ttb->curr; j++) {
    fprintf(fp1,"\n%s %s ",ttb->s_cols[1][j],ttb->s_cols[0][j]);
    for (i=2; i<pr_cols; i++) {
      if(&da1[cp[i]][0] != NULL) {
      fprintf(fp1," %f",da1[cp[i]][j]);
      }
    }
  }
  return;
}
int pro_correct_filter(int iplane, double sigcut)
{
  int    ic, im, ip, icnt;
  struct id_mic *m, *c;  /* access to tables for monitors and correctors */
  struct table *ttb;
  static char  pl[2] = "xy";
  double **da1;
  double bx_m,by_m;
  double xsig, ysig;
  double xmea, ymea;
  double xlim, ylim;
  double xn, yn;
  double *dmat;
  ttb = twiss_table;
  da1 = ttb->d_cols;
  ic = 0; im = 0; icnt = 0;
  ip = iplane - 1;
  printf("A (normalized) cut of %-2.2f is requested\n",sigcut);
      m = correct_orbit->mon_table;
      xmea= 0.0; ymea = 0.0;
      while(m) {
        if (get_option("debug")) {
      printf("monitor flag: %d\n",m->enable);
        }
      if(m->enable == 1) {
          if(ip == 0) {
          bx_m = da1[3][m->id_ttb];
          } else if(ip == 1) {
          bx_m = da1[6][m->id_ttb];
          }
          xn = m->val.before[ip]/sqrt(bx_m);
          xmea += xn;
        if (get_option("debug")) {
          printf("==> %s %-4.3f %-4.3f \n",m->p_node->name,bx_m,m->val.before[ip]);
          printf("==> %-4.3f %-4.3f\n",xn,yn);
        }
        im++;
      }
      m = m->next;
      };
          xmea = xmea/im;
        if (get_option("debug")) {
          printf("Mean values: %-4.3f \n",xmea);
        }
      m = correct_orbit->mon_table;
      im = 0;
      xsig= 0.0;
      while(m) {
      if(m->enable == 1) {
          if(ip == 0) {
          bx_m = da1[3][m->id_ttb];
          } else if(ip == 1) {
          bx_m = da1[6][m->id_ttb];
          }
          xn = m->val.before[ip]/sqrt(bx_m);
          xsig += (xmea - xn)*(xmea - xn);
        im++;
      }
      m = m->next;
      };
          xsig = sqrt(xsig/im);
        if (get_option("debug")) {
          printf("Sigma values: %-4.3f \n",xsig);
        }
      m = correct_orbit->mon_table;
      while(m) {
      if(m->enable == 1) {
          if(ip == 0) {
          bx_m = da1[3][m->id_ttb];
          } else if(ip == 1) {
          bx_m = da1[6][m->id_ttb];
          }
          xn = (m->val.before[ip]/sqrt(bx_m)) - xmea;
          if(fabs(xn) > (sigcut*xsig)) {
             printf("disabled %s %c = %-4.3f (%-4.3f), limit is %-2.2f*%-4.3f\n",
                     m->p_node->name,pl[ip],xn,m->val.before[ip],sigcut,xsig);
             m->enable = 0;
             icnt++;
          }
      }
      m = m->next;
      };
  return(icnt);
}
double* pro_correct_response_ring(int ip, int nc, int nm)
{
  int    ic, im;
  struct id_mic *m, *c;  /* access to tables for monitors and correctors */
  struct table *ttb;
  double **da1;
  double bx_c,by_c,pix_c,piy_c;
  double bx_m,by_m,pix_m,piy_m;
  double qx0, qy0;
  double respx1, respy1;
  double respx, respy;
  double *dmat;
  ttb = twiss_table;
  da1 = ttb->d_cols;
  ic = 0; im = 0;
  dmat = (double *)mycalloc("pro_correct_response_ring",nc*nm,sizeof(double));
  correct_orbit->qx0 = da1[5][ttb->curr-1];
  correct_orbit->qy0 = da1[8][ttb->curr-1];
  qx0 = correct_orbit->qx0;
  qy0 = correct_orbit->qy0;
  c = correct_orbit->cor_table;
  ic = 0;
  while(c) {
    if (get_option("debug")) {
    printf("corrector flag: %d\n",c->enable);
    }
    if(c->enable == 1) {
      bx_c = da1[3][c->id_ttb];
      by_c = da1[6][c->id_ttb];
      pix_c = da1[5][c->id_ttb];
      piy_c = da1[8][c->id_ttb];
      m = correct_orbit->mon_table;
      im = 0;
      while(m) {
        if (get_option("debug")) {
      printf("monitor flag: %d\n",m->enable);
        }
      if(m->enable == 1) {
        bx_m = da1[3][m->id_ttb];
        by_m = da1[6][m->id_ttb];
        pix_m = da1[5][m->id_ttb];
        piy_m = da1[8][m->id_ttb];
        respx = 0.0;
        respy = 0.0;
        if(ip == 1) {
            respx1 = cos((fabs(pix_m - pix_c)*twopi) - qx0*pi);
            respx = respx1*sqrt(bx_m*bx_c)/(2.0*sin(pi*qx0));
            setup_(&respx, dmat, &im, &ic, &nm, &nc);
        }  else if (ip == 2) {
            respy1 = cos((fabs(piy_m - piy_c)*twopi) - qy0*pi);
            respy = respy1*sqrt(by_m*by_c)/(2.0*sin(pi*qy0));
            setup_(&respy, dmat, &im, &ic, &nm, &nc);
        }
        im++;
      }
      m = m->next;
      };
      ic++;
    }
    c = c->next;
  };
  return(dmat);
}
double* pro_correct_response_line(int ip, int nc, int nm)
{
  int    ic, im;
  struct id_mic *m, *c;  /* access to tables for monitors and correctors */
  struct table *ttb;
  double **da1;
  double bx_c,by_c,pix_c,piy_c;
  double bx_m,by_m,pix_m,piy_m;
  double qx0, qy0;
  double respx1, respy1;
  double respx, respy;
  double *dmat;
  ttb = twiss_table;
  da1 = ttb->d_cols;
  ic = 0; im = 0;
  dmat = (double *)mycalloc("pro_correct_response_ring",nc*nm,sizeof(double));
  correct_orbit->qx0 = da1[5][ttb->curr-1];
  correct_orbit->qy0 = da1[8][ttb->curr-1];
  qx0 = correct_orbit->qx0;
  qy0 = correct_orbit->qy0;
  c = correct_orbit->cor_table;
  ic = 0;
  while(c) {
    if(c->enable ==1) {
      bx_c = da1[3][c->id_ttb];
      by_c = da1[6][c->id_ttb];
      pix_c = da1[5][c->id_ttb];
      piy_c = da1[8][c->id_ttb];
      m = correct_orbit->mon_table;
      im = 0;
      while(m) {
      if(m->enable ==1) {
        bx_m = da1[3][m->id_ttb];
        by_m = da1[6][m->id_ttb];
        pix_m = da1[5][m->id_ttb];
        piy_m = da1[8][m->id_ttb];
        respx = 0.0;
        respy = 0.0;
        if(ip == 1) {
            if(pix_m > pix_c) {
              respx1 = sin((pix_m - pix_c)*twopi);
              respx = respx1*sqrt(bx_m*bx_c);
            } else {
              respx = 0.0;
            }
          /*          printf("resp: %d %d %f\n",im,ic,respx);
      printf("++ %s %s %le   \n",m->p_node->name,c->p_node->name,respx);
            */
            setup_(&respx, dmat, &im, &ic, &nm, &nc);
        }  else if (ip == 2) {
            if(piy_m > piy_c) {
              respy1 = sin((piy_m - piy_c)*twopi);
              respy = respy1*sqrt(by_m*by_c);
            } else {
              respy = 0.0;
            }
            setup_(&respy, dmat, &im, &ic, &nm, &nc);
        }
        im++;
      }
      m = m->next;
      };
      ic++;
    }
    c = c->next;
  };
  return(dmat);
}
void pro_correct_make_corr_table()
{
  struct table *ttb;
  int j;
  static char atm[5][4] = {"hmon","vmon","hkic","vkic","kick"};
  ttb = twiss_table;
  for (j=0; j < ttb->curr; j++) {
    if((strncmp(atm[2],ttb->p_nodes[j]->base_name,4) == 0) ||
       (strncmp(atm[3],ttb->p_nodes[j]->base_name,4) == 0) ||
       (strncmp(atm[4],ttb->p_nodes[j]->base_name,4) == 0))  {
      string_to_table("corr","name",ttb->p_nodes[j]->name);
      augment_count("corr");
    }
  }
}
void pro_correct_make_mon_table()
{
  struct table *ttb;
  int j;
  static char atm[3][4] = {"hmon","vmon","moni"};
  ttb = twiss_table;
  for (j=0; j < ttb->curr; j++) {
    if((strncmp(atm[0],ttb->p_nodes[j]->base_name,4) == 0) ||
       (strncmp(atm[1],ttb->p_nodes[j]->base_name,4) == 0) ||
       (strncmp(atm[2],ttb->p_nodes[j]->base_name,4) == 0))  {
      string_to_table("mon","name",ttb->p_nodes[j]->name);
      augment_count("mon");
    }
  }
}
void pro_correct_fill_corr_table(int ip ,char *name, double old, double new)
{
  struct table *cor;
  int j;
  cor =  corr_table;
  for (j=0; j < cor->curr; j++) {
    if(strcmp(name,cor->s_cols[0][j]) == 0) {
      cor->d_cols[ip][j] = old;
      cor->d_cols[ip+2][j] = new;
    }
  }
}
void pro_correct_fill_mon_table(int ip ,char *name, double old, double new)
{
  struct table *mon;
  int j;
  mon =  mon_table;
  for (j=0; j < mon->curr; j++) {
    if(strcmp(name,mon->s_cols[0][j]) == 0) {
      mon->d_cols[ip][j] = old;
      mon->d_cols[ip+2][j] = new;
    }
  }
}
void pro_correct_write_results(double *monvec, double *resvec, double *corvec, int *nx, int *nc, int *nm, int imon, int icor, int ip)
{
  int i;
  int rst;
  double corrm;
  struct id_mic *m, *c;      /* access to tables for monitors and correctors */
  m = correct_orbit->mon_table;
  c = correct_orbit->cor_table;
  if(fddata != NULL) {
     rst = get_variable("n");
     fprintf(fddata,"%d %d %e %e %e %e %e %e\n",ip,rst,cprp(monvec,imon),cprp(resvec,imon),crms(monvec,imon),crms(resvec,imon),copk(monvec,imon),copk(resvec,imon));
  }
  if(print_correct_opt > 0) {
     printf("CORRECTION SUMMARY:   \n\n");
     printf("rms before correction: %f mm\nrms after correction:  %f mm\n\n",crms(monvec,imon),crms(resvec,imon));
     printf("ptp before correction: %f mm\nptp after correction:  %f mm\n\n",cprp(monvec,imon),cprp(resvec,imon));
  }
  if(print_correct_opt > 1) {
  printf("Monitor:  Before:     After:    Difference:\n");
  printf("           (mm)        (mm)         (mm)   \n");
  }
  for(i=0;i<imon;i++) {
    if(print_correct_opt > 1) {
      printf("%s   %-4.3f     %-4.3f     %-4.3f\n",m[nm[i]].p_node->name,monvec[i],resvec[i],resvec[i]-monvec[i]);
    }
    m[nm[i]].val.after[ip-1] = resvec[i];
    pro_correct_fill_mon_table(ip,m[nm[i]].p_node->name,monvec[i],resvec[i]);
  }
  corrm = copk(corvec,icor);
  printf("Max strength: %e should be less than %e\n",corrm,corrl);
  if(corrm > corrl) {
     printf("++++++ warning: maximum corrector strength larger than limit\n");
  }
  set_variable("corrmax",&corrm);
  if(print_correct_opt > 1) {
     printf("Max strength: %e\n",copk(corvec,icor));
     printf("Corrector:  Before:     After:    Difference:\n");
     printf("             (mrad)     (mrad)       (mrad)  \n");
  }
  for(i=0;i<icor;i++) {
    if(fddata != NULL) {
       fprintf(fcdata,"%s %e\n",c[nc[i]].p_node->name,corvec[nx[i]-1]);
    }
    if(print_correct_opt > 1) {
      printf("%s %-3.6f %-3.6f %-3.6f\n",c[nc[i]].p_node->name,c[nc[i]].val.before[ip-1],corvec[nx[i]-1],corvec[nx[i]-1]-c[nc[i]].val.before[ip-1]);
    }
    c[nc[i]].val.after[ip-1] = corvec[nx[i]-1];
    if(ip == 1) {
      c[nc[i]].p_node->chkick += 0.001*corvec[nx[i]-1];
    } else if (ip == 2) {
      c[nc[i]].p_node->cvkick += 0.001*corvec[nx[i]-1];
    }
    pro_correct_fill_corr_table(ip,c[nc[i]].p_node->name,c[nc[i]].val.before[ip-1],corvec[nx[i]-1]);
  }
}
int pro_correct_getactive(int ip, int *nm, int *nx, int *nc, double *corvec, double *monvec,char *conm)
{
  int  imon, icor;
  int  imona, icora;
  struct id_mic *m, *c;
  m = correct_orbit->mon_table;
  imon = 0;
  imona = 0;
  while(m) {
    if (get_option("debug")) {
      printf("from list: %d %s %s\n",m->id_ttb,m->p_node->name,m->p_node->base_name);
      printf("orbit readings: %d %f %f\n",ip,m->val.before[0],m->val.before[1]);
    }
    if(m->enable == 1) {
      monvec[imon] = m->val.before[ip-1];
      nm[imon] = imona;
      imon++;
    }
    imona++;
    m = m->next;
  };
  c = correct_orbit->cor_table;
  icor = 0;
  icora = 0;
  while(c) {
    if (get_option("debug")) {
      printf("from list: %d %d %s %s\n",c->enable,c->id_ttb,c->p_node->name,c->p_node->base_name);
      printf("kicker readings: %f %f\n",c->val.before[0],c->val.before[1]);
    }
    if(c->enable == 1) {
      corvec[icor] = c->val.before[ip-1];
      nx[icor] = icora;
      nc[icor] = icora;
      strcpy(conm,c->p_node->name);
      conm+=16;
      /*          printf("nc: %d %d \n",icor,nc[icor]); */
      icor++;
    }
    icora++;
    c = c->next;
  };
  return(10000*imon + icor);
}
void correct_option(struct in_cmd* cmd)
{
  struct name_list* nl = cmd->clone->par_names;
  int i, debug;
  int val, pos, seed;
  if ((debug=get_option("debug")))  {
     fprintf(prt_file, "in coption routine\n");
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
 val = command_par_value("print", cmd->clone);
 if(val == 0) {
      if (debug)  fprintf(prt_file, "print option not set\n");
      print_correct_opt = 0;
 }else {
      if (debug)  fprintf(prt_file, "print option set\n");
      print_correct_opt = val;
 }
 val = command_par_value("debug", cmd->clone);
 if(val == 0) {
      if (debug)  fprintf(prt_file, "debug option not set\n");
      debug_correct_opt = 0;
 }else {
      if (debug)  fprintf(prt_file, "debug option set\n");
      debug_correct_opt = val;
 }
}
void correct_getorbit(struct in_cmd* cmd)
{
}
void correct_putorbit(struct in_cmd* cmd)
{
  int i;
  struct name_list* nl;
  char* filename = command_par_string("file", cmd->clone);
  char* table_name;
  current_twiss = clone_command(find_command("twiss", defined_commands));
  nl = current_twiss->par_names;
  for (i = 0; i < nl->curr; i++) nl->inform[i] = 0;
  pro_twiss();
  table_name = permbuff("orbit");
  orbit_table = make_table(table_name, "orbit", orbit_table_cols,
            orbit_table_types, current_sequ->n_nodes);
  add_to_table_list(orbit_table, table_register);
  fill_orbit_table(orbit_table, twiss_table);
  out_table("orbit", orbit_table, filename);
  current_twiss = delete_command(current_twiss);
}
void correct_usekick(struct in_cmd* cmd)
{
  char temp[12];
  int count = set_enable("kicker", cmd);
  sprintf(temp, "%d", count);
  put_info(temp, "corrector(s) affected");
}
void correct_usemonitor(struct in_cmd* cmd)
{
  char temp[12];
  int count = set_enable("monitor", cmd);
  sprintf(temp, "%d", count);
  put_info(temp, "monitor(s) affected");
}
double crms(double *r, int m)
{
  double xave = {0.0};
  double xrms = {0.0};
  int    i;
  for(i=0; i<m; i++) {
    xave = xave + r[i];
  }
  xave = xave/m;
  for(i=0; i<m; i++) {
    xrms = xrms + (xave - r[i])*(xave - r[i]);
  }
  xrms = sqrt(xrms/m);
  return(xrms);
}
double cprp(double *r, int m)
{
  double xhi  = {-9999.};
  double xlo  = {9999.};
  double xptp = {0.0};
  int    i;
  for(i=0; i<m; i++) {
     if(r[i] < xlo) xlo = r[i];
     if(r[i] > xhi) xhi = r[i];
  }
     xptp = xhi - xlo;
  return(xptp);
}
double copk(double *r, int m)
{
  double xpk  = {-9999.};
  int    i;
  for(i=0; i<m; i++) {
     if(fabs(r[i]) > xpk) xpk = fabs(r[i]);
  }
  return(xpk);
}
unsigned int locf(iadr)            /* changed from locf_ by JMJ, 8/4/2003 */
#define NADUPW 4   /* Number of ADdress Units Per Word */
#define LADUPW 2   /* Logarithm base 2 of ADdress Units Per Word */
     char *iadr;
{
  return( ((unsigned) iadr) >> LADUPW );
}
int c_micit(double *dmat,char *conm, double *monvec,double *corvec,double *resvec,int *nx,float rms,int imon,int icor,int niter)
{
  int *ny;
  int ifail;
  float *ax,*cinx,*xinx,*resx;
  float *rho,*ptop,*rmss,*xrms,*xptp,*xiter;

  /* allocate auxiliary vectors used by correction algorithms */
  ny   =(int *)mycalloc("c_micit_ny",icor,sizeof(int));
  ax   =(float *)mycalloc("c_micit_ax",imon*icor,sizeof(float));
  cinx =(float *)mycalloc("c_micit_cinx",icor,sizeof(float));
  xinx =(float *)mycalloc("c_micit_xinx",imon,sizeof(float));
  resx =(float *)mycalloc("c_micit_resx",imon,sizeof(float));
  rho  =(float *)mycalloc("c_micit_rho",3*icor,sizeof(float));
  ptop =(float *)mycalloc("c_micit_ptop",icor,sizeof(float));
  rmss =(float *)mycalloc("c_micit_rmss",icor,sizeof(float));
  xrms =(float *)mycalloc("c_micit_xrms",icor,sizeof(float));
  xptp =(float *)mycalloc("c_micit_xptp",icor,sizeof(float));
  xiter=(float *)mycalloc("c_micit_xiter",icor,sizeof(float));

  micit_(dmat,conm,monvec,corvec,resvec,nx,&rms,&imon,&icor,&niter,ny,ax,cinx,xinx,resx,rho,ptop,rmss,xrms,xptp,xiter,&ifail);

  free(ny); free(ax); free(cinx); free(xinx); free(resx); free(rho);
  free(ptop); free(rmss); free(xrms); free(xptp); free(xiter);
/*
*/
  return(ifail);
}
void c_haveit(double *dmat,double *monvec,double *corvec,double *resvec,int *nx,int imon,int icor)
{
  double *cb,*xmeas,*xres,*y,*z,*xd;
  cb   =(double *)mycalloc("c_haveit_cb",icor,sizeof(double));
  xmeas=(double *)mycalloc("c_haveit_xmeas",imon,sizeof(double));
  xres =(double *)mycalloc("c_haveit_xres",imon,sizeof(double));
  y    =(double *)mycalloc("c_haveit_y",icor*imon,sizeof(double));
  z    =(double *)mycalloc("c_haveit_z",icor*icor,sizeof(double));
  xd   =(double *)mycalloc("c_haveit_xd",icor,sizeof(double));
  haveit_(dmat,monvec,corvec,resvec,nx,&imon,&icor,cb,xmeas,xres,y,z,xd);

  free(cb); free(xmeas); free(xres); free(y);  free(z); free(xd);

  return;
}
int  c_svddec(double *dmat, int imon, int icor, int *sing)
{
  int    i;
  int    flag;
  int    dbg;
  double *s, *u, *v, *w, *ut, *vt, *wt;
  double *ws, *wv;
  int    *sw;
  s   =(double *)mycalloc("c_svddec_s",icor*imon,sizeof(double));
  u   =(double *)mycalloc("c_svddec_u",icor*imon,sizeof(double));
  v   =(double *)mycalloc("c_svddec_v",icor*imon,sizeof(double));
  w   =(double *)mycalloc("c_svddec_w",icor*imon,sizeof(double));
  ut  =(double *)mycalloc("c_svddec_ut",icor*imon,sizeof(double));
  vt  =(double *)mycalloc("c_svddec_vt",icor*imon,sizeof(double));
  wt  =(double *)mycalloc("c_svddec_wt",icor*imon,sizeof(double));
  ws  =(double *)mycalloc("c_svddec_ws",icor,sizeof(double));
  wv  =(double *)mycalloc("c_svddec_wv",icor,sizeof(double));
  sw  =(int *)mycalloc("c_svddec_sw",icor,sizeof(int));
  dbg = debug_correct_opt;
#ifdef _ORBDBG_
  if(imon >= icor ) {
      svddec_m_(dmat,s,u,v,w,ut,vt,wt,ws,wv,sw,&imon,&icor,&flag,sing,&dbg);
  } else {
      svddec_c_(dmat,s,u,v,w,ut,vt,wt,ws,wv,sw,&imon,&icor,&flag,sing,&dbg);
  }
#endif
  free(s); free(u); free(v); free(w); free(ut); free(vt); free(wt);
  free(ws); free(wv); free(sw);
  return(flag);
}
int  c_svdcorr(double *dmat, double *xin, double *cor, double *res, int *nx, int imon, int icor)
{
  int    i;
  int    flag;
  int    dbg;
  double *s, *u, *v, *w, *ut, *vt, *wt;
  double *xa, *xb, *xp, *wv, *ws;
  int    *sw;
  s   =(double *)mycalloc("c_svdcorr_s",icor*imon,sizeof(double));
  u   =(double *)mycalloc("c_svdcorr_u",icor*imon,sizeof(double));
  v   =(double *)mycalloc("c_svdcorr_v",icor*imon,sizeof(double));
  w   =(double *)mycalloc("c_svdcorr_w",icor*imon,sizeof(double));
  ut  =(double *)mycalloc("c_svdcorr_ut",icor*imon,sizeof(double));
  vt  =(double *)mycalloc("c_svdcorr_vt",icor*imon,sizeof(double));
  wt  =(double *)mycalloc("c_svdcorr_wt",icor*imon,sizeof(double));
  xa  =(double *)mycalloc("c_svdcorr_xa",imon,sizeof(double));
  xb  =(double *)mycalloc("c_svdcorr_xb",imon,sizeof(double));
  xp  =(double *)mycalloc("c_svdcorr_xp",imon,sizeof(double));
  ws  =(double *)mycalloc("c_svdcorr_xp",icor,sizeof(double));
  wv  =(double *)mycalloc("c_svdcorr_xp",icor,sizeof(double));
  sw  =(int *)mycalloc("c_svdcorr_sw",icor,sizeof(int));
  dbg = debug_correct_opt;
#ifdef _ORBDBG_
  if(imon >= icor ) {
      svdcorr_m_(dmat,s,u,v,w,ut,vt,wt,xin,cor,res,
                 xa,xb,xp,ws,wv,sw,
                 nx,&imon,&icor,&flag,&dbg);
  } else {
      svdcorr_c_(dmat,s,u,v,w,ut,vt,wt,xin,cor,res,
                 xa,xb,xp,ws,wv,sw,
                 nx,&imon,&icor,&flag,&dbg);
  }
#endif
  free(s); free(u); free(v); free(w); free(ut); free(vt); free(wt);
  free(sw); free(xa); free(xb); free(xp); free(ws); free(wv);
  return(flag);
}
void   f_ctof(int *j, char *string, int *nel)
{
    long i, flg;
    flg = 0L;
    for(i=0;i<*nel;i++)    {
        if(flg == 1L)    {
            string[i] = ' ';
            continue;
        }
    if(string[i] == '\0')     {
        string[i] = ' ';
        flg = 1L;
        continue;
    }
    }
        *j = i;
}

#ifdef _WIN32
#include <time.h> /* for gettimeofday */
#else
#include <sys/time.h> /* for gettimeofday */
#endif
float fextim()
{
   #if ( __GNUC__==2 && __GNUC_MINOR__ > 94 ) || __GNUC__ > 2 //hbu  gettimeofday available
     float mytime;
     struct timeval tp;
     struct timezone tzp;
     gettimeofday(&tp,&tzp);
     mytime = (float)(tp.tv_sec%10000) + 1.e-6 * tp.tv_usec; // seconds from epoch, modulo 10 000
   #else // use old ftime
     struct timeb tp;
     float mytime;
     ftime(&tp);
     mytime = (float)(tp.time%10000) + 0.001*tp.millitm;
   #endif
     /* printf("Time now:  %-6.3f\n",mytime);    */
     return(mytime);
}
