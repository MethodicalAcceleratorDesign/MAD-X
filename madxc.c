void correct_correct(struct in_cmd* cmd)
{
  int ix, im, ip, it;
  int j,err,nnnseq;
  int imon, icor;
  int ncorr, nmon;
  int niter;
  int resout;
  int twism;
  float  rms;
  double tmp1, tmp2, tmp3, tmp4;
  char    *clist, *mlist;    /* file names for monitor and corrector output */ 
  double  *dmat;             /* response matrix, double precision */
  double  *corvec, *monvec;  /* vectors to hold measured orbit and correctors */
  double  *resvec;           /* vector to hold corrected orbit */
  char    *conm;             /* vector to hold corrector names (for MICADO) */
  static int     *nm, *nx, *nc;

  ip = pro_correct_getcommands(cmd);
  im = pro_correct_gettables(ip);
  ncorr = im%10000; nmon  = im/10000;
  printf("%d monitors and %d correctors found in input\n",nmon,ncorr);
  setbuf(stdout,NULL);

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
  nx = (int *)mycalloc("correct_correct_nx",ncorr,sizeof(int));
  nc = (int *)mycalloc("correct_correct_nc",ncorr,sizeof(int));
  nm = (int *)mycalloc("correct_correct_nm",nmon,sizeof(int));
  corvec = (double *)mycalloc("correct_correct_corvec",ncorr,sizeof(double));
  monvec = (double *)mycalloc("correct_correct_monvec",nmon,sizeof(double));
  resvec = (double *)mycalloc("correct_correct_resvec",nmon,sizeof(double));
  conm = (char *)mycalloc("correct_correct_conm",ncorr*16,sizeof(char));

  /* get input orbit, default is from input Twiss-table */
  it = pro_correct_getorbit(cmd);

  ix = pro_correct_getactive(ip, nm, nx, nc, corvec, monvec, conm);
  icor = ix%10000; imon  = ix/10000;
  printf("%d monitors and %d correctors enabled\n",imon,icor);

  ncorr = icor;
  nmon = imon;
  /* set up response matrix for ring or line */
  if(strcmp("ring",command_par_string("flag",cmd->clone)) == 0) {
    dmat = (double *)pro_correct_response_ring(ip,ncorr,nmon); } 
  else if(strcmp("line",command_par_string("flag",cmd->clone)) == 0) {
    dmat = (double *)pro_correct_response_line(ip,ncorr,nmon); } 
  else { printf("INVALID MACHINE TYPE\n"); exit(-1);}

  if (get_option("debug")) {
    pro_correct_prtwiss();
    pro_correct_write_cocu_table();
  }

/*
  ix = pro_correct_getactive(ip, nm, nx, nc, corvec, monvec, conm);
  icor = ix%10000; imon  = ix/10000;
  printf("%d monitors and %d correctors enabled\n",imon,icor);
*/

  /* SVD correction, use all available correctors */
  /* if NCORR not set, use all available correctors */
  if(strcmp("svd",command_par_string("mode",cmd->clone)) == 0) {
    /*frs haveit_(dmat,monvec,corvec,resvec,nx,&imon,&icor); */
    c_haveit(dmat,monvec,corvec,resvec,nx,imon,icor);
    pro_correct_write_results(monvec, resvec, corvec, nx, nc, nm, imon, icor, ip);
  }

  /* MICADO correction, get desired number of correctors from command */
  if(strcmp("micado",command_par_string("mode",cmd->clone)) == 0) {
    if((niter =command_par_value("ncorr",cmd->clone)) == 0) niter = icor;
    if((niter =command_par_value("ncorr",cmd->clone)) < 0) niter = 0;     
    if(niter > icor) {
          printf("Fewer correctors available than requested by ncorr\n");
          printf("ncorr reset to %d\n",icor);
          niter = icor;
    }
    
    rms  = 1000.0*command_par_value("error",cmd->clone);
    /*frs       micit_(dmat,monvec,corvec,resvec,nx,&rms,&imon,&icor,&niter); */
    c_micit(dmat,conm,monvec,corvec,resvec,nx,rms,imon,icor,niter);
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

  free(dmat);
  free(nx); free(nc); free(nm);
  free(corvec); free(monvec); free(resvec);
  free(conm);
  return;
}

int pro_correct_getcommands(struct in_cmd* cmd)
{

  static char att[10][8] = {"iterate","plane","ncorr","error","clist",
                            "mlist", "flag","mode","",""};

  int n_iter, ncorr;
  static int iplane = 1;
  char plane[20];
 
  if (get_option("debug")) printf("enter CORRECT module\n");
  if (current_sequ == NULL || current_sequ->ex_start == NULL)
    {
      warning("CORRECT, but no active sequence:", "ignored");
      return(-1);
    }
  setbuf(stdout,(char *)0);

  n_iter = command_par_value(att[0],cmd->clone);
  ncorr  = command_par_value(att[2],cmd->clone);
  printf("mode is: %s\n",command_par_string(att[7],cmd->clone));
  printf("flag is: %s\n",command_par_string(att[6],cmd->clone));
  strcpy(plane,command_par_string(att[1],cmd->clone));
  if(strcmp("x",plane) == 0) {
    iplane = 1; 
  } else if (strcmp("y",plane) == 0) {
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
    corr_table = make_table("corr", corr_table_cols,                     
			    corr_table_types, 5000);
    add_to_table_list(corr_table, table_register);
    pro_correct_make_corr_table();
  }
  if(mon_table == NULL) {
    mon_table = make_table("mon", mon_table_cols,                     
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

  setbuf(stdout,(char *)0);

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
  setbuf(stdout,(char *)0);

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
  setbuf(stdout,(char *)0);

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
  setbuf(stdout,(char *)0);

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
  setbuf(stdout,(char *)0);

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
  struct id_mic *m, *c;      /* access to tables for monitors and correctors */

  m = correct_orbit->mon_table;
  c = correct_orbit->cor_table;

  printf("CORRECTION SUMMARY:   \n\n");                                    
  printf("rms before correction: %f mm\nrms after correction:  %f mm\n\n",crms(monvec,imon),crms(resvec,imon));
  printf("ptp before correction: %f mm\nptp after correction:  %f mm\n\n",cprp(monvec,imon),cprp(resvec,imon));
  if(fddata != NULL) {
     rst = get_variable("n");
     fprintf(fddata,"%d %d %e %e %e %e %e %e\n",ip,rst,cprp(monvec,imon),cprp(resvec,imon),crms(monvec,imon),crms(resvec,imon),copk(monvec,imon),copk(resvec,imon));
  }
  printf("Monitor:  Before:     After:    Difference:\n");                                    
  printf("           (mm)        (mm)         (mm)   \n");                                    
  for(i=0;i<imon;i++) {
    printf("%s   %-4.3f     %-4.3f     %-4.3f\n",m[nm[i]].p_node->name,monvec[i],resvec[i],resvec[i]-monvec[i]);
    m[nm[i]].val.after[ip-1] = resvec[i];
    pro_correct_fill_mon_table(ip,m[nm[i]].p_node->name,monvec[i],resvec[i]);
  }
  printf("Max strength: %e\n",copk(corvec,icor));
  printf("Corrector:  Before:     After:    Difference:\n");                                    
  printf("             (mrad)     (mrad)       (mrad)  \n");                                    
  for(i=0;i<icor;i++) {
  if(fddata != NULL) {
     fprintf(fcdata,"%s %e\n",c[nc[i]].p_node->name,corvec[nx[i]-1]);
  }
    printf("%s %-3.6f %-3.6f %-3.6f\n",c[nc[i]].p_node->name,c[nc[i]].val.before[ip-1],corvec[nx[i]-1],corvec[nx[i]-1]-c[nc[i]].val.before[ip-1]);
    c[nc[i]].val.after[ip-1] = corvec[nx[i]-1];
    if(ip == 1) {
      c[nc[i]].p_node->chkick = 0.001*corvec[nx[i]-1];
    } else if (ip == 2) {
      c[nc[i]].p_node->cvkick = 0.001*corvec[nx[i]-1];
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
unsigned int locf_(iadr)
#define NADUPW 4   /* Number of ADdress Units Per Word */
#define LADUPW 2   /* Logarithm base 2 of ADdress Units Per Word */
     char *iadr;
{
  return( ((unsigned) iadr) >> LADUPW );
}

void c_micit(double *dmat,char *conm, double *monvec,double *corvec,double *resvec,int *nx,float rms,int imon,int icor,int niter)
{
  int *ny;
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
  
  micit_(dmat,conm,monvec,corvec,resvec,nx,&rms,&imon,&icor,&niter,ny,ax,cinx,xinx,resx,rho,ptop,rmss,xrms,xptp,xiter);

  free(ny);
  free(ax);
  free(cinx);
  free(xinx);
  free(resx);
  free(rho);
  free(ptop);
  free(rmss);
  free(xrms);
  free(xptp);
  free(xiter);
/*
*/

  return;
}

void c_haveit(double *dmat,double *monvec,double *corvec,double *resvec,int *nx,int imon,int icor)
{
  double *cb,*xmeas,*xres,*y,*z,*xd;

  cb   =(double *)mycalloc("c_haveit_cb",icor,sizeof(float));
  xmeas=(double *)mycalloc("c_haveit_xmeas",imon,sizeof(float));
  xres =(double *)mycalloc("c_haveit_xres",imon,sizeof(float));
  y    =(double *)mycalloc("c_haveit_y",icor*imon,sizeof(float));
  z    =(double *)mycalloc("c_haveit_z",icor*icor,sizeof(float));
  xd   =(double *)mycalloc("c_haveit_xd",icor,sizeof(float));

  haveit_(dmat,monvec,corvec,resvec,nx,&imon,&icor,cb,xmeas,xres,y,z,xd);

  free(cb);   
  free(xmeas);
  free(xres); 
  free(y);    
  free(z);    
  free(xd);   

  return;
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
