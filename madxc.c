#ifdef _WRAP_FORTRAN_CALLS
#include "fortran_wrappers.h"
#endif
#ifdef _WRAP_C_CALLS
#include "c_wrappers.h"
#endif

#ifdef _WIN32
#ifndef _UINTPTR_T_
#define _UINTPTR_T_
#define uintptr_t unsigned int	/* 32 bytes-long (should be 64 on WIN64) */
#endif
#else
#include <stdint.h>		/* uintptr_t, to fit pointers into integers of correct size */
#endif
/* frs 21.10.2008 revert to old version after Thys Risselada's fix of Micado */ 

void setupi_(int*, int*, int*, int*, int*, int*);
void primat_(int*, int*, int*);
void prdmat_(double*, int*, int*);
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
  else if (strcmp(cmd->tok_list->p[0], "readmytable") == 0)
    {
     read_my_table(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "readcorr") == 0)
    {
     correct_readcorr(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "setcorr") == 0)
    {
     correct_setcorr(cmd);
    }
  else if (strcmp(cmd->tok_list->p[0], "coption") == 0)
    {
     correct_option(cmd);
    }
}


void correct_correct(struct in_cmd* cmd)
/* Steering routine for orbit corrections */
{
/*
  char rout_name[] = "correct_correct";
*/
  char  *orbtab1, *orbtab2;
  

  /* Call for one or two ring orbit correction */
  if(command_par_value("tworing",cmd->clone)) {
           if ((orbtab1 = command_par_string("beam1tab",cmd->clone)) == NULL) {
              fatal_error("Two beam correction requested but no table supplied for beam 1",orbtab1);
           }
           if ((orbtab2 = command_par_string("beam2tab",cmd->clone)) == NULL) {
              fatal_error("Two beam correction requested but no table supplied for beam 2",orbtab2);
           }
           printf("Want to use orbits from: %s and : %s\n",orbtab1,orbtab2);
           correct_correct2(cmd);
  }  else {
           printf("Want to correct orbit of a single ring\n");                    
           if ((orbtab1 = command_par_string("beam1tab",cmd->clone)) != NULL) {
              warning(" "," ");
              warning("Single beam correction requested but beam 1 table supplied:",orbtab1);
              warning("Requested table ignored:",orbtab1);
              warning(" "," ");
           }
           if ((orbtab2 = command_par_string("beam2tab",cmd->clone)) != NULL) {
              warning(" "," ");
              warning("Single beam correction requested but beam 2 table supplied:",orbtab2);
              warning("Requested table ignored:",orbtab2);
              warning(" "," ");
           }
           correct_correct1(cmd);
  }
}


void correct_correct2(struct in_cmd* cmd)
/* Steering routine for orbit corrections of two beams */
{
  char rout_name[] = "correct_correct2";

/*
  struct name_list* spos = sequences->list;
  struct table *twb1;
  struct table *twb2;
  int idrop;
  int pos;
  struct timeb tp;
  int sflag, svdflg;
  double  sigcut;         
*/

  int ix, im, ip, it;
  int i,j,err,nnnseq;
  int imon, icor;
  int ncorr, nmon;
  int niter;
  int resout;
  int twism;
  int ifail;
  float  rms;
  double rrms;
  double tmp1, tmp2, tmp3, tmp4;
  char    *clist, *mlist;    /* file names for monitor and corrector output */
  char    clist1[100], clist2[100];   /* file names for corrector output ring 1 and ring 2 */
  double  *dmat = {NULL};    /* response matrix, double precision */
  double  *corvec, *monvec;  /* vectors to hold measured orbit and correctors */
  double  *resvec;           /* vector to hold corrected orbit */
  char    *conm;             /* vector to hold corrector names (for MICADO) */
  int     *sing;             /* array to store pointer to singular correctors */
  static int     *nm, *nx, *nc;
  struct id_mic2  *c;
  /*
  struct id_mic2   *m;
  */

  strcpy(clist1,"\0");
  strcpy(clist2,"\0");

  printf("for two beam orbit corrections ...\n");
  ip = pro_correct_getcommands(cmd);
  im = pro_correct2_gettables(ip,cmd);
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
     if(fgdata == NULL) {
        if((fgdata = fopen("plot.orb","w")) == NULL)
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
  nx  = (int *)mycalloc("correct_correct2_nx",ncorr,sizeof(int));
  nc  = (int *)mycalloc("correct_correct2_nc",ncorr,sizeof(int));
  nm  = (int *)mycalloc("correct_correct2_nm",nmon,sizeof(int));
  sing = (int *)mycalloc("correct_correct2_sing",ncorr*2,sizeof(int));
  corvec = (double *)mycalloc("correct_correct2_corvec",ncorr,sizeof(double));
  monvec = (double *)mycalloc("correct_correct2_monvec",nmon,sizeof(double));
  resvec = (double *)mycalloc("correct_correct2_resvec",nmon,sizeof(double));
  conm = (char *)mycalloc("correct_correct2_conm",ncorr*16,sizeof(char));

  /* get original settings of correctors from input Twiss-table */
  it = pro_correct2_getcorrs(cmd);
  /* get input orbit, default is from input Twiss-table */
  it = pro_correct2_getorbit(cmd);


  /* find and prepare enabled correctors and monitors, may be repeated */
  ix = pro_correct2_getactive(ip, nm, nx, nc, corvec, monvec, conm);
  icor = ix%10000; imon  = ix/10000;
  printf("%d monitors and %d correctors enabled\n",imon,icor);


    if (get_option("debug")) {
  for (i=0;i<icor;i++) {
    printf("C: %d %d \n",nx[i],nc[i]);
  }
  for (i=0;i<imon;i++) {
    printf("M: %d %e \n",nm[i],monvec[i]);
  }
  }

  if(strcmp("ring",command_par_string("flag",cmd->clone)) == 0) {
    if(dmat != NULL) myfree(rout_name,dmat);
    /* icor and imon used to set up correct matrix size !! */
    dmat = (double *)pro_correct2_response_ring(ip,icor,imon);
  }
  else { printf("INVALID MACHINE TYPE\n"); exit(-1);
  }

  /* MICADO correction, get desired number of correctors from command */
  corrl = command_par_value("corrlim",cmd->clone);
  set_variable("corrlim",&corrl);
  if(strcmp("micado",command_par_string("mode",cmd->clone)) == 0) {
    printf("enter MICADO correction ...\n");
    if((niter = command_par_value("ncorr",cmd->clone)) == 0) {
          printf("Requested %d correctors (\?\?\?) set to %d\n",niter,icor);
          niter = icor;
    }
    else if((niter = command_par_value("ncorr",cmd->clone)) < 0) {
          printf("Requested %d correctors (\?\?\?) set to 0\n",niter);
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
    printf("Back from micado %d\n",ifail);
    if(ifail != 0) {
       printf("MICADO correction completed with error code %d\n\n",ifail);
       warning("MICADO back with error",", no correction done");
    }
    rrms = crms(monvec,imon);
    printf("RMS before %e\n",rrms);
    rrms = crms(resvec,imon);
    printf("RMS after  %e\n",rrms);
    if(fgdata != NULL) {
    for (i=0; i<nmon; i++) {
       fprintf(fgdata,"%e %e \n",monvec[i],resvec[i]);
    }
    }
/*
    for (i=0; i<nmon; i++) {
       printf("monvec: %d %e \n",i,monvec[i]);
    }
    printf("\n");
    for (i=0; i<nmon; i++) {
       printf("resvec: %d %e \n",i,resvec[i]);
    }
    m = correct_orbit12->mon_table;
    for (i=0; i<nmon; i++) {
       printf("resvec: %s %e \n",m[nm[i]].p_node->name,resvec[i]);
    }
    printf("\n");
    for (i=0; i<ncorr; i++) {
       printf("corvec: %d %e \n",i,corvec[i]);
    }
    printf("\n");
*/

    c = correct_orbit12->cor_table;
    for(i=0;i<icor;i++) {
       printf("%s %e\n",c[nc[i]].p_node->name,corvec[nx[i]-1]);
    }
    printf("\n");
 
    /* printf("Time after micado:  %-6.3f\n",fextim());   */
    if(ifail != 0) {
       printf("MICADO correction completed with error code %d\n\n",ifail);
       warning("MICADO back with error",", no correction done");
    }
    if(ifail == 0) {
       pro_correct2_write_results(monvec, resvec, corvec, nx, nc, nm, imon, icor, ip);
    }
  }
 
  /* write corrector output to tfs table */
  if ((clist = command_par_string("clist",cmd->clone)) != NULL) {
    strcat(clist1,clist);
    strcat(clist1,"_1");
    strcat(clist2,clist);
    strcat(clist2,"_2");
    out_table("corr1",corr_table1,clist1);
    out_table("corr2",corr_table2,clist2);
  }

  /* write monitor output to tfs table */
  if ((mlist = command_par_string("mlist",cmd->clone)) != NULL) {
    out_table("mon",mon_table,mlist);
  }


  /* Clean up at the end of the module */
   
  myfree(rout_name,nm);myfree(rout_name,dmat);myfree(rout_name,nx);
  myfree(rout_name,nc);myfree(rout_name,corvec);
  myfree(rout_name,monvec);myfree(rout_name,resvec); myfree(rout_name,conm);
   
  return;

}

int  pro_correct2_gettables(int iplane, struct in_cmd* cmd)
{

  char rout_name[] = "pro_correct2_gettables";

  struct id_mic2 *cor_l1,  *cor_l2;
  struct id_mic2 *mon_l1,  *mon_l2;
  struct id_mic2 *cor_l12, *mon_l12;
  struct id_mic2 *prt;

  struct table *ttb;

  struct table *b1 = NULL;
  struct table *b2 = NULL;

  char* orbtab1;
  char* orbtab2;

  int t1, t2;   

  int ebl1, ebl2;

  int j,k;
  int set0;
  int cntm1 = {0};
  int cntc1 = {0};
  int cntm2 = {0};
  int cntc2 = {0};
  int cntm12 = {0};
  int cntc12 = {0};
  
  double ounits;

/*
  static char atm[6][4] = {"hmon","vmon","moni","hkic","vkic","kick"};
*/

/* Get access to tables, for orbit and model the default is twiss_table */

  if ((orbtab1 = command_par_string("beam1tab",cmd->clone)) != NULL) {
      printf("Want to use orbit from: %s\n",orbtab1);
    if ((t1 = name_list_pos(orbtab1, table_register->names)) > -1) {
       b1 = table_register->tables[t1];
    } else {
       fatal_error("Beam 1 ORBIT table requested, but not provided:",orbtab1);
    }
  } else {

  }
  if ((orbtab2 = command_par_string("beam2tab",cmd->clone)) != NULL) {
      printf("Want to use orbit from: %s\n",orbtab2);
    if ((t2 = name_list_pos(orbtab2, table_register->names)) > -1) {
       b2 = table_register->tables[t2];
    } else {
       fatal_error("Beam 2 ORBIT table requested, but not provided:",orbtab2);
    }
  } else {

  }

/* store as globals for later use */
  if((b1 != NULL) && (b2 != NULL)) {
     twiss_table_beam1 = b1;
     twiss_table_beam2 = b2;
  } else {
       fatal_error("Beam 1 and 2 orbit tables not found:",orbtab1);
  }
 
/* reserve space for orbit correction structures */
  if(correct_orbit12 == NULL) {
    correct_orbit12 = (struct orb_cor2*)mycalloc("pro_correct2_gettables",1, 
                     sizeof(struct orb_cor2));
  }



  if(correct_orbit12->cor_table != NULL) myfree(rout_name,correct_orbit12->cor_table);
  if(correct_orbit12->mon_table != NULL) myfree(rout_name,correct_orbit12->mon_table);

  correct_orbit12->cor_table = (struct id_mic2 *)mycalloc("pro_correct2_gettables_cor",5200, sizeof(struct id_mic2));
  correct_orbit12->mon_table = (struct id_mic2 *)mycalloc("pro_correct2_gettables_mon",5200, sizeof(struct id_mic2));

/* orbit table available, get units, if defined */
   if((ounits = command_par_value("units",cmd->clone)) > 0) {
          correct_orbit12->units=ounits;                      
   } else {
          correct_orbit12->units=1.0;      
   }

  ttb = model_table;
/* no more need, we have b1 and b2 as pointers .. */

  correct_orbit12->mon_table->previous = NULL;
  correct_orbit12->mon_table->next = NULL;
  correct_orbit12->cor_table->previous = NULL;
  correct_orbit12->cor_table->next = NULL;

  mon_l1 = correct_orbit12->mon_table;
  cor_l1 = correct_orbit12->cor_table;
 
  for (j=0; j < b1->curr; j++) {
   if((strncmp(atm[iplane-1],b1->p_nodes[j]->base_name,4) == 0) ||
      (strncmp(atm[2],      b1->p_nodes[j]->base_name,4) == 0))  {
/*    printf("1m: %s %ld\n", b1->p_nodes[j]->name, strstr(".b2", b1->p_nodes[j]->name)); */
      if(strstr(b1->p_nodes[j]->name,".b1") != NULL) {
        mon_l1->id_ttb[0] = j;
        mon_l1->id_ttb[1] = -1;
        mon_l1->enable = b1->p_nodes[j]->enable;
        mon_l1->p_node = b1->p_nodes[j];
        mon_l1->next = mon_l1;
        mon_l1->next++; mon_l1++;
        cntm1++;
      } else {
/*      printf("Removed: %s\n",b1->p_nodes[j]->name); */
      }
    }
    if((strncmp(atc[iplane-1],b1->p_nodes[j]->base_name,4) == 0) ||
       (strncmp(atc[2],       b1->p_nodes[j]->base_name,4) == 0))  {
/*    printf("1c: %s %ld\n", b1->p_nodes[j]->name, b1->p_nodes[j]->name); */
      if(strstr(b1->p_nodes[j]->name,".b1") != NULL) {
        cor_l1->id_ttb[0] = j;
        cor_l1->id_ttb[1] = -1;
        cor_l1->enable = b1->p_nodes[j]->enable;
        cor_l1->p_node = b1->p_nodes[j];
        cor_l1->p_node_s1 = b1->p_nodes[j];
        cor_l1->p_node_s2 = NULL;            
        if((set0 = command_par_value("corzero",cmd->clone)) > 0) {
          if(iplane == 1) cor_l1->p_node_s1->chkick = 0.0;
          if(iplane == 2) cor_l1->p_node_s1->cvkick = 0.0;
        }
        cor_l1->next = cor_l1;
        cor_l1->next++; cor_l1++;
        cntc1++;
      } else {
/*      printf("Removed: %s\n",b1->p_nodes[j]->name); */
      }
    }
  }

  mon_l2 = mon_l1;
  cor_l2 = cor_l1;
  for (j=0; j < b2->curr; j++) {
   if((strncmp(atm[iplane-1],b2->p_nodes[j]->base_name,4) == 0) ||
      (strncmp(atm[2],      b2->p_nodes[j]->base_name,4) == 0))  {
/*    printf("2m: %s %ld\n", b2->p_nodes[j]->name, b2->p_nodes[j]->name); */
      if(strstr(b2->p_nodes[j]->name,".b2") != NULL) {
        mon_l2->id_ttb[0] = -1;
        mon_l2->id_ttb[1] = j;
        mon_l2->enable = b2->p_nodes[j]->enable;
        mon_l2->p_node = b2->p_nodes[j];
        mon_l2->next = mon_l2;
        mon_l2->next++; mon_l2++;
        cntm2++;
      } else {
/*      printf("Removed: %s\n",b2->p_nodes[j]->name); */
      }
    }
    if((strncmp(atc[iplane-1],b2->p_nodes[j]->base_name,4) == 0) ||
       (strncmp(atc[2],       b2->p_nodes[j]->base_name,4) == 0))  {
/*    printf("2c: %s %ld\n", b2->p_nodes[j]->name, b2->p_nodes[j]->name); */
      if(strstr(b2->p_nodes[j]->name,".b2") != NULL) {
        cor_l2->id_ttb[0] = -1;
        cor_l2->id_ttb[1] = j;
        cor_l2->enable = b2->p_nodes[j]->enable;
        cor_l2->p_node = b2->p_nodes[j];
        cor_l2->p_node_s2 = b2->p_nodes[j];
        cor_l2->p_node_s1 = NULL;            
        if((set0 = command_par_value("corzero",cmd->clone)) > 0) {
          if(iplane == 1) cor_l2->p_node_s2->chkick = 0.0;
          if(iplane == 2) cor_l2->p_node_s2->cvkick = 0.0;
        }
        cor_l2->next = cor_l2;
        cor_l2->next++; cor_l2++;
        cntc2++;
      } else {
/*      printf("Removed: %s\n",b2->p_nodes[j]->name); */
      }
    }
  }

  mon_l12 = mon_l2;
  cor_l12 = cor_l2;
  for (j=0; j < b1->curr; j++) {
   if((strncmp(atm[iplane-1],b1->p_nodes[j]->base_name,4) == 0) ||
      (strncmp(atm[2],      b1->p_nodes[j]->base_name,4) == 0))  {
/*    printf("12m: %s \n", b1->p_nodes[j]->name); */
      if((strstr(b1->p_nodes[j]->name,".b1") == NULL)  &&
         (strstr(b1->p_nodes[j]->name,".b2") == NULL)) {
        mon_l12->id_ttb[0] = j;
         for (k=0; k < b2->curr; k++) {
           if(strcmp(b2->p_nodes[k]->name,b1->p_nodes[j]->name) == 0) {
            mon_l12->id_ttb[1] = k;
           }
         }
        mon_l12->enable = b1->p_nodes[j]->enable;
        mon_l12->p_node = b1->p_nodes[j];
        mon_l12->next = mon_l12;
        mon_l12->next++; mon_l12++;
        cntm12++;
      } else {
/*      printf("Removed: %s\n",b1->p_nodes[j]->name); */
      }
    }
    if((strncmp(atc[iplane-1],b1->p_nodes[j]->base_name,4) == 0) ||
       (strncmp(atc[2],       b1->p_nodes[j]->base_name,4) == 0))  {
/*    printf("12c: %s \n", b1->p_nodes[j]->name);     */
      if((strstr(b1->p_nodes[j]->name,".b1") == NULL)  &&
         (strstr(b1->p_nodes[j]->name,".b2") == NULL)) {
         cor_l12->id_ttb[0] = j;
         for (k=0; k < b2->curr; k++) {
           if(strcmp(b2->p_nodes[k]->name,b1->p_nodes[j]->name) == 0) {
            cor_l12->id_ttb[1] = k;
           }
         }
        cor_l12->p_node = b1->p_nodes[j];
        cor_l12->p_node_s1 = b1->p_nodes[cor_l12->id_ttb[0]];
        cor_l12->p_node_s2 = b2->p_nodes[cor_l12->id_ttb[1]];
        ebl1 = b1->p_nodes[cor_l12->id_ttb[0]]->enable;
        ebl2 = b2->p_nodes[cor_l12->id_ttb[1]]->enable;
        cor_l12->enable = ebl1*ebl2;                  
        if((set0 = command_par_value("corzero",cmd->clone)) > 0) {
          if(iplane == 1) cor_l12->p_node_s1->chkick = 0.0;
          if(iplane == 2) cor_l12->p_node_s1->cvkick = 0.0;
          if(iplane == 1) cor_l12->p_node_s2->chkick = 0.0;
          if(iplane == 2) cor_l12->p_node_s2->cvkick = 0.0;
        }
        cor_l12->next = cor_l12;
        cor_l12->next++; cor_l12++;
        cntc12++;
      } else {
/*      printf("Removed: %s\n",b1->p_nodes[j]->name); */
      }
    }
  }
  /* terminate linked list   */
  mon_l12--; mon_l12->next = NULL;
  cor_l12--; cor_l12->next = NULL;

  printf("mons and corrs (beam 1)   : %ld %ld\n",(long int)cntm1, (long int)cntc1);
  printf("mons and corrs (beam 2)   : %ld %ld\n",(long int)cntm2, (long int)cntc2);
  printf("mons and corrs (beam 1+2) : %ld %ld\n",(long int)cntm12, (long int)cntc12);

    if (get_option("debug")) {
     prt = correct_orbit12->mon_table;
     while(prt != NULL) {
       printf("Monitors beam12: %s %ld %ld\n",prt->p_node->name,(long int)prt->id_ttb[0],(long int)prt->id_ttb[1]);
       prt = prt->next;
     }

     prt = correct_orbit12->cor_table;
     while(prt != NULL) {
       printf("Correctors beam12: %s %ld %ld\n",prt->p_node->name,(long int)prt->id_ttb[0],(long int)prt->id_ttb[1]);
       prt = prt->next;
     }
    }

/*
     prt = correct_orbit12->cor_table;
     while(prt != NULL) {
       printf("Correctors beam12: %s %ld %ld\n",prt->p_node->name,prt->id_ttb[0],prt->id_ttb[1]);
       for (j=0; j < b2->curr; j++) {
         if(strcmp(b2->p_nodes[j]->name,prt->p_node->name) == 0) {
            prt->id_ttb[1] = j;
            printf("matched correctors beam12: %s %ld %ld\n",prt->p_node->name,prt->id_ttb[0],prt->id_ttb[1]);
         }
       }
       prt = prt->next;
       
     }
*/

  if(corr_table1 == NULL) {
    corr_table1 = make_table("corr1", "corr1", corr_table_cols,
            corr_table_types, 15000);
    add_to_table_list(corr_table1, table_register);
  }
  if(corr_table2 == NULL) {
    corr_table2 = make_table("corr2", "corr2", corr_table_cols,
            corr_table_types, 15000);
    add_to_table_list(corr_table2, table_register);
  }
  pro_correct2_make_corr_table();

  if(mon_table == NULL) {
    mon_table = make_table("mon", "mon", mon_table_cols,
           mon_table_types, 15000);
    add_to_table_list(mon_table, table_register);
    pro_correct2_make_mon_table();
  }

  return(10000*(cntm1+ cntm2+ cntm12) + (cntc1 + cntc2 + cntc12));


}

int pro_correct2_getorbit(struct in_cmd* cmd)
{
  struct name_list* nl;
  int i;
  int pos;
/*
  int pps, ppt;
*/
  struct id_mic2 *m;  /* access to tables for monitors and correctors */
  double **da1;
  double **da2;
  double xlimit;

  char   strx[40];
  char   stry[40];

  int    posx, posy, pospx, pospy;

  da1 = twiss_table_beam1->d_cols;
  da2 = twiss_table_beam2->d_cols;

  nl = cmd->clone->par_names;

  m = correct_orbit12->mon_table;

  strcpy(strx,"x");
  strcpy(stry,"y");

  if((posx = name_list_pos(strx,twiss_table_beam1->columns)) < 0) { 
      fatal_error("orbit x not found in input table",", MAD-X terminates ");
  }
  if((posy = name_list_pos(stry,twiss_table_beam1->columns)) < 0) { 
      fatal_error("orbit y not found in input table",", MAD-X terminates ");
  }
  if (get_option("debug")) {
    if((pospx = name_list_pos("px",twiss_table_beam1->columns)) < 0) { 
        warning("orbit px not found in input table",", MAD-X continues ");
    }
    if((pospy = name_list_pos("py",twiss_table_beam1->columns)) < 0) { 
        warning("orbit py not found in input table",", MAD-X continues ");
    }
    printf("====c1===>  %d %d %d %d \n",posx,posy,pospx,pospy);
  }


  while(m) {

/* If correction to target orbit, subtract the wanted orbit ... */
    if(m->id_ttb[0] > 0) {
    m->val.before[0] = m->p_node->other_bv*da1[ 9][m->id_ttb[0]];
    m->val.before[1] = m->p_node->other_bv*da1[11][m->id_ttb[0]];
    m->val.before[0] = m->p_node->other_bv*da1[ 9][m->id_ttb[0]]*1000.;
    m->val.before[1] = m->p_node->other_bv*da1[11][m->id_ttb[0]]*1000.;
    } else if (m->id_ttb[1] > 0) {
    m->val.before[0] = m->p_node->other_bv*da2[ 9][m->id_ttb[1]];
    m->val.before[1] = m->p_node->other_bv*da2[11][m->id_ttb[1]];
    m->val.before[0] = m->p_node->other_bv*da2[ 9][m->id_ttb[1]]*1000.;
    m->val.before[1] = m->p_node->other_bv*da2[11][m->id_ttb[1]]*1000.;
    } else {
      printf("BIG SHIT .... \n");
      exit(-10);
    }

    pos = name_list_pos("monon", nl);
       if(nl->inform[pos] > 0) {
          xlimit = command_par_value("monon",cmd->clone);
          if(frndm() > xlimit) {
             m->enable = 0;
             printf("Monitor %s disabled\n",m->p_node->name);
          }
       }
    if (get_option("debug")) {
        printf("m-list: %d %d %s %s\n",m->id_ttb[0],m->id_ttb[1],m->p_node->name,m->p_node->base_name);
        printf("initial reading: %e %e\n\n",m->val.before[0],m->val.before[1]);
    }
       /*
       */
    m = m->next;
  };
  i = 0;
  return(i);
}

int pro_correct2_getcorrs(struct in_cmd* cmd)
{
  int i;
  struct id_mic2 *c;  /* access to tables for monitors and correctors */
  double **da1;
  double **da2;

/*
  double xlimit;
*/
  da1 = twiss_table_beam1->d_cols;
  da2 = twiss_table_beam2->d_cols;

  c = correct_orbit12->cor_table;
  while(c) {
    if(c->id_ttb[0] > 0) {
/*     c->val.before[0] = da1[59][c->id_ttb[0]]*1000.; */
/*     c->val.before[1] = da1[60][c->id_ttb[0]]*1000.; */
       c->val.before[0] = c->p_node_s1->chkick*1000.;
       c->val.before[1] = c->p_node_s1->cvkick*1000.;
    } else if(c->id_ttb[1] > 0) {
/*     c->val.before[0] = da2[59][c->id_ttb[1]]*1000.; */
/*     c->val.before[1] = da2[60][c->id_ttb[1]]*1000.; */
       c->val.before[0] = c->p_node_s2->chkick*1000.;
       c->val.before[1] = c->p_node_s2->cvkick*1000.;
    }
    if (get_option("debug")) {
      printf("c-list: %d %d %s %s\n",c->id_ttb[0],c->id_ttb[1],c->p_node->name,c->p_node->base_name);
      printf("initial strengths: %e %e\n",c->val.before[0],c->val.before[1]);
    }
    /*
    */

    c = c->next;
  };
  i = 0;
  return(i);
}


int pro_correct2_getactive(int ip, int *nm, int *nx, int *nc, double *corvec, double *monvec,char *conm)
{
  int  imon, icor;
  int  imona, icora;
  struct id_mic2 *m, *c;

  m = correct_orbit12->mon_table;
  imon = 0;
  imona = 0;
  while(m) {
    if (get_option("debug")) {
      printf("from list: %d %d %s %s\n",m->id_ttb[0],m->id_ttb[1],m->p_node->name,m->p_node->base_name);
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

  c = correct_orbit12->cor_table;
  icor = 0;
  icora = 0;
  while(c) {
    if (get_option("debug")) {
      printf("from list: %d %d %d %s %s\n",c->enable,c->id_ttb[0],c->id_ttb[1],c->p_node->name,c->p_node->base_name);
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

double* pro_correct2_response_ring(int ip, int nc, int nm)
{
  int    ic, im;
  struct id_mic2 *m, *c;  /* access to tables for monitors and correctors */

  double **da1;
  double **da2;
  double bx_c,by_c,pix_c,piy_c;
  double bx_m,by_m,pix_m,piy_m;
  double qx0, qy0;
  double respx1, respy1;
  double respx, respy;
  double *dmat;
  int  *imat;
  int    mp;
  int    i_zero, i_one;
  int icb;
  int    i, j;

  setbuf(stdout,(char *)0);

  ic = 0; im = 0;
  i_zero = 0; i_one = 1;

  da1 = twiss_table_beam1->d_cols;
  da2 = twiss_table_beam2->d_cols;

  dmat = (double *)mycalloc("pro_correct2_response_ring",nc*nm,sizeof(double));
  imat = (int *)mycalloc("pro_correct2_response_ring",nc*nm,sizeof(int));

/* initialize imat: */
  for (i=0; i<nc; i++) {
     for (j=0; j<nm; j++) {
         setupi_(&i_zero, imat, &j, &i, &nm, &nc);
     }
  }

  c = correct_orbit12->cor_table;
  ic = 0;

  while(c) {
    if (get_option("debug")) {
    printf("corrector flag: %d\n",c->enable);
    }
    if(c->enable == 1) {

      for(icb = 0; icb < 2; icb++) {
      if(c->id_ttb[icb] > 0) {

         if(icb == 0) {
           correct_orbit12->qx0 = da1[5][twiss_table_beam1->curr-1];
           correct_orbit12->qy0 = da1[8][twiss_table_beam1->curr-1];
           qx0 = correct_orbit12->qx0;
           qy0 = correct_orbit12->qy0;
           if(c->id_ttb[icb] > 0) {
              bx_c = da1[3][c->id_ttb[icb]];
              by_c = da1[6][c->id_ttb[icb]];
              pix_c = da1[5][c->id_ttb[icb]];
              piy_c = da1[8][c->id_ttb[icb]];
           } else {
              bx_c = 0.0;                    
              by_c = 0.0;                     
              pix_c = 0.0;                     
              piy_c = 0.0;                         
           }
         } else {
           correct_orbit12->qx0 = da2[5][twiss_table_beam2->curr-1];
           correct_orbit12->qy0 = da2[8][twiss_table_beam2->curr-1];
           qx0 = correct_orbit12->qx0;
           qy0 = correct_orbit12->qy0;
           if(c->id_ttb[icb] > 0) {
              bx_c = da2[3][c->id_ttb[icb]];
              by_c = da2[6][c->id_ttb[icb]];
              pix_c = da2[5][c->id_ttb[icb]];
              piy_c = da2[8][c->id_ttb[icb]];
           } else {
              bx_c = 0.0;                    
              by_c = 0.0;                     
              pix_c = 0.0;                     
              piy_c = 0.0;                         
           }
         }

         m = correct_orbit12->mon_table;
         im = 0;
         while(m) {
              if (get_option("debug")) {
            printf("monitor flag: %d\n",m->enable);
              }
            if(m->enable == 1) {
             if((m->id_ttb[icb] > 0) && (c->id_ttb[icb] > 0)) { 
              if(m->id_ttb[icb] > 0) { 
                if(icb == 0) {
                   mp = m->id_ttb[icb];
                   bx_m = da1[3][mp];
                   by_m = da1[6][mp];
                   pix_m = da1[5][mp];
                   piy_m = da1[8][mp];
                } else {
                   mp = m->id_ttb[icb];
                   bx_m = da2[3][mp];
                   by_m = da2[6][mp];
                   pix_m = da2[5][mp];
                   piy_m = da2[8][mp];
                }
              } else {
                   bx_m = 0.0;        
                   by_m = 0.0;        
                   pix_m = 0.0;       
                   piy_m = 0.0;         
              }
      
              respx = 0.0;
              respy = 0.0;

/*  print Twiss parameters ... */
              if (get_option("debug"))  {
                printf("%s %d %e %e %e %e -- %s %e %e %e %e\n",
                c->p_node->name,icb,bx_c,by_c,pix_c,piy_c,
                m->p_node->name,bx_m,by_m,pix_m,piy_m);
              }

              if(ip == 1) {
                  respx1 = cos((fabs(pix_m - pix_c)*twopi) - qx0*pi);
                  respx = respx1*sqrt(bx_m*bx_c)/(2.0*sin(pi*qx0));
                  if(icb != 0) { respx =  respx; }
                  setup_(&respx, dmat, &im, &ic, &nm, &nc);
              }  else if (ip == 2) {
                  respy1 = cos((fabs(piy_m - piy_c)*twopi) - qy0*pi);
                  respy = respy1*sqrt(by_m*by_c)/(2.0*sin(pi*qy0));
                  if(icb != 0) { respy =  respy; }
                  setup_(&respy, dmat, &im, &ic, &nm, &nc);
              }
              if((fabs(respy) > 0.000006) || (fabs(respx) > 0.000006)) { 
                 if (get_option("debug"))  {
                    printf("true %d %d",ic,im); 
                 }
                    setupi_(&i_one, imat, &im, &ic, &nm, &nc);
              } else {
                 if (get_option("debug"))  {
                    printf("false "); 
                 }
                    setupi_(&i_zero, imat, &im, &ic, &nm, &nc);
              }
              if (get_option("debug"))  {
                printf("Response:  %d %d %e %e %e \n",ic,im,respx, respy,fabs(respy));
              }
         }
              im++;
            }
            m = m->next;
         };
      }
      }
      ic++;
    }
    c = c->next;
  };
              if (get_option("debug"))  {
                 primat_(imat,&nm, &nc);
                 prdmat_(dmat,&nm, &nc);
                 printf("\n");
                 printf("\n");
              }
  return(dmat);
}


void pro_correct2_write_results(double *monvec, double *resvec, double *corvec, int *nx, int *nc, int *nm, int imon, int icor, int ip)
{
/*                                              */
/* Writes a summary of the correction           */
/* Writes correctors strengths into sequences   */
/* Fills TFS tables for correctors and monitors */
/* Fills the 'stren.out' output                 */
/* Makes various prints on request              */
/*                                              */
  int i;
  int rst;
  double corrm;
  struct id_mic2 *m, *c;      /* access to tables for monitors and correctors */

  m = correct_orbit12->mon_table;
  c = correct_orbit12->cor_table;

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
    pro_correct2_fill_mon_table(ip,m[nm[i]].p_node->name,monvec[i],resvec[i]);
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

  for(i=0;i<icor;i++) {   /* loop over all correctors */

    c[nc[i]].val.after[ip-1] = corvec[nx[i]-1];
    if(print_correct_opt > 1) {
      printf("%s %-3.6f %-3.6f %-3.6f\n",c[nc[i]].p_node->name,
                                         c[nc[i]].val.before[ip-1],
                                         corvec[nx[i]-1]+c[nc[i]].val.before[ip-1], 
                                         corvec[nx[i]-1]);
    }

    if(ip == 1) {
      /* Fill horizontal corrections for beam 1  */
      if(c[nc[i]].id_ttb[0] > 0) {
         c[nc[i]].p_node_s1->chkick += c[nc[i]].p_node_s1->other_bv*0.001*corvec[nx[i]-1];
         pro_correct2_fill_corr_table(0,
                                      ip,
                                      c[nc[i]].p_node->name,
                                      c[nc[i]].val.before[ip-1]*0.001,
                                      c[nc[i]].p_node_s1->chkick);
/*                                    c[nc[i]].p_node_s1->other_bv*0.001*corvec[nx[i]-1]); */
         if(fcdata != NULL) {
           fprintf(fcdata,"[1] %s = %e;\n",c[nc[i]].p_node->name,c[nc[i]].p_node_s1->other_bv*0.001*corvec[nx[i]-1]);
         }
      }
      /* Fill horizontal corrections for beam 2  */
      if(c[nc[i]].id_ttb[1] > 0) {
         c[nc[i]].p_node_s2->chkick += 0.001*corvec[nx[i]-1];
         pro_correct2_fill_corr_table(1,
                                      ip,
                                      c[nc[i]].p_node->name,
                                      c[nc[i]].val.before[ip-1]*0.001,
                                      c[nc[i]].p_node_s2->chkick);
/*                                    c[nc[i]].p_node_s2->other_bv*0.001*corvec[nx[i]-1]); */
         if(fcdata != NULL) {
           fprintf(fcdata,"[2] %s = %e;\n",c[nc[i]].p_node->name,0.001*corvec[nx[i]-1]);
         }
      }

    } else if (ip == 2) {
      /* Fill vertical corrections for beam 1  */
      if(c[nc[i]].id_ttb[0] > 0) {
         c[nc[i]].p_node_s1->cvkick += c[nc[i]].p_node_s1->other_bv*0.001*corvec[nx[i]-1];
         pro_correct2_fill_corr_table(0,
                                      ip,
                                      c[nc[i]].p_node->name,
                                      c[nc[i]].val.before[ip-1]*0.001,
                                      c[nc[i]].p_node_s1->cvkick);
/*                                    c[nc[i]].p_node_s1->other_bv*0.001*corvec[nx[i]-1]); */
         if(fcdata != NULL) {
           fprintf(fcdata,"[1] %s = %e;\n",c[nc[i]].p_node->name,c[nc[i]].p_node_s1->other_bv*0.001*corvec[nx[i]-1]);
         }
      }
      if(c[nc[i]].id_ttb[1] > 0) {
      /* Fill vertical corrections for beam 2  */
         c[nc[i]].p_node_s2->cvkick += 0.001*corvec[nx[i]-1];
         pro_correct2_fill_corr_table(1,
                                      ip,
                                      c[nc[i]].p_node->name,
                                      c[nc[i]].val.before[ip-1]*0.001,
                                      c[nc[i]].p_node_s2->cvkick);
/*                                    c[nc[i]].p_node_s2->other_bv*0.001*corvec[nx[i]-1]); */
         if(fcdata != NULL) {
           fprintf(fcdata,"[2] %s = %e;\n",c[nc[i]].p_node->name,0.001*corvec[nx[i]-1]);
         }
      }
    }

  }  /* end loop over correctors */
}


void correct_correct1(struct in_cmd* cmd)
/* Steering routine for orbit corrections of one beam */
{
  char rout_name[] = "correct_correct";
  int ix, im, ip, it, idrop;
  int j,err,nnnseq;
  int imon, icor;
  int ncorr, nmon;
  int niter;
  int resout;
  int twism;
  int dbg = 0;
  int ifail, sflag, svdflg;
  float  rms;
  double sngcut, sngval;
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
  im = pro_correct_gettables(ip, cmd);
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
  /* if flag "extern" is true: can be from external table */
  if(command_par_value("extern",cmd->clone)) {
     it = pro_correct_getorbit_ext(cmd);
  } else {
     it = pro_correct_getorbit(cmd);
  }


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
    if(dmat != NULL) myfree(rout_name,dmat);
    /* icor and imon used to set up correct matrix size !! */
    dmat = (double *)pro_correct_response_ring(ip,icor,imon);
    if((svdflg = command_par_value("cond",cmd->clone)) == 1) {
       sngcut = command_par_value("sngcut", cmd->clone);
       sngval = command_par_value("sngval", cmd->clone);
       printf("SVD conditioning requested ...\n");
       if (get_option("debug")) printf("Conditioning parameters: %e %e\n",sngcut, sngval);

       /* printf("Time before svd-comd:  %-6.3f\n",fextim());    */
       sflag=c_svddec(dmat,imon,icor,sing,&sngcut,&sngval);
       printf("Initially found %d singular values\n",sflag);
       /* printf("Time after svd-cond:  %-6.3f\n",fextim());     */
       /* printf("sflag: %d\n",sflag); */
        printf("sflag: %d\n",sflag); 
       for(ix=0;ix<sflag;ix++) {
         corl[nx[sing[2*ix+0]]].enable = 0;
           printf("Removed:   %d %s\n",nx[sing[2*ix+0]],
                 corl[nx[sing[2*ix+0]]].p_node->name);
         if(dbg == 1) {
           printf("Removed:   %d %s\n",nx[sing[2*ix+0]],
                 corl[nx[sing[2*ix+0]]].p_node->name);
         }
       }
       ix = pro_correct_getactive(ip, nm, nx, nc, corvec, monvec, conm);
       icor = ix%10000; imon  = ix/10000;
       printf("After SVD conditioning:             \n");
       printf("%d monitors and %d correctors enabled\n\n",imon,icor);
       if(dmat != NULL) myfree(rout_name,dmat);
       /* icor and imon used to set up correct matrix size !! */
       dmat = (double *)pro_correct_response_ring(ip,icor,imon);
       sflag=c_svddec(dmat,imon,icor,sing,&sngcut,&sngval);
       printf("Finally found %d singular values\n",sflag);
    }
  }
  else if(strcmp("line",command_par_string("flag",cmd->clone)) == 0) {
    if(dmat != NULL) myfree(rout_name,dmat);
          printf("make response for line\n");
    dmat = (double *)pro_correct_response_line(ip,icor,imon); 

    if((svdflg = command_par_value("cond",cmd->clone)) == 1) {
       sngcut = command_par_value("sngcut", cmd->clone);
       sngval = command_par_value("sngval", cmd->clone);
       printf("SVD conditioning requested ...\n");
       if (get_option("debug")) printf("Conditioning parameters: %e %e\n",sngcut, sngval);

       /* printf("Time before svd-comd:  %-6.3f\n",fextim());    */
       sflag=c_svddec(dmat,imon,icor,sing,&sngcut,&sngval);
       printf("Initially found %d singular values\n",sflag);
       /* printf("Time after svd-cond:  %-6.3f\n",fextim());     */
       /* printf("sflag: %d\n",sflag); */
       for(ix=0;ix<sflag;ix++) {
         corl[nx[sing[2*ix+0]]].enable = 0;
         if(dbg == 1) {
           printf("Removed:   %d %s\n",nx[sing[2*ix+0]],
                 corl[nx[sing[2*ix+0]]].p_node->name);
         }
       }
       ix = pro_correct_getactive(ip, nm, nx, nc, corvec, monvec, conm);
       icor = ix%10000; imon  = ix/10000;
       printf("After SVD conditioning:             \n");
       printf("%d monitors and %d correctors enabled\n\n",imon,icor);
       if(dmat != NULL) myfree(rout_name,dmat);
       /* icor and imon used to set up correct matrix size !! */
       dmat = (double *)pro_correct_response_ring(ip,icor,imon);
       sflag=c_svddec(dmat,imon,icor,sing,&sngcut,&sngval);
       printf("Finally found %d singular values\n",sflag);
    }
  }
 
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
          printf("Requested %d correctors (\?\?\?) set to %d\n",niter,icor);
          niter = icor;
    }
    else if((niter = command_par_value("ncorr",cmd->clone)) < 0) {
          printf("Requested %d correctors (\?\?\?) set to 0\n",niter);
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
       warning("MICADO back with error",", no correction done");
    }
if(ifail == 0) {
    pro_correct_write_results(monvec, resvec, corvec, nx, nc, nm, imon, icor, ip);
 }
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
  myfree(rout_name,nm);myfree(rout_name,dmat);myfree(rout_name,nx);
  myfree(rout_name,nc);myfree(rout_name,corvec);
  myfree(rout_name,monvec);myfree(rout_name,resvec); myfree(rout_name,conm);
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

int  pro_correct_gettables(int iplane, struct in_cmd* cmd)
{

  char rout_name[] = "pro_correct_gettables";

  struct id_mic *cor_l;
  struct id_mic *mon_l;

  struct table *ttb;

  char* orbtab;
  char* tartab;
  char* modtab;

  int j;
  int pps, ppt;
  int set0;
  int cntm = {0};
  int cntc = {0};

  double ounits;

/*
  static char atm[6][4] = {"hmon","vmon","moni","hkic","vkic","kick"};
*/

/* Get access to tables, for orbit and model the default is twiss_table */

  if ((orbtab = command_par_string("orbit",cmd->clone)) != NULL) {
      printf("Want to use orbit from: %s\n",orbtab);
    if ((pps = name_list_pos(orbtab, table_register->names)) > -1) {
       orbin_table = table_register->tables[pps];
    } else {
       fatal_error("ORBIT table for correction requested, but not provided:",orbtab);
    }
  } else {
       if((orbin_table = twiss_table) == NULL) {
         printf("FATAL ERROR:\n");
         printf("You request the ORBIT from a non-existing TWISS table\n");
         printf("You MUST run TWISS before trying to correct the orbit\n");
         printf("MAD-X stops\n");
         exit(81);
       } else {
         if (get_option("debug")) {
            printf("TWISS table: %ld\n",(long int)twiss_table);
         }
       }
       pps = -1;
  }


  if ((tartab = command_par_string("target",cmd->clone)) != NULL) {
       printf("Want to use target orbit from: %s\n",tartab);
    if ((ppt = name_list_pos(tartab, table_register->names)) > -1) {
       target_table = table_register->tables[ppt];
    } else {
       fatal_error("TARGET table for correction requested, but not provided:",tartab);
    }
  } else {
       if (get_option("debug")) {
         printf("No target orbit requested\n");
       }
       ppt = -1;
  }

  if ((modtab = command_par_string("model",cmd->clone)) != NULL) {
       printf("Want to use model orbit from: %s\n",modtab);
    if ((ppt = name_list_pos(modtab, table_register->names)) > -1) {
       model_table = table_register->tables[ppt];
    } else {
       fatal_error("MODEL table for correction requested, but not provided:",modtab);
    }
  } else {
       if((model_table = twiss_table) == NULL) {
         printf("FATAL ERROR:\n");
         printf("You request the MODEL from a non-existing TWISS table\n");
         printf("You MUST run TWISS before trying to correct the orbit\n");
         printf("MAD-X stops\n");
         exit(81);
       } else {
         if (get_option("debug")) {
            printf("TWISS table: %ld\n",(long int)twiss_table);
         }
       }
       ppt = -1;
  }



       if (get_option("debug")) {
            printf("The tables are: %ld %ld %ld %ld\n",
                   (long int) orbin_table,(long int)twiss_table,(long int)target_table,(long int)model_table);
       }
       if (get_option("debug")) {
       }

  if(correct_orbit == NULL) {
    correct_orbit = (struct orb_cor*)mycalloc("pro_correct_gettables",1, sizeof(struct orb_cor));
  }
       if (get_option("debug")) {
           printf("-0-\n");
       }

/*    if(corr_table == NULL) {   */
    corr_table = make_table("corr", "corr", corr_table_cols,
            corr_table_types, 5000);
       if (get_option("debug")) {
           printf("-01-\n");
       }
    add_to_table_list(corr_table, table_register);
       if (get_option("debug")) {
           printf("-02-\n");
       }
    pro_correct_make_corr_table();
       if (get_option("debug")) {
           printf("-03-\n");
       }
/*    }                          */
       if (get_option("debug")) {
           printf("-1-\n");
       }

/*    if(mon_table == NULL) { */
    mon_table = make_table("mon", "mon", mon_table_cols,
           mon_table_types, 5000);
       if (get_option("debug")) {
           printf("-11-\n");
       }
    add_to_table_list(mon_table, table_register);
       if (get_option("debug")) {
           printf("-12-\n");
       }
    pro_correct_make_mon_table();
       if (get_option("debug")) {
           printf("-13-\n");
       }
/*    }                       */
       if (get_option("debug")) {
           printf("-2-\n");
       }


  if(correct_orbit->cor_table != NULL) myfree(rout_name,correct_orbit->cor_table);
  if(correct_orbit->mon_table != NULL) myfree(rout_name,correct_orbit->mon_table);
  correct_orbit->cor_table = (struct id_mic *)mycalloc("pro_correct_gettables_cor",5200, sizeof(struct id_mic));
  correct_orbit->mon_table = (struct id_mic *)mycalloc("pro_correct_gettables_mon",5200, sizeof(struct id_mic));

/* orbit table available, get units, if defined */
   if((ounits = command_par_value("units",cmd->clone)) > 0) {
          correct_orbit->units=ounits;                      
   } else {
          correct_orbit->units=1.0;      
   }

  ttb = model_table;
  correct_orbit->mon_table->previous = NULL;
  correct_orbit->mon_table->next = NULL;
  correct_orbit->cor_table->previous = NULL;
  correct_orbit->cor_table->next = NULL;
       if (get_option("debug")) {
           printf("-3-\n");
       }

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
    if((strncmp(atc[iplane-1],ttb->p_nodes[j]->base_name,4) == 0) ||
       (strncmp(atc[2],       ttb->p_nodes[j]->base_name,4) == 0))  {
      cor_l->id_ttb = j;
      cor_l->enable = ttb->p_nodes[j]->enable;
      cor_l->p_node = ttb->p_nodes[j];
      if((set0 = command_par_value("corzero",cmd->clone)) > 0) {
        if(iplane == 1) cor_l->p_node->chkick = 0.0;
        if(iplane == 2) cor_l->p_node->cvkick = 0.0;
      }
      cor_l->next = cor_l;
      cor_l->next++; cor_l++;
      cntc++;
    }
  }
       if (get_option("debug")) {
           printf("-4-\n");
       }
  mon_l--; mon_l->next = NULL;
  cor_l--; cor_l->next = NULL;
       if (get_option("debug")) {
         printf("done: %d %d\n",cntm,cntc);
       }

  return(10000*cntm + cntc);
}


int pro_correct_getorbit(struct in_cmd* cmd)
{
  struct name_list* nl;
  int i;
  char *tartab;
  int pos;
  struct id_mic *m;  /* access to tables for monitors and correctors */
  struct table *ttb;
  struct table *tar = NULL;
  double **da1;
  double **da2 = NULL;
  double xlimit;
  double rx, ry, dpsi;

  char   strx[40];
  char   stry[40];

  int    posx, posy, pospx, pospy;
  int    tosx = -1;
  int    tosy = -1;
  int    tospx, tospy;

  ttb = orbin_table;
  da1 = ttb->d_cols;

  if(target_table != NULL) {
     tar = target_table;
     da2 = tar->d_cols;
  }

  nl = cmd->clone->par_names;

  m = correct_orbit->mon_table;

  strcpy(strx,"x");
  strcpy(stry,"y");
                                                                                                          
  if((posx = name_list_pos(strx,ttb->columns)) < 0) { 
      fatal_error("orbit x not found in input table",", MAD-X terminates ");
  }
  if((posy = name_list_pos(stry,ttb->columns)) < 0) { 
      fatal_error("orbit y not found in input table",", MAD-X terminates ");
  }
  if (get_option("debug")) {
    if((pospx = name_list_pos("px",ttb->columns)) < 0) { 
        fatal_error("orbit px not found in input table",", MAD-X terminates ");
    }
    if((pospy = name_list_pos("py",ttb->columns)) < 0) { 
        fatal_error("orbit py not found in input table",", MAD-X terminates ");
    }
    printf("====c1===>  %d %d %d %d \n",posx,posy,pospx,pospy);
  }

  if ((tartab = command_par_string("target",cmd->clone)) != NULL) {
    if((tosx = name_list_pos("x",tar->columns)) < 0) { 
        fatal_error("target orbit x not found in table",", MAD-X terminates ");
    }
    if((tosy = name_list_pos("y",tar->columns)) < 0) { 
        fatal_error("target orbit y not found in table",", MAD-X terminates ");
    }
    if (get_option("debug")) {
      if((tospx = name_list_pos("px",tar->columns)) < 0) { 
          fatal_error("target orbit px not found in table",", MAD-X terminates ");
      }
      if((tospy = name_list_pos("py",tar->columns)) < 0) { 
          fatal_error("target orbit px not found in table",", MAD-X terminates ");
      }
      printf("====c1===>  %d %d %d %d \n",tosx,tosy,tospx,tospy);
    }
  }


  while(m) {

/* If correction to target orbit, subtract the wanted orbit ... */
  if ((tartab = command_par_string("target",cmd->clone)) != NULL) {
    m->val.before[0] =  da1[posx][m->id_ttb] - da2[tosx][m->id_ttb];
    m->val.before[1] =  da1[posy][m->id_ttb] - da2[tosy][m->id_ttb];
    m->val.before[0] = (da1[posx][m->id_ttb] - da2[tosx][m->id_ttb])*1000.*correct_orbit->units;
    m->val.before[1] = (da1[posy][m->id_ttb] - da2[tosy][m->id_ttb])*1000.*correct_orbit->units;
  } else {
    m->val.before[0] = da1[posx][m->id_ttb];
    m->val.before[1] = da1[posy][m->id_ttb];
    m->val.before[0] = da1[posx][m->id_ttb]*1000.*correct_orbit->units;
    m->val.before[1] = da1[posy][m->id_ttb]*1000.*correct_orbit->units;
  }

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
           printf("scales: %e %e\n",m->p_node->p_al_err->a[12],m->p_node->p_al_err->a[13]);
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
          dpsi = m->p_node->p_al_err->a[5];
          rx = m->val.before[0];
          ry = m->val.before[1];
          printf("\nA: %e %e %e\n",m->val.before[0],m->val.before[1],dpsi);
          m->val.before[0] =  rx * cos(dpsi) +  ry * sin(dpsi);
          m->val.before[1] = -rx * sin(dpsi) +  ry * cos(dpsi);
          printf("B: %e %e %e\n",m->val.before[0],m->val.before[1],dpsi);
          m->val.before[0] += m->p_node->p_al_err->a[6]*1000.;
          m->val.before[1] += m->p_node->p_al_err->a[7]*1000.;
          printf("C: %e %e %e\n",m->val.before[0],m->val.before[1],dpsi);
        }
      }
    }
    m = m->next;
  };
  i = 0;
  return(i);
}

int pro_correct_getorbit_ext(struct in_cmd* cmd)
{
  struct name_list* nl;
  int i;
  int j;
  char *tartab;
  int pos;
  struct id_mic *m;  /* access to tables for monitors and correctors */
  struct table *ttb;
  struct table *tar = NULL;
  double **da1;
  double **da2 = NULL;
  double xlimit;
  char     name[NAME_L];
  char   l1name[NAME_L];
  char   l2name[NAME_L];
  char   l3name[NAME_L];
  char   l4name[NAME_L];
  double rx, ry, dpsi;

  char   *nam_col;
  char   *x_col;
  char   *y_col;

  char   strx[40];
  char   stry[40];
  char   strn[40];

  int    posx, posy, pospx, pospy;
  int    tosx = -1;
  int    tosy = -1;
  int    tospx, tospy;

  int dbk = 0;
  int yok;

  int    jjx, jjy, jj;

  ttb = orbin_table;
  da1 = ttb->d_cols;

  if(target_table != NULL) {
     tar = target_table;
     da2 = tar->d_cols;
  }

  nl = cmd->clone->par_names;

  m = correct_orbit->mon_table;

  if ((x_col = command_par_string("x_col",cmd->clone)) != NULL) {
         printf("X orbit in column: %s\n",x_col);
         strcpy(strx,x_col);
  } else {
         strcpy(strx,"x");
  }
  if ((y_col = command_par_string("y_col",cmd->clone)) != NULL) {
         printf("y orbit in column: %s\n",y_col);
         strcpy(stry,y_col);
  } else {
         strcpy(stry,"y");
  }
  if ((nam_col = command_par_string("name_col",cmd->clone)) != NULL) {
         printf("names in column: %s\n",nam_col);
         strcpy(strn,"name");
  } else {
         strcpy(strn,"name");
  }
  
  if((posx = name_list_pos(strx,ttb->columns)) < 0) { 
      fatal_error("orbit x not found in input table",", MAD-X terminates ");
  }
  if((posy = name_list_pos(stry,ttb->columns)) < 0) { 
      fatal_error("orbit y not found in input table",", MAD-X terminates ");
  }
  if (get_option("debug")) {
    if((pospx = name_list_pos("px",ttb->columns)) < 0) { 
        warning("orbit px not found in input table",", MAD-X continues ");
    }
    if((pospy = name_list_pos("py",ttb->columns)) < 0) { 
        warning("orbit py not found in input table",", MAD-X continues ");
    }
    printf("====c1===>  %d %d %d %d \n",posx,posy,pospx,pospy);
  }

  if ((tartab = command_par_string("target",cmd->clone)) != NULL) {
    if((tosx = name_list_pos("x",tar->columns)) < 0) { 
        fatal_error("target orbit x not found in table",", MAD-X terminates ");
    }
    if((tosy = name_list_pos("y",tar->columns)) < 0) { 
        fatal_error("target orbit y not found in table",", MAD-X terminates ");
    }
    if (get_option("debug")) {
      if((tospx = name_list_pos("px",tar->columns)) < 0) { 
          warning("target orbit px not found in table",", MAD-X continues ");
      }
      if((tospy = name_list_pos("py",tar->columns)) < 0) { 
          warning("target orbit px not found in table",", MAD-X continues ");
      }
      printf("====c1===>  %d %d %d %d \n",tosx,tosy,tospx,tospy);
    }
  }


  if (get_option("debug")) {
    printf("Number in table: %d\n",ttb->curr);
    for (j=1; j < (ttb->curr)+1; j++) {
      i = str_from_tablet(ttb, "name", &j, name);
    }
  }

  jj = 0;
  while(m) {


    strcpy(l1name,m->p_node->name);
    stolower(l1name);
    strcpy(l2name,strip(l1name));
    supp_tb(l2name);

     if (dbk ==1) printf("monitor name: %s\n",l2name);

    jjx = -1;
    jjy = -1;
    jj++;
    yok = 0;

    for (j=1; j < (ttb->curr)+1; j++) {
      i = str_from_tablet(ttb, "name", &j, name);
      strcpy(l3name,name);
      stolower(l3name);
      strcpy(l4name,strip(l3name));
      supp_tb(l4name);
      if(strlen(l4name) == strlen(l2name)) {
         if(strncmp(l4name,l2name,strlen(l2name)) == 0) {
              jjx = j-1;
              jjy = jj-1;
              yok = 1;
           if(dbk ==1) printf("monitor names found: %s %s %d %d\n",l2name,l4name,strlen(l2name),yok);
             }
      }
    }
     if(dbk ==1) printf("jjx,jjy %d %d\n",jjx,jjy);

/* If correction to target orbit, subtract the wanted orbit ... */

    if((jjy >= 0) && (yok == 1)) { 
/*  if(jjx >= 0)  {  */
     if ((tartab = command_par_string("target",cmd->clone)) != NULL) {
if(dbk == 1) {
       printf("x ==> %d %d %e %e\n",jjx,m->id_ttb,da1[posx][jjx],da2[tosx][jjy]);
       printf("y ==> %e %e\n",da1[posy][jjx],da2[tosy][jjy]);
}
       m->val.before[0] =  da1[posx][jjx] - da2[tosx][jjy];
       m->val.before[1] =  da1[posy][jjx] - da2[tosy][jjy];
       m->val.before[0] = (da1[posx][jjx] - da2[tosx][jjy])*1000.*correct_orbit->units;
       m->val.before[1] = (da1[posy][jjx] - da2[tosy][jjy])*1000.*correct_orbit->units;
if(dbk == 1) {
       printf("bxy ==> %s %d %e %e\n",m->p_node->name,jjx,m->val.before[0],m->val.before[1]);
}
     } else {
if(dbk == 1) {
       printf("x ==> %e %e\n",da1[posx][jjx],da2[tosx][jjx]);
       printf("y ==> %e %e\n",da1[posy][jjx],da2[tosy][jjx]);
}
       m->val.before[0] = da1[posx][jjx];
       m->val.before[1] = da1[posy][jjx];
       m->val.before[0] = da1[posx][jjx]*1000.*correct_orbit->units;
       m->val.before[1] = da1[posy][jjx]*1000.*correct_orbit->units;
if(dbk == 1) {
       printf("bxy ==> %s %d %e %e\n",m->p_node->name,jjx,m->val.before[0],m->val.before[1]);
}
     }

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
           printf("scales: %e %e\n",m->p_node->p_al_err->a[12],m->p_node->p_al_err->a[13]);
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
          dpsi = m->p_node->p_al_err->a[5];
          rx = m->val.before[0];
          ry = m->val.before[1];
          printf("\nA: %e %e %e\n",m->val.before[0],m->val.before[1],dpsi);
          m->val.before[0] =  rx * cos(dpsi) +  ry * sin(dpsi);
          m->val.before[1] = -rx * sin(dpsi) +  ry * cos(dpsi);
          printf("B: %e %e %e\n",m->val.before[0],m->val.before[1],dpsi);
          m->val.before[0] += m->p_node->p_al_err->a[6]*1000.;
          m->val.before[1] += m->p_node->p_al_err->a[7]*1000.;
          printf("C: %e %e %e\n",m->val.before[0],m->val.before[1],dpsi);
        }
      }
     }
    } else { m->enable = 0;} /* Only enable monitors found in input */
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

  ttb = model_table;

  da1 = ttb->d_cols;

  c = correct_orbit->cor_table;
  while(c) {
/*  c->val.before[0] = da1[59][c->id_ttb]*1000.;   */
/*  c->val.before[1] = da1[60][c->id_ttb]*1000.;   */
    c->val.before[0] = c->p_node->chkick*1000.;
    c->val.before[1] = c->p_node->cvkick*1000.;

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

  ttb = model_table;

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
  ttb = model_table;

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
  struct id_mic *m;  /* access to tables for monitors */

  struct table *ttb;
  static char  pl[2] = "xy";
  double **da1;
  double bx_m = -9999.;
  double xsig;
  double xmea, ymea;
  double xn;

  ttb = model_table;
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
          printf("==> %-4.3f \n",xn);
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

  ttb = model_table;
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

  ttb = model_table;
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

/*
  static char atm[5][4] = {"hmon","vmon","hkic","vkic","kick"};
*/

/*
  ttb = orbin_table;
*/
  ttb = model_table;

  for (j=0; j < ttb->curr; j++) {
    if((strncmp(atc[0],ttb->p_nodes[j]->base_name,4) == 0) ||
       (strncmp(atc[1],ttb->p_nodes[j]->base_name,4) == 0) ||
       (strncmp(atc[2],ttb->p_nodes[j]->base_name,4) == 0))  {
      string_to_table("corr","name",ttb->p_nodes[j]->name);
      augment_count("corr");
    }
  }
}

void pro_correct2_make_corr_table()
{
  struct id_mic2 *ttb;

/*
  static char atm[5][4] = {"hmon","vmon","hkic","vkic","kick"};
*/

  ttb = correct_orbit12->cor_table;

  while(ttb != NULL) {
    if((strncmp(atc[0],ttb->p_node->base_name,4) == 0) ||
       (strncmp(atc[1],ttb->p_node->base_name,4) == 0) ||
       (strncmp(atc[2],ttb->p_node->base_name,4) == 0))  {
       if(ttb->id_ttb[0] > 0) {
         string_to_table("corr1","name",ttb->p_node->name);
         augment_count("corr1");
       }

       if(ttb->id_ttb[1] > 0) {
         string_to_table("corr2","name",ttb->p_node->name);
         augment_count("corr2");
       }
    }
    ttb = ttb->next;
  }
}

void pro_correct_make_mon_table()
{
  struct table *ttb;
  int j;

/*
  static char atm[3][4] = {"hmon","vmon","moni"};
*/

  ttb = model_table;

  for (j=0; j < ttb->curr; j++) {
    if((strncmp(atm[0],ttb->p_nodes[j]->base_name,4) == 0) ||
       (strncmp(atm[1],ttb->p_nodes[j]->base_name,4) == 0) ||
       (strncmp(atm[2],ttb->p_nodes[j]->base_name,4) == 0))  {
      string_to_table("mon","name",ttb->p_nodes[j]->name);
      augment_count("mon");
    }
  }
}

void pro_correct2_make_mon_table()
{
  struct id_mic2 *ttb;
/*
  static char atm[3][4] = {"hmon","vmon","moni"};
*/

  ttb = correct_orbit12->mon_table;

  while(ttb != NULL) {
    if((strncmp(atm[0],ttb->p_node->base_name,4) == 0) ||
       (strncmp(atm[1],ttb->p_node->base_name,4) == 0) ||
       (strncmp(atm[2],ttb->p_node->base_name,4) == 0))  {
      string_to_table("mon","name",ttb->p_node->name);
      augment_count("mon");
    }
    ttb = ttb->next;
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

void pro_correct2_fill_corr_table(int b, int ip ,char *name, double old, double new)
{

  struct table *cor = NULL;

  int j;
  
  long longB = b;
  
  if((b != 1) && (b != 0)) {
      fatal_error("Invalid beam requested:",(char *) longB);
  }

  if(b == 0) cor =  corr_table1;
  if(b == 1) cor =  corr_table2;

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
      mon->d_cols[ip][j] = old*0.001;
      mon->d_cols[ip+2][j] = new*0.001;
    }
  }
}

void pro_correct2_fill_mon_table(int ip ,char *name, double old, double new)
{
  struct table *mon;

  int j;

  mon =  mon_table;

  for (j=0; j < mon->curr; j++) {
    if(strcmp(name,mon->s_cols[0][j]) == 0) {
      mon->d_cols[ip][j] = old*0.001;
      mon->d_cols[ip+2][j] = new*0.001;
    }
  }
}

void pro_correct_write_results(double *monvec, double *resvec, double *corvec, int *nx, int *nc, int *nm, int imon, int icor, int ip)
/*                                              */
/* Writes a summary of the correction           */
/* Writes correctors strengths into sequences   */
/* Fills TFS tables for correctors and monitors */
/* Fills the 'stren.out' output                 */
/* Makes various prints on request              */
/*                                              */
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

  for(i=0;i<icor;i++) {  /* loop over all correctors */

    if(print_correct_opt > 1) {
      printf("%s %-3.6f %-3.6f %-3.6f\n",c[nc[i]].p_node->name,
                                         c[nc[i]].val.before[ip-1],
                                         corvec[nx[i]-1]+c[nc[i]].val.before[ip-1], 
                                         corvec[nx[i]-1]);
    }

    c[nc[i]].val.after[ip-1] = corvec[nx[i]-1];
    if(ip == 1) {
      c[nc[i]].p_node->chkick += c[nc[i]].p_node->other_bv*0.001*corvec[nx[i]-1];
      pro_correct_fill_corr_table(ip,
                                  c[nc[i]].p_node->name,
                                  c[nc[i]].val.before[ip-1]*0.001,
                                  c[nc[i]].p_node->chkick);
/*                                c[nc[i]].p_node->other_bv*0.001*corvec[nx[i]-1]); */
      if(fcdata != NULL) {
         fprintf(fcdata,"%s = %e;\n",c[nc[i]].p_node->name,c[nc[i]].p_node->other_bv*0.001*corvec[nx[i]-1]);
      }
    } else if (ip == 2) {
      c[nc[i]].p_node->cvkick += c[nc[i]].p_node->other_bv*0.001*corvec[nx[i]-1];
      pro_correct_fill_corr_table(ip,
                                  c[nc[i]].p_node->name,
                                  c[nc[i]].val.before[ip-1]*0.001,
                                  c[nc[i]].p_node->cvkick);
/*                                c[nc[i]].p_node->other_bv*0.001*corvec[nx[i]-1]); */
      if(fcdata != NULL) {
         fprintf(fcdata,"%s = %e;\n",c[nc[i]].p_node->name,c[nc[i]].p_node->other_bv*0.001*corvec[nx[i]-1]);
      }
    }
  } /* end of loop ove correctors */
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
  fill_orbit_table(orbit_table, orbin_table);
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

uintptr_t locf_(char *iadr)    
#define NADUPW 4   /* Number of ADdress Units Per Word */
#define LADUPW 2   /* Logarithm base 2 of ADdress Units Per Word */
{
	return( (uintptr_t) iadr >> LADUPW );
}

int c_micit(double *dmat,char *conm, double *monvec,double *corvec,double *resvec,int *nx,float rms,int imon,int icor,int niter)
{
  char rout_name[] = "c_micit";
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

  myfree(rout_name,ny); myfree(rout_name,ax); myfree(rout_name,cinx);
  myfree(rout_name,xinx); myfree(rout_name,resx); myfree(rout_name,rho);
  myfree(rout_name,ptop); myfree(rout_name,rmss); myfree(rout_name,xrms);
  myfree(rout_name,xptp); myfree(rout_name,xiter);
/*
*/

  return(ifail);
}

void c_haveit(double *dmat,double *monvec,double *corvec,double *resvec,int *nx,int imon,int icor)
{
  char rout_name[] = "c_haveit";
  double *cb,*xmeas,*xres,*y,*z,*xd;

  cb   =(double *)mycalloc("c_haveit_cb",icor,sizeof(double));
  xmeas=(double *)mycalloc("c_haveit_xmeas",imon,sizeof(double));
  xres =(double *)mycalloc("c_haveit_xres",imon,sizeof(double));
  y    =(double *)mycalloc("c_haveit_y",icor*imon,sizeof(double));
  z    =(double *)mycalloc("c_haveit_z",icor*icor,sizeof(double));
  xd   =(double *)mycalloc("c_haveit_xd",icor,sizeof(double));

  haveit_(dmat,monvec,corvec,resvec,nx,&imon,&icor,cb,xmeas,xres,y,z,xd);

  myfree(rout_name,cb); myfree(rout_name,xmeas); myfree(rout_name,xres);
  myfree(rout_name,y);  myfree(rout_name,z); myfree(rout_name,xd);

  return;
}

int  c_svddec(double *dmat, int imon, int icor, int *sing, double *sngcut, double *sngval)

{
  char rout_name[] = "c_svddev";
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

  if(imon >= icor ) {
      svddec_m_(dmat,s,u,v,w,ut,vt,wt,ws,wv,sw,sngcut,sngval,&imon,&icor,&flag,sing,&dbg);
  } else {
      svddec_c_(dmat,s,u,v,w,ut,vt,wt,ws,wv,sw,sngcut,sngval,&imon,&icor,&flag,sing,&dbg);
  }
  myfree(rout_name,s); myfree(rout_name,u);
  myfree(rout_name,v); myfree(rout_name,w);
  myfree(rout_name,ut); myfree(rout_name,vt);
  myfree(rout_name,wt); myfree(rout_name,ws);
  myfree(rout_name,wv); myfree(rout_name,sw);

  return(flag);
}

int  c_svdcorr(double *dmat, double *xin, double *cor, double *res, int *nx, int imon, int icor)
{
  char rout_name[] = "c_svdcorr";
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

  if(imon >= icor ) {
      svdcorr_m_(dmat,s,u,v,w,ut,vt,wt,xin,cor,res,
                 xa,xb,xp,ws,wv,sw,
                 nx,&imon,&icor,&flag,&dbg);
  } else {
      svdcorr_c_(dmat,s,u,v,w,ut,vt,wt,xin,cor,res,
                 xa,xb,xp,ws,wv,sw,
                 nx,&imon,&icor,&flag,&dbg);
  }
  myfree(rout_name,s); myfree(rout_name,u); myfree(rout_name,v);
  myfree(rout_name,w); myfree(rout_name,ut); myfree(rout_name,vt);
  myfree(rout_name,wt); myfree(rout_name,sw); myfree(rout_name,xa);
  myfree(rout_name,xb); myfree(rout_name,xp); myfree(rout_name,ws);
  myfree(rout_name,wv);

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
   #if ( __GNUC__==2 && __GNUC_MINOR__ > 94 ) || __GNUC__ > 2 /*hbu  gettimeofday available */
     float mytime;
     struct timeval tp;
     struct timezone tzp;
     gettimeofday(&tp,&tzp);
     mytime = (float)(tp.tv_sec%10000) + 1.e-6 * tp.tv_usec; /* seconds from epoch, modulo 10 000 */
   #else /* use old ftime */
     struct timeb tp;
     float mytime;

     ftime(&tp);

     mytime = (float)(tp.time%10000) + 0.001*tp.millitm;
   #endif
     /* printf("Time now:  %-6.3f\n",mytime);    */

     return(mytime);
}

int str_from_table(char* table, char* name, int* row, char* val)
     /* WH 22.06.2004, corrected from: char_from_table */
     /* returns val at position row in column with name "name".
        function value return:
        0  OK
        -1 table  does not exist
        -2 column does not exist
        -3 row    does not exist
     */
{
  int pos;
  struct table* t;

  strcpy(val,"No-Name");
  mycpy(c_dum->c, table);
  if ((pos = name_list_pos(c_dum->c, table_register->names)) > -1)
    t = table_register->tables[pos];
  else return -1;
  mycpy(c_dum->c, name);
  if ((pos = name_list_pos(c_dum->c, t->columns)) < 0) return -2;
  if (*row > t->curr)  return -3;
   strncpy(val,t->s_cols[pos][*row-1],NAME_L);
  while (strlen(val)<NAME_L) val[strlen(val)]=' ';
  val[NAME_L-1] = '\0';
  return 0;
}


int str_from_tablet(struct table *t, char* name, int* row, char* val)
     /* WH 22.06.2004, corrected from: char_from_table */
     /* returns val at position row in column with name "name".
        function value return:
        0  OK
        -1 table  does not exist
        -2 column does not exist
        -3 row    does not exist
     */
{
  int pos;

  strcpy(val,"No-Name");
  mycpy(c_dum->c, name);
  if ((pos = name_list_pos(c_dum->c, t->columns)) < 0) return -2;
  if (*row > t->curr)  return -3;
   strncpy(val,t->s_cols[pos][*row-1],NAME_L);
  while (strlen(val)<NAME_L) val[strlen(val)]=' ';
  val[NAME_L-1] = '\0';
  return 0;
}


struct table* read_my_table(struct in_cmd* cmd)
     /* reads and stores TFS table */
{
  struct table* t = NULL;
  struct char_p_array* tcpa = NULL;
  struct name_list* tnl = NULL;
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int pos = name_list_pos("file", nl);
  int i, k, error = 0;
  short  sk;
  char *cc, *filename, *type = NULL, *tmp, *name;

  char* namtab;

  if ((namtab = command_par_string("table",cmd->clone)) != NULL) {
       printf("Want to make named table: %s\n",namtab);
  } else {
       if (get_option("debug")) {
         printf("No table name requested\n");
         printf("Use default name (i.e. name from file) \n");
       }
       namtab = NULL;
  }

  if(nl->inform[pos] && (filename = pl->parameters[pos]->string) != NULL)
    {
     if ((tab_file = fopen(filename, "r")) == NULL)
       {
         fatal_error("cannot open file:", filename); return NULL; /* frs: to avoid unwanted results */
       }
    }
  else
    {
     warning("no filename,","ignored"); return NULL;
    }
  while (fgets(aux_buff->c, aux_buff->max, tab_file))
    {
     cc = strtok(aux_buff->c, " \"\n");
     if (*cc == '@')
       {
       if ((tmp = strtok(NULL, " \"\n")) != NULL
              && strcmp(tmp, "TYPE") == 0)
        {
         if ((name = strtok(NULL, " \"\n")) != NULL) /* skip format */
           {
            if ((name = strtok(NULL, " \"\n")) != NULL)
                  type = permbuff(stolower(name));
           }
        }
       }
     else if (*cc == '*' && tnl == NULL)
       {
      tnl = new_name_list("table_names", 20);
        while ((tmp = strtok(NULL, " \"\n")) != NULL)
            add_to_name_list(permbuff(stolower(tmp)), 0, tnl);
       }
     else if (*cc == '$' && tcpa == NULL)
       {
      if (tnl == NULL)
        {
         warning("formats before names","skipped"); return NULL;
        }
      tcpa = new_char_p_array(20);
        while ((tmp = strtok(NULL, " \"\n")) != NULL)
        {
         if (tcpa->curr == tcpa->max) grow_char_p_array(tcpa);
           if (strcmp(tmp, "%s") == 0)       tnl->inform[tcpa->curr] = 3;
           else if (strcmp(tmp, "%hd") == 0) tnl->inform[tcpa->curr] = 1;
           else if (strcmp(tmp, "%d") == 0)  tnl->inform[tcpa->curr] = 1;
           else                              tnl->inform[tcpa->curr] = 2;
           tcpa->p[tcpa->curr++] = permbuff(tmp);
        }
       }
     else
       {
        if(t == NULL)
          {
         if (type == NULL)
           {
            warning("TFS table without type,","skipped"); error = 1;
           }
         else if (tcpa == NULL)
           {
            warning("TFS table without formats,","skipped"); error = 1;
           }
         else if (tnl == NULL)
           {
            warning("TFS table without column names,","skipped"); error = 1;
           }
         else if (tnl->curr == 0)
           {
            warning("TFS table: empty column name list,","skipped");
              error = 1;
           }
         else if (tnl->curr != tcpa->curr)
           {
            warning("TFS table: number of names and formats differ,",
                       "skipped");
              error = 1;
           }
           if (error)
           {
            delete_name_list(tnl); return NULL;
           }
           if(namtab != NULL) {
             t = new_table(namtab, type,    500, tnl);
           } else {
             t = new_table(type, type,    500, tnl);
           }
        }
      for (i = 0; i < tnl->curr; i++)
        {
         if (t->curr == t->max) grow_table(t);
         tmp = tcpa->p[i];
           if (strcmp(tmp,"%s") == 0) t->s_cols[i][t->curr] = stolower(tmpbuff(cc));
           else if (strcmp(tmp,"%d") == 0 )
           {
            sscanf(cc, tmp, &k); t->d_cols[i][t->curr] = k;
           }
           else if (strcmp(tmp,"%hd") == 0 )
           {
            sscanf(cc, tmp, &sk); t->d_cols[i][t->curr] = sk;
           }
           else sscanf(cc, tmp, &t->d_cols[i][t->curr]);
           if (i+1 < tnl->curr)
           {
              if ((cc =strtok(NULL, " \"\n")) == NULL)
              {
               warning("incomplete table line starting with:", aux_buff->c);
                 return NULL;
              }
           }
        }
        t->curr++;
       }
    }
  fclose(tab_file);
  t->origin = 1;
  add_to_table_list(t, table_register);
  return NULL;
}

void correct_setcorr(struct in_cmd* cmd)
{

/* read the correctors from named table  and stores
   them in the nodes of the sequence at 
   "->chkick" and "->cvkick". Subsequent Twiss will
   use them correctly.
   ===> Must be preceded by a call to "read_table"
   ===> (unless table exists in memory !)
   ===> Watch out, does not yet take care of existing corrector
   ===> settings already present in sequence
   ===> Uses table with name specified in parameter: table=
*/

  int i, ix;

  struct node *ndexe;
  struct node *nextnode;

  char     name[NAME_L];
  char   slname[NAME_L];

  char     nname[NAME_L];
  char   slnname[NAME_L];

  char*    namtab;
  int      t1;

  double   xnew, ynew;

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
          /* fatal_error("Corrector table requested, but not existing:",namtab); */
          /* exit(-77); */ 
          printf("No such corrector table in memory: %s\n",namtab);
       }

  } else {
       if (get_option("debug")) {
         printf("No table name requested\n");
         printf("Use default name\n");
       }
       strcpy(namtab,"corr");
  }


  i = 1; ix=0;
  while(ix == 0) {
      ix =   str_from_table(namtab, "name", &i, name);
      ix = double_from_table(namtab, "px.correction", &i, &xnew);
      ix = double_from_table(namtab, "py.correction", &i, &ynew);
      if(ix == 0) {
             stolower(name);
             strcpy(slname,strip(name));
             supp_tb(slname);

          /* printf("corrs: %s %d %e %e %e %e\n",name,ix,xold,yold,xnew,ynew); */   
          nextnode = mysequ->ex_start;
          while (nextnode != ndexe) {
             stolower(name);
             strcpy(slname,strip(name));
             supp_tb(slname);
            
             strcpy(nname,nextnode->name);
             stolower(nname);
             strcpy(slnname,strip(nname));
             supp_tb(slnname);
          
             /* printf("seq and input (0): %s %d %s %d\n", nname,strlen(nname),  name,strlen(name));
             printf("seq d in (2): %s %d %s %d\n",slnname,strlen(slnname),slname,strlen(slname)); */
       
             if(strcmp(slname,slnname) == 0) {
                /*
                printf("Corrector selection found: %s, %s %d\n",lname,nextnode->name,nextnode->sel_err);
                printf("corrs: %s %d %e %e %e %e\n",name,ix,xold,yold,xnew,ynew);    
                printf("corrs in sequence: %s %e %e\n",nextnode->name,nextnode->chkick,nextnode->cvkick);
                */
                nextnode->chkick += xnew;
                nextnode->cvkick += ynew;
                /*
                printf("corrs in sequence: %s %e %e\n",nextnode->name,nextnode->chkick,nextnode->cvkick);
                */
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


void correct_readcorr(struct in_cmd* cmd)
{

/* read the correctors from table "corr" and stores
   them in the nodes of the sequence at 
   "->chkick" and "->cvkick". Subsequent Twiss will
   use them correctly.
   ===> Must be preceded by a call to "read_table"
   ===> Watch out, does not yet take care of existing corrector
   ===> settings already present in sequence
   ===> Always uses table with name "corr", will change ...
*/

  int i, ix;

  struct node *ndexe;
  struct node *nextnode;

  char     name[NAME_L];
  char    lname[NAME_L];
  char   slname[NAME_L];
  char* uslname;

  char     nname[NAME_L];
  char    lnname[NAME_L];
  char   slnname[NAME_L];
  char* uslnname;

  double   xnew, ynew;

/* set up pointers to current sequence for later use */
  struct sequence* mysequ = current_sequ;
  nextnode = mysequ->ex_start;
  ndexe = mysequ->ex_end;

/* printf("Pointers: %d %d %d\n",mysequ,nextnode,ndexe); */


  i = 1; ix=0;
  while(ix == 0) {
      ix =   str_from_table("corr", "name", &i, name);
      ix = double_from_table("corr", "px.correction", &i, &xnew);
      ix = double_from_table("corr", "py.correction", &i, &ynew);
      if(ix == 0) {
          /* printf("corrs: %s %d %e %e %e %e\n",name,ix,xold,yold,xnew,ynew); */   
          nextnode = mysequ->ex_start;
          while (nextnode != ndexe) {
             strcpy(lname,name);
             stolower(lname);
             strcpy(slname,strip(lname));
             uslname = supp_tb(slname);
            
             strcpy(nname,nextnode->name);
             strcpy(lnname,nname);
             stolower(lnname);
             strcpy(slnname,strip(lnname));
             uslnname = supp_tb(slnname);
          
             /* printf("seq and input (0): %s %d %s %d\n", nname,strlen(nname),  name,strlen(name));
             printf("seq d in (1): %s %d %s %d\n",lnname,strlen(lnname),lname,strlen(lname));
             printf("seq d in (2): %s %d %s %d\n",slnname,strlen(slnname),slname,strlen(slname));
             printf("seq d in (3): %s %d %s %d\n",uslnname,strlen(uslnname),uslname,strlen(uslname)); 
             printf("compare: %s %d %s %d \n",uslname,strlen(uslname),uslnname,strlen(uslnname)); */
       
             if(strcmp(uslname,uslnname) == 0) {
                /*
                printf("Corrector selection found: %s, %s %d\n",lname,nextnode->name,nextnode->sel_err);
                printf("corrs: %s %d %e %e %e %e\n",name,ix,xold,yold,xnew,ynew);    
                printf("corrs in sequence: %s %e %e\n",nextnode->name,nextnode->chkick,nextnode->cvkick);
                */
                nextnode->chkick += xnew;
                nextnode->cvkick += ynew;
                /*
                printf("corrs in sequence: %s %e %e\n",nextnode->name,nextnode->chkick,nextnode->cvkick);
                */
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


