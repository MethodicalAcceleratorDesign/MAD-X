#include "madx.h"

#ifdef _ONLINE

#include "SDDS.h"

static int
head_split(char* buf, struct char_p_array* list)
{
  /* splits header information into tokens */
  int j;
  char* p;
  if ((p =strtok(buf, " \n")) == NULL) return 0;
  list->curr = 0;
  list->p[list->curr++] = p;
  while ((p = strtok(NULL, " \n")) != NULL)
  {
    if (list->curr == list->max) grow_char_p_array(list);
    list->p[list->curr++] = p;
  }
  /* remove '@' in strings */
  for (j = 0; j < list->curr; j++)
    if(*list->p[j] == '\"' || *list->p[j] == '\'') /* quote */
      replace(list->p[j], '@', ' ');
  return list->curr;
}
                                                                                                     
static int
sdds_get_parm(SDDS_TABLE *SDDS_table, struct table *tfs_table)
{
  PARAMETER_DEFINITION *pardef;
  void    *parval;
  long    *parvall;
  short   *parvals;
  double  *parvald;
  float   *parvalf;
  char    **parvalstr;

  int     h_length;

  char    **cpar;

  int   i2, npar;

      /* access parameter data  SDDS_GetParameter(...) */
      /*   */
      cpar = SDDS_GetParameterNames(SDDS_table,(int32_t *)&npar);
      h_length = npar;
       if (tfs_table->header == NULL)  tfs_table->header = new_char_p_array(h_length);
      for(i2=0;i2<npar;i2++) {
         pardef = SDDS_GetParameterDefinition(SDDS_table,cpar[i2]);

         parval = NULL;
         parval = SDDS_GetParameter(SDDS_table,cpar[i2],NULL);

         if(pardef->type == SDDS_LONG) {
           parvall = (long *)parval;
           if (get_option("debug")) printf("Parameter: %s, value: %ld\n",cpar[i2],*parvall);
           table_add_header(tfs_table, "@ %-16s %%ld  %ld", cpar[i2],*parvall);
         }
         if(pardef->type == SDDS_SHORT) {
           parvals = (short *)parval;
           if (get_option("debug")) printf("Parameter: %s, value: %ld\n",cpar[i2],(long)*parvals);
           table_add_header(tfs_table, "@ %-16s %%ld  %ld", cpar[i2],(long)*parvals);
         }
         if(pardef->type == SDDS_FLOAT) {
           parvalf = (float *)parval;
           if (get_option("debug")) printf("Parameter: %s, value: %e\n",cpar[i2],(double)*parvalf);
           table_add_header(tfs_table, "@ %-16s %%e  %e", cpar[i2],(double)*parvalf);
         }
         if(pardef->type == SDDS_DOUBLE) {
           parvald = (double *)parval;
           if (get_option("debug")) printf("Parameter: %s, value: %e\n",cpar[i2],*parvald);
           table_add_header(tfs_table, "@ %-16s %%le  %le", cpar[i2],*parvald);
         }
         if(pardef->type == SDDS_STRING) {
           parvalstr = SDDS_GetParameter(SDDS_table,cpar[i2],NULL);
           if (get_option("debug")) printf("Parameter: %s, value: %s \n",cpar[i2],*parvalstr);
           table_add_header(tfs_table, "@ %-16s %%%02ds \"%s\"", cpar[i2],strlen(*parvalstr),*parvalstr);
         }

      }

  return(npar);
}

static int
sdds_readt(char *filename, char *tfsname)
{
  int     i1, i2, i3; // i5; not used
  int     j1, j2;

//  long    lb; // not used
  long    narr, nall;
  long    arrdim[MAX_TFS_COL];
  long    arrtyp[MAX_TFS_COL];
  long    arrele[MAX_TFS_COL];

  double  *datmp;
  float   *datmpf;
  long    *datmpl;
  short   *datmps;
  double  *datd[MAX_TFS_COL];
  char    **datstr[MAX_TFS_COL];

  char    **c0;
  char    *c1[MAX_TFS_COL];

  SDDS_TABLE SDDS_table;

  SDDS_ARRAY *arr=0;

  char*   tfs_table_cols[MAX_TFS_COL];

  int     tfs_table_types[MAX_TFS_COL];

  struct  table  *tfstab;

  /* convert SDDS to TFS types */
  int tfs_sdds_types[] =
  {
    -1, MADX_DOUBLE, MADX_DOUBLE, MADX_LONG, MADX_LONG, MADX_LONG, MADX_LONG, MADX_STRING
  };

  if(sdds_pat != NULL) {
       for(i2=0; i2<sdds_pat->curr; i2++) {
          if (get_option("debug")) printf("my selected pattern: %d %s\n",i2,sdds_pat->p[i2]);
       }
  }

  /* open the file and read the SDDS header */ 
  if(SDDS_InitializeInput(&SDDS_table, filename) != 1) {
     SDDS_PrintErrors(stdout,SDDS_VERBOSE_PrintErrors);
  }

  /* read and process each data table in the data set */ 
  while (SDDS_ReadTable(&SDDS_table) > 0) { // (lb = not used
      /* set all rows and all columns to initially be "of interest" */ 
      SDDS_SetColumnFlags(&SDDS_table, 1); 
      SDDS_SetRowFlags(&SDDS_table, 1);

  }

      /* access array data  SDDS_GetArray(...) */
      /* */
      {
        int32_t nall32 = 0;
        c0 = SDDS_GetArrayNames(&SDDS_table,&nall32);
        nall = nall32;
      }

      if (get_option("debug")) printf("Found %ld arrays in total\n",nall);

      narr = 0;

        for(i1=0;i1<nall;i1++){ 
           if(sdds_pat != NULL) {
             if (get_option("debug")) printf("Check %d patterns\n",sdds_pat->curr);
             for(i2=0; i2 < sdds_pat->curr; i2++) {
               if (get_option("debug")) printf("Check pattern %s %s\n",sdds_pat->p[i2],c0[i1]);
               if(myregex(sdds_pat->p[i2],c0[i1]) == 0) {
                  if (get_option("debug")) printf("found now %s %s %ld\n",sdds_pat->p[i2],c0[i1],narr);
                  c1[narr++] = c0[i1];
                  if (get_option("debug")) printf("found now %s %s %ld\n",sdds_pat->p[i2],c0[i1],narr);
               }
             }
           } else {
                  if (get_option("debug")) printf("no check, use %s %ld\n",c0[i1],narr);
                  c1[narr++] = c0[i1];
           }
        }

      for(i1=0;i1<narr;i1++){ 
        if ((SDDS_CheckArray(&SDDS_table, c1[i1], NULL, 
                      0, stderr)) != SDDS_CHECK_OKAY) {
           fprintf(stderr, "array %s is not in the data file",c1[i1]); 
           exit(1); 
        } 

        arr = SDDS_GetArray(&SDDS_table, c1[i1], NULL);
        arrdim[i1] = (long)arr->definition->dimensions;
        arrtyp[i1] = (long)arr->definition->type;
        arrele[i1] = (long)arr->elements;          
                                                                                                            
  /* check whether type is valid and compatible with TFS */
        if((arrtyp[i1] < 1) || (arrtyp[i1] > SDDS_STRING))     {
           fatal_error("Type not valid or compatible for: ",c1[i1]);
        }
                                                                                                            
  /* check whether array dimension and size is compatible with TFS */
        if(arrdim[i1] != 1 )     {
           fatal_error("Array is 2-Dimensional: ",c1[i1]);
        }
        if(i1 > 0) {
           if (get_option("debug")) printf("===> %d %ld %ld\n",i1,arrele[i1],arrele[i1-1]);
           if(arrele[i1] != arrele[i1-1]) {
              printf("found two arrays with different lengths: %ld %ld\n",arrele[i1-1],arrele[i1]);
              fatal_error("Two arrays with different length, ", "conversion aborted");
           }
        }

        datd[i1] = mycalloc_atomic("double ptr buffer for SDDS", arrele[i1]+8, sizeof *datd[0]);
        if(arrtyp[i1] == SDDS_STRING) datstr[i1] = (char **)arr->data;

        for(i3=0;i3<arr->elements;i3++){
           datmp = (double *)arr->data;   
           datmpl = (long *)arr->data;   
           datmps = (short *)arr->data;   
           datmpf = (float *)arr->data;   
           if (get_option("debug")) {
             if(arrtyp[i1] != SDDS_STRING) printf("data: %e %e \n",(double)datmp[i3],datmp[i3]);
             if(arrtyp[i1] == SDDS_STRING) printf("data: %s \n",datstr[i1][i3]);
           }
           if(arrtyp[i1] == SDDS_DOUBLE) datd[i1][i3] = (double)datmp[i3];     
           /* for float data: convert to double float */
           if(arrtyp[i1] == SDDS_FLOAT)  datd[i1][i3] = datmpf[i3];     
           /* for integer data: convert to double float */
           if(arrtyp[i1] == SDDS_LONG)   datd[i1][i3] = datmpl[i3];     
           if(arrtyp[i1] == SDDS_SHORT)  datd[i1][i3] = datmps[i3];     
        }
        if (get_option("debug")) {
           printf("For TFS table: %s %ld %ld %d\n",c1[i1],arrdim[i1],arrtyp[i1],tfs_sdds_types[arrtyp[i1]]);
        }
        if(arrdim[i1] != 1 )     {
           warning("Array is 2-Dimensional",c1[i1]);
        }
  /*
  */
        tfs_table_cols[i1] = stolower(c1[i1]);
        tfs_table_types[i1] = tfs_sdds_types[arrtyp[i1]];
       }

       char c1_tmp[] = " ";
       c1[narr] = c1_tmp;
       tfs_table_cols[narr] = c1[narr];
  /*
       c1[3] = " ";
       tfs_table_cols[3] = c1[3];
       strcpy(&tfs_table_cols[narr]," ");
  */
        
       tfstab = make_table(tfsname, tfsname, (const char* const*)tfs_table_cols, tfs_table_types, 5000);
       add_to_table_list(tfstab, table_register);

       sdds_get_parm(&SDDS_table, tfstab); // i5 = not used

       if (get_option("debug")) {
           for(j1=0;j1<narr; j1++){
              for(j2=0;j2<arr->elements; j2++){
                 printf("data: %d %d %e\n",j1,j2,datd[j1][j2]);
              }
           }
       }

           for(i2=0;i2<arr->elements; i2++) {
                   for(i1=0;i1<narr;i1++) {
                     if((arrtyp[i1] >= 1) && (arrtyp[i1] < SDDS_STRING)) {
                        double_to_table_curr(tfsname,c1[i1],&datd[i1][i2]);
                     } else if(arrtyp[i1] == SDDS_STRING)                {
                        string_to_table_curr(tfsname,c1[i1],datstr[i1][i2]);
                     } else {
                        fatal_error("Type not valid or compatible for: ",c1[i1]);
                     }
                   }
  /*
  */
                   augment_count(tfsname);
           }
                                                                                                            


       if (get_option("debug")) printf("--> %d %d\n",tfstab->curr,tfstab->num_cols);

       if (get_option("debug")) out_table(tfsname,tfstab,"outtfs.1");

  /* free all allocated space ....  */
       for(i1=0; i1 < narr; i1++) {
         if(datd[i1] != NULL) myfree("free double array",datd[i1]);
            datd[i1] = NULL;
       }

  return(narr);
}

static int
treat_tfs_header_set(SDDS_TABLE *SDDS_table, struct table* t)
{
  struct char_p_array* head_buf;
  int j, k; // i, not used
  char  dumc[1000];

  double  dbuf;
  long    lbuf;

  head_buf = new_char_p_array(1000);

  printf("number of headers: %d\n",t->header->curr);
  for(j=0; j < t->header->curr; j++) {
    if (get_option("debug")) printf("for set header: %s\n", t->header->p[j]);
    if (get_option("debug")) printf("header: %s\n", t->header->p[j]);
    pre_split(t->header->p[j], l_wrk, 0);
    head_split(l_wrk->c,head_buf); // i = not used
    if (get_option("debug")) printf("for set curr: %d\n",head_buf->curr);
    if (get_option("debug")) printf("curr: %d\n",head_buf->curr);

    if(head_buf->curr > 0) {
      for(k=0; k < head_buf->curr; k++) {
        if (get_option("debug")) printf("for set header: %d %s ", k,  head_buf->p[k]);

        if(strcmp(head_buf->p[2],"%ld") == 0) {
          sscanf(head_buf->p[3],"%ld",&lbuf);
          if (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                               head_buf->p[1], lbuf, NULL)) {
             SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
             exit(1);
          }
        } else if(strcmp(head_buf->p[2],"%le") == 0) {
          sscanf(head_buf->p[3],"%le",&dbuf);
          if (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                               head_buf->p[1], dbuf, NULL)) {
             SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
             exit(1);
          }
        } else  {
          strcpy(dumc,head_buf->p[3]);
          replace(dumc, '\"', ' ');
          if (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,     
                               head_buf->p[1], dumc, NULL)) {
             SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
             exit(1);
          }
        }

      }
    }
  }

  return(head_buf->curr);
}

static int
treat_tfs_header_define(SDDS_TABLE *SDDS_table, struct table* t)
{
  struct char_p_array* head_buf;
  int j, k; // i, not used

  head_buf = new_char_p_array(1000);

  if (get_option("debug")) printf("number of headers: %d\n",t->header->curr);
  for(j=0; j < t->header->curr; j++) {
    if (get_option("debug")) printf("header: %s\n", t->header->p[j]);
    pre_split(t->header->p[j], l_wrk, 0);
    head_split(l_wrk->c,head_buf); // i = not used 
    if (get_option("debug")) printf("curr: %d\n",head_buf->curr);

    if(head_buf->curr > 0) {
      for(k=0; k < head_buf->curr; k++) {
  /*
        SDDS_DefineParameter(SDDS_table, head_buf->p[1], NULL, NULL, NULL, NULL, SDDS_STRING, head_buf->p[3]);
  */
        if(strcmp(head_buf->p[2],"%le") == 0) {
        SDDS_DefineParameter(SDDS_table, head_buf->p[1], NULL, NULL, NULL, NULL, SDDS_DOUBLE, NULL);
        } else if(strcmp(head_buf->p[2],"%ld") == 0) {
        SDDS_DefineParameter(SDDS_table, head_buf->p[1], NULL, NULL, NULL, NULL, SDDS_LONG, NULL);
        } else {
        SDDS_DefineParameter(SDDS_table, head_buf->p[1], NULL, NULL, NULL, NULL, SDDS_STRING, NULL);
        }
      }
    }
  }

  return(head_buf->curr);
}
                                                                                                    
static int
sdds_writet_sel(char *filename, struct table *tfstab)
{
//  int     i1;
  int     j1, j2;

  double  **da1;
  char    ***sa1;

  int     pos[1000];

  long    *pl;

  SDDS_TABLE SDDS_table;

  struct int_array* col = tfstab->col_out;
  struct int_array* row = tfstab->row_out;

  /* convert TFS to SDDS types */
  int sdds_tfs_types[] =
  {
    -1, SDDS_LONG, SDDS_DOUBLE, SDDS_STRING
  };
                                                                                                            
  /* set up to put data into file "atest.out" */
    if (!SDDS_InitializeOutput(&SDDS_table, SDDS_ASCII, 1, NULL, NULL, filename)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
     

  /*     */ 
  if (get_option("debug")) {
    printf("tfs table with %d columns \n",tfstab->num_cols);
    printf("tfs table with %d rows    \n",tfstab->curr    );
    printf("tfs table with %d selected cols    \n",col->curr    );
    printf("tfs table with %d selected rows    \n",row->curr    );
  }

    da1 = tfstab->d_cols;
    sa1 = tfstab->s_cols;

    for(j1=0; j1<col->curr; j1++) {
       /* define data arrays from TFS columns*/
       if (SDDS_DefineArray(&SDDS_table, tfstab->columns->names[col->i[j1]], 
                            NULL, NULL, NULL, NULL, 
                            sdds_tfs_types[tfstab->columns->inform[col->i[j1]]], 0, 1, NULL)<0)   { 
         SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
         exit(1);
       }
      if (get_option("debug")) printf("column: %s  type: %d\n",tfstab->columns->names[col->i[j1]],
                                                               tfstab->columns->inform[col->i[j1]]);
      pos[j1] = name_list_pos(tfstab->columns->names[col->i[j1]],tfstab->columns);
      if (get_option("debug")) printf("position: %d %d\n",j1,pos[j1]);
    }


    /* define parameters from TFS table header */
    if(tfstab->header != NULL) {
      if (get_option("debug")) printf("c ==> %d\n",tfstab->header->curr);
     treat_tfs_header_define(&SDDS_table, tfstab); // i1 = not used
    }

    if (!SDDS_SaveLayout(&SDDS_table)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    if (!SDDS_WriteLayout(&SDDS_table)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    if (!SDDS_StartTable(&SDDS_table,1000)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }

  /* set parameters from TFS table header */
    if(tfstab->header != NULL) {
      if (get_option("debug")) printf("cc ==> %d\n",tfstab->header->curr);
     treat_tfs_header_set(&SDDS_table, tfstab); // i1 = not used
    }

  /* fill SDDS table with data */
    for(j1=0; j1<col->curr; j1++) {
      if (get_option("debug")) printf("column: %s  type: %d\n",tfstab->columns->names[col->i[j1]],
                                                               tfstab->columns->inform[col->i[j1]]);
      pos[j1] = name_list_pos(tfstab->columns->names[col->i[j1]],tfstab->columns);
      if (get_option("debug")) printf("position: %d %d\n",j1,pos[j1]);

      if(tfstab->columns->inform[col->i[j1]] == 1) {
          /* need a long buffer for double to long conversion */
          pl = mycalloc_atomic("long buffer for SDDS", tfstab->curr, sizeof *pl);
          for(j2=0;j2<tfstab->curr; j2++) {
            /* convert from double to long */
            pl[j2] = da1[pos[j1]][j2];
            if (get_option("debug")) printf(" %ld\n",pl[j2]);
          }
          if (!SDDS_SetArrayVararg(&SDDS_table, tfstab->columns->names[col->i[j1]], 
               SDDS_POINTER_ARRAY, pl,tfstab->curr)) {
               SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
               exit(1);
          }
          if(pl != NULL) myfree("write sdds",pl);
      }

      if(tfstab->columns->inform[col->i[j1]] == 2) {
          /* pd = (double *)da1[pos[j1]]; */
          for(j2=0;j2<tfstab->curr; j2++) {
            if (get_option("debug")) {
                 printf(" %e\n",da1[pos[j1]][j2]);
                 printf("FILLING ? %s %d %e\n",sa1[pos[0]][j2],row->i[j2],da1[pos[j1]][j2]);
            }
            /* for row->i[j2] == 1  ==> row is selected, not implemented */
          }
          if (!SDDS_SetArrayVararg(&SDDS_table, tfstab->columns->names[col->i[j1]], 
               SDDS_POINTER_ARRAY, &da1[pos[j1]][0],tfstab->curr)) {
               SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
               exit(1);
          }
      }

      if(tfstab->columns->inform[col->i[j1]] == 3) {
          for(j2=0;j2<tfstab->curr; j2++) {
            if (get_option("debug")) printf(" %s\n",tfstab->s_cols[pos[j1]][j2]);
          }
          if (!SDDS_SetArrayVararg(&SDDS_table, tfstab->columns->names[col->i[j1]], 
               SDDS_POINTER_ARRAY, &sa1[pos[j1]][0],tfstab->curr)) {
               SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
               exit(1);
          }
      }

    }
  /*                           */
    if (!SDDS_WriteTable(&SDDS_table)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    if (!SDDS_Terminate(&SDDS_table)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    
  return(0);
}

static int
pass_select_tab(char* name, struct command* sc)
  /* checks name against class (if element) and pattern that may
     (but need not) be contained in command sc;
     0: does not pass, 1: passes */
{
  struct name_list* nl = sc->par_names;
  struct command_parameter_list* pl = sc->par;
  struct element* el = find_element(strip(name), element_list);
  int pos, in = 0, any = 0;
  char *class, *pattern;
  pos = name_list_pos("class", nl);
  if (pos > -1 && nl->inform[pos])  /* parameter has been read */
  {
    el = find_element(strip(name), element_list);
    if (el != NULL)
    {
      class = pl->parameters[pos]->string;
      in = belongs_to_class(el, class);
      if (in == 0) return 0;
    }
  }
  any = in = 0;
  pos = name_list_pos("pattern", nl);
  if (pos > -1 && nl->inform[pos])  /* parameter has been read */
  {
    any = 1;
    pattern = stolower(pl->parameters[pos]->string);
    if(myregex(pattern, strip(name)) == 0)  in = 1;
  }
  if (any == 0) return 1;
  else return in;
}

static void
set_selected_rows_tab(struct table* t, struct command_list* select, struct command_list* deselect)
{
  int i, j, n = 0;
  if (select != 0)
  {
    for (j = 0; j < t->curr; j++)  t->row_out->i[j] = 0;
       t->row_out->curr = 0;
    for (i = 0; i < select->curr; i++)
    {
      for (j = 0; j < t->curr; j++)
      {
        if (t->row_out->i[j] == 0) t->row_out->i[j]
                                     = pass_select_tab(t->s_cols[0][j], select->commands[i]);
        if (t->row_out->i[j] == 1) n++;
      }
    }
  }
  if (deselect != NULL)
  {
    for (i = 0; i < deselect->curr; i++)
    {
      for (j = 0; j < t->curr; j++)
      {
        if (t->row_out->i[j] == 1) t->row_out->i[j]
                                     = 1 - pass_select_tab(t->s_cols[0][j], deselect->commands[i]);
        if (t->row_out->i[j] == 1) n++;
      }
    }
  }
  t->row_out->curr = n;
}

static void
sel_table(char* tname, struct table* t)
  /* output of a table */
{
  int j;
  struct command_list* scl = find_command_list(tname, table_select);
  struct command_list* dscl = find_command_list(tname, table_deselect);
  while (t->num_cols > t->col_out->max)
    grow_int_array(t->col_out);
  while (t->curr > t->row_out->max)
    grow_int_array(t->row_out);
  t->row_out->curr = t->curr;
  if (par_present("full", NULL, scl))
    put_info("obsolete option 'full'"," ignored on 'select'");
  for (j = 0; j < t->curr; j++) t->row_out->i[j] = 1;
  for (j = 0; j < t->num_cols; j++) t->col_out->i[j] = j;
  t->col_out->curr = t->num_cols;
  if ((scl != NULL && scl->curr > 0) || (dscl != NULL && dscl->curr > 0))
  {
    set_selected_columns(t, scl);
    set_selected_rows_tab(t, scl, dscl);
  }
}
                                                                                                          
static int
sdds_ior(struct in_cmd* cmd)
{
    char *sdds_table_file;
    char *tfs_table_name;
    int   i;

    if((sdds_table_file = command_par_string("file",cmd->clone)) == NULL) {
         fatal_error("No file name to read SDDS table ","\n");
    }
    if((tfs_table_name = command_par_string("table",cmd->clone)) == NULL) {
         fatal_error("No table name to read SDDS table ","\n");
    }
    printf("access SDDS table: %s %s\n",sdds_table_file,tfs_table_name);
    i = sdds_readt(sdds_table_file,tfs_table_name);
    return(i);
}

static int 
sdds_iow(struct in_cmd* cmd)
{
    char *sdds_table_file;
    char *tfs_table_name;
    struct table *tfs_table=0;
    int   i, pos;
    if((sdds_table_file = command_par_string("file",cmd->clone)) == NULL) {
         fatal_error("No file name to write SDDS table ","\n");
    }
    if((tfs_table_name = command_par_string("table",cmd->clone)) == NULL) {
         fatal_error("No table name to write SDDS table ","\n");
    }

    mycpy(c_dum->c, tfs_table_name);
    if ((pos = name_list_pos(c_dum->c, table_register->names)) > -1) {
      tfs_table = table_register->tables[pos];
    }

    printf("create SDDS table: %s %d\n",sdds_table_file,pos);
    sel_table(tfs_table_name, tfs_table);
    i = sdds_writet_sel(sdds_table_file,tfs_table);
    return(i);
}

// public interface

void
pro_sdds(struct in_cmd* cmd)
{
  if (strcmp(cmd->tok_list->p[0], "sddsin") == 0)
    sdds_ior(cmd);
  else if (strcmp(cmd->tok_list->p[0], "sddsout") == 0)
    sdds_iow(cmd);
}

#endif // _ONLINE

