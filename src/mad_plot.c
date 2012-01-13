#include "madx.h"

#if 0 // not used...
static int
embedded_plot(void)
  /* returns the embedded_flag */
{
  int ret;
  ret = embedded_flag;
  return ret;
}
#endif

/* Append the gnuplot ps file to the main ps file */
static void
gnuplot_append(char *gplfilename, char *psfilename)
{
  char line[1000];
  FILE *newpsfile;
  FILE *oldpsfile;
  FILE *gplpsfile;
  int np=0;
  int page=0;
  /* if psfilename does not exist rename it as psfilename and exit*/
  newpsfile=fopen(psfilename,"r");
  if( newpsfile==NULL) {
    rename(gplfilename,psfilename);
    return;
  }
  else
  {
    /* the file has to be closed it is going to change */
    fclose(newpsfile);
  };

  /* else append the gnuplot ps file to psfilename
     Save old value */
  rename(psfilename,"tmpoldplot.ps");
  newpsfile=fopen(psfilename,"w");
  oldpsfile=fopen("tmpoldplot.ps","r");
  gplpsfile=fopen(gplfilename,"r");
  /* read old ps file and copy on the new ps file */
  while(fgets(line,1000,oldpsfile)!=NULL){
    /* don't print after %%Trailer */
    if (strncmp("%%Trailer",line,9)==0)  np=1;
    /* Count the pages and rewrite the line */
    if (strncmp("%%Page:",line,7)==0) {
      page++;
      fprintf(newpsfile,"%%%%Page: %d %d\n",page,page);
    }
    else
    {
      /* write the lines */
      if(np==0) {
        fprintf(newpsfile,"%s",line);
      }
    }
  }
  fclose(oldpsfile);
  remove("tmpoldplot.ps");
  /* read gnuplot ps file and append on the final file */
  while(fgets(line,1000,gplpsfile)!=NULL){
    /* don't print after %%Trailer */
    if (strncmp("%%Trailer",line,9)==0)  np=1;
    /* Count the pages and rewrite the line */
    if (strncmp("%%Page:",line,7)==0) {
      page++;
      fprintf(newpsfile,"%%%%Page: %d %d\n",page,page);
    }
    else
    {
      if(np==0) {
        /* write the lines */
        fprintf(newpsfile,"%s",line);
      }
    }
    /* Print after prologue */
    if (strncmp("%%EndProlog",line,11)==0) np=0;
  }
  fclose(gplpsfile);
  remove("tmpplot.ps");
  /* Print the trailer */
  fprintf(newpsfile,"%%%%Trailer\n");
  fprintf(newpsfile,"%%%%DocumentFonts: Times-Roman\n");
  fprintf(newpsfile,"%%%%Pages: %d\n",page);
  fprintf(newpsfile,"%%%%EOF\n");
  fclose(newpsfile);
}

// public interface

void
get_title(char* tlt, int* l)
  /* copies title from buffer into tl without trailing '\0' */
{
  *l = 0;
  if (title != NULL) {
    *l = strlen(title);
    strncpy(tlt, title, *l);
  }
}

void
get_version(char* tlt, int* l)
  /* returns version number */
{
  time_t tmp;
  struct tm* tm;
  int n = strlen(version_name);
  time(&tmp);
  tm = localtime(&tmp); /* split system time */
  strcpy(tlt, version_name);
  tlt += n;
  sprintf(tlt, "  %02d/%02d/%02d %02d.%02d.%02d\n",
          tm->tm_mday, tm->tm_mon+1, tm->tm_year%100,
          tm->tm_hour, tm->tm_min, tm->tm_sec);
  *l = n + 19;
}

double
plot_option(char* name)
  /* returns the value of setplot parameters */
{
  double val = zero;
  int i;
  mycpy(c_dum->c, name);
  if (plot_options != NULL
      && (i = name_list_pos(c_dum->c, plot_options->par_names)) > -1)
    val = plot_options->par->parameters[i]->double_value;
  return val;
}

void
exec_plot(struct in_cmd* cmd)
{
  int i, j, k, ierr, pos, nt = strcmp(title,"no-title") == 0 ? 1 : 0;
  int nointerp = 0, multiple = 0, noversion = 0, nolegend = 0, s_haxis = 1, track_flag = 0;
  int tsm1 = TITLE_SIZE - 1, tsm2 = TITLE_SIZE - 2;
  int part_idx[100], curr, track_cols_length, haxis_idx = 0, vaxis_idx = 0;
  int size_plot_title = tsm1, size_version = tsm1;
  int *title_length = &size_plot_title, *version_length = &size_version;
  char* pt = title, *haxis_name = NULL, *vaxis_name = NULL, *file_name = NULL;
  char* particle_list;
  struct name_list* nl_plot = NULL;
  struct command_parameter_list* pl_plot = NULL;
  struct table* p_table = NULL;
  char *table_name, *last_twiss_table, *trackfile;
  char track_file_name[NAME_L], ps_file_name[NAME_L];
  char plot_title[TITLE_SIZE], version[TITLE_SIZE];
  FILE *gpu;

  /* use correct beam for sequence to be plotted - HG 031127 */
  struct command* keep_beam = current_beam;
  if (attach_beam(current_sequ) == 0)
    fatal_error("TWISS - sequence without beam:", current_sequ->name);
  /* end part1 of HG 031127 */

  /* Check table name is the same as in the last twiss command */
  if (this_cmd != NULL && this_cmd->clone != NULL)
  {
    nl_plot = this_cmd->clone->par_names;
    pl_plot = this_cmd->clone->par;

    /* get vaxis_name */

    pos = name_list_pos("vaxis", nl_plot);

    vaxis_name = pl_plot->parameters[pos]->m_string->p[0];

    /* get interpolation */

    pos = name_list_pos("interpolation", nl_plot);
    nointerp = 1 - nl_plot->inform[pos];

    /* get haxis_name & s_haxis flag */

    pos = name_list_pos("haxis", nl_plot);
    if(nl_plot->inform[pos])
    {
      if ((haxis_name = pl_plot->parameters[pos]->string) == NULL)
        haxis_name = pl_plot->parameters[pos]->call_def->string;
      s_haxis = strcmp(haxis_name,"s");
    }

    /* get table_name & track_flag */

    pos = name_list_pos("table", nl_plot);
    if(nl_plot->inform[pos]) /* table name specified */
    {

      if ((table_name = pl_plot->parameters[pos]->string) == NULL)
        table_name = pl_plot->parameters[pos]->call_def->string;

      if(strcmp(table_name,"track") == 0)
        track_flag = 1;
    }
    else
    {
      printf("Plot - default table plotted: twiss\n");
      table_name = "twiss";
    }
    /* check if table name is the same of the last twiss call if haxis is "s" and no interpolation */

    if(nointerp == 0 && s_haxis == 0)
    {
      last_twiss_table = current_sequ->tw_table->name;
      if (strcmp(table_name,"aperture") != 0 )
      {
        if(strcmp(table_name,last_twiss_table) != 0)
        {
          printf("Only allowed table attribute in plot command is \"aperture\". Else, table name is automatically changed to %s \n",last_twiss_table );
          if ((pl_plot->parameters[pos]->string = last_twiss_table) == NULL)
            pl_plot->parameters[pos]->call_def->string =last_twiss_table ;
        }
      }
    }

    /* HG 21.10.09 allow plot from external table, part1 */
    if ((pos = name_list_pos(table_name, table_register->names)) > -1)
        p_table = table_register->tables[pos];
    else
      {
/*       fatal_error("Plot - non-existing table:", table_name); return; */
       warning("Plot - potentially non-existing table:", table_name);
       /*p_table = table_register->tables[pos];*/

      }
    /* HG 21.10.09 allow plot from external table, end part1 */

    /* get file_name */

    pos = name_list_pos("file", nl_plot);
    if(nl_plot->inform[pos]) /* file name specified */
    {
      if ((file_name = pl_plot->parameters[pos]->string) == NULL)
        file_name = pl_plot->parameters[pos]->call_def->string;
    }
    else
    {
      if (track_flag)
        file_name = "madx_track";
      else
        file_name = "madx";
    }
  }
  else
    fatal_error("Plot "," - non existing command");

  /* If table name is "track" use the gnuplot package */

  if (track_flag)
  {
    /* get track file name */

    trackfile = command_par_string("trackfile", this_cmd->clone);

    /* get particle */

    pos = name_list_pos("particle", nl_plot);
    curr = pl_plot->parameters[pos]->m_string->curr;
    for (i = 0; i < curr; i++)
    {
      particle_list = pl_plot->parameters[pos]->m_string->p[i];
      part_idx[i] = atoi(particle_list);
    }

    /* get multiple */

    pos = name_list_pos("multiple", nl_plot);
    multiple = nl_plot->inform[pos];

    /* get noversion */

    pos = name_list_pos("noversion", nl_plot);
    noversion = nl_plot->inform[pos];

    /* get nolegend */

    pos = name_list_pos("nolegend", nl_plot);
    nolegend = nl_plot->inform[pos];

    /* find the column numbers corresponding to haxis_name & vaxis_name */

    track_cols_length = track_table_cols_len-1; // sizeof(track_table_cols)/sizeof(uintptr_t) - 1;
    for (j = 0; j < track_cols_length; j++)
    {
      if(strcmp(track_table_cols[j],haxis_name) == 0 && haxis_idx == 0)
        haxis_idx = j + 1;
      if(strcmp(track_table_cols[j],vaxis_name) == 0 && vaxis_idx == 0)
        vaxis_idx = j + 1;
    }

    /* build-up the title */

    for (j = 0; j < tsm1; j++)
    {
      plot_title[j] = ' ';
      version[j] = ' ';
    }
    plot_title[tsm1] = '\0';
    version[tsm1] = '\0';
    get_title(plot_title,title_length);
    for (k = *title_length + 1; k > 0; k--)
    {
      plot_title[k] = plot_title[k - 1];
    }
    plot_title[0]= '\"';
    if (noversion)
    {
      plot_title[*title_length+1] =  '\"';
      plot_title[*title_length+2] =  '\0';
    }
    else
    {
      plot_title[tsm2] =  '\"';
      get_version(version,version_length);
      k = tsm2 - *version_length;
      for (j = k; j < tsm2; j +=1)
      {
        plot_title[j] = version[j - k];
      }
    }

    /* build-up the gnuplot command file */
    mycpy(track_plot_filename,file_name);
    sprintf(ps_file_name,"%s",track_plot_filename);
    strcat(ps_file_name,".ps");

    gpu = fopen("gnu_plot.cmd","w");
    fprintf(gpu,"set terminal postscript color\n");
    fprintf(gpu,"set pointsize 0.48\n");
    fprintf(gpu,"set output '%s'\n","tmpplot.ps");

    fprintf(gpu,"set title %s\n",plot_title);
    fprintf(gpu,"set xlabel '%s'\n",haxis_name);
    fprintf(gpu,"set ylabel '%s'\n",vaxis_name);
    for (j = 0; j < curr; j++)
    {
      sprintf(track_file_name, "%s.obs%04d.p%04d", trackfile, 1, part_idx[j]);
      if (fopen(track_file_name,"r") == NULL)
        printf("file %s does not exist \n",track_file_name);
      else
      {
        if (j == 0) fprintf(gpu,"plot ");
        else
        {
          if (multiple == 0)
            fprintf(gpu,"\nplot ");
          else
            fprintf(gpu,", \\\n     ");
        }
        fprintf(gpu,"'%s' using %d:%d ",track_file_name,haxis_idx,vaxis_idx);
        if (nolegend)
          fprintf(gpu,"notitle with points %d ",part_idx[j]);
        else
          fprintf(gpu,"title 'particle %d' with points %d ",part_idx[j],part_idx[j]);

      }
    }
    fclose(gpu);
    /* gnuplot command file ready. it produces the file "tmpplot.ps"*/
    system("gnuplot 'gnu_plot.cmd'");
    /* Copy or append the gnuplot ps file in the target ps_file */
    gnuplot_append("tmpplot.ps",ps_file_name);
    /* Remove the gnuplot command */
/*    remove("gnu_plot.cmd");*/
  }
  else

    /* normal plot */

  {
    embedded_twiss_cmd = cmd;
    /* HG 21.10.09 allow plot from external table, part 2 */
    if (p_table->origin) title = p_table->name;
    else if (nt && current_sequ != NULL) title = current_sequ->name;
    pesopt_(&ierr);
    if (ierr == 0)
    {
      if (p_table->origin == 0)
      {
       adjust_beam();
       probe_beam = clone_command(current_beam);
       adjust_probe(twiss_deltas->a[0]); /* sets correct gamma, beta, etc. */
       adjust_rfc(); /* sets freq in rf-cavities from probe */
      }
      pefill_(&ierr);
      pemima_();
      plotit_(&plots_made);
      plots_made = 1;
    }
    /* HG 21.10.09 allow plot from external table, end part 2 */
    if (nt) title = pt;
  }

  /* part 2 of HG 031127 */
  current_beam = keep_beam;
  /* end of part 2 of HG 031127 */
}

int
interp_node(int *nint)
{

  /* Creates interpolating nodes for the plotting routine */

  struct node *first_node, *clone;
  struct element* el;
  struct command_parameter* cp;
  int i;
  int j, number_nodes;
  double bvk, angle, e1, e2, h1, h2, fint, hgap;
  double zero = 0.0, minus_one = -1.0, length, step, numint;
  char *elem_name;
  int bend_flag = 0;

  numint = (double)*nint;
  number_nodes = *nint - 1;


  /* Set up length, angle and e2 of the first slice
     (first node in the original sequence) */

  first_node = current_node;
  el = first_node->p_elem;
  elem_name = el->base_type->name;
  rbend = (strcmp(elem_name, "rbend") == 0);
  bend_flag = (strcmp(elem_name, "sbend")*(rbend-1) == 0);

/*  bv = node_value("dipole_bv");*/
  bvk = node_value("other_bv");

  if (bend_flag)
  {
    angle = command_par_value("angle", el->def);
    e1 = command_par_value("e1", el->def);
    e2 = command_par_value("e2", el->def);
    h1 = command_par_value("h1", el->def);
    h2 = command_par_value("h2", el->def);
    fint = command_par_value("fint", el->def);
    fintx_plot = command_par_value("fintx", el->def);
    hgap = command_par_value("hgap", el->def);

    if (rbend)
    {
      e1 = e1 + bvk * angle / two;
      e2 = e2 + bvk * angle / two;
      strcpy(elem_name,"sbend");
      el->def->mad8_type = 3;
    }
    angle = angle/numint;
    /*    store_node_value("angle",&angle); */
    i = name_list_pos("angle", el->def->par_names);
    cp = el->def->par->parameters[i];
    if(cp->expr != NULL) backup_expr = cp->expr;
    backup_type = cp->type;
    cp->type = 2;
    cp->expr = NULL;
    cp->double_value = angle;
    store_node_value("e1",&e1);
    store_node_value("e2",&zero);
    store_node_value("h1",&h1);
    store_node_value("h2",&zero);
    store_node_value("fint",&fint);
    store_node_value("fintx",&zero);
    store_node_value("hgap",&hgap);
  }
  length = first_node->length;
  step = length/numint;
  first_node->length = step;

  /* Set first_node in range_start of the sequence */

  current_sequ->range_start = first_node;

  /* clone the current node */

  clone = clone_node(first_node,0);
  if (bend_flag)
  {
    clone->p_elem = clone_element(first_node->p_elem);
    clone->p_elem->def = clone_command(first_node->p_elem->def);
  }

  /* Reset to first node */

  current_node = first_node;

  /* advance to next node */

  current_node = current_node->next;

  /* set last node in the range to the current node */

  current_sequ->range_end = current_node;


  /* insert nint - 1 nodes in between the two main nodes */

  for (j = 1; j <= number_nodes; j++)
  {
    link_in_front(clone,current_node);
    current_node = current_node->previous;
    current_node->previous->next = current_node;
    store_node_value("angle",&angle);
/*    store_node_value("dipole_bv",&bv); */
    store_node_value("other_bv",&bvk);
    if (bend_flag)
    {
      if (j == 1)
      {
        store_node_value("e2",&e2);
        store_node_value("h2",&h2);
        store_node_value("hgap",&hgap);
        if (fintx_plot < zero)
          store_node_value("fintx",&fint);
        else
          store_node_value("fintx",&fintx_plot);
        store_node_value("fint",&zero);
      }
      else
      {
        store_node_value("e2",&zero);
        store_node_value("h2",&zero);
        store_node_value("fint",&zero);
        store_node_value("fintx",&minus_one);
        store_node_value("hgap",&zero);
      }
      store_node_value("e1",&zero);
      store_node_value("h1",&zero);
    }
    clone = clone_node(first_node,0);
    if (bend_flag)
    {
      clone->p_elem = clone_element(first_node->p_elem);
      clone->p_elem->def = clone_command(first_node->p_elem->def);
    }
  }

  current_node = current_node->previous;

  return 0;
}

int
reset_interpolation(int *nint)
{
  struct node *c_node, *second_node;
  struct command_parameter* cp;
  int i;
  int j,bend_flag = 0;
  double angle=0,length,e1,e2,numint, h1, h2, fint, fintx, hgap, bvk;

  /* Deletes the interpolating nodes expanded by the routine interp_node */

  numint = (double)*nint;

  /* reset first and last node in the sequence range */

  current_sequ->range_start = current_sequ->ex_start;
  current_sequ->range_end = current_sequ->ex_end;

  /* reset current_node at first node */

  for (j = 1; j <= *nint ; j++)
    current_node = current_node->previous;

  /* reset length of first node */

  length = numint*current_node->length;
  current_node->length = length;

  /* resets angle and saves e1 if the element is a bending magnet */

  bend_flag = strcmp(current_node->p_elem->base_type->name, "sbend") == 0 || rbend;
  if (bend_flag)
  {
    angle = numint*node_value("angle");
    i = name_list_pos("angle", current_node->p_elem->def->par_names);
    cp = current_node->p_elem->def->par->parameters[i];
    cp->expr = backup_expr;
    cp->type = backup_type;
    cp->double_value = angle;
       /*    store_node_value("angle",&angle); */
    e1 = node_value("e1");
    h1 = node_value("h1");
    fint = node_value("fint");
    fintx = fintx_plot;
    hgap = node_value("hgap");
  }

  /* advance to nint-th  node (second node in original sequence) */

  for (j = 1; j <= *nint; j++)
    advance_node();
  second_node = current_node;

  /* back to the last interpolated node */

  retreat_node();

  /* saves e2 if the element is a bending magnet */

  if (bend_flag)
  {
    e2 = node_value("e2");
    h2 = node_value("h2");
  }

  /* delete the interpolating nodes */

  for (j = 2; j <= *nint; j++)
  {
    c_node = current_node;

    retreat_node();
    if (bend_flag)
    {
      c_node->p_elem->def = delete_command(c_node->p_elem->def);
      c_node->p_elem = delete_element(c_node->p_elem);
    }
    delete_node(c_node);
  }

  /* current_node points now to the first node of the original sequence */
  /* sets next pointer of first node to second node of original sequence */

  current_node->next = second_node;

  /* sets pointer of second node to first node of original sequence */

  current_node->next->previous = current_node;

  /* Updates the values of e1 and e2 and stores them in first node */

/*  bv = node_value("dipole_bv");*/
  bvk = node_value("other_bv");

  if (bend_flag)
  {
    if (rbend)
    {
      strcpy(current_node->p_elem->base_type->name,"rbend");
      current_node->p_elem->def->mad8_type = 2;
      e1 = e1 - bvk * angle / two;
      e2 = e2 - bvk * angle / two;
    }
    store_node_value("e1",&e1);
    store_node_value("e2",&e2);
    store_node_value("h1",&h1);
    store_node_value("h2",&h2);
    store_node_value("fint",&fint);
    store_node_value("fintx",&fintx_plot);
    store_node_value("hgap",&hgap);
  }

  return 0;
}

