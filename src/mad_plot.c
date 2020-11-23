#include "madx.h"

static void
gnuplot_append(const char *gplfilename, char *psfilename)
/* Append the gnuplot ps file to the main ps file */
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
  else {
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
    else {
      /* write the lines */
      if(np==0) fprintf(newpsfile,"%s",line);
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
    else {
      /* write the lines */
      if(np==0) fprintf(newpsfile,"%s",line);
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
    *l = strlen(title)+1;
    strfcpy(tlt, title, *l);
  }
}

void
get_version(char* tlt, int* l)
  /* returns version number */
{
  time_t tmp;
  struct tm* tm;
  time(&tmp);
  tm = localtime(&tmp);
  strcpy(tlt, "MAD-X ");
  strcat(tlt, version_name);
  sprintf(tlt+strlen(tlt),
          "  %02d/%02d/%02d %02d.%02d.%02d ",
          tm->tm_mday, tm->tm_mon+1, tm->tm_year%100,
          tm->tm_hour, tm->tm_min, tm->tm_sec);
  *l = strlen(tlt);
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
  int i, j, k, ierr, notitle = !strcmp(title,"no-title");
  int nointerp = 0, multiple = 0, noversion = 0, nolegend = 0, s_haxis = 0, track_flag = 0;
  int tsm1 = TITLE_SIZE - 1, tsm2 = TITLE_SIZE - 2;
  int part_idx[100], curr, track_cols_length, haxis_idx = 0, vaxis_idx = 0;
  int size_plot_title = tsm1, size_version = tsm1;
  int *title_length = &size_plot_title, *version_length = &size_version;
  char* pt = title, *haxis_name = NULL, *vaxis_name = NULL;
  char* particle_list;
  struct table* p_table = NULL;
  const char *table_name = NULL, *file_name = NULL;
  char *last_twiss_table, *trackfile;
  char track_file_name[NAME_L], ps_file_name[NAME_L];
  char plot_title[TITLE_SIZE], version[TITLE_SIZE];
  FILE *gpu;
  struct command_parameter* cp;

  int debug = get_option("debug");

  /* use correct beam for sequence to be plotted - HG 031127 */
  struct command* keep_beam = current_beam;
  if (attach_beam(current_sequ) == 0)
    fatal_error("PLOT - sequence without beam:", current_sequ->name);
  /* end part1 of HG 031127 */


  if (this_cmd == NULL || this_cmd->clone == NULL)
    fatal_error("Plot "," - non existing command");

  /* Check table name is the same as in the last twiss command */

  /* get vaxis_name */
  command_par("vaxis", this_cmd->clone, &cp);
  vaxis_name = cp->m_string->p[0];

  /* get interpolation */
  nointerp = !par_present("interpolation", this_cmd->clone); // TG: should be "interpolate"?!

  /* get haxis_name & s_haxis flag */
  haxis_name = command_par_string_user("haxis", this_cmd->clone);
  if (haxis_name) {
    s_haxis = !strcmp(haxis_name,"s");
  }

  /* get table_name & track_flag */
  table_name = command_par_string_user("table", this_cmd->clone);
  if(table_name) {
    if(strncmp(table_name,"track",5) == 0)
      track_flag = 1;

    if (debug)
      printf("tablename and trackflag: %s %d\n", table_name, track_flag);
  }
  else {
    printf("Plot - default table plotted: twiss\n");
    table_name = "twiss";
  }

  /* get file_name */
  file_name = command_par_string_user("file", this_cmd->clone);
  if(!file_name) {
    if (track_flag) file_name = "madx_track";
    else file_name = "madx";
  }


  /* If table name begins with "track" use the gnuplot package, and access files not tables*/
  if (track_flag) {
    /* get track file name */
    trackfile = command_par_string("trackfile", this_cmd->clone);
    printf("trackfile is: %s\n", trackfile);

    // now trackfile contains the basis name from which we will get the real file names

    /* get particle */
    command_par("particle", this_cmd->clone, &cp);
    curr = cp->m_string->curr;
    for (i = 0; i < curr; i++) {
      particle_list = cp->m_string->p[i];
      part_idx[i] = atoi(particle_list);
    }

    multiple  = par_present("multiple",  this_cmd->clone);
    noversion = par_present("noversion", this_cmd->clone);
    nolegend  = par_present("nolegend",  this_cmd->clone);

    /* find the column numbers corresponding to haxis_name & vaxis_name */
    track_cols_length = track_table_cols_len-1;
    for (j = 0; j < track_cols_length; j++) {
      if(strcmp(track_table_cols[j],haxis_name) == 0 && haxis_idx == 0)
	haxis_idx = j + 1;
      if(strcmp(track_table_cols[j],vaxis_name) == 0 && vaxis_idx == 0)
	vaxis_idx = j + 1;
    }

    /* build-up the title */
    for (j = 0; j < tsm1; j++) {
      plot_title[j] = ' ';
      version[j] = ' ';
    }
    plot_title[tsm1] = '\0';
    version[tsm1] = '\0';
    get_title(plot_title,title_length);
    for (k = *title_length + 1; k > 0; k--) {
      plot_title[k] = plot_title[k - 1];
    }
    plot_title[0]= '\"';
    if (noversion) {
      plot_title[*title_length+1] =  '\"';
      plot_title[*title_length+2] =  '\0';
    }
    else {
      plot_title[tsm2] =  '\"';
      get_version(version,version_length);
      k = tsm2 - *version_length;
      for (j = k; j < tsm2; j +=1) {
	plot_title[j] = version[j - k];
      }
    }

    /* build-up the gnuplot command file */
    mycpy(track_plot_filename,file_name);
    sprintf(ps_file_name,"%s",track_plot_filename);
    strcat(ps_file_name,".ps");

    gpu = fopen("gnu_plot.gp","w");
    fprintf(gpu,"set terminal postscript color\n");
    fprintf(gpu,"set pointsize 0.48\n");
    fprintf(gpu,"set output '%s'\n","tmpplot.ps");
    fprintf(gpu,"set title %s\n",plot_title);
    fprintf(gpu,"set xlabel '%s'\n",haxis_name);
    fprintf(gpu,"set ylabel '%s'\n",vaxis_name);

    if (strcmp(table_name,"trackone") == 0) {  // case of a single file with all particle coordinates
      sprintf(track_file_name, "%sone", trackfile);
      printf("looking for file %s \n", track_file_name);
      if (fopen(track_file_name,"r") == NULL)
	printf("file %s does not exist \n",track_file_name);
      else
	for (j = 0; j < curr; j++) {
	  if (j == 0) fprintf(gpu,"plot "); // first plot command
	  else if (multiple == 0) fprintf(gpu,"\nplot "); // one plot per particle
	  else fprintf(gpu,", \\\n     "); // all particles on one plot
	  fprintf(gpu,"'%s' using %d:($1==%d ? $%d : NaN) ",track_file_name,haxis_idx,part_idx[j],vaxis_idx);
	  if (nolegend)
	    fprintf(gpu,"notitle with points pointtype %d ",part_idx[j]);
	  else
	    fprintf(gpu,"title 'particle %d' with points pointtype %d ",part_idx[j],part_idx[j]);
	}
    }
    else { // case of multiple files each containing single particle coordinates
      for (j = 0; j < curr; j++) {
	sprintf(track_file_name, "%s.obs%04d.p%04d", trackfile, 1, part_idx[j]);
	printf("looking for file %s \n", track_file_name);
	if (fopen(track_file_name,"r") == NULL)
	  printf("file %s does not exist \n",track_file_name);
	else {
	  if (j == 0) fprintf(gpu,"plot ");
	  else if (multiple == 0) fprintf(gpu,"\nplot ");
	  else fprintf(gpu,", \\\n     ");
	  fprintf(gpu,"'%s' using %d:%d ",track_file_name,haxis_idx,vaxis_idx);
	  if (nolegend)
	    fprintf(gpu,"notitle with points pointtype %d ",part_idx[j]);
	  else
	    fprintf(gpu,"title 'particle %d' with points pointtype %d ",part_idx[j],part_idx[j]);
	}
      }
    }

    fclose(gpu);

#if 0
    {  // backup gnuplot command file for debugging purpose
      static int i = 0;
      char cmd[1000];
#ifdef _WIN32
      const char *cp = "copy /Y";
#else
      const char *cp = "cp -f";
#endif
      sprintf(cmd, "%s gnu_plot.gp gnu_plot_%d.gp", cp, ++i);
      if (system(cmd) == -1)
	warning("Plot - system cannot run the command: ", cmd);
    }
#endif

    /* gnuplot command file ready. it produces the file "tmpplot.ps"*/
    if (system("gnuplot gnu_plot.gp") == -1)
      warning("Plot - system cannot run the command: ", "gnuplot gnu_plot.gp");
      // and we leave the .gp file behind...
    else { // Copy or append the gnuplot ps file in the target ps_file */
      gnuplot_append("tmpplot.ps",ps_file_name);
      remove("gnu_plot.gp"); // remove the .gp file
    }

  } // end of plotting of tracking data : if (track_flag)...

  else { /* normal plot */

    /* if haxis is "s" and no interpolation, check if table name is the same of the last twiss call */
    if (nointerp == 0 && s_haxis == 1) {

      if (!current_sequ->tw_table) {
	warning("PLOT - no TWISS table present", "PLOT command ignored");
	return;
      }

      last_twiss_table = current_sequ->tw_table->name;

      if (strcmp(table_name,"aperture") != 0 ) {
	if(strcmp(table_name,last_twiss_table) != 0) {
	  printf("Only allowed table attribute in plot command is \"aperture\"."
		 " Else, table name is automatically changed to %s \n", last_twiss_table);
	  if ((cp->string = last_twiss_table) == NULL)
	    cp->call_def->string =last_twiss_table ;
	}
      }
    }

    //if (debug)
    //  printf("name_list_pos(table_name, table_register->names) = %d \n",
    //        name_list_pos(table_name, table_register->names));

    /* HG 21.10.09 allow plot from external table, part1 */
    p_table = find_table(table_name);
    if (!p_table) {
      /* fatal_error("Plot - non-existing table:", table_name); return; */
      warning("Plot - potentially non-existing table:", table_name);
      return;
    }
    /* HG 21.10.09 allow plot from external table, end part1 */

    embedded_twiss_cmd = cmd;
    /* HG 21.10.09 allow plot from external table, part 2 */
    //    printf("P_TABLE CHECK %p\n", p_table);
    if (p_table->origin) title = p_table->name;
    else if (notitle && current_sequ != NULL) title = current_sequ->name;

    pesopt_(&ierr);

    if (ierr == 0) {
      if (p_table->origin == 0) {
        // LD 2016.04.19
	adjust_beam();
	probe_beam = clone_command(current_beam);
	adjust_probe_fp(twiss_deltas->a[0]); /* sets correct gamma, beta, etc. */
      }
      pefill_(&ierr);
      pemima_();
      plotit_(&plots_made);
      plots_made = 1;
    }

    /* HG 21.10.09 allow plot from external table, end part 2 */
    if (notitle) title = pt;
  }

  /* part 2 of HG 031127 */
  current_beam = keep_beam;
  /* end of part 2 of HG 031127 */
}

