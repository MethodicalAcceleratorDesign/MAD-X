/* extract SXF file from mad-X */
/* #define _call_tree_ */


void pro_sxf(struct in_cmd* cmd)
     /* controls reading and writing of SXF format files */
{
  struct name_list* nl = cmd->clone->par_names;
  struct command_parameter_list* pl = cmd->clone->par;
  int pos = name_list_pos("file", nl);
  char* filename = NULL;
  FILE* inout;

  if (nl->inform[pos])
    {
     if ((filename = pl->parameters[pos]->string) == NULL)
       {
        if (pl->parameters[pos]->call_def != NULL)
        filename = pl->parameters[pos]->call_def->string;
       }
    }
  else
    {
     if (pl->parameters[pos]->call_def != NULL)
     filename = pl->parameters[pos]->call_def->string;
    }
  if (filename == NULL) filename = permbuff("dummy");
  if (strcmp(cmd->tok_list->p[0], "sxfread") == 0)
    {
     if ((inout = fopen(filename, "r")) == NULL)
       {
        warning("cannot open input file: ", filename);
        return;
       }
     sxf_read(cmd->clone, inout);
    }
  else if (strcmp(cmd->tok_list->p[0], "sxfwrite") == 0)
    {
     if ((inout = fopen(filename, "w")) == NULL)
       {
        warning("cannot open output file: ", filename);
        return;
       }
     sxf_write(cmd->clone, inout);
    }
}

void sxf_read(struct command* comm, FILE* in)
     /* reads an expanded sequence including errors from an SXF file */
{
  puts("entered sxfread");
}

void sxf_write(struct command* comm, FILE* out)
     /* writes the currently USEd sequence in SXF format to a file */
{
  if (current_sequ == NULL || current_sequ->ex_start == NULL)
    {
     warning("SXFWRITE, but no active sequence:", "ignored");
     return;
    }
  sxf_rtag();
  put_line(out, "// SXF version 2.0");
  sprintf(c_dummy, "%s sequence", current_sequ->name); 
  put_line(out, c_dummy);
  s_indent(add_indent[0]); put_line(out, "{");
  current_node = current_sequ->range_start;
  sequ_start = current_node->position;
  while (current_node)
    {
     if (strchr(current_node->name, '$') == NULL)  pro_elem_sxf(out);
     if ((current_node = current_node->next) == current_sequ->range_end) break;
    }
  sequ_end = current_node->position;
  sequ_length = sequ_end - sequ_start;
  sprintf(c_dummy, "endsequence at = %.12g", sequ_length);
  put_line(out, c_dummy);
  indent = b_indent[--b_level];
  put_line(out, "}");
  put_line(out, "// SXF end");
  printf("\nSXF_ex terminated - total number of elements: %d\n", sxf_elem_cnt);
  printf("              elements with alignment errors: %d\n", sxf_align_cnt);
  printf("              elements with field     errors: %d\n", sxf_field_cnt);
}

void sxf_rtag()
{
  int j;
  char type[4][12] = {"kicker", "RBEND", "monitor", "vmonitor"};
  char code[4][12] = {"kick", "rb", "mon", "vmon"};

  tag_cnt = 4;
  for (j = 0; j < tag_cnt; j++)
    {
     strcpy(tag_type[j], type[j]);
     strcpy(tag_code[j], code[j]);
     lower(tag_type[j]); lower(tag_code[j]); 
    }
  if (tag_cnt > 0)  tag_flag = 2;
}

void put_line(FILE* out, char* s)
{
  char tline[2*LINE_MAX];
  int i;
  if (s != line) reset_line(out);
  for (i = 0; i < indent; i++) tline[i] = ' ';
  strcpy(&tline[indent], s);
  fprintf(out, "%s\n",tline);
}

void reset_line(FILE* out)
{
  if (all_blank(line) == 0)  put_line(out,line);
  line[0] = '\0';
}

void s_indent(int i)
{
  b_indent[b_level] = indent;
  indent += i;
  b_level++;
}

int all_blank(char* s)
{
  int i;

  for (i = 0; i < strlen(s); i++)  if(s[i] != ' ') return 0;
  return 1;
}

void pro_elem_sxf(FILE* out)
{
  if (strcmp(current_node->base_name, "drift") == 0)  return; /* skip drifts */
  write_elstart(out);
  write_body(out); sxf_elem_cnt++;
  if (current_node->p_fd_err && current_node->p_fd_err->curr > 0)
    {
     write_field(out, current_node->p_fd_err); sxf_field_cnt++;
    }
  if (current_node->p_al_err && current_node->p_al_err->curr > 0)
    {
     write_align(out, current_node->p_al_err); sxf_align_cnt++;
    }
  write_elend(out);
}

void write_elstart(FILE* out)
{
  char name[NAME_L];
  char* pc;
  s_indent(add_indent[1]);
  if (current_node->occ_cnt > 1)  /* add count to name, print warning */
    {
     if (occnt_add++ == 0)
     printf("+=+=+= SXF_ex warning - making names unique\n\n");
     strcpy(name, current_node->name); replace(name, ':', '.');
    }
  else strcpy(name, current_node->p_elem->name);
  put_line(out, name); s_indent(add_indent[2]);
  sprintf(c_dummy, "%s {", current_node->base_name); put_line(out, c_dummy);
  s_indent(add_indent[3]);
  if (tag_flag == 1 &&  current_node->p_elem != current_node->p_elem->parent)
    {
     sprintf(c_dummy, "tag = %s", current_node->p_elem->parent->name); 
     put_line(out, c_dummy);
    }
  else if (tag_flag == 2 && (pc = tag_spec(current_node->base_name)) != NULL)
    {
     sprintf(c_dummy, "tag = %s", pc); 
     put_line(out,c_dummy);
    }
  if (current_node->length > zero)
      {
       if (strstr(current_node->base_name, "bend") == NULL)
          sprintf(c_dummy, "l = %.12g", current_node->length);
       else sprintf(c_dummy, "arc = %.12g", current_node->length);
       put_line(out,c_dummy);
      }
  sprintf(c_dummy, "at = %.12g", current_node->position);
  put_line(out,c_dummy);
}

char* tag_spec(char* intype)
{
  int j;

  for (j = 0; j < tag_cnt; j++)
    {
     if (strcmp(intype, tag_type[j]) == 0)  return tag_code[j];
    }
  return NULL;
}

void write_field(FILE* out, struct double_array* fd)
{
  int j;
  int k = -1;

  for (j = 0; j < fd->curr; j++)  if (fd->a[j] != zero) k = j;
  if (++k > 0)
    {
     put_line(out, "body.dev = {");  s_indent(add_indent[4]);
     fill_dump(out, 1, "kl ", fd->a, k, 2);
     fill_dump(out, 1, "kls", &fd->a[1], k, 2);
     put_line(out, "}"); r_indent();
    }
}

void write_align(FILE* out, struct double_array* al)
{
  int j;
  int k = -1;

  for (j = 0; j < al->curr; j++)  if (al->a[j] != zero) k = j;
  if (++k > 0)
    {
     put_line(out, "align.dev = {");  s_indent(add_indent[4]);
     fill_dump(out, 1, "al", al->a, k, 1);
     put_line(out, "}"); r_indent();
    }
}

void write_elend(FILE* out)
{
  put_line(out, "};"); 
  r_indent(); r_indent(); r_indent();
}

void accu_line(FILE* out, char* s)
{
  if (strlen(line) + strlen(s) + indent > LINE_MAX)  reset_line(out);
  strcpy(&line[strlen(line)], s);
}

void fill_dump(FILE* out, int flag, char* label, double* values, int count, 
               int inc)
{
  int j;
  double temp;

  if (flag == 0) sprintf(c_dummy, " %s = ", label);
  else           sprintf(c_dummy, " %s = [", label);
  accu_line(out, c_dummy);
  for (j = 0; j < count; j += inc)
    {
     sprintf(c_dummy, " %.12g", values[j]); accu_line(out, c_dummy);
    }
  if (flag != 0) 
    {
     accu_line(out, "]"); reset_line(out);
    }
}

void write_body(FILE* out)
{
  int i, flag, nval, pos, set = 0;
  double val[100];
  char out_name[NAME_L];
  struct command* eldef;
  char npart[] = "npart";
  eldef = current_node->p_elem->def;
  for (i = 0; i < eldef->par_names->curr; i++)
    {
     if (strcmp(eldef->par_names->names[i], "l") != 0
         && (pos = name_list_pos(eldef->par_names->names[i], sxf_list)) > -1)
       {
  	if ((nval = kl_trans(sxf_list->names[pos], out_name, val, &flag)) > 0)
  	  {
           if (set++ == 0) 
  	     {
              put_line(out, "body = {"); s_indent(add_indent[4]);
  	     }
           fill_dump(out, flag, out_name, val, nval, 1);
  	  }
       }
    }
  if (strcmp(current_node->base_name, "beambeam") == 0)
    {
     if (set++ == 0) 
        {
         put_line(out, "body = {"); s_indent(add_indent[4]);
  	}
     val[0] = command_par_value(npart,current_beam);
     fill_dump(out, flag, npart, val, 1, 1);
    }
  if(set > 0) 
    {
     put_line(out, "}"); r_indent(); 
    }
}

int kl_trans(char* name, char* out_name, double* val, int* flag)
{
  int j, length, rep;
  double corr;
  *flag = 0;
  if (strstr(name, "kick") != NULL)
    {
     if (strchr(name, 'v') == NULL)  
       {
        corr = current_node->chkick;
        strcpy(out_name, "kl");
       }
     else
       {
        corr = current_node->cvkick;
        strcpy(out_name, "kls");
       }
     val[0] = node_value(name) + corr;
     if (val[0] != zero)  return 1;
     else return 0;
    }
  else if (*name == 'k' && strcmp(name, "ks") != 0)
    {
     *flag = 1;
     if (strcmp(name, "knl") == 0 || strcmp(name, "ksl") == 0)
       {
        if (strchr(name, 's') == NULL)  strcpy(out_name, "kl");
        else                            strcpy(out_name, "kls");
        length = element_vector(current_node->p_elem, name, val);
	if (length > 1 || val[0] != zero)   return length;
        else return 0;
       }
     sscanf(&name[1], "%1d", &rep);
     for (j = 0; j < rep; j++)  val[j] = zero;
     if (strchr(name, 's') == NULL)  strcpy(out_name, "kl");
     else                            strcpy(out_name, "kls");
     if ((val[rep] = current_node->length * node_value(name)) != zero)  
          return (rep+1);
     else return 0;
    }
  else
    {
     strcpy(out_name, name);
     val[0] = node_value(name);
     if (strcmp(name, "e1") == 0 || strcmp(name, "e2") == 0)
       val[0] = command_par_value(name, current_node->p_elem->def);
     else val[0] = node_value(name);
     if (val[0] != zero)  return 1;
     else return 0;
    }
}

void r_indent()
{
  if (b_level == 0)
    {
     printf("+=+=+= SXF_ex fatal - too many closing '}'\n");
     exit(1);
    }
  indent = b_indent[--b_level];
}
