#include "madx.h"

/* extract SXF file from mad-X, or read SXF file into mad-X */

static double
  sequ_length,         /* length of  sequence */
  sequ_start,
  sequ_end;

// forward declarations

static void sxf_fill_command(struct command*, int ntok, char** toks);

// private functions

static void put_line(FILE* out, const char* s);

static void
r_indent(void)
{
  if (b_level == 0)
  {
    printf("+=+=+= SXF_ex fatal - too many closing '}'\n");
    exit(1);
  }
  indent = b_indent[--b_level];
}

#if 0 // not used
static char*
tag_spec(char* intype)
{
  int j;

  for (j = 0; j < tag_cnt; j++)
  {
    if (strcmp(intype, tag_type[j]) == 0)  return tag_code[j];
  }
  return NULL;
}
#endif

static int
join_prefix(const char* prefix, int ntok, char** toks)
  /* joins two tokens into one if the first is the "prefix" */
{
  int j, k;
  for (j = 0; j < ntok; j++)
  {
    if (strcmp(toks[j], prefix) == 0 && j+1 < ntok)
    {
      strcat(toks[j], toks[j+1]);
      for (k = j+1; k < ntok-1; k++)  toks[k] = toks[k+1];
      ntok--;
    }
  }
  return ntok;
}

#if 0 // not used
static int
get_token_list                 /* returns no. of tokens */
(char* text,    /* input: blank separated tokens */
 char** tokens, /* output: token pointer list */
 int max_tok)   /* length of tokens */
{
  char* p;
  int count = 0;

  p = strtok(text," =;\n");
  while (p != NULL && count < max_tok)
  {
    tokens[count++] = p;
    p = strtok(NULL," =;\n");
  }
  return count;
}
#endif

#if 0 // not used
static int
version_header(char* line)            /* processes and checks header */
{
  int j = 0;
  char* header[] = {"//", "sxf", "version", "1.0"};
  char* token_list[10];
  if(get_token_list(line, token_list, 10) != 4) return 0;
  for (j = 0; j < 3; j++)
  {
    if (strcmp(header[j], stolower(token_list[j])) != 0) return 0;
  }
  return 1;
}
#endif

static void
reset_line(FILE* out)
{
  if (all_blank(line) == 0)  put_line(out,line);
  line[0] = '\0';
}

static void
put_line(FILE* out, const char* s)
{
  char tline[2*MADX_LINE_MAX];
  int i;
  if (s != line) reset_line(out);
  for (i = 0; i < indent; i++) tline[i] = ' ';
  strcpy(&tline[indent], s);
  fprintf(out, "%s\n",tline);
}

static void
s_indent(int i)
{
  b_indent[b_level] = indent;
  indent += i;
  b_level++;
}

static void
accu_line(FILE* out, const char* s)
{
  if (strlen(line) + strlen(s) + indent > MADX_LINE_MAX)  reset_line(out);
  strcpy(&line[strlen(line)], s);
}

static void
fill_dump(FILE* out, int flag, const char* label, double* values, int count, int inc)
{
  int j;

  if (flag == 0) sprintf(c_dum->c, " %s = ", label);
  else           sprintf(c_dum->c, " %s = [", label);
  accu_line(out, c_dum->c);
  for (j = 0; j < count; j += inc)
  {
    sprintf(c_dum->c, " %.12g", values[j]); accu_line(out, c_dum->c);
  }
  if (flag != 0)
  {
    accu_line(out, "]"); reset_line(out);
  }
}

static int
kl_trans(const char* name, char* out_name, double* val, int* flag)
{
  int j, length, rep;
  double corr;
  *flag = 0;

  if (strstr(name, "angle") != NULL){
    *flag = 1;
    strcpy(out_name, "kl");
    val[0] = node_value(name);
    return 1;
  }


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

static void
write_elstart(FILE* out)
{
  char name[NAME_L];
  s_indent(add_indent[1]);
  if (current_node->occ_cnt > 1)  /* add count to name, print warning */
  {
    if (occnt_add++ == 0)
      printf("+=+=+= SXF_ex warning - making names unique\n\n");
    strcpy(name, current_node->name); replace(name, ':', '.');
  }
  else strcpy(name, current_node->p_elem->name);
  put_line(out, name); s_indent(add_indent[2]);
  sprintf(c_dum->c, "%s {", current_node->base_name); put_line(out, c_dum->c);
  s_indent(add_indent[3]);
  /*
    if (tag_flag == 1 &&  current_node->p_elem != current_node->p_elem->parent)
    {
    sprintf(c_dum->c, "tag = %s", current_node->p_elem->parent->name);
    put_line(out, c_dum->c);
    }
    else if (tag_flag == 2 && (pc = tag_spec(current_node->base_name)) != NULL)
    {
    sprintf(c_dum->c, "tag = %s", pc);
    put_line(out,c_dum->c);
    }
  */

  sprintf(c_dum->c, "tag = %s", current_node->p_elem->name);
  put_line(out, c_dum->c);

  if (current_node->length > zero)
  {
    if (strstr(current_node->base_name, "bend") == NULL)
      sprintf(c_dum->c, "l = %.12g", current_node->length);
    else sprintf(c_dum->c, "arc = %.12g", current_node->length);
    put_line(out,c_dum->c);
  }

  /*
    sprintf(c_dum->c, "at = %.12g", current_node->position);
    nm printf("%s l = %.12g at = %.12g, prev at= %.12g\n", current_node->name, current_node->length,
    nm current_node->position, global_tmp_at);
    if(current_node->position < global_tmp_at) {
    printf("error: %s position (%.12g) < prev position (%.12g)\n",  current_node->name,
    current_node->position, global_tmp_at);
    sprintf(c_dum->c, "at = %.12g", global_tmp_at);
    }
    else {
    sprintf(c_dum->c, "at = %.12g", current_node->position);
    global_tmp_at = current_node->position;
    }
  */
  sprintf(c_dum->c, "at = %.12g", current_node->position);
  put_line(out,c_dum->c);
}

static void
write_field(FILE* out, struct double_array* fd)
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

static void
write_align(FILE* out, struct double_array* al)
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

static void
write_elend(FILE* out)
{
  put_line(out, "};");
  r_indent(); r_indent(); r_indent();
}

static void
write_body(FILE* out)
{
  int i, flag = 0, nval, pos, set = 0;
  double val[100];
  char out_name[NAME_L];
  struct command* eldef;
  char npart[] = "npart";
  eldef = current_node->p_elem->def;

  /* printf(" %s, , angle = %e \n", current_node->name,  angle); */
  for (i = 0; i < eldef->par_names->curr; i++)
  {
    /* nm printf("  i = %d, %s, value = %e\n", i, eldef->par_names->names[i], node_value(eldef->par_names->names[i])); */
    double v = node_value(eldef->par_names->names[i]);
    if(v == 0.0 && strcmp(eldef->par_names->names[i], "knl") != 0) continue;
    /* nm printf("  i = %d, %s, angle = %e\n", i, eldef->par_names->names[i], angle); */
    if (strcmp(eldef->par_names->names[i], "l") != 0
        && (pos = name_list_pos(eldef->par_names->names[i], sxf_list)) > -1)
    {
      /*
        nm printf("kl_trans\n");
        nm: ad hoc solution to avoid the SXF conceptual bug
        if((nval = angle_kl_trans(angle, sxf_list->names[pos], out_name, val, &flag, &angle_flag)) > 0)
      */
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

static void
pro_elem_sxf(FILE* out)
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

static void
sxf_rtag(void)
{
  int j;
  char type[4][12] = {"kicker", "RBEND", "monitor", "vmonitor"};
  char code[4][12] = {"kick", "rb", "mon", "vmon"};

  tag_cnt = 4;
  for (j = 0; j < tag_cnt; j++)
  {
    strcpy(tag_type[j], type[j]);
    strcpy(tag_code[j], code[j]);
    stolower(tag_type[j]); stolower(tag_code[j]);
  }
  if (tag_cnt > 0)  tag_flag = 2;
}

static void
sxf_write(struct command* cmd, FILE* out)
  /* writes the currently USEd sequence in SXF format to a file */
{
  (void)cmd;
  if (current_sequ == NULL || current_sequ->ex_start == NULL)
  {
    warning("SXFWRITE, but no active sequence:", "ignored");
    return;
  }
  sxf_rtag();
  put_line(out, "// SXF version 2.0");
  sprintf(c_dum->c, "%s sequence", current_sequ->name);
  put_line(out, c_dum->c);
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
  sprintf(c_dum->c, "endsequence at = %.12g", sequ_length);
  put_line(out, c_dum->c);
  indent = b_indent[--b_level];
  put_line(out, "}");
  put_line(out, "// SXF end");
  printf("\nSXF_ex terminated - total number of elements: %d\n", sxf_elem_cnt);
  printf("              elements with alignment errors: %d\n", sxf_align_cnt);
  printf("              elements with field     errors: %d\n", sxf_field_cnt);
}

static double
find_value(const char* name, int ntok, char** toks)
  /* returns value found in construct "name = value", or INVALID */
{
  double val = INVALID;
  int j;
  for (j = 0; j < ntok; j++)
  {
    if (strcmp(toks[j], name) == 0)
    {
      if (j+2 < ntok && *toks[j+1] == '=')
      {
        sscanf(toks[j+2], "%lf", &val);
        break;
      }
    }
  }
  return val;
}

static int
sxf_align_fill(int start, int end, int ntok, char** toks, double* vec)
{
  int cnt = 0, j, rss, res;

  (void)end;
  if (strcmp(toks[start+1], "al") == 0)
  {
    get_bracket_t_range(toks, '[', ']', start+3, ntok, &rss, &res);
    if (rss == start+3)  /* square bracket found */
    {
      rss++; res--;
      for (j = rss; j <= res; j++) sscanf(toks[j], "%lf", &vec[cnt++]);
    }
  }
  return cnt;
}

static int
sxf_field_fill(int start, int end, int ntok, char** toks, double* vec)
{
  int cnt[2], i, j, k, rss, res;

  for (k = 0; k < 2; k++)
  {
    cnt[k] = 0;
    for (i = start+1; i < end; i++)
    {
      if ((k == 0 && strcmp(toks[i], "kl") == 0) ||
          (k == 1 && strcmp(toks[i], "kls") == 0))
      {
        get_bracket_t_range(toks, '[', ']', i+2, ntok, &rss, &res);
        if (rss == i+2)  /* square bracket found */
        {
          rss++; res--; cnt[k] = k;
          for (j = rss; j <= res; j++)
          {
            sscanf(toks[j], "%lf", &vec[cnt[k]]);
            cnt[k] += 2;
          }
        }
      }
    }
  }
  return (cnt[0] < cnt[1] ? cnt[1]-1 : cnt[0]-1);
}

static void
sxf_body_fill(struct command* comm, int start, int end, int ntok, char** toks, double length)
{
  struct name_list* nl = comm->par_names;
  struct command_parameter_list* pl = comm->par;
  int cnt, j, k, pos, rss, res, skew;
  double vec[20];

  /* scan for magnet strengths first - special format */
  for (skew = 0; skew < 2; skew++)
  {
    for (k = start+1; k < end; k++)
    {
      if      (skew == 0 && strcmp("kl", toks[k]) == 0)  break;
      else if (skew == 1 && strcmp("kls", toks[k]) == 0) break;
    }
    if (k < end)
    {
      get_bracket_t_range(toks, '[', ']', k, ntok, &rss, &res);
      if (rss == k + 2)  /* square bracket found */
      {
        rss++; res--;
      }
      else rss = res = k+2;
      cnt = 0;
      for (j = rss; j <= res; j++) sscanf(toks[j], "%lf", &vec[cnt++]);
      if (strstr(toks[1], "bend"))
      {
        if (length == zero) fatal_error("bend without length:", toks[0]);
        /*if (skew == 0) pos = name_list_pos("k0", nl);
          else           pos = name_list_pos("k0s", nl);*/
        /* nm printf("bend v[0] = %e; \n", vec[0]); */
        pos = name_list_pos("angle", nl);
        pl->parameters[pos]->double_value = vec[0];
        /* pl->parameters[pos]->double_value = vec[0] / length; */
      }
      else if (strstr(toks[1], "quad"))
      {
        if (length == zero) fatal_error("quad without length:", toks[0]);
        if (skew == 0)  pos = name_list_pos("k1", nl);
        else            pos = name_list_pos("k1s", nl);
        pl->parameters[pos]->double_value = vec[1] / length;
      }
      else if (strstr(toks[1], "sext"))
      {
        if (length == zero) fatal_error("sextupole without length:",
                                        toks[0]);
        if (skew == 0)  pos = name_list_pos("k2", nl);
        else            pos = name_list_pos("k2s", nl);
        pl->parameters[pos]->double_value = vec[2] / length;
      }
      else if (strstr(toks[1], "oct"))
      {
        if (length == zero) fatal_error("octupole without length:",
                                        toks[0]);
        if (skew == 0)  pos = name_list_pos("k3", nl);
        else            pos = name_list_pos("k3s", nl);
        pl->parameters[pos]->double_value = vec[3] / length;
      }
      else if (strstr(toks[1], "kick"))
      {
        if (skew == 0)  pos = name_list_pos("hkick", nl);
        else            pos = name_list_pos("vkick", nl);
        pl->parameters[pos]->double_value = vec[0];
      }
      else if (strstr(toks[1], "multi"))
      {
        if (skew == 0)  pos = name_list_pos("knl", nl);
        else            pos = name_list_pos("ksl", nl);
        while(pl->parameters[pos]->double_array->max < cnt)
          grow_double_array(pl->parameters[pos]->double_array);
        for (j = 0; j < cnt; j++)
        {
          pl->parameters[pos]->double_array->a[j] = vec[j];
        }
        pl->parameters[pos]->double_array->curr = cnt;
        pl->parameters[pos]->expr_list
          = delete_expr_list(pl->parameters[pos]->expr_list);
      }
    }
  }
  /* now all the other parameters */
  for (k = start+1; k < end; k++)
  {
    if (isalpha(*toks[k]) && strstr(toks[k], "kl") == NULL)
    {
      if ((pos = name_list_pos(toks[k], nl)) < 0)
        warning("unknown input parameter skipped: ", toks[k]);
      else if (k+2 < end && *toks[k+1] == '=' && *toks[k+2] != '[')
        sscanf(toks[k+2], "%lf", &pl->parameters[pos]->double_value);
    }
  }
}

static void
sxf_fill_command(struct command* comm, int ntok, char** toks)
{
  struct name_list* nl = comm->par_names;
  struct command_parameter_list* pl = comm->par;
  int n, pos, rs, re;
  double length = zero;

  if ((pos = name_list_pos("l", nl)) > -1)
  {
    if (strstr(toks[1], "bend"))
    {
      if ((length = find_value("arc", ntok, toks)) == INVALID)
        length = find_value("l", ntok, toks);
    }
    else length = find_value("l", ntok, toks);
    if (length == INVALID) length = zero;
    pl->parameters[pos]->double_value = length;
  }
  for (n = 0; n < ntok; n++)  if(strcmp("body", toks[n]) == 0) break;
  if (n < ntok)
  {
    get_bracket_t_range(toks, '{', '}', n, ntok, &rs, &re);
    if (rs < n) fatal_error("element body empty:", toks[0]);
    sxf_body_fill(comm, rs, re, ntok, toks, length);
  }
}

static int
sxf_decin(char* p, int count) /* decode one SXF input item, store */
{
  int j, n, ntok, pos, rs, re;
  char** toks = tmp_p_array->p;
  struct command* clone;
  struct element* el;
  double at, vec[FIELD_MAX];

  tmp_p_array->curr = 0;
  pre_split(p, aux_buff, 0);
  ntok = mysplit(aux_buff->c, tmp_p_array);
  ntok = join_prefix("-", ntok, toks);
  if (count == 0)
  {
    if (ntok < 2 || strcmp(toks[1], "sequence") != 0) return -1;
    else current_sequ = new_sequence(toks[0], 0);
    if (occ_list == NULL)
      occ_list = new_name_list("occ_list", 10000);  /* for occurrence count */
    else occ_list->curr = 0;
    current_sequ->cavities = new_el_list(100);
    current_sequ->crabcavities = new_el_list(100);
    pos = name_list_pos("marker", defined_commands->list);
    clone = clone_command(defined_commands->commands[pos]);
    sprintf(c_dum->c, "%s$start", current_sequ->name);
    el = make_element(c_dum->c, "marker", clone, 0);
    make_elem_node(el, 1);
    current_sequ->start = current_node;
    for (j = 3; j < ntok; j++) /* push first element down */
    {
      toks[j-3] = toks[j];
    }
    ntok -= 3;
  }
  else if (strcmp(toks[0], "endsequence") == 0)
  {
    current_sequ->length = find_value("at", ntok, toks);
    pos = name_list_pos("marker", defined_commands->list);
    clone = clone_command(defined_commands->commands[pos]);
    sprintf(c_dum->c, "%s$end", current_sequ->name);
    el = make_element(c_dum->c, "marker", clone, 0);
    make_elem_node(el, 1);
    current_node->at_value = current_sequ->length;
    current_sequ->end = current_node;
    current_sequ->start->previous = current_sequ->end;
    current_sequ->end->next = current_sequ->start;
    return 1;
  }
  if ((pos = name_list_pos(toks[1], defined_commands->list)) < 0)
    fatal_error("element type not found:", toks[1]);
  clone = clone_command(defined_commands->commands[pos]);
  sxf_fill_command(clone, ntok, toks);
  el = make_element(toks[0], toks[1], clone, sequ_is_on+1);
  if (strcmp(el->base_type->name, "rfcavity") == 0 &&
      find_element(el->name, current_sequ->cavities) == NULL)
    add_to_el_list(&el, 0, current_sequ->cavities, 0);
  if (strcmp(el->base_type->name, "crabcavity") == 0 &&
      find_element(el->name, current_sequ->crabcavities) == NULL)
    add_to_el_list(&el, 0, current_sequ->crabcavities, 0);
  add_to_name_list(el->name, 1, occ_list);
  make_elem_node(el, 1);
  sxf_suml += el->length / 2;
  if ((at = find_value("at", ntok, toks)) == INVALID) at = sxf_suml;
  else sxf_suml = at;
  sxf_suml += el->length / 2;
  current_node->at_value = at;
  for (n = 0; n < ntok; n++)  if(strcmp("align.dev", toks[n]) == 0) break;
  if (n < ntok)
  {
    get_bracket_t_range(toks, '{', '}', n, ntok, &rs, &re);
    if (rs < n) fatal_error("alignment errors empty:", toks[0]);
    n = sxf_align_fill(rs, re, ntok, toks, vec);
    current_node->p_al_err = new_double_array(n);
    current_node->p_al_err->curr = n;
    for (j = 0; j < n; j++) current_node->p_al_err->a[j] = vec[j];
  }
  for (n = 0; n < ntok; n++)  if(strcmp("body.dev", toks[n]) == 0) break;
  if (n < ntok)
  {
    get_bracket_t_range(toks, '{', '}', n, ntok, &rs, &re);
    if (rs < n) fatal_error("field errors empty:", toks[0]);
    for (j = 0; j < FIELD_MAX; j++) vec[j] = zero;
    n = sxf_field_fill(rs, re, ntok, toks, vec);
    current_node->p_fd_err = new_double_array(n);
    current_node->p_fd_err->curr = n;
    for (j = 0; j < n; j++) current_node->p_fd_err->a[j] = vec[j];
  }
  return 0;
}

static void
sxf_read(struct command* cmd)
  /* reads an expanded sequence including errors from an SXF file */
{
  struct sequence* keep_sequ = current_sequ;
  int echo, err, izero = 0, count = 0; // n, not used
  FILE* in_file = in->input_files[in->curr];
  char *p, *pp;

  (void)cmd;
  sxf_suml = zero;
  if (fgets(aux_buff->c, aux_buff->max, in_file) == NULL)
  {
    warning("SXF input file empty,"," ignored");
    return;
  }
  /*
    if ((rcode = version_header(aux_buff->c)) == 0)
    {
    warning("SXF header missing or wrong,"," ignored");
    return;
    }
  */
  sequ_is_on = 1;
  echo = get_option("echo");
  set_option("echo", &izero);
  get_stmt(in_file, 1); // n = not used
  replace(in->buffers[in->curr]->c_a->c, ',', ' ');
  replace(in->buffers[in->curr]->c_a->c, '\n', ' ');
  p = strtok(in->buffers[in->curr]->c_a->c, ";");
  while(p)
  {
    pp = &p[strlen(p)+1];
    if ((err = sxf_decin(p, count++)) == -1)
    {
      warning("No sequence name found, ", "ignored");
      goto term;
    }
    else if (err == 1)  break;
    p = strtok(pp, ";");
  }
  if (current_sequ->length == zero)
  {
    warning("No endsequence with length found, ", "ignored");
    current_sequ = keep_sequ;
    goto term;
  }
  printf("SXF -- sequence %s: declared length = %e; element l_sum = %e\n",
         current_sequ->name, current_sequ->length, sxf_suml);
  add_to_sequ_list(current_sequ, sequences);
  if (attach_beam(current_sequ) == 0)
    fatal_error("USE - sequence without beam:", current_sequ->name);
  current_sequ->beam = current_beam;
  current_range = tmpbuff("#s/#e");
  expand_curr_sequ(1);
  term:
  set_option("echo", &echo);
  sequ_is_on = 0;
  return;
}

// public functions

void
get_sxf_names(void)
  /* reads and stores names for SXF I/O from madxl.h */
{
  for (int i=0; sxf_table_names[i][0] != ' '; i++)
    add_to_name_list(sxf_table_names[i], 0, sxf_list);
}

void
pro_sxf(struct in_cmd* cmd)
  /* controls reading and writing of SXF format files */
{
  struct command_parameter* cp;
  char* filename = NULL;
  FILE* inout;

  if (command_par("file", cmd->clone, &cp))
  {
    if ((filename = cp->string) == NULL)
    {
      if (cp->call_def != NULL)
        filename = cp->call_def->string;
    }
  }
  else
  {
    if (cp->call_def != NULL)
      filename = cp->call_def->string;
  }
  if (filename == NULL) filename = permbuff("dummy");
  if (strcmp(cmd->tok_list->p[0], "sxfread") == 0)
  {
    if (down_unit(filename) == 0)  return;
    sxf_read(cmd->clone);
    fclose(in->input_files[in->curr--]);
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

