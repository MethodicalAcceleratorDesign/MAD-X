 #include "madx.h"

static struct element*
new_element(const char* name)
{
  const char *rout_name = "new_element";
  struct element* el = mycalloc(rout_name, 1, sizeof *el);
  el->aper =           mycalloc(rout_name, 1, sizeof *el->aper);
  el->aper->aperture = mycalloc(rout_name, 4, sizeof *el->aper->aperture);
  el->aper->aper_offset = mycalloc(rout_name, 2, sizeof *el->aper->aper_offset);
  strcpy(el->name, name);
  el->stamp = 123456;
  el->def = 0x0;
  el->parent = 0x0;
  el->base_type = 0x0;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", el->name);
  return el;
}

void
grow_el_list(struct el_list* p)
{
  const char *rout_name = "grow_el_list";
  p->max *= 2;
  p->elem = myrecalloc(rout_name, p->elem, p->curr * sizeof *p->elem, p->max * sizeof *p->elem);
}

#if 0 // not used...
static void
dump_el_list(struct el_list* ell)
{
  for (int i = 0; i < ell->curr; i++)
    dump_element(ell->elem[i]);
}
#endif

static void
export_element(struct element* el, struct el_list* ell, FILE* file, int noexpr)
  /* recursive to have parents always in front for MAD-8 */
{
  int pos = name_list_pos(el->name, ell->list);
  char out[AUX_LG];
  if (pos >= 0)
  {
    if (ell->list->inform[pos] == 0)  /* not yet written */
    {
      export_element(el->parent, ell, file, noexpr);
      strcpy(out, el->name);

      strcat(out, ": ");
      strcat(out, el->parent->name);
      export_el_def(el, out, noexpr);
      write_nice(out, file);
      ell->list->inform[pos] = 1;
    }
  }
}

static void
export_elem_8(struct element* el, struct el_list* ell, FILE* file)
  /* exports an element in mad-8 format */
  /* recursive to have parents always in front for MAD-8 */
{
  int pos = name_list_pos(el->name, ell->list);
  char out[AUX_LG];
  if (pos >= 0)
  {
    if (ell->list->inform[pos] == 0)  /* not yet written */
    {
      export_elem_8(el->parent, ell, file);
      strcpy(out, el->name);
      strcat(out, ": ");
      strcat(out, el->parent->name);
      export_el_def_8(el, out);
      write_nice_8(out, file);
      ell->list->inform[pos] = 1;
    }
  }
}

static void
export_el_par_8(struct command_parameter* par, char* string)
  /* exports an element parameter in mad-8 format */
{
  int i, k, last, vtilt = 0;
  char num[2*NAME_L], tmp[16], tmpt[16];
  switch(par->type)
  {
    case 0:
      strcat(string, ",");
      strcat(string, par->name);
      strcat(string, " =");
      if (par->double_value == zero) strcat(string, "false");
      else                           strcat(string, "true");
      break;
    case 1:
    case 2:
      strcat(string, ",");
      strcat(string, par->name);
      strcat(string, "=");
      if (par->expr != NULL && strcmp(par->name, "harmon") != 0)
        strcat(string, par->expr->string);
      else
      {
        if (par->type == 1)
        {
          k = par->double_value; sprintf(num, v_format("%I"), k);
        }
        else sprintf(num, v_format("%F"), par->double_value);
        strcat(string, supp_tb(num));
      }
      break;
    case 3:
      if (par->string)
      {
        strcat(string, ",");
        strcat(string, par->name);
        strcat(string, "=");
        strcat(string, par->string);
      }
      break;
    case 11:
    case 12:
      vtilt = strcmp(par->name, "ks") == 0 ? 1 : 0;
      for (last = par->double_array->curr-1; last > 0; last--)
      {
        if (par->expr_list->list[last] != NULL)
        {
          if (zero_string(par->expr_list->list[last]->string) == 0) break;
        }
        else if (par->double_array->a[last] != zero) break;
      }
      for (i = 0; i <= last; i++)
      {
        if (par->expr_list->list[i] != NULL
            && !zero_string(par->expr_list->list[i]->string))
        {
          strcat(string, ",");
          sprintf(tmp, " k%dl =", i);
          sprintf(tmpt, ", t%d", i);
          strcat(string, tmp);
          strcat(string, par->expr_list->list[i]->string);
          if (vtilt) strcat(string, tmpt);
        }
        else if (par->double_array->a[i] != zero)
        {
          strcat(string, ",");
          sprintf(tmp, " k%dl =", i);
          sprintf(tmpt, ", t%d", i);
          if (par->type == 11)
          {
            k = par->double_array->a[i]; sprintf(num, "%d", k);
          }
          else sprintf(num, v_format("%F"), par->double_array->a[i]);
          strcat(string, tmp);
          strcat(string, supp_tb(num));
          if (vtilt) strcat(string, tmpt);
        }
      }
  }
}

static void
enter_elm_reference(struct in_cmd* cmd, struct element* el, int flag, int isupdating)
  /* enters an element in a sequence */
{
  int i, nupdates=1, k = 1;
  double at;

  if (strcmp(el->base_type->name, "rfcavity") == 0 &&
      find_element(el->name, current_sequ->cavities) == NULL)
    add_to_el_list(&el, 0, current_sequ->cavities, 0);
  if (!par_present("at", cmd->clone))
    fatal_error("element reference without 'at':",
                join(cmd->tok_list->p, cmd->tok_list->curr));
  at = command_par_value("at", cmd->clone);
  if ((i = name_list_pos(el->name, occ_list)) < 0)
    i = add_to_name_list(el->name, k, occ_list);
  else if (flag)
    fatal_error("multiple element definition inside sequence:", el->name);
  else k = ++occ_list->inform[i];
  make_elem_node(el, k);
  current_node->at_value = at;
  current_node->at_expr = command_par_expr("at", cmd->clone);
  const char* from = command_par_string_user("from", cmd->clone);
  if (from){
    current_node->from_name = permbuff(from);
    nupdates = 2;
  }
    if (isupdating==0) check_for_update_in_seq(el, cmd->clone, nupdates);
}

static int
par_out_flag(char* base_name, char* par_name)
{
  /* marks the element parameters that are to be written on "save" */
  if (strcmp(par_name,"at") == 0 || strcmp(par_name,"from") == 0) return 0;
  if (strcmp(base_name, "multipole") == 0
      && strcmp(par_name,"l") == 0) return 0;

  return 1;
}

// public interface

char*
compound(char* e_name, int occ)
  /* makes node name from element name and occurrence count */
{
  sprintf(c_dum->c,"%s:%d", e_name, occ);
  return c_dum->c;
}

struct node*
new_elem_node(struct element* el, int occ_cnt)
{
  struct node* p;
  p = new_node(compound(el->name, occ_cnt));
  p->p_elem = el;
  p->length = el->length;
  p->base_name = el->base_type->name;
  p->occ_cnt = occ_cnt;
  return p;
}

struct el_list*
new_el_list(int length)
{
  const char *rout_name = "new_el_list";
  struct el_list* ell = mycalloc(rout_name, 1, sizeof *ell);
  strcpy(ell->name, "el_list");
  ell->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", ell->name);
  ell->list = new_name_list(ell->name, length);
  ell->elem = mycalloc(rout_name, length, sizeof *ell->elem);
  ell->max = length;
  return ell;
}

struct element*
clone_element(struct element* el)
{
  struct element* clone = new_element(el->name);
  clone->length = el->length;
  clone->bv = el->bv;
  clone->def = el->def;
  clone->parent = el;
  clone->base_type = el->base_type;
  return clone;
}

struct element*
make_element(const char* name, const char* parent, struct command* def, int flag)
  /* makes a new element from declaration, stores in list */
{
/*  double length; */
  struct element* el = new_element(name);
  el->def = def;
  if (strcmp(name, parent) == 0)  /* basic element type like drift etc. */
  {
    add_to_el_list(&el, def->mad8_type, base_type_list, 1);
    el->parent = el->base_type = el;
  }
  else
  {
    if((el->parent = find_element(parent, element_list)) == NULL)
      fatal_error("unknown class type:", parent);
    el->base_type = el->parent->base_type;
    if(command_par_value("l",def) !=0 && belongs_to_class(el,"multipole"))
      warning("Multipole defined with non-zero length:", el->name);
    el->length = el_par_value("l", el);
    set_aperture_element(el, def);
  }
  
  add_to_el_list(&el, def->mad8_type, element_list, flag);
  return el;
}
void set_aperture_element(struct element *el, struct command* def){
  char *type;
//enum en_apertype{circle, ellipse, rectangle, lhcscreen, rectcircle, rectellipse, racetrack, octagon};
  type = command_par_string("apertype", def);
  el->aper->custom_inter = 0; 
  if(type!=NULL){
    if(strcmp(type,"circle")==0){
      
      double vector [4]; 
      element_vector(el,"aperture", vector);
      if(vector[0] > ten_m_12)
        el->aper->apertype = circle;
      else
        el->aper->apertype = notdefined;
    }
    else if(strcmp(type,"ellipse")==0)
      el->aper->apertype = ellipse;
    else if(strcmp(type,"rectangle")==0)
      el->aper->apertype = rectangle;
    else if(strcmp(type,"lhcscreen")==0)
      el->aper->apertype = lhcscreen;
    else if(strcmp(type,"rectcircle")==0)
      el->aper->apertype = rectcircle;
    else if(strcmp(type,"rectellipse")==0)
      el->aper->apertype = rectellipse;
    else if(strcmp(type,"racetrack")==0)
      el->aper->apertype = racetrack;
    else if(strcmp(type,"octagon")==0)
      el->aper->apertype = octagon;
    else{
      el->aper->apertype = custom;
      int lines=0, ch;
      FILE *fp = fopen(type,"r");
      if(fp==NULL){
          fatal_error("Aperture File is not existing ",type);
        }

      while(!feof(fp))
      {
        ch = fgetc(fp);
        if(ch == '\n'){
          lines++;
        }
      }
      
      el->aper->xlist = mycalloc("aperlist", lines+1, sizeof *el->aper->xlist);
      el->aper->ylist = mycalloc("aperlist", lines+1, sizeof *el->aper->ylist);
      rewind(fp);
      int i=0;
      while (2==fscanf(fp, "%lf %lf", &el->aper->xlist[i], &el->aper->ylist[i])) i++;
      /* closing the shape: a last point is inserted in table
     with coordinates equal to those of the first point */
    el->aper->length = i; // this minus 1 has to be there because of how the algorithm is done.  
    el->aper->xlist[i]=el->aper->xlist[0];
    el->aper->ylist[i]=el->aper->ylist[0];   
    fclose(fp);
    }
  }

  element_vector(el, "aperture" ,el->aper->aperture);
  element_vector(el, "aper_offset",el->aper->aper_offset);


  double tmpx [MAXARRAY];
  double tmpy [MAXARRAY];
  for(int i=0;i<MAXARRAY;i++){
    tmpx[i] = -999;
    tmpy[i] = -999;
  }
  
  int lx = element_vector(el, "aper_vx", tmpx);
  int ly = element_vector(el, "aper_vy", tmpy);
  int tmp_l=MAXARRAY+1;
  if(tmpx[0]!=-1 && ly > 1 && lx >1){
    for(int i=0;i<MAXARRAY;i++){
      if(tmpx[i]==-999 && tmpy[i]==-999){
        tmp_l = i;
        break; 
      }
    }

    if(tmp_l > MAXARRAY){
      mad_error("Different length of aper_vx and aper_vy for element:",el->name);
    }
    else{
      el->aper->custom_inter = 1;//sets the flagg that it should be used
      el->aper->xlist = mycalloc("aperlist", tmp_l+1, sizeof *el->aper->xlist);
      el->aper->ylist = mycalloc("aperlist", tmp_l+1, sizeof *el->aper->ylist);

      for(int i=0;i<tmp_l;i++){
        el->aper->xlist[i] = tmpx[i];
        el->aper->ylist[i] = tmpy[i];
      }
      //printf("2nd last %f, and last %f %d", el->aper->xlist[tmp_l-2], el->aper->xlist[tmp_l-1], tmp_l);


      el->aper->length = tmp_l; // minus 1 or not ?? has to be there because of how the algorithm is done.  
      el->aper->xlist[tmp_l]=el->aper->xlist[0];
      el->aper->ylist[tmp_l]=el->aper->ylist[0];
      if(el->aper->apertype==notdefined){ //If no other aperture is defined then a 10 meter rectangle is set! 
        el->aper->apertype=custom_inter; // sets it to a rcircle so the check is still done
      }
    }
  }

}

void update_node_aperture(void){
  char *type;
//enum en_apertype{circle, ellipse, rectangle, lhcscreen, rectcircle, rectellipse, racetrack, octagon};
  type = command_par_string("apertype", current_node->p_elem->def);
  if(type!=NULL && current_node->p_elem->aper->apertype!=custom_inter){
    if(strcmp(type,"circle")==0){
      
      double vector [4]; 
      element_vector(current_node->p_elem,"aperture", vector);
      if(vector[0] > ten_m_12)
        current_node->p_elem->aper->apertype = circle;
      else
        current_node->p_elem->aper->apertype = notdefined;
    }
    else if(strcmp(type,"ellipse")==0)
      current_node->p_elem->aper->apertype = ellipse;
    else if(strcmp(type,"rectangle")==0)
      current_node->p_elem->aper->apertype = rectangle;
    else if(strcmp(type,"lhcscreen")==0)
      current_node->p_elem->aper->apertype = lhcscreen;
    else if(strcmp(type,"rectcircle")==0)
      current_node->p_elem->aper->apertype = rectcircle;
    else if(strcmp(type,"rectellipse")==0)
      current_node->p_elem->aper->apertype = rectellipse;
    else if(strcmp(type,"racetrack")==0)
      current_node->p_elem->aper->apertype = racetrack;
    else if(strcmp(type,"octagon")==0)
      current_node->p_elem->aper->apertype = octagon;
  }

  element_vector(current_node->p_elem, "aperture", current_node->p_elem->aper->aperture);
  element_vector(current_node->p_elem, "aper_offset",current_node->p_elem->aper->aper_offset);

  if(current_node->p_elem->aper->custom_inter ==1){

    element_vector(current_node->p_elem, "aper_vx", current_node->p_elem->aper->xlist);
    element_vector(current_node->p_elem, "aper_vy", current_node->p_elem->aper->ylist);
    
  }
}

int is_custom_set(void){

  return current_node->p_elem->aper->custom_inter;
}

void
make_elem_node(struct element* el, int occ_cnt)
  /* makes + links a new node at the end of the current sequence */
{
  prev_node = current_node;
  current_node = new_elem_node(el, occ_cnt);
  current_node->occ_cnt = occ_cnt;
  current_node->chkick = el_par_value("chkick", el);
  current_node->cvkick = el_par_value("cvkick", el);
  add_to_node_list(current_node, 0, current_sequ->nodes);

  if (prev_node != NULL) prev_node->next = current_node;
  current_node->previous = prev_node;
  current_node->next = NULL;

}

struct element*
delete_element(struct element* el)
{
  const char *rout_name = "delete_element";
  if (el == NULL)  return NULL;
  if (stamp_flag && el->stamp != 123456)
    fprintf(stamp_file, "d_e double delete --> %s\n", el->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", el->name);
  myfree(rout_name, el->aper);
  myfree(rout_name, el);
  return NULL;
}

struct el_list*
delete_el_list(struct el_list* ell)
{
  const char *rout_name = "delete_el_list";
  if (ell->list == NULL) return NULL;
  if (stamp_flag && ell->stamp != 123456)
    fprintf(stamp_file, "d_e_l double delete --> %s\n", ell->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", ell->name);
  delete_name_list(ell->list);
  if (ell->elem != NULL) myfree(rout_name, ell->elem);
  myfree(rout_name, ell);
  return NULL;
}

void
dump_element(struct element* el)
{
  fprintf(prt_file, v_format("+++ dumping element %S  parent %S\n"),
          el->name, el->parent->name);
  dump_command(el->def);
}

void
write_elems(struct el_list* ell, struct command_list* cl, FILE* file, int noexpr)
{
  int i;
  for (i = 0; i < ell->curr; i++)
  {
    if (pass_select_list_el(ell->elem[i], cl))
      export_element(ell->elem[i], ell, file, noexpr);
  }
}

void
write_elems_8(struct el_list* ell, struct command_list* cl, FILE* file)
{
  int i;
  for (i = 0; i < ell->curr; i++)
  {
    if (pass_select_list_el(ell->elem[i], cl))
      export_elem_8(ell->elem[i], ell, file);
  }
}

void
export_el_def(struct element* el, char* string, int noexpr)
  /* exports an element definition in mad-X format */
{
  int i;
  struct command* def = el->def;
  struct command_parameter* par;
  for (i = 0; i < def->par->curr; i++) {
    par = def->par->parameters[i];
    if (def->par_names->inform[i] && par_out_flag(el->base_type->name, par->name))
      export_comm_par(par, string, noexpr);
  }
}

void
export_el_def_8(struct element* el, char* string)
  /* exports an element definition in mad-8 format */
{
  struct command* def = el->def;
  struct command_parameter* par;
  int i, div = 1;
  double val[] = {0, 0, 0};
  char* base_name = el->base_type->name;
  char num[2*NAME_L];
  /* special treatment for tilt */
  for (i = 0; i < def->par->curr; i++)
  {
    par = def->par->parameters[i];
    if (def->par_names->inform[i])
    {
      if (strcmp(base_name, "quadrupole") == 0)
      {
        div = 2;
        if (strcmp(par->name,"k1") == 0) val[0] = command_par_special("k1", el);
        else if (strcmp(par->name,"k1s") == 0)
          val[1] = command_par_special("k1s", el);
        else if (strcmp(par->name,"tilt") == 0)
          val[2] = command_par_special("tilt", el);
        else if(par_out_flag(el->base_type->name, par->name))
          export_el_par_8(par, string);
      }
      else if (strcmp(base_name, "sextupole") == 0)
      {
        div = 3;
        if (strcmp(par->name,"k2") == 0) val[0] = command_par_special("k2", el);
        else if (strcmp(par->name,"k2s") == 0)
          val[1] = command_par_special("k2s", el);
        else if (strcmp(par->name,"tilt") == 0)
          val[2] = command_par_special("tilt", el);
        else if(par_out_flag(el->base_type->name, par->name))
          export_el_par_8(par, string);
      }
      else if (strcmp(base_name, "octupole") == 0)
      {
        div = 4;
        if (strcmp(par->name,"k3") == 0) val[0] = command_par_special("k3", el);
        else if (strcmp(par->name,"k3s") == 0)
          val[1] = command_par_special("k3s", el);
        else if (strcmp(par->name,"tilt") == 0)
          val[2] = command_par_special("tilt", el);
        else if(par_out_flag(el->base_type->name, par->name))
          export_el_par_8(par, string);
      }
      else if (strcmp(base_name, "elseparator") == 0)
      {
        if (strcmp(par->name,"ex") == 0) val[0] = command_par_special("ex", el);
        else if (strcmp(par->name,"ey") == 0)
          val[1] = command_par_special("ey", el);
        else if(par_out_flag(el->base_type->name, par->name))
          export_el_par_8(par, string);
      }
      else if(par_out_flag(el->base_type->name, par->name))
        export_el_par_8(par, string);
    }
  }
  if (val[1] != zero)
    val[2] = atan2(val[1], val[0]) / div;
  if (val[2] != zero) val[0] = sqrt(val[0]*val[0]+val[1]*val[1]);
  if (val[0] != zero)
  {
    strcat(string, ",");
    if (strcmp(base_name, "quadrupole") == 0)
    {
      strcat(string, "k1 =");
    }
    else if (strcmp(base_name, "sextupole") == 0)
    {
      strcat(string, "k2 =");
    }
    else if (strcmp(base_name, "octupole") == 0)
    {
      strcat(string, "k3 =");
    }
    else if (strcmp(base_name, "elseparator") == 0)
    {
      strcat(string, "e =");
    }
    sprintf(num, v_format("%F"), val[0]);
    strcat(string, supp_tb(num));
    if (val[2] != zero)
    {
      strcat(string, ",tilt =");
      sprintf(num, v_format("%F"), val[2]);
      strcat(string, supp_tb(num));
    }
  }
}

int
belongs_to_class(struct element* el, const char* class)
  /* returns 1 if an element belongs to a class, else 0 */
{
  int in = 0;
  if (strcmp(el->name, class) == 0) in = 1;
  else
  {
    while (el->parent != el)
    {
      if (strcmp(el->parent->name, class) == 0)
      {
        in = 1; break;
      }
      el = el->parent;
    }
  }
  return in;
}

void
element_name(char* name, int* l)
  /* returns current node element name in Fortran format */
  /* l is max. allowed length in name */
{
  strfcpy(name, current_node->p_elem->name, *l);
}

double
element_value(const struct node* node, const char* par)
  /* all element parameter values except vectors are modified here
     resp. in el_par_value if any */
{
  double e_val;

  if (node == 0) {
     mad_error("element_value","node parameter is NULL.");
     return 0.0;
   }

  const struct element* el = node->p_elem;

   if (el == 0) {
     mad_error("element_value","node has NULL element pointer.");
     return 0.0;
   }

   if (strcmp(el->name,"in_cmd") == 0) {
     mad_error("element_value","node '%.47s' refers to invalid element (improper (re)definition?).", node->name);
     return 0.0;
   }

   const struct command* def = el->def;

   if (def == 0) {
     mad_error("element_value","element has NULL defintion pointer.");
     return 0.0;
   }

  if (strcmp(par, "mad8_type") == 0)
    e_val = el->def->mad8_type;
  else
    e_val = el_par_value(par, el);
  return e_val;
}

int
element_vector(const struct element* el, const char* par, double* vector)
  /* returns length + vector of parameter par for element el */
{
  int i, l = 0;
  struct double_array* da;
  struct expr_list* ell;
  if ((i = name_list_pos(par, el->def->par_names)) > -1)
  {
    if ((da = el->def->par->parameters[i]->double_array) != NULL)
    {
      if ((ell = el->def->par->parameters[i]->expr_list) != NULL)
        update_vector(ell, da);
      l = da->curr;
      copy_double(da->a, vector, l);
    }
  }
  return l;
}

void
get_node_vector(const char* par, int* length, double* vector)
  /* returns vector for parameter par of current element */
{
  char lpar[NAME_L];
  mycpy(lpar, par);

  if (strcmp(lpar, "orbit0") == 0) copy_double(orbit0, vector, 6);
  else if (strcmp(lpar, "obs_orbit") == 0)
  {
    if (current_node->obs_orbit)
    {
      *length = current_node->obs_orbit->curr;
      copy_double(current_node->obs_orbit->a, vector, *length);
    }
    else *length = 0;
  }
  else if (strcmp(lpar, "orbit_ref") == 0)
  {
    if (current_node->orbit_ref)
    {
      *length = current_node->orbit_ref->curr;
      copy_double(current_node->orbit_ref->a, vector, *length);
    }
  }
  else if (strcmp(lpar, "surv_data") == 0)
    {
     copy_double(current_node->surv_data, vector, 7);
     *length = 7;
    }
  else
   {
//     printf("get_node_vector: going to element_vector\n");
     *length = element_vector(current_node->p_elem, lpar, vector);
/*     printf("get_node_vector: got %d elements\n",*length);
       int i;
     for (i =0; i < *length; i++)
      {
        printf("%f ",vector[i]);
      }
     printf("\n"); */
   }
}

double
el_par_value(const char* par, const struct element* el)
  /* returns an element parameter value */
{
  int k = 0; // , n; not used
  char tmp[8];
  double val = zero, angle = zero, l, vec[100];
  double fact = strcmp(el->base_type->name, "rbend") == 0 ? one : zero;
  int mult = strcmp(el->base_type->name, "multipole") == 0 ? 1 : 0;
  int mark = strcmp(el->base_type->name, "marker") == 0 ? 1 : 0;
  if (fact != zero || !strcmp(el->base_type->name, "sbend")) /* bend */
  {
    if ((l = command_par_value("l", el->def)) == zero) {
//      printf("*** el_par_value(mad_elem): %s, %s, %s, %g\n", el->name, el->def->name, el->base_type->name, l);
      fatal_error("bend with zero length:",el->name);
    }

    angle = command_par_value("angle", el->def);
         if (!strcmp(par, "angle")) val = angle;
    else if (!strcmp(par, "tilt") ) val = command_par_value("tilt", el->def);
    else if (!strcmp(par, "k0")   ) val = command_par_value("k0"  , el->def);
    else if (!strcmp(par, "k0s")  ) val = command_par_value("k0s" , el->def);
    else if (!strcmp(par, "l")    )
    {
      if (fact != zero && get_option("rbarc") && fabs(angle) > 1e-8)
        val = l * angle / (two * sin(angle/two)); // l_arc = l_cord / sinc(angle/2)
      else val = l;
    }
    else if (strcmp(par, "e1") == 0)
      val = command_par_value("e1", el->def);/* + fact * angle / two; dipole_bv kill initiative SF TR FS */
    else if (strcmp(par, "e2") == 0)
      val = command_par_value("e2", el->def);/* + fact * angle / two; dipole_bv kill initiative SF TR FS */
    else if (strcmp(par, "rhoinv") == 0) val = angle / l;
    else if (strcmp(par, "blen") == 0) val = l;
    else val = command_par_value(par, el->def);
  }
  /* all elements except bends */
  else if (strcmp(par, "rhoinv") == 0) val = zero;
  else if (strcmp(par, "blen") == 0) val = zero;
  else if (mark) /* marker */
  {
    if ((l = command_par_value("l", el->def)) != zero)
      fatal_error("marker with nonzero length:",el->name);
    val = command_par_value(par, el->def);
  }
  else if (mult)  /* multipole */
  {
    if (strcmp(par, "l") == 0) val = zero;
    else if (par[0] == 'k' && isdigit(par[1]) && par[strlen(par)-1] == 'l')
      /* single component requested for multipole */
    {
      if (strchr(par, 's')) strcpy(tmp, "ksl");
      else                  strcpy(tmp, "knl");
      sscanf(&par[1], "%d", &k);
      if (element_vector(el, tmp, vec) > k)  val = vec[k]; // n = not used
    }
    else val = command_par_value(par, el->def);
  }
  else val = command_par_special(par, el);
  /* extra code for kickers */
  if (val == zero && strcmp(el->base_type->name, "hkicker") == 0)
  {
    if (strcmp(par,"hkick") == 0)
      val = command_par_value("kick", el->def);
    else if (strcmp(par,"kick") == 0)
      val = command_par_value("hkick", el->def);
  }
  else if (val == zero && strcmp(el->base_type->name, "vkicker") == 0)
  {
    if (strcmp(par,"vkick") == 0)
      val = command_par_value("kick", el->def);
    else if (strcmp(par,"kick") == 0)
      val = command_par_value("vkick", el->def);
  }
  return val;
}

int
el_par_vector(int* total, double* vect)
/* returns a complete vector of element parameters for current node */
/* the integer return value is the length */
{
  struct command* elc = current_node->p_elem->def;
  struct command_parameter_list* parl = elc->par;
  struct command_parameter* cp;
  int i, len = 0;
  double val;
  for (i = 0; i < *total; i++) {
     if (i < elc->par->curr) {
        cp = parl->parameters[i];
        if (cp->type < 3) {
          if (cp->expr == NULL) val = cp->double_value;
          else val = expression_value(cp->expr, 2);
          vect[len++] = val;
	      }
        else val = 0;
     }
  }
  return len;
}

/* returns parameter if it has been modified, otherwise NULL */
struct command_parameter*
return_param(const char* par, const struct element* elem)
{
  /* don't return base type definitions */
  if (elem==elem->parent) return NULL;

  struct command_parameter* cp;
  return command_par(par, elem->def, &cp) ? cp : NULL;
}

/* returns parameter if it has been modified, otherwise NULL  - recursively */
struct command_parameter*
return_param_recurse(const char* par, const struct element* elem)
{
  struct command_parameter* param;
  param = return_param(par,elem);

  if (param) return param;
  if (elem!=elem->parent)
    return return_param_recurse(par,elem->parent);
  return NULL;
}

/* returns first parameter value found and recusively checks sub_elements */
double
el_par_value_recurse(const char* par, const struct element* elem)
{
  if (return_param(par, elem)) return el_par_value(par,elem);
  if (elem != elem->parent)
    return el_par_value_recurse(par,elem->parent);
  return 0;
}

void
enter_element(struct in_cmd* cmd)
  /* enters an element in the list (and the sequence if applicable) */
{
//  struct command_parameter_list* pl; // not used
  char** toks = cmd->tok_list->p;
  struct element *el, *parent;
  struct command* comm;
  int flag = 0, k = cmd->type == 1 ? 2 : 0;
  if ((parent = find_element(toks[k], element_list)) == NULL)
  {
    fatal_error("unknown class type:", toks[k]);
  }
  else
  {
    cmd->cmd_def = parent->def;
    cmd->clone = clone_command(cmd->cmd_def);
    strcpy(cmd->clone->name, toks[0]);
    scan_in_cmd(cmd);
    if (k == 0 || strcmp(toks[0], toks[2]) == 0) {
      el = parent;
    }
    else
    {
      if ((el = make_element(toks[0], parent->name,
                             cmd->clone, 1+sequ_is_on)) == NULL) return;
      el->def_type = sequ_is_on;
      flag = 1; /* new element - definition only once in sequence allowed */
    }
    // pl = cmd->clone->par; // not used
    if (el != parent)
    {
      if (par_present("bv", cmd->clone)) el->bv = command_par_value("bv", cmd->clone);
      else if ((comm = find_command(el->parent->name, defined_commands))
              && par_present("bv", comm))
        el->bv = command_par_value("bv", comm);
      else el->bv = parent->bv;
    }
    if (sequ_is_on) enter_elm_reference(cmd, el, flag, k);
  }
}

void
fill_elem_var_list(struct element* el, struct el_list* ell, struct var_list* varl)
  /* puts all variables an element depends on, in a list */
{
  struct command* cmd = el->def;
  int i;
  for (i = 0; i < cmd->par->curr; i++)
    fill_par_var_list(ell, cmd->par->parameters[i], varl);
}

struct element*
find_element(const char* name, struct el_list* ell)
{
  int pos;
  if ((pos = name_list_pos(name, ell->list)) < 0)
    return NULL;
  return ell->elem[pos];
}

void
check_for_update_in_seq(struct element* el, struct command* update, int nupdates)
  /* checks if someone tries to update the element in sequence creation */
{
  struct command_parameter_list* e_pl = el->def->par;
  int pos, cupdate=0;
  for (pos = 0; pos < e_pl->curr; pos++)
  {
    if (update->par_names->inform[pos])  /* parameter has been read */
    {
      cupdate++;
      if(cupdate>nupdates)
        fatal_error("Not possible to update attribute for element in sequence definition: ", el->name );
    }
  }
}

void
update_element(struct element* el, struct command* update)
  /* updates the parameters of el from those read into update */
{
  struct command_parameter_list* e_pl = el->def->par;
  struct command_parameter_list* pl = update->par;
  struct command_parameter *e_par, *par;
  int pos;
  for (pos = 0; pos < e_pl->curr; pos++)
  {
    if (update->par_names->inform[pos])  /* parameter has been read */
    {
      el->def->par_names->inform[pos]=update->par_names->inform[pos]; /*hbu activate this parameter in the element */
      e_par = e_pl->parameters[pos];
      par = pl->parameters[pos];
      switch (par->type)
      {
        case 0:
        case 1:
        case 2:
          e_par->double_value = par->double_value;
          e_par->expr = clone_expression(par->expr);
          /* fix for bv flag start */
          if (strcmp(e_par->name, "bv") == 0)
            el->bv = e_par->double_value;
          /* fix for bv flag end */
          break;
        case 3:
          e_par->string = permbuff(par->string);
          break;
        case 11:
        case 12:
          e_par->double_array = clone_double_array(par->double_array);
          e_par->expr_list = clone_expr_list(par->expr_list);
      }
    }
  }
  set_aperture_element(el, update); //updates contains all the info of the element
}

void
update_element_children(struct element* el, struct command* update)
  /* updates the parameters of the children to el. Note that it is only updating one layer (not recursive) */
{

  for(int i=0; i<element_list->max;i++){
    if(element_list->elem[i]==NULL) break;

    if(strcmp(el->name,element_list->elem[i]->parent->name)==0)
      update_element(element_list->elem[i], update);
  }
}
void
add_to_el_list( /* adds element to alphabetic element list */
  struct element** el, int inf, struct el_list* ell, int flag)
  /* inf is entered in the namelist */
  /*  flag < 0: do not delete if already present, do not warn */
  /*       = 0: delete, but do not warn */
  /*       = 1: delete & warn */
  /*       = 2: warn and ignore if already present - resets *el to old */
{
  int pos, j;
  struct node* p_node;
  if ((pos = name_list_pos((*el)->name, ell->list)) > -1)
  {
    if (flag > 1)
    {
      warning("implicit element re-definition ignored:", (*el)->name);
      *el = ell->elem[pos];
    }
    else
    {
      if (flag > 0)
      {
        put_info("element redefined:", (*el)->name);
/*
 *         printf("File %s line %d\n",filenames[in->curr], currentline[in->curr] );
 *         printf("Old Definition:\n");
 *         dump_element(ell->elem[pos]);
 *         printf("New Definition:\n");
 *         dump_element(*el);
 */

      }
      if (flag >= 0 && ell == element_list)
      {
        for (j = 0; j < ell->curr; j++) /* set existing pointers to new */
        {
          if (ell->elem[j] != ell->elem[pos]
              && ell->elem[j]->parent == ell->elem[pos])
            ell->elem[j]->parent = *el;
        }
        for (j = 0; j < sequences->curr; j++)
        {
          p_node = sequences->sequs[j]->start;
          while (p_node && p_node != sequences->sequs[j]->end)
          {
            if (p_node->p_elem == ell->elem[pos]) p_node->p_elem = *el;
            p_node = p_node->next;
          }
          if (strcmp((*el)->base_type->name, "rfcavity") == 0 &&
	      find_element((*el)->name, sequences->sequs[j]->cavities) != NULL)
	    sequences->sequs[j]->cavities->elem[name_list_pos((*el)->name,
                                                              sequences->sequs[j]->cavities->list)] = *el;
        }
        delete_element(ell->elem[pos]);
      }
      ell->elem[pos] = *el;
    }
  }
  else
  {
    if (ell->curr == ell->max) grow_el_list(ell);
    j = add_to_name_list(permbuff((*el)->name), inf, ell->list);
    ell->elem[ell->curr++] = *el;
  }
}

