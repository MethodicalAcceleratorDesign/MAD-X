#include "madx.h"
#include "mad_mac.h"

static void
remove_from_node_list(struct node* node, struct node_list* nodes)
{
  int i;
  if ((i = remove_from_name_list(node->name, nodes->list)) > -1)
    nodes->nodes[i] = nodes->nodes[--nodes->curr];
}

// public interface

struct node*
new_node(char* name)
{
  const char *rout_name = "new_node";
  struct node* p = mycalloc(rout_name, 1, sizeof *p);
  strcpy(p->name, name);
  p->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", p->name);
  return p;
}

struct node*
clone_node(struct node* p, int flag)
{
  /* Transfers errors from original nodes if flag != 0;
     this is needed for SXF input  */
  struct node* clone = new_node(p->name);
  strcpy(clone->name,p->name);

  /*clone->base_name = p->base_name;*/
  clone->base_name = permbuff(p->base_name);

  clone->occ_cnt = p->occ_cnt;
  clone->sel_err = p->sel_err;
  clone->position = p->position;
  clone->at_value = p->at_value;
  clone->length = p->length;
  clone->at_expr = p->at_expr;
  clone->from_name = p->from_name;
  clone->p_elem = p->p_elem;
  clone->p_sequ = p->p_sequ;
  clone->savebeta = p->savebeta;
  if (flag) {
    clone->p_al_err = p->p_al_err;
    clone->p_fd_err = p->p_fd_err;
    // AL: RF-Multipole errors (EFCOMP)
    clone->p_ph_err = p->p_ph_err;
  }
  // AL: RF-Multipole
  //clone->pnl = p->pnl;
  //clone->psl = p->psl;
  clone->rfm_lag = p->rfm_lag;
  clone->rfm_freq = p->rfm_freq;
  clone->rfm_volt = p->rfm_volt;
  clone->rfm_harmon = p->rfm_harmon;

  clone->chkick = p->chkick;
  clone->cvkick = p->cvkick;
  return clone;
}

struct node*
delete_node(struct node* p)
{
  const char *rout_name = "delete_node";
  if (p == NULL) return NULL;
  if (stamp_flag && p->stamp != 123456)
    fprintf(stamp_file, "d_n double delete --> %s\n", p->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", p->name);
  if (p->p_al_err) p->p_al_err = delete_double_array(p->p_al_err);
  if (p->p_fd_err) p->p_fd_err = delete_double_array(p->p_fd_err);
  // AL: RF-Multipole errors (EFCOMP)
  if (p->p_ph_err) p->p_ph_err = delete_double_array(p->p_ph_err);
  // AL: RF-Multipole
  // if (p->pnl) p->pnl = delete_double_array(p->pnl);
  // if (p->psl) p->psl = delete_double_array(p->psl);
  //
  myfree(rout_name, p);
  return NULL;
}

struct node_list*
new_node_list(int length)
{
  const char *rout_name = "new_node_list";
  struct node_list* nll = mycalloc(rout_name, 1, sizeof *nll);
  strcpy(nll->name, "node_list");
  nll->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", nll->name);
  nll->list = new_name_list(nll->name, length);
  nll->nodes = mycalloc(rout_name, length, sizeof *nll->nodes);
  nll->max = length;
  return nll;
}

struct node_list*
delete_node_list(struct node_list* l)
{
  const char *rout_name = "delete_node_list";
  if (l == NULL)  return NULL;
  if (stamp_flag && l->stamp != 123456)
    fprintf(stamp_file, "d_no_l double delete --> %s\n", l->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", l->name);
  if (l->nodes != NULL)  myfree(rout_name, l->nodes);
  if (l->list != NULL)  delete_name_list(l->list);
  myfree(rout_name, l);
  return NULL;
}

struct node*
delete_node_ring(struct node* start)
{
  struct node *p, *q;
  if (start == NULL) return NULL;
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", "node_ring");
  q = start->next;
  while (q != NULL && q != start)
  {
    p = q; q = q->next;
    /*printf("delete_node_ring deleting node %#x %s\n",p,p->base_name);*/
    delete_node(p);
  }

  /*printf("delete_node_ring deleting start node %#x %s\n",start,start->base_name);*/
  start = delete_node(start);

  return NULL;
}

static void
grow_node_list(struct node_list* p)
{
  const char *rout_name = "grow_node_list";
  p->max *= 2;
  if (p->max == 0) p->max++;
  p->nodes = myrecalloc(rout_name, p->nodes, p->curr * sizeof *p->nodes, p->max * sizeof *p->nodes);
}

void
dump_node(struct node* node)
{
  int i;
  char pname[NAME_L] = "NULL", nname[NAME_L] = "NULL",
    from_name[NAME_L] = "NULL";
  if (node->previous != NULL) strcpy(pname, node->previous->name);
  if (node->next != NULL) strcpy(nname, node->next->name);
  if (node->from_name != NULL) strcpy(from_name, node->from_name);
  fprintf(prt_file,
          v_format("name: %S  occ: %I base: %S  from_name: %S at_value: %F  position: %F\n"),
          node->name, node->occ_cnt, node->base_name, from_name,
          node->at_value, node->position);
  fprintf(prt_file, v_format("  names of - previous: %S  next: %S\n"),
          pname, nname);
  if (node->cl != NULL)  for (i = 0; i < node->cl->curr; i++)
    dump_constraint(node->cl->constraints[i]);
}

struct node*
expand_node(struct node* node, struct sequence* top_sequ, struct sequence* sequ, double position)
  /* replaces a (sequence) node by a sequence of nodes - recursive */
{
  struct sequence* nodesequ = node->p_sequ;
  struct node *p, *q = nodesequ->start;
  int i;
  (void)sequ;

  double node_pos, refpos;

  int debug;
  debug=get_option("debug");

  p = clone_node(q, 0);

  if (debug)
    printf("\n Entering expand_node with node->name: %s and position %e\n",node->name,position);

  if ((i = name_list_pos(p->p_elem->name, occ_list)) < 0)
    i = add_to_name_list(p->p_elem->name, 1, occ_list);
  else strcpy(p->name, compound(p->p_elem->name, ++occ_list->inform[i]));

  add_to_node_list(p, 0, top_sequ->ex_nodes);
  p->previous = node->previous;
  p->previous->next = p;

  if (debug) {
    node_pos = get_node_pos(p, nodesequ);
    refpos = get_refpos(nodesequ);
    p->position = position + node_pos - refpos;
    printf("  elem_name %s at position %e + get_node_pos %e - get_refpos %e = p->position %e \n",
	   p->p_elem->name, position, node_pos, refpos, p->position);
  } else
    p->position = position + get_node_pos(p, nodesequ) - get_refpos(nodesequ);

  while (p != NULL) {
    if (q == nodesequ->end) break;

    if (strcmp(p->base_name, "rfcavity") == 0 &&
        find_element(p->p_elem->name, top_sequ->cavities) == NULL)
      add_to_el_list(&p->p_elem, 0, top_sequ->cavities, 0);

    p->next = clone_node(q->next, 0);
    p->next->previous = p;
    p = p->next;
    q = q->next;

    if (debug) {
      node_pos = get_node_pos(p, nodesequ);
      refpos = get_refpos(nodesequ);
      p->position = position + node_pos - refpos;
      printf("  elem name %s at position %e + get_node_pos %e - get_refpos %e = p->position %e \n",
	     p->p_elem->name, position, node_pos, refpos, p->position);
    } else
      p->position = position + get_node_pos(p, nodesequ) - get_refpos(nodesequ);


    if (p->p_sequ == NULL){ // simple element, not a subsequence
      add_to_node_list(p, 0, top_sequ->ex_nodes);
    }

    else { // subsequence
      if (p->p_sequ->refpos != NULL) { // REFPOS given for subsequence, ignore REFER of current sequence
	if (debug)
	  printf("\n\n Recursively expand sub-sequence %s with initial position %e, final position %e, length %e, ref_flag %d, refpos '%s'\n",
		 p->name, p->position, p->position - nodesequ->ref_flag*p->p_sequ->length/2.,
		 p->length, nodesequ->ref_flag, p->p_sequ->refpos);
	p = expand_node(p, top_sequ, p->p_sequ, p->position - nodesequ->ref_flag*p->p_sequ->length/2. );
	if (debug) printf("\n\n");
      }
      else { // no REFPOS given
	if (debug)
	  printf("\n\n Recursively expand sub_sequence  %s with position %e, length %e, ref_flag %d\n",
		 p->name, p->position, p->length, sequ->ref_flag);
	p = expand_node(p, top_sequ, p->p_sequ, p->position);
	if (debug) printf("\n\n");
      }
    }
  }

  delete_node(node);
  return p;
}

double
get_node_pos(struct node* node, struct sequence* sequ)
  /* recursive via call of hidden_node_pos */
  /* returns node position from declaration for expansion */
{
  double fact = 0.5 * sequ->ref_flag; /* element half-length offset */
  double pos, from = 0;

  if (node->at_expr == NULL) pos = node->at_value;
  else                       pos = expression_value(node->at_expr, 2);

  if (get_option("debug"))
    printf("   in get_node_pos: name: %s, pos: %e, fact: %e, length: %e, from_name: %s\n",
	   node->p_elem->name, pos, fact, node->length, node->from_name);

  pos += fact * node->length; /* make centre position */

  if (node->from_name != NULL) {
    if ((from = hidden_node_pos(node->from_name, sequ)) == INVALID)
      fatal_error("'from' reference to unknown element:", node->from_name);
  }

  pos += from;

  if (get_option("debug"))
    printf("\t in get_node_pos: name: %s, from: %e\t\t\t  ---> final pos: %e \n",
	   node->p_elem->name, from, pos);

  return pos;
}

double
get_refpos(struct sequence* sequ)
  /* returns the position of a refpos element of a sequence, or the sequence half-length */
{
  int i;
  if (sequ != NULL && sequ->refpos != NULL) {
    sprintf(c_dum->c, "%s:1", sequ->refpos);
    if ((i = name_list_pos(c_dum->c, sequ->nodes->list)) < 0)
      fatal_error("'refpos' reference to unknown element:", sequ->refpos);
    return get_node_pos(sequ->nodes->nodes[i], sequ);
  }

  /* 2013-Jul-19  20:54:49  ghislain: if no element was specified for the refpos,
     the ref position has to be the middle of the sequence and get_refpos should
     therefore return half the sequence length. (TRAC ticket #206) */
  // else return zero;
  else return sequ->length/2.;
}
double get_length_(void){
  return current_node->length;
}

double
node_value(const char* par)
  /* returns value for parameter par of current element */
{
  double value;
  char lpar[NAME_L];

  if (current_node == 0x0)
   {
     warning("node_value","current_node pointer in NULL. Sequence not set?");
     return 0.0;
   }

  mycpy(lpar, par);
  if (strcmp(lpar, "l") == 0) value = current_node->length;
/*  else if (strcmp(lpar, "dipole_bv") == 0) value = current_node->dipole_bv;*/
  else if (strcmp(lpar, "other_bv") == 0) value = current_node->other_bv;
  else if (strcmp(lpar, "chkick") == 0) value = current_node->chkick;
  else if (strcmp(lpar, "cvkick") == 0) value = current_node->cvkick;
  else if (strcmp(lpar, "obs_point") == 0) value = current_node->obs_point;
  else if (strcmp(lpar, "sel_sector") == 0) value = current_node->sel_sector;
  else if (strcmp(lpar, "enable") == 0) value = current_node->enable;
  else if (strcmp(lpar, "occ_cnt") == 0) value = current_node->occ_cnt;
  else if (strcmp(lpar, "pass_flag") == 0) value = current_node->pass_flag;
/* AL: added by A. Latina 16 Feb 2012 */
  else if (strcmp(lpar, "rfm_freq") == 0) value = current_node->rfm_freq;
  else if (strcmp(lpar, "rfm_volt") == 0) value = current_node->rfm_volt;
  else if (strcmp(lpar, "rfm_lag") == 0) value = current_node->rfm_lag;
  else if (strcmp(lpar, "rfm_harmon") == 0) value = current_node->rfm_harmon;
//
  else value =  element_value(current_node, lpar);
  return value;
}
double node_obs_point(void){
  return current_node->obs_point;
}


void set_tt_multipoles(int *maxmul){
  int tmp_n, tmp_s;
  double tmp_nv[*maxmul] ;
  double tmp_sv[*maxmul] ;
  current_node->p_elem->multip = mycalloc("alloc mult struct", 1, sizeof (*current_node->p_elem->multip));
  current_node->p_elem->multip->knl = mycalloc("alloc multip normal", *maxmul, sizeof (*current_node->p_elem->multip->knl));
  current_node->p_elem->multip->ksl = mycalloc("alloc multip skew"  , *maxmul, sizeof (*current_node->p_elem->multip->ksl));

  get_node_vector("knl", &tmp_n, tmp_nv);
  get_node_vector("ksl", &tmp_s, tmp_sv);
  current_node->p_elem->multip->nn = tmp_n;
  current_node->p_elem->multip->ns = tmp_s;
  
  for(int i=0;i<tmp_n;i++){
    current_node->p_elem->multip->knl[i] = tmp_nv[i];
  }
  for(int i=0;i<tmp_s;i++){
    current_node->p_elem->multip->ksl[i] = tmp_sv[i];
  }


}

void get_tt_multipoles(int *nn, double *knl, int *ns, double *ksl){
    nn[0]=current_node->p_elem->multip->nn;
    ns[0]=current_node->p_elem->multip->ns;
    for(int i=0;i<*nn;i++){
      knl[i] = current_node->p_elem->multip->knl[i];
    }
    for(int i=0;i<*ns;i++){
      ksl[i] = current_node->p_elem->multip->ksl[i];
    }


}
void alloc_tt_attrib(int *length){
  current_node->p_elem->tt_attrib = mycalloc("tmp_array_tt", (*length+1), sizeof (*current_node->p_elem->tt_attrib));
}

void set_tt_attrib(int *index, double *value){
  current_node->p_elem->tt_attrib[*index] = *value;
}

double get_tt_attrib(int *index){
  return current_node->p_elem->tt_attrib[*index];
}

void
link_in_front(struct node* new, struct node* el)
  /* links a node 'new' in front of a node 'el' */
{
  el->previous->next = new;
  new->previous = el->previous; new->next = el;
  el->previous = new;
}

double
hidden_node_pos(char* name, struct sequence* sequ) /*recursive */
  /* (for 'from' calculation:) returns the position of a node
     in the current or a sub-sequence (hence hidden) */
{
  double pos;
  struct node* c_node;
  char tmp[2*NAME_L];
  int i;

  strcpy(tmp, name);
  square_to_colon(tmp);

  if ((i = name_list_pos(tmp, sequ->nodes->list)) > -1)
    return get_node_pos(sequ->nodes->nodes[i], sequ);

  else {
    c_node = sequ->start;
    while (c_node != NULL) {
      if (c_node->p_sequ != NULL) {
        if ((pos = hidden_node_pos(name, c_node->p_sequ)) != INVALID) {
          pos += get_node_pos(c_node, sequ) - get_refpos(c_node->p_sequ);
          return pos;
        }
      }
      if (c_node == sequ->end) break;
      c_node = c_node->next;
    }
    return INVALID;
  }
}

void
resequence_nodes(struct sequence* sequ)
  /* resequences occurrence list */
{
  struct node* c_node = sequ->start;
  int i, cnt;

  while (c_node != NULL) {
    if (c_node->p_elem != NULL) {
      if ((i = name_list_pos(c_node->p_elem->name, occ_list)) < 0) {
        i = add_to_name_list(c_node->p_elem->name, 1, occ_list);
        cnt = 1;
      }
      else cnt = ++occ_list->inform[i];
      strcpy(c_node->name, compound(c_node->p_elem->name, cnt));
    }
    if (c_node == sequ->end) break;
    c_node = c_node->next;
  }
}

int
retreat_node(void)
  /* replaces current node by previous node; 0 = already at start, else 1 */
{
  if (current_node == current_sequ->range_start)  return 0;
  current_node = current_node->previous;
  return 1;
}
void store_orbit_correctors(void){


    restart_sequ();
  while(1){
    set_command_par_value("chkick",current_node->p_elem->def,current_node->chkick);
    set_command_par_value("cvkick",current_node->p_elem->def, current_node->cvkick);
    if (advance_node()==0) break;
      
  }

}

void
store_node_value(const char* par, double* value)
  /* stores value for parameter par at current node */
{
  char lpar[NAME_L];
  struct element* el = current_node->p_elem;

  mycpy(lpar, par);
  if (strcmp(lpar, "chkick") == 0) current_node->chkick = *value;
  else if (strcmp(lpar, "cvkick") == 0)current_node->cvkick = *value;
/*  else if (strcmp(lpar, "dipole_bv") == 0) current_node->dipole_bv = *value;*/
  else if (strcmp(lpar, "other_bv") == 0) current_node->other_bv = *value;
  else if (strcmp(lpar, "obs_point") == 0) current_node->obs_point = *value;
  else if (strcmp(lpar, "sel_sector") == 0) current_node->sel_sector = *value;
  else if (strcmp(lpar, "enable") == 0) current_node->enable = *value;

  else if (strcmp(lpar, "k0") == 0) {
    store_comm_par_value("k0",*value,el->def);
    el->def->par_names->inform[8] = 1;
  }
  
  else if (strcmp(lpar, "k1") == 0) store_comm_par_value("k1",*value,el->def);
  else if (strcmp(lpar, "k2") == 0) store_comm_par_value("k2",*value,el->def);
  // The inform is to make sure they are written out to a new sequence. 
  else if (strcmp(lpar, "k1tap") == 0) {
    store_comm_par_value("k1tap",*value,el->def);
    el->def->par_names->inform[9] = 1;
  }
  else if (strcmp(lpar, "k1stap") == 0) {
    store_comm_par_value("k1stap",*value,el->def);
    el->def->par_names->inform[10] = 1;
  }
  else if (strcmp(lpar, "k2tap") == 0){
    store_comm_par_value("k2tap",*value,el->def);
    el->def->par_names->inform[9] = 1;
  }
  else if (strcmp(lpar, "k2stap") == 0){
    store_comm_par_value("k2stap",*value,el->def);
    el->def->par_names->inform[10] = 1;
  }
  else if (strcmp(lpar, "lagtap") == 0) {
    store_comm_par_value("lagtap",*value,el->def);
    el->def->par_names->inform[9] = 1;
  }
  else if (strcmp(lpar, "lag") == 0) store_comm_par_value("lag",*value,el->def);

  /* added by E. T. d'Amico 27 feb 2004 */

  else if (strcmp(lpar, "e1") == 0) store_comm_par_value("e1",*value,el->def);
  else if (strcmp(lpar, "e2") == 0) store_comm_par_value("e2",*value,el->def);
  else if (strcmp(lpar, "angle") == 0) store_comm_par_value("angle",*value,el->def);

  /* added by E. T. d'Amico 12 may 2004 */

  else if (strcmp(lpar, "h1") == 0) store_comm_par_value("h1",*value,el->def);
  else if (strcmp(lpar, "h2") == 0) store_comm_par_value("h2",*value,el->def);
  else if (strcmp(lpar, "fint") == 0) store_comm_par_value("fint",*value,el->def);
  else if (strcmp(lpar, "fintx") == 0) store_comm_par_value("fintx",*value,el->def);
  else if (strcmp(lpar, "hgap") == 0) store_comm_par_value("hgap",*value,el->def);
  else if (strcmp(lpar, "pass_flag") == 0) current_node->pass_flag = *value;

  /* AL: added by A. Latina 16 Feb 2012 */

  else if (strcmp(lpar, "rfm_freq") == 0) store_comm_par_value("rfm_freq", *value, el->def);
  else if (strcmp(lpar, "rfm_volt") == 0) store_comm_par_value("rfm_volt", *value, el->def);
  else if (strcmp(lpar, "rfm_lag") == 0) store_comm_par_value("rfm_lag", *value, el->def);
  else if (strcmp(lpar, "rfm_harmon") == 0) store_comm_par_value("rfm_harmon", *value, el->def);

 /* added by FRS 29 August 2012 */

  else if (strcmp(lpar, "rm11") == 0) store_comm_par_value("rm11",*value,el->def);
  else if (strcmp(lpar, "rm12") == 0) store_comm_par_value("rm12",*value,el->def);
  else if (strcmp(lpar, "rm13") == 0) store_comm_par_value("rm13",*value,el->def);
  else if (strcmp(lpar, "rm14") == 0) store_comm_par_value("rm14",*value,el->def);
  else if (strcmp(lpar, "rm15") == 0) store_comm_par_value("rm15",*value,el->def);
  else if (strcmp(lpar, "rm16") == 0) store_comm_par_value("rm16",*value,el->def);
  else if (strcmp(lpar, "rm21") == 0) store_comm_par_value("rm21",*value,el->def);
  else if (strcmp(lpar, "rm22") == 0) store_comm_par_value("rm22",*value,el->def);
  else if (strcmp(lpar, "rm23") == 0) store_comm_par_value("rm23",*value,el->def);
  else if (strcmp(lpar, "rm24") == 0) store_comm_par_value("rm24",*value,el->def);
  else if (strcmp(lpar, "rm25") == 0) store_comm_par_value("rm25",*value,el->def);
  else if (strcmp(lpar, "rm26") == 0) store_comm_par_value("rm26",*value,el->def);
  else if (strcmp(lpar, "rm31") == 0) store_comm_par_value("rm31",*value,el->def);
  else if (strcmp(lpar, "rm32") == 0) store_comm_par_value("rm32",*value,el->def);
  else if (strcmp(lpar, "rm33") == 0) store_comm_par_value("rm33",*value,el->def);
  else if (strcmp(lpar, "rm34") == 0) store_comm_par_value("rm34",*value,el->def);
  else if (strcmp(lpar, "rm35") == 0) store_comm_par_value("rm35",*value,el->def);
  else if (strcmp(lpar, "rm36") == 0) store_comm_par_value("rm36",*value,el->def);
  else if (strcmp(lpar, "rm41") == 0) store_comm_par_value("rm41",*value,el->def);
  else if (strcmp(lpar, "rm42") == 0) store_comm_par_value("rm42",*value,el->def);
  else if (strcmp(lpar, "rm43") == 0) store_comm_par_value("rm43",*value,el->def);
  else if (strcmp(lpar, "rm44") == 0) store_comm_par_value("rm44",*value,el->def);
  else if (strcmp(lpar, "rm45") == 0) store_comm_par_value("rm45",*value,el->def);
  else if (strcmp(lpar, "rm46") == 0) store_comm_par_value("rm46",*value,el->def);
  else if (strcmp(lpar, "rm51") == 0) store_comm_par_value("rm51",*value,el->def);
  else if (strcmp(lpar, "rm52") == 0) store_comm_par_value("rm52",*value,el->def);
  else if (strcmp(lpar, "rm53") == 0) store_comm_par_value("rm53",*value,el->def);
  else if (strcmp(lpar, "rm54") == 0) store_comm_par_value("rm54",*value,el->def);
  else if (strcmp(lpar, "rm55") == 0) store_comm_par_value("rm55",*value,el->def);
  else if (strcmp(lpar, "rm56") == 0) store_comm_par_value("rm56",*value,el->def);
  else if (strcmp(lpar, "rm61") == 0) store_comm_par_value("rm61",*value,el->def);
  else if (strcmp(lpar, "rm62") == 0) store_comm_par_value("rm62",*value,el->def);
  else if (strcmp(lpar, "rm63") == 0) store_comm_par_value("rm63",*value,el->def);
  else if (strcmp(lpar, "rm64") == 0) store_comm_par_value("rm64",*value,el->def);
  else if (strcmp(lpar, "rm65") == 0) store_comm_par_value("rm65",*value,el->def);
  else if (strcmp(lpar, "rm66") == 0) store_comm_par_value("rm66",*value,el->def);
  // This needs to be cleaned up.
  
  

  /* end of additions */
}

void
set_new_position(struct sequence* sequ)
  /* sets a new node position for all nodes */
{
  struct node* c_node = sequ->start;
  double zero_pos = c_node->position;
  int flag = 0;
  while (c_node != NULL)
  {
    if (c_node->from_name == NULL)
    {
      c_node->position -= zero_pos;
      if (c_node->position < zero || (flag && c_node->position == zero))
        c_node->position += sequence_length(sequ);
      if (c_node->position > zero) flag = 1;
      c_node->at_value = c_node->position;
      c_node->at_expr = NULL;
    }
    if (c_node == sequ->end) break;
    c_node = c_node->next;
  }
  c_node->position = c_node->at_value = sequence_length(sequ);
}

void
set_node_bv(struct sequence* sequ)
  /* sets bv flag for all nodes */
{
  struct node* c_node = sequ->ex_start;
  double beam_bv;
  beam_bv = command_par_value("bv", current_beam);
  while (c_node != NULL)
  {
    c_node->other_bv = beam_bv;
    c_node->dipole_bv = beam_bv;
/* dipole_bv kill initiative SF TR FS */
/*      if (c_node->p_elem->bv) c_node->dipole_bv = beam_bv;
        else                    c_node->dipole_bv = 1;*/
    if (c_node == sequ->ex_end) break;
    c_node = c_node->next;
  }
}

void
add_to_node_list( /* adds node to alphabetic node list */
  struct node* node, int inf, struct node_list* nll)
{
  // int j, pos; // not used
  if (name_list_pos(node->name, nll->list) < 0) // (pos = not used
  {
    if (nll->curr == nll->max) grow_node_list(nll);
    add_to_name_list(node->name, inf, nll->list); // j = not used
    nll->nodes[nll->curr++] = node;
  }
}

int
count_nodes(struct sequence* sequ)
{
  int count = 1;
  struct node* c_node = sequ->start;
  while(c_node != sequ->end)
  {
    c_node = c_node->next;
    count++;
  }
  return count;
}

void
current_node_name(char* name, int* l)
/* returns node_name blank padded up to l */
{
  strfcpy(name, current_node->name, *l);
}

void
node_name(char* name, int* l)
  /* returns current node name in Fortran format */
  /* l is max. allowed length in name */
{
  strfcpy(name, current_node->name, *l);
  stoupper(name);
}

int
get_node_count(struct node* node)
  /* finds the count of a node in the current expanded sequence */
{
  int cnt = 0;
  current_node = current_sequ->ex_start;
  while (current_node != NULL)
  {
    if (current_node == node) return cnt;
    cnt++;
    if (current_node == current_sequ->ex_end) break;
    current_node = current_node->next;
  }
  return -1;
}

void
store_node_vector(char* par, int* length, double* vector)
  /* stores vector at node */
{
  char lpar[NAME_L];
  mycpy(lpar, par);
  if (strcmp(lpar, "orbit0") == 0)  copy_double(vector, orbit0, 6);
  else if (strcmp(lpar, "orbit_ref") == 0)
  {
    if (current_node->orbit_ref)
    {
      while (*length > current_node->orbit_ref->max)
        grow_double_array(current_node->orbit_ref);
    }
    else current_node->orbit_ref = new_double_array(*length);
    copy_double(vector, current_node->orbit_ref->a, *length);
    current_node->orbit_ref->curr = *length;
  }
  else if (strcmp(lpar, "surv_data") == 0)
        copy_double(vector, current_node->surv_data, *length);

}

int
store_no_fd_err(double* errors, int* curr)
{
  if (current_node->p_fd_err == NULL) {
    current_node->p_fd_err = new_double_array(FIELD_MAX);
    current_node->p_fd_err->curr = FIELD_MAX;
  }
  else {
    if(current_node->p_fd_err->curr < *curr)
      grow_double_array(current_node->p_fd_err);
  }
  current_node->p_fd_err->curr = *curr;
  copy_double(errors, current_node->p_fd_err->a,current_node->p_fd_err->curr);
  return current_node->p_fd_err->curr;
}

int
advance_node(void)
  /* advances to next node in expanded sequence;
     returns 0 if end of range, else 1 */
{

  if (current_node == current_sequ->range_end)  return 0;
  current_node = current_node->next;

  return 1;
}

int
advance_to_pos(char* table, int* t_pos)
  /* advances current_node to node at t_pos in table */
{
  struct table* t;
  int cnt = 0, ret = 0;
  mycpy(c_dum->c, table);
  if ((t = find_table(c_dum->c)))
  {
    ret = 1;
    if (t->origin == 1)  return 1; /* table is read, has no node pointers */
    while (current_node)
    {
      if (current_node == t->p_nodes[*t_pos-1]) break;
      if ((current_node = current_node->next)
          == current_sequ->ex_start) cnt++;
      if (cnt > 1) return 0;
    }
  }
  return ret;
}

double
line_nodes(struct char_p_array* flat)
  /* creates a linked node list from a flat element list of a line */
{
  int i, j, k;
  double pos = zero, val;
  struct element* el;
  
  for (j = 0; j < flat->curr; j++)
   {
    if ((el = find_element(flat->p[j], element_list)) == NULL)
     {
      
      fatal_error("line contains unknown element:", flat->p[j]);
     }
     
    if (strcmp(el->base_type->name, "rfcavity") == 0 &&
        find_element(el->name, current_sequ->cavities) == NULL)
     {   
      add_to_el_list(&el, 0, current_sequ->cavities, 0);
     }
     
     
    val = el_par_value("l", el);
    pos += val / 2;
    k = 1;
    if ((i = name_list_pos(el->name, occ_list)) < 0)
      i = add_to_name_list(el->name, k, occ_list);
    else k = ++occ_list->inform[i];
    make_elem_node(el, k);
    current_node->at_value = current_node->position = pos;
    pos += val / 2;
  }
  return pos;
}

void
node_string(const char* key, char* string, int* l)
  /* returns current node string value for "key" in Fortran format */
  /* l is max. allowed length in string */
{
  char tmp[2*NAME_L]; mycpy(tmp, key);
  char* p;
  if ((p = command_par_string(tmp, current_node->p_elem->def)))
    strfcpy(string, p, *l);
  else
    memset(string, ' ', *l);
}

int node_apertype(void){
  return current_node->p_elem->aper->apertype;
}
void node_aperture_vector(double *vec){
  for(int i=0;i<4; i++){
  vec[i] = current_node->p_elem->aper->aperture[i];
  }
}
void node_aperture_offset(double *vec){
  for(int i=0;i<2; i++){
  vec[i] = current_node->p_elem->aper->aper_offset[i];
  }
}

int inside_userdefined_geometry(double* x, double *y){
  return aper_chk_inside(*x, *y, current_node->p_elem->aper->xlist, 
    current_node->p_elem->aper->ylist, current_node->p_elem->aper->length );
}

/*
Returns aperture defined as arbitrary polygon
x -> x coordinates
y -> y coordinates
maxlen -> length of the x and y arrays 
returns: total length of the aperture (can be bigger than maxlen)
*/
int get_userdefined_geometry(double* x, double *y, int* maxlen)
{
  double* xi = current_node->p_elem->aper->xlist;
  double* yi = current_node->p_elem->aper->ylist;
  int mx = current_node->p_elem->aper->length;
  
  if (*maxlen < mx) mx = *maxlen;
  
  for(int i=0; i<mx; i++)
   {
     x[i] =  xi[i];
     y[i] =  yi[i];
   }

  return current_node->p_elem->aper->length;

}

int get_userdefined_geometry_len()
{
  return current_node->p_elem->aper->length;
}

int
remove_one(struct node* node)
{
  int pos;
  /* removes a node from a sequence being edited */
  if ((pos = name_list_pos(node->p_elem->name, occ_list)) < 0)  return 0;
  if (node->previous != NULL) node->previous->next = node->next;
  if (node->next != NULL) node->next->previous = node->previous;
  if (occ_list->inform[pos] == 1)
  {
    remove_from_node_list(node, edit_sequ->nodes);
    remove_from_name_list(node->p_elem->name, occ_list);
  }
  else --occ_list->inform[pos];
  /* myfree(rout_name, node); */
  return 1;
}

void
replace_one(struct node* node, struct element* el)
  /* replaces an existing node by a new one made from el */
{
  int i, k = 1;
  remove_from_node_list(node, edit_sequ->nodes);
  if ((i = name_list_pos(el->name, occ_list)) < 0)
    i = add_to_name_list(el->name, 1, occ_list);
  else k = ++occ_list->inform[i];
  strcpy(node->name, compound(el->name, k));
  add_to_node_list(node, 0, edit_sequ->nodes);
  node->p_elem = el;

  node->base_name = el->base_type->name;
  //strcpy(node->base_name, el->base_type->name);
  node->length = el->length;
  if (strcmp(el->base_type->name, "rfcavity") == 0 &&
      find_element(el->name, edit_sequ->cavities) == NULL)
    add_to_el_list(&el, 0, edit_sequ->cavities, 0);
}

// return node index into nl. must be within [start, stop].
struct node*
find_node_by_name(const char* name, struct node_list* nl, struct node* start, struct node* stop)
{
  if (*name == '#') {
    // using `strncmp` to allow '#start' etc (TG: this is broken IMO…)
    // e.g. '#STARTAD' (test-sequence-3)
    if (strncmp(name, "#s", 2) == 0) return start;
    if (strncmp(name, "#e", 2) == 0) return stop;
    return NULL;
  }

  char tmp[2*NAME_L];
  strcpy(tmp, name);
  if (square_to_colon(tmp) == 0)
    return NULL;

  // TG: not checking for boundaries, because:
  // - sometimes a node outside the range is requested (test-track-11).
  // - start/stop may not even be registered into the node list (test-dynap).
  int pos = name_list_pos(tmp, nl->list);
  if (pos >= 0)
      return nl->nodes[pos];

  // This fallback is required for selecting implicit DRIFTs which are
  // currently not registered in the node list:
  struct node* n = start;
  while (n && strcmp(n->name, tmp)) {
    // TG: note that `stop->next=start` is possible, therefore, have to compare
    // against `stop` itself - in the loop body:
    if (n == stop)
      return NULL;
    n = n->next;
  }
  return n;
}
