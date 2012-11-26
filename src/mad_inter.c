#include "madx.h"

static struct backup {
  struct node*       node;
  int                rbend;
  double             length;
  double             e1, e2;
  int                angle_type;
  struct expression* angle_expr;
  double             angle_value;
} backup;

// Creates interpolating nodes for the plotting routine
int
interpolate_node(int *nint)
{
  struct node *first_node, *clone;
  struct element* el;
  struct command_parameter* cp;
  int i, j, number_nodes;
  double bvk, angle, e1, e2, h1, h2, fint, hgap;
  double zero = 0.0, minus_one = -1.0, numint;
  char *elem_name;
  int rbend_flag, bend_flag = 0;

  numint = *nint;
  number_nodes = *nint - 1;

  /* Set up length, angle and e2 of the first slice
     (first node in the original sequence) */
  
  if (backup.node)
    warning("interpolate_node: node interpolation ongoing, undefined behavior will follow", "");

  first_node = current_node;
  backup.node = first_node;

  el = first_node->p_elem;
  elem_name = el->base_type->name;
  rbend_flag = strcmp(elem_name, "rbend") == 0;
   bend_flag = strcmp(elem_name, "sbend") == 0 || rbend_flag;
  backup.rbend = rbend_flag;

//  bv = node_value("dipole_bv");
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

    if (rbend_flag)
    {
      backup.e1 = e1;
      backup.e2 = e2;

      e1 += bvk * angle / 2.0;
      e2 += bvk * angle / 2.0;
      strcpy(elem_name,"sbend");
      el->def->mad8_type = 3;
    }

    angle /= numint;
    // store_node_value("angle",&angle);
    i = name_list_pos("angle", el->def->par_names);
    cp = el->def->par->parameters[i];

    backup.angle_type  = cp->type;
    backup.angle_expr  = cp->expr;
    backup.angle_value = cp->double_value;
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
  backup.length = first_node->length; 
  first_node->length /= numint;

  // set first_node in range_start of the sequence
  current_sequ->range_start = first_node;

  // clone the current node
  clone = clone_node(first_node,0);
  if (bend_flag) {
    clone->p_elem = clone_element(first_node->p_elem);
    clone->p_elem->def = clone_command(first_node->p_elem->def);
  }

  // Reset to first node
  current_node = first_node;

  // advance to next node
  current_node = current_node->next;

  // set last node in the range to the current node
  current_sequ->range_end = current_node;

  // insert nint - 1 nodes in between the two main nodes
  for (j = 1; j <= number_nodes; j++) {
    link_in_front(clone,current_node);
    current_node = current_node->previous;
    current_node->previous->next = current_node;
    store_node_value("angle",&angle);
/*    store_node_value("dipole_bv",&bv); */
    store_node_value("other_bv",&bvk);
    if (bend_flag) {
      if (j == 1) {
        store_node_value("e2",&e2);
        store_node_value("h2",&h2);
        store_node_value("hgap",&hgap);
        if (fintx_plot < zero)
          store_node_value("fintx",&fint);
        else
          store_node_value("fintx",&fintx_plot);
        store_node_value("fint",&zero);
      }
      else {
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
    if (bend_flag) {
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
  int i, j, rbend_flag, bend_flag = 0;
  double angle=0,e1,e2,numint, h1, h2, fint, hgap, bvk;

  // Deletes the interpolating nodes expanded by the routine interp_node
  numint = *nint;

  // reset first and last node in the sequence range
  current_sequ->range_start = current_sequ->ex_start;
  current_sequ->range_end = current_sequ->ex_end;

  // reset current_node at first node
  for (j = 1; j <= *nint ; j++) {
    if (!current_node)
      error("reset_interpolation: missing current node (deleted?) since last interpolation, undefined behavior will follow", "");
    current_node = current_node->previous;
  }

  if (!backup.node || backup.node != current_node)
    warning("reset_interpolation: current node changed since last interpolation, undefined behavior will follow", "");

  // reset length of first node
  current_node->length = backup.length;

  // resets angle and saves e1 if the element is a bending magnet
  rbend_flag = strcmp(current_node->p_elem->base_type->name, "rbend") == 0 || backup.rbend;
  bend_flag = strcmp(current_node->p_elem->base_type->name, "sbend") == 0 || rbend_flag;

  if (bend_flag) {
//    angle = numint*node_value("angle");
    i = name_list_pos("angle", current_node->p_elem->def->par_names);
    cp = current_node->p_elem->def->par->parameters[i];
    cp->expr = backup.angle_expr;
    cp->type = backup.angle_type;
    cp->double_value = backup.angle_value;

    // store_node_value("angle",&angle);
    e1 = node_value("e1");
    h1 = node_value("h1");
    fint = node_value("fint");
    hgap = node_value("hgap");
  }

  // advance to nint-th  node (second node in original sequence)
  for (j = 1; j <= *nint; j++) advance_node();
  second_node = current_node;

  // back to the last interpolated node
  retreat_node();

  // saves e2 if the element is a bending magnet
  if (bend_flag) {
    e2 = node_value("e2");
    h2 = node_value("h2");
  }

  // delete the interpolating nodes
  for (j = 2; j <= *nint; j++) {
    c_node = current_node;

    retreat_node();
    if (bend_flag) {
      c_node->p_elem->def = delete_command(c_node->p_elem->def);
      c_node->p_elem = delete_element(c_node->p_elem);
    }
    delete_node(c_node);
  }

  /* current_node points now to the first node of the original sequence
     sets next pointer of first node to second node of original sequence */
  current_node->next = second_node;

  // sets pointer of second node to first node of original sequence
  current_node->next->previous = current_node;

  // Updates the values of e1 and e2 and stores them in first node
  //  bv = node_value("dipole_bv");
  bvk = node_value("other_bv");

  if (bend_flag) {
    if (rbend_flag) {
      strcpy(current_node->p_elem->base_type->name, "rbend");
      current_node->p_elem->def->mad8_type = 2;
      e1 = backup.e1; // e1 = e1 - bvk * angle / 2.0;
      e2 = backup.e2; // e2 = e2 - bvk * angle / 2.0;
    }

    store_node_value("e1",&e1);
    store_node_value("e2",&e2);
    store_node_value("h1",&h1);
    store_node_value("h2",&h2);
    store_node_value("fint",&fint);
    store_node_value("fintx",&fintx_plot);
    store_node_value("hgap",&hgap);
  }

  // reset backup
  backup.node = 0;

  return 0;
}

