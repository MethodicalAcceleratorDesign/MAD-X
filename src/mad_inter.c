#include "madx.h"

static struct backup {
  struct node *current_node,
              *first_node, *last_node,
              *range_start, *range_end;
  int bend_flag, rbend_flag;
  int nint;
} backup;

// Creates interpolating nodes for the plotting routine
int
interpolate_node(int *nint)
{
  struct node *clone;
  struct element* el;
  int j;
  double bvk, angle, e1, e2, h1, h2, fint, fintx, hgap;
  double zero = 0.0;
  int rbend_flag, bend_flag = 0;

  /* Set up length, angle and e2 of the first slice
     (first node in the original sequence) */

  if (backup.current_node)
    fatal_error("interpolate_node: node interpolation ongoing, undefined behavior will follow", "");

  backup.current_node = current_node;
  backup.range_start = current_sequ->range_start;
  backup.range_end = current_sequ->range_end;
  backup.nint = *nint;


  el = current_node->p_elem;
  rbend_flag = strcmp(el->base_type->name, "rbend") == 0;
  bend_flag  = strcmp(el->base_type->name, "sbend") == 0 || rbend_flag;

  backup.bend_flag = bend_flag;
  backup.rbend_flag = rbend_flag;

  // create template node
  clone = clone_node(current_node,0);
  if (bend_flag) {
    clone->p_elem = el = clone_element(clone->p_elem);
    clone->p_elem->def = clone_command(clone->p_elem->def);
  }
  clone->length /= *nint;
  clone->other_bv = bvk = backup.current_node->other_bv;

  backup.first_node = clone;
  current_node      = clone;

  if (bend_flag) {
    angle = command_par_value("angle", el->def);
    e1    = command_par_value("e1",    el->def);
    e2    = command_par_value("e2",    el->def);
    h1    = command_par_value("h1",    el->def);
    h2    = command_par_value("h2",    el->def);
    fint  = command_par_value("fint",  el->def);
    fintx = command_par_value("fintx", el->def);
    hgap  = command_par_value("hgap",  el->def);

    if (rbend_flag) {
      e1 += bvk * angle / 2.0;
      e2 += bvk * angle / 2.0;
      el->base_type = find_element("sbend", base_type_list);
      el->def->mad8_type = 3;
    }

    angle /= *nint;
    store_node_value("angle", &angle);
    store_node_value("e1",    &zero);
    store_node_value("e2",    &zero);
    store_node_value("h1",    &zero);
    store_node_value("h2",    &zero);
    store_node_value("fint",  &zero);
    store_node_value("fintx", &zero);
    store_node_value("hgap",  &zero);
  }

  for (j = 2; j <= *nint; j++) {
    clone = clone_node(current_node, 0);
    current_node->next = clone;
    clone->previous = current_node;
    current_node = clone;
    current_node->other_bv = bvk;
  }

  backup.last_node = current_node;

  // link it into the sequence... (needed?)
  backup.current_node->previous->next = backup.first_node;
  backup.current_node->next->previous = backup.last_node;
  backup.first_node->previous         = backup.current_node->previous;
  backup.last_node->next              = backup.current_node->next;

  if (bend_flag) {
    // handle entry node
    current_node = backup.first_node;
    if (*nint >= 2) {
      current_node->p_elem      = clone_element(clone->p_elem);
      current_node->p_elem->def = clone_command(clone->p_elem->def);
    }
    store_node_value("e1",   &e1);
    store_node_value("h1",   &h1);
    store_node_value("hgap", &hgap);
    store_node_value("fint", &fint);

    // handle exit node
    current_node = backup.last_node;
    if (*nint >= 3) {
      current_node->p_elem = clone_element(current_node->p_elem);
      current_node->p_elem->def = clone_command(current_node->p_elem->def);
    }
    store_node_value("e2",    &e2);
    store_node_value("h2",    &h2);
    store_node_value("hgap",  &hgap);
    store_node_value("fintx", fintx < zero ? &fint : &fintx);
  }

  current_node              = backup.first_node;
  current_sequ->range_start = backup.first_node;
  current_sequ->range_end   = backup.last_node;

  return 0;
}

// Deletes the interpolating nodes expanded by the routine interp_node
int
reset_interpolation(void)
{
  struct node *curr, *next;

  if (!backup.current_node)
    fatal_error("reset_interpolation: current node changed since last interpolation, undefined behavior will follow", "");

  current_sequ->range_start = backup.range_start;
  current_sequ->range_end   = backup.range_end;
  current_node              = backup.current_node;

  if (backup.bend_flag) {
    if (backup.nint >= 1) {
        curr = backup.first_node;
        delete_command(curr->p_elem->def);
        delete_element(curr->p_elem);
    }
    if (backup.nint >= 2) {
      curr = backup.first_node->next;
      delete_command(curr->p_elem->def);
      delete_element(curr->p_elem);
    }
    if (backup.nint >= 3) {
      curr = backup.last_node;
      delete_command(curr->p_elem->def);
      delete_element(curr->p_elem);
    }
  }

  // delete the interpolating nodes
  backup.last_node->next = NULL;
  curr = backup.first_node;
  while (curr) {
    next = curr->next;
    delete_node(curr);
    curr = next;
  }

  current_node = backup.current_node;
  current_node->next->previous = current_node;
  current_node->previous->next = current_node;

  backup.current_node = 0;
  return 0;
}


// interpolation without creating intermediate nodes:

int
start_interp_node(int* i)
{
  *i = 0;
  return current_node->interp_at && current_node->interp_at->curr > 0;
}

int
fetch_interp_node(int* i, double* dl)
{
  *dl = current_node->interp_at->a[*i] * current_node->length;
  return ++*i <= current_node->interp_at->curr;
}

void
select_interp(struct command* cmd)
{
  struct select_iter* iter = start_iter_select(cmd, NULL, NULL);
  struct sequence* sequ;

  int clear = log_val("clear", cmd) || strcmp(cmd->name, "deselect") == 0;
  int nint = command_par_value("slice", cmd);
  double step = command_par_value("step", cmd);
  struct double_array* at = command_par_array("at", cmd);

  at = clear || !at || at->curr == 0 ? NULL : clone_double_array(at);
  int fixed_at = clear || at;

  struct node* node;
  while (fetch_node_select(iter, &node, &sequ)) {
    if (node->length == 0)
      continue;
    if (step > 0)
      nint = node->length / step;
    // allocate new `at` only if needed:
    if (!fixed_at && nint > 1 && (!at || at->curr != nint)) {
      at = new_double_array(nint);
      at->curr = nint;
      for (int i = 0; i < nint; ++i)
        at->a[i] = (double) (i+1) / nint;
    }
    if (fixed_at || nint > 1) {
      if (node->interp_at) sequ->num_interp -= node->interp_at->curr-1;
      if (at)              sequ->num_interp +=              at->curr-1;
      node->interp_at = at;
    }
  }
}
