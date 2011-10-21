#include <math.h>

#include "madx.h"

struct element*
get_drift(double length)
  /* makes a drift space with the required length */
{
  struct element *p, *bt;
  struct command* clone;
  char key[NAME_L];
  int i;
  for (i = 0; i < drift_list->curr; i++)
  {
    p = drift_list->elem[i];
    if (fabs(p->length - length) < ten_m_12) return p;
  }
  sprintf(key, "drift_%d", drift_list->curr);
  bt = find_element("drift", base_type_list);
  clone = clone_command(bt->def);
  store_comm_par_value("l", length, clone);
  p = make_element(key, "drift", clone, 0);
  add_to_el_list(&p, 1, drift_list, 0);
  return p;
}

int
add_drifts(struct node* c_node, struct node* end)
{
  struct node *d1;
  struct element *drift;
  double pos, dl, el2;
  int cnt = 0;
  pos = c_node->position
    - c_node->length / two;
  while (c_node != NULL)
  {
    cnt++;
    el2 = c_node->length / two;
    dl = c_node->position - el2 - pos;
    if (dl + ten_m_6 < zero)
    {
      sprintf(c_dum->c, " %s and %s, length %e", c_node->previous->name,
              c_node->name, dl);
      fatal_error("negative drift between elements", c_dum->c);
    }
    else if (dl > ten_m_6)
    {
      cnt++;
      drift = get_drift(dl);
      d1 = new_elem_node(drift, 0);
      link_in_front(d1, c_node);
      d1->position = pos + dl / two;
    }
    pos = c_node->position + el2;
    if (c_node == end) break;
    c_node = c_node->next;
  }
  return cnt;
}


