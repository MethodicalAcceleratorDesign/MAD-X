#include "madx.h"

void
print_rfc(void)
  /* prints the rf cavities present */
{
  double freq0, harmon, freq;
  int i, n = current_sequ->cavities->curr;
  struct element* el;
  if (n == 0)  return;
  freq0 = command_par_value("freq0", probe_beam);
  printf("\n RF system: \n");
  printf(v_format(" %S %NFs %NFs %NFs %NFs %NFs\n"),
         "Cavity","length[m]","voltage[MV]","lag","freq[MHz]","harmon");
  for (i = 0; i < n; i++)
  {
    el = current_sequ->cavities->elem[i];
    if ((harmon = el_par_value("harmon", el)) > zero)
    {
      freq = freq0 * harmon;
      printf(v_format(" %S %F %F %F %F %F\n"),
             el->name, el->length, el_par_value("volt", el),
             el_par_value("lag", el), freq, harmon);
    }
  }
}

void
adjust_rfc(void)
{
  /* adjusts rfc frequency to given harmon number */
  double freq0, harmon, freq;
  int i;
  struct element* el;
  freq0 = command_par_value("freq0", probe_beam);
  if (current_sequ == 0x0)
   {
     printf("adjust_rfc: mad sequence not initialized, can not calculate frequencies from harmonic\n");
     return;
   }
  for (i = 0; i < current_sequ->cavities->curr; i++)
  {
    el = current_sequ->cavities->elem[i];
    if ((harmon = command_par_value("harmon", el->def)) > zero)
    {
      freq = freq0 * harmon;
      store_comm_par_value("freq", freq, el->def);
    }
  }
}

double
rfc_slope(void)
  /* calculates the accumulated "slope" of all cavities */
{
  double slope = zero, harmon, charge, pc;
  struct node* c_node = current_sequ->range_start;
  struct element* el;
  charge = command_par_value("charge", current_beam);
  pc = command_par_value("pc", current_beam);
  do
  {
    el = c_node->p_elem;
    if (strcmp(el->base_type->name, "rfcavity") == 0 &&
        (harmon = command_par_value("harmon", el->def)) > zero)
    {
      double volt = command_par_value("volt", el->def);
      double lag = command_par_value("lag", el->def);
      slope += ten_m_3 * charge * volt * harmon * cos(twopi * lag) / pc;
    }
    if (c_node == current_sequ->range_end) break;
    c_node = c_node->next;
  }
  while (c_node != NULL);
  return slope;
}


