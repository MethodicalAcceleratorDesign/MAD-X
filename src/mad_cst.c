#include "madx.h"

void
get_defined_constants(void)
{
  /* reads + stores the constants defined in madxdict.h */
  supp_char('\n', constant_def);
  pro_input(constant_def);
  start_var = variable_list->curr;
}


