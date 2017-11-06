#include "madx.h"

static struct table* bb6d_ixy;

void
make_bb6d_ixy(int* bb6d_ixy_max_rows)
{
  int k;
  struct table* table = find_table("bb6d_ixy");

  if (table)
  {
    delete_table(table);
    k = remove_from_name_list(table->name, table_register->names);
    table_register->tables[k] = table_register->tables[--table_register->curr];
  }
  /* initialise table */
  bb6d_ixy = make_table("bb6d_ixy", "bb6d_ixy", bb6d_ixy_cols,
                         bb6d_ixy_types, *bb6d_ixy_max_rows);
  add_to_table_list(bb6d_ixy, table_register);
  bb6d_ixy->dynamic = 1;
  reset_count("bb6d_ixy");
}

