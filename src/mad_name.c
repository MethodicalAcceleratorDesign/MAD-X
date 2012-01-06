#include "madx.h"

char*
get_new_name(void)
  /* makes a new internal element or variable name */
{
  char name[NAME_L] = "__";
  sprintf(&name[2], "%d", new_name_count++);
  strcat(name, "__");
  return permbuff(name);
}

struct name_list*
clone_name_list(struct name_list* p)
{
  int i, l = p->curr > 0 ? p->curr : 1;
  char name[2*NAME_L];
  struct name_list* clone;
  strcpy(name, p->name); strcat(name, "_clone");
  clone = new_name_list(name, l);
  for (i = 0; i < p->curr; i++) clone->index[i] = p->index[i];
  for (i = 0; i < p->curr; i++) clone->inform[i] = p->inform[i];
  for (i = 0; i < p->curr; i++) clone->names[i] = p->names[i];
  clone->curr = p->curr;
  return clone;
}

struct name_list*
new_name_list(char* list_name, int length)
{
  char rout_name[] = "new_name_list";
  struct name_list* il =
    (struct name_list*) mycalloc(rout_name,1, sizeof(struct name_list));
  strcpy(il->name, list_name);
  il->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", il->name);
  il->names = (char**) mycalloc(rout_name,length, sizeof(char*));
  il->index = (int*) mycalloc(rout_name,length, sizeof(int));
  il->inform = (int*) mycalloc(rout_name,length, sizeof(int));
  il->max = length;
  return il;
}

struct name_list*
delete_name_list(struct name_list* l)
{
  char rout_name[] = "delete_name_list";
  if (l == NULL) return NULL;
  if (stamp_flag && l->stamp != 123456)
    fprintf(stamp_file, "d_n_l double delete --> %s\n", l->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", l->name);
  if (l->index != NULL)  myfree(rout_name, l->index);
  if (l->inform != NULL)  myfree(rout_name, l->inform);
  if (l->names != NULL)   myfree(rout_name, l->names);
  myfree(rout_name, l);
  return NULL;
}

struct vector_list*
new_vector_list(int length)
  /* creates a name list and pointer list
     for double arrays with initial length "length".
  */
{
  char rout_name[] = "new_vector_list";
  struct vector_list* vector = mycalloc(rout_name,1, sizeof(struct vector_list));
  vector->max = length;
  vector->names = new_name_list("vector_list", length);
  vector->vectors = mycalloc(rout_name, length, sizeof(struct double_array*));
  return vector;
}

struct vector_list*
delete_vector_list(struct vector_list* vector)
{
  char rout_name[] = "delete_vector_list";
  int j;
  if (vector == NULL) return NULL;
  if (vector->names != NULL)
  {
    for (j = 0; j < vector->names->curr; j++)
      if (vector->vectors[j]) delete_double_array(vector->vectors[j]);
    delete_name_list(vector->names);
  }
  if (vector->vectors != NULL) myfree(rout_name, vector->vectors);
  myfree(rout_name, vector);
  return NULL;
}

void
dump_name_list(struct name_list* nl)
{
  int i;
  puts(" ");
  for (i = 0; i < nl->curr; i++)
  {
    fprintf(prt_file, v_format("%S %I\n"),
            nl->names[nl->index[i]], nl->inform[nl->index[i]]);
  }
}

void
copy_name_list(struct name_list* out, struct name_list* in)
  /* copies namelist in to namelist out */
{
  int i, l = in->curr > 0 ? in->curr : 1;
  while (out->max < l) grow_name_list(out);
  for (i = 0; i < in->curr; i++) out->index[i]  = in->index[i];
  for (i = 0; i < in->curr; i++) out->inform[i] = in->inform[i];
  for (i = 0; i < in->curr; i++) out->names[i]  = in->names[i];
  out->curr = in->curr;
}

void
grow_name_list(struct name_list* p)
{
  char rout_name[] = "grow_name_list";
  char** n_loc = p->names;
  int* l_ind = p->index;
  int* l_inf = p->inform;
  int j, new = 2*p->max;

  p->max = new;
  p->names = (char**) mycalloc(rout_name,new, sizeof(char*));
  p->index = (int*) mycalloc(rout_name,new, sizeof(int));
  p->inform = (int*) mycalloc(rout_name,new, sizeof(int));
  for (j = 0; j < p->curr; j++)
  {
    p->names[j] = n_loc[j];
    p->index[j] = l_ind[j];
    p->inform[j] = l_inf[j];
  }
  myfree(rout_name, n_loc);
  myfree(rout_name, l_ind);
  myfree(rout_name, l_inf);
}

void
grow_vector_list(struct vector_list* p)
{
  char rout_name[] = "grow_vector_list";
  struct double_array** v_loc = p->vectors;
  int j, new = 2*p->max;

  p->max = new;
  p->vectors
    = (struct double_array**) mycalloc(rout_name,new,
                                       sizeof(struct double_array*));
  for (j = 0; j < p->curr; j++) p->vectors[j] = v_loc[j];
  myfree(rout_name, v_loc);
}

int
add_to_name_list(char* name, int inf, struct name_list* vlist)
  /* adds name to alphabetic name list vlist */
  /* inf is an integer kept with name */
{
  int j, num, low = 0, mid, high = vlist->curr - 1, pos = 0, ret;

  if (name == NULL) return -1;

  ret = name_list_pos(name, vlist);
  if ( ret < 0)
  {
    while (low <= high)
    {
      mid = (low + high) / 2;
      if ((num = strcmp(name, vlist->names[vlist->index[mid]])) < 0)
      {
        high = mid - 1; pos = mid;
      }
      else if (num > 0) {
        low  = mid + 1; pos = low;
      }
    }
    ret = vlist->curr;
    if (vlist->curr == vlist->max) grow_name_list(vlist);
    for (j = vlist->curr; j > pos; j--) vlist->index[j] = vlist->index[j-1];
    vlist->index[pos] = vlist->curr;
    vlist->inform[vlist->curr] = inf;
    vlist->names[vlist->curr++] = name;
  }
  else  vlist->inform[ret] = inf;
  return ret;
}

int
name_list_pos(char* p, struct name_list* vlist)
{
  int num, mid, low = 0, high = vlist->curr - 1;
  while (low <= high)
  {
    mid = (low + high) / 2;
    if ((num=strcmp(p, vlist->names[vlist->index[mid]])) < 0)  high = mid - 1;
    else if ( num > 0) low  = mid + 1;
    else               return vlist->index[mid];
  }
  return -1;
}


