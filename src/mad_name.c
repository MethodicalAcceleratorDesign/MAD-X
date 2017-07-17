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
new_name_list(const char* list_name, int length)
{
  const char *rout_name = "new_name_list";
  struct name_list* il = mycalloc(rout_name, 1, sizeof *il);
  strcpy(il->name, list_name);
  il->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", il->name);
  il->names  = mycalloc       (rout_name, length, sizeof *il->names);
  il->index  = mycalloc_atomic(rout_name, length, sizeof *il->index);
  il->inform = mycalloc_atomic(rout_name, length, sizeof *il->inform);
  il->max = length;
  return il;
}

struct name_list*
delete_name_list(struct name_list* l)
{
  const char *rout_name = "delete_name_list";
  if (l == NULL) return NULL;
  if (stamp_flag && l->stamp != 123456)
    fprintf(stamp_file, "d_n_l double delete --> %s\n", l->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", l->name);
  if (l->index != NULL)   myfree(rout_name, l->index);
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
  const char *rout_name = "new_vector_list";
  struct vector_list* vector = mycalloc(rout_name, 1, sizeof *vector);
  vector->max = length;
  vector->names = new_name_list("vector_list", length);
  vector->vectors = mycalloc(rout_name, length, sizeof *vector->vectors);
  return vector;
}

struct vector_list*
delete_vector_list(struct vector_list* vector)
{
  const char *rout_name = "delete_vector_list";
  if (vector == NULL) return NULL;
  if (vector->names != NULL) {
    for (int j = 0; j < vector->names->curr; j++)
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
  puts(" ");
  for (int i = 0; i < nl->curr; i++)
    fprintf(prt_file, v_format("%S %I\n"), nl->names[nl->index[i]], nl->inform[nl->index[i]]);
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
  const char *rout_name = "grow_name_list";
  p->max *= 2;
  if (p->max == 0) p->max++;
  p->names  = myrecalloc(rout_name, p->names,  p->curr * sizeof *p->names,  p->max * sizeof *p->names);
  p->index  = myrecalloc(rout_name, p->index,  p->curr * sizeof *p->index,  p->max * sizeof *p->index);
  p->inform = myrecalloc(rout_name, p->inform, p->curr * sizeof *p->inform, p->max * sizeof *p->inform);
}

void
grow_vector_list(struct vector_list* p)
{
  const char *rout_name = "grow_vector_list";
  struct double_array** v_loc = p->vectors;
  int new = 2*p->max;
  if (new == 0) new++;

  p->max = new;
//  p->vectors = myrealloc(rout_name, p->vectors, new * sizeof *p->vectors);
  p->vectors = mycalloc(rout_name, new, sizeof *p->vectors);
  for (int j = 0; j < p->curr; j++) p->vectors[j] = v_loc[j];
  myfree(rout_name, v_loc);
}


#if 0 // kept for debugging
static void
show_name_list(const char *from, struct name_list* nl)
{
  fprintf(stderr, "**** show name list %s from %s with %d elements\n", nl->name, from, nl->curr);
  for (int i = 0; i < nl->curr; i++) {
    fprintf(stderr, "**** [%d]: name[%d]=%-20s, name[%d]=%-20s\n",
       i, i, nl->names[i], nl->index[i], nl->names[nl->index[i]]);
  }
}
#endif

int add_to_name_list(const char* name, int inf, struct name_list* vlist)
  /* adds name to alphabetic name list vlist */
  /* inf is an integer kept with name */
{
  assert(name);
  assert(vlist);

//  if (!strcmp(name, "map_table")) show_name_list(__func__,vlist);

  int ret = name_list_pos(name, vlist);
  if (ret < 0) {
    int low=0, high = vlist->curr - 1, pos = 0;
    while (low <= high) {
      int mid = (low + high) / 2;
      int num = strcmp(name, vlist->names[vlist->index[mid]]);
           if (num < 0) high = mid - 1, pos = mid;
      else if (num > 0) low  = mid + 1, pos = low;
    }
    ret = vlist->curr;
    if (vlist->curr == vlist->max) grow_name_list(vlist);
    for (int j = vlist->curr; j > pos; j--) vlist->index[j] = vlist->index[j-1];
    vlist->index[pos] = vlist->curr;
    vlist->inform[vlist->curr] = inf;
    vlist->names[vlist->curr++] = name; // only string literal accepted!
  }
  else  vlist->inform[ret] = inf;
  return ret;
}

int name_list_pos(const char* name, struct name_list* vlist)
{
  assert(name);

  if (!vlist) return -1; // empty vlist
  
//  if (!strcmp(name, "map_table")) show_name_list(__func__,vlist);

  int low = 0, high = vlist->curr - 1;
  
  while (low <= high) {
    int mid = (low + high) / 2;
    int num = strcmp(name, vlist->names[vlist->index[mid]]);
         if (num < 0) high = mid - 1;
    else if (num > 0) low  = mid + 1;
    else return vlist->index[mid];
  }
  
  return -1;
}

int
remove_from_name_list(const char* name, struct name_list* nl)
{
  
  char buf[5*NAME_L];
  assert(name);
  assert(nl);

  int j, i, k = -1;
  for (i = 0; i < nl->curr; i++)
    if (strcmp(nl->names[nl->index[i]], name) == 0) break;
  
  
  if (i < nl->curr)
  {
   
    k = nl->index[i];
    for (j = i+1; j < nl->curr; j++) nl->index[j-1] = nl->index[j];
    
    for (j = 0; j < nl->curr-1; j++)
      if(nl->index[j] == nl->curr-1) break;
      
    nl->index[j] = k;
    nl->inform[k] = nl->inform[nl->curr-1];
    nl->names[k] = nl->names[--nl->curr];
  }
  else
  {
    warning("remove_from_name_list",(sprintf(buf,"Table %s not found",name),buf));
  }
  
  return k;
}

