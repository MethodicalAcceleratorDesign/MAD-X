#include "madx.h"

void
copy_double(double* source, double* target, int n)
  /* copies n double precision values from source to target */
{
  for (int j = 0; j < n; j++)
    target[j] = source[j];
}

int
char_p_pos(char* name, struct char_p_array* p)
  /* returns the position of name in character pointer array p,
     or -1 if not found */
{
  for (int i = 0; i < p->curr; i++)
    if (strcmp(name, p->p[i]) == 0) return i;

  return -1;
}

struct char_array*
new_char_array(int length)
{
  const char *rout_name = "new_char_array";
  struct char_array* il = mycalloc(rout_name, 1, sizeof *il);
  il->stamp = 123456;
  il->curr = 0;
  il->max = length;
  il->c = mycalloc_atomic(rout_name, length, sizeof *il->c);
  return il;
}

struct char_p_array*
new_char_p_array(int length)
{
  const char *rout_name = "new_char_p_array";
  struct char_p_array* il = mycalloc(rout_name, 1, sizeof *il);
  strcpy(il->name, "char_p_array");
  il->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", il->name);
  il->curr = 0;
  il->max = length;
  il->p = mycalloc(rout_name, length, sizeof *il->p);
  return il;
}

int
addto_char_p_array(struct char_p_array* ch_p_arr, struct char_array* ch_arr)
{

  if (ch_p_arr == 0x0)
   {
     printf("addto_char_p_array: argument ch_p_arr is NULL\n");
     return 1;
   }

  if (ch_arr == 0x0)
   {
     printf("addto_char_p_array: argument ch_arr is NULL\n");
     return 2;
   }

  if( ch_p_arr->curr >= ch_p_arr->max )
   {
     printf("addto_char_p_array: No more space in the array.\n");
     return 3;
   }

  ch_p_arr->p[ch_p_arr->curr++] = tmpbuff(ch_arr->c);

  return 0;

}

struct int_array*
new_int_array(int length)
{
  const char *rout_name = "new_int_array";
  struct int_array* il = mycalloc(rout_name, 1, sizeof *il);
  strcpy(il->name, "int_array");
  il->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", il->name);
  il->curr = 0;
  il->max = length;
  il->i = mycalloc_atomic(rout_name, length, sizeof *il->i);
  return il;
}

struct double_array*
new_double_array(int length)
{
  const char *rout_name = "new_double_array";
  struct double_array* il = mycalloc(rout_name,1, sizeof *il);
  il->stamp = 123456;
  il->curr = 0;
  il->max = length;
  il->a = mycalloc_atomic(rout_name, length, sizeof *il->a);
  return il;
}

struct char_array_list*
new_char_array_list(int size)
{
  const char *rout_name = "new_char_array_list";
  struct char_array_list* tl = mycalloc(rout_name, 1, sizeof *tl);
  strcpy(tl->name, "char_array_list");
  tl->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", tl->name);
  tl->max = size;
  tl->ca  = mycalloc(rout_name, size, sizeof *tl->ca);
  return tl;
}

struct char_p_array*
clone_char_p_array(struct char_p_array* p)
{
  int i;
  struct char_p_array* clone = new_char_p_array(p->max);
  for (i = 0; i < p->curr; i++) clone->p[i] = permbuff(p->p[i]);
  clone->curr = p->curr;
  return clone;
}

struct double_array*
clone_double_array(struct double_array* p)
{
  struct double_array* clone = new_double_array(p->curr);
  clone->curr = p->curr;
  for (int i = 0; i < p->curr; i++) clone->a[i] = p->a[i];
  return clone;
}

struct int_array*
clone_int_array(struct int_array* p)
{
  struct int_array* clone = new_int_array(p->curr);
  clone->curr = p->curr;
  for (int i = 0; i < p->curr; i++) clone->i[i] = p->i[i];
  return clone;
}

struct char_array*
delete_char_array(struct char_array* pa)
{
  const char *rout_name = "delete_char_array";
  if (pa == NULL)  return NULL;
  if (pa->c != NULL)  myfree(rout_name, pa->c);
  myfree(rout_name, pa);
  return NULL;
}

struct char_p_array*
delete_char_p_array(struct char_p_array* pa, int flag)
  /* flag = 0: delete only pointer array, = 1: delete char. buffers, too */
{
  const char *rout_name = "delete_char_p_array";
  if (pa == NULL)  return NULL;
  if (stamp_flag && pa->stamp != 123456)
    fprintf(stamp_file, "d_c_p_a double delete --> %s\n", pa->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", pa->name);
  if (flag)
    for (int i = 0; i < pa->curr; i++)
     {
       myfree(rout_name, pa->p[i]);
     }

  if (pa->p != NULL)  myfree(rout_name, pa->p);
  myfree(rout_name, pa);
  return NULL;
}

struct int_array*
delete_int_array(struct int_array* i)
{
  const char *rout_name = "delete_int_array";
  if (i == NULL)  return NULL;
  if (stamp_flag && i->stamp != 123456)
    fprintf(stamp_file, "d_i_a double delete --> %s\n", i->name);
  if (watch_flag) fprintf(debug_file, "deleting --> %s\n", i->name);
  if (i->i != NULL) myfree(rout_name, i->i);
  myfree(rout_name, i);
  return NULL;
}

struct double_array*
delete_double_array(struct double_array* a)
{
  const char *rout_name = "delete_double_array";
  if (a != NULL)
  {
    if (a->a != NULL) myfree(rout_name, a->a);
    myfree(rout_name, a);
  }
  return NULL;
}

void
dump_char_array(struct char_array* a)
{
  char* c = a->c;
  int n = 0, l_cnt = 60, k;
  while (n < a->curr)
  {
    k = a->curr - n; if (k > l_cnt) k = l_cnt;
    strncpy(c_dum->c, c, k);
    c += k; n += k;
    c_dum->c[k] = '\0';
    fprintf(prt_file, "%s\n", c_dum->c);
  }
}

void
dump_char_p_array(struct char_p_array* p)
{
  for (int i = 0; i < p->curr; i++)
    fprintf(prt_file, "%s\n", p->p[i]);
}

void
dump_int_array(struct int_array* ia)
{
  fprintf(prt_file, "dump integer array, length: %d\n", ia->curr);
  for (int i = 0; i < ia->curr; i++) {
    fprintf(prt_file, v_format("%d "), ia->i[i]);
    if ((i+1)%10 == 0) fprintf(prt_file, "\n");
  }
  if (ia->curr%10 != 0) fprintf(prt_file, "\n");
}

void
grow_char_array(struct char_array* p)
{
  const char *rout_name = "grow_char_array";
  p->max *= 2;
  if (p->max == 0) p->max++;
  p->c = myrecalloc(rout_name, p->c, p->curr * sizeof *p->c, p->max * sizeof *p->c);
}

void
grow_char_p_array(struct char_p_array* p)
{
  const char *rout_name = "grow_char_p_array";
  p->max *= 2;
  if (p->max == 0) p->max++;
  p->p = myrecalloc(rout_name, p->p, p->curr * sizeof *p->p, p->max * sizeof *p->p);
}

void
grow_int_array(struct int_array* p)
{
  const char *rout_name = "grow_int_array";
  p->max *= 2;
  if (p->max == 0) p->max++;
  p->i = myrecalloc(rout_name, p->i, p->curr * sizeof *p->i, p->max * sizeof *p->i);
}

void
grow_double_array(struct double_array* p)
{
  const char *rout_name = "grow_double_array";
  p->max *= 2;
  if (p->max == 0) p->max++;
  p->a = myrecalloc(rout_name, p->a, p->curr * sizeof *p->a, p->max * sizeof *p->a);
}

void
grow_char_array_list(struct char_array_list* p)
{
  const char *rout_name = "grow_char_array_list";
  p->max *= 2;
  if (p->max == 0) p->max++;
  p->ca = myrecalloc(rout_name, p->ca, p->curr * sizeof *p->ca, p->max * sizeof *p->ca);
}

void
ftoi_array(struct double_array* da, struct int_array* ia)
  /* converts and copies double array into integer array */
{
  int l = da->curr;
  while (l >= ia->max) grow_int_array(ia);
  for (int i = 0; i < l; i++)
    ia->i[i] = da->a[i];
  ia->curr = l;
}

int
int_in_array(int k, int n, int* array)
  /* returns 1 if k in first n elements of array, else 0 */
{
  for (int j = 0; j < n; j++)
    if (k == array[j]) return 1;
  return 0;
}
