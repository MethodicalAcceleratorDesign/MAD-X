#include "madx.h"

int
down_unit(char* file_name)
  /* makes a called file the current input unit */
{
  FILE* new;
  if ((new = fopen(file_name, "r")) == NULL)
  {
    if (interactive) warning("cannot open input file:", file_name);
    else             fatal_error("cannot open input file:", file_name);
    return 0;
  }
  if (in->curr+1 == in->max) grow_in_buff_list(in);
  in->input_files[++in->curr] = new;
  strcpy(filenames[in->curr],file_name);
  currentline[in->curr] = 0;
  return 1;
}

void
grow_in_buff_list(struct in_buff_list* p)
{
  char rout_name[] = "grow_in_buff_list";
  struct in_buffer** e_loc = p->buffers;
  FILE** f_loc = p->input_files;
  int j, new = 2*p->max;
  p->max = new;
  p->buffers
    = (struct in_buffer**) mycalloc(rout_name,new, sizeof(struct in_buffer*));
  for (j = 0; j < p->curr; j++) p->buffers[j] = e_loc[j];
  myfree(rout_name, e_loc);
  p->input_files = mycalloc(rout_name, new, sizeof(FILE*));
  for (j = 0; j < p->curr; j++) p->input_files[j] = f_loc[j];
  myfree(rout_name, f_loc);
}

struct in_buffer*
new_in_buffer(int length)
{
  char rout_name[] = "new_in_buffer";
  struct in_buffer* new =
    (struct in_buffer*) mycalloc(rout_name,1, sizeof(struct in_buffer));
  strcpy(new->name, "in_buffer");
  new->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", new->name);
  new->c_a = new_char_array(length);
  new->flag = -1;
  return new;
}

struct in_buff_list*
new_in_buff_list(int length)
{
  char rout_name[] = "new_inbuf_list";
  struct in_buff_list* bll =
    (struct in_buff_list*) mycalloc(rout_name,1, sizeof(struct in_buff_list));
  strcpy(bll->name, "in_buff_list");
  bll->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", bll->name);
  bll->buffers =
    (struct in_buffer**) mycalloc(rout_name,length, sizeof(struct in_buffer*));
  bll->input_files =
    (FILE**) mycalloc(rout_name,length, sizeof(FILE*));
  bll->max = length;
  return bll;
}


