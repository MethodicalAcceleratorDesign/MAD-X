#include "madx.h"

struct in_buffer*
new_in_buffer(int length)
{
  const char *rout_name = "new_in_buffer";
  struct in_buffer* new = mycalloc(rout_name, 1, sizeof *new);
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
  const char *rout_name = "new_inbuf_list";
  struct in_buff_list* bll = mycalloc(rout_name, 1, sizeof *bll);
  strcpy(bll->name, "in_buff_list");
  bll->stamp = 123456;
  if (watch_flag) fprintf(debug_file, "creating ++> %s\n", bll->name);
  bll->buffers = mycalloc(rout_name, length, sizeof *bll->buffers);
  bll->input_files = mycalloc(rout_name, length, sizeof *bll->input_files);
  bll->max = length;
  return bll;
}

void
grow_in_buff_list(struct in_buff_list* p)
{
  const char *rout_name = "grow_in_buff_list";
  struct in_buffer** e_loc = p->buffers;
  FILE** f_loc = p->input_files;
  int j, new = 2*p->max;
  if (new == 0) new++;
  p->max = new;
  p->buffers = mycalloc(rout_name, new, sizeof *p->buffers);
  for (j = 0; j < p->curr; j++) p->buffers[j] = e_loc[j];
  myfree(rout_name, e_loc);
  p->input_files = mycalloc(rout_name, new, sizeof *p->input_files);
  for (j = 0; j < p->curr; j++) p->input_files[j] = f_loc[j];
  myfree(rout_name, f_loc);
}

int
down_unit(char* file_name)
  /* makes a called file the current input unit */
{
  FILE* new;
  file_name = str2path(file_name);
  
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

