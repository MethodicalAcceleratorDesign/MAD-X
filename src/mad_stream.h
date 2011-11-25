#ifndef MAD_STREAM_H
#define MAD_STREAM_H

// types

struct char_array;

struct in_buffer
{
  char name[NAME_L];
  int flag;                     /* flag for logical tests */
  struct char_array* c_a;
  int stamp;
};

struct in_buff_list
{
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation = call level */
  FILE** input_files;           /* input file pointers */
  int stamp;
  struct in_buffer** buffers;   /* in_buff pointer list */
};

// interface

int   down_unit(char* file_name);
void  grow_in_buff_list(struct in_buff_list* p);
struct in_buffer*     new_in_buffer(int length);
struct in_buff_list*  new_in_buff_list(int length);

#endif // MAD_STREAM_H

