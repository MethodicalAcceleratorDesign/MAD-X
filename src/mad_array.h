#ifndef MAD_ARRAY_H
#define MAD_ARRAY_H

// types

struct char_array               /* dynamic array of char */
{
  int stamp;
  int max,                      /* max. array size */
      curr;                     /* current occupation */
  char* c;
};

struct char_p_array             /* dynamic array of char pointers */
{
  char name[NAME_L];
  int  max,                     /* max. array size */
       curr,                    /* current occupation */
       flag;                    /* ancillary flag */
  int stamp;
  char** p;
};

struct int_array                /* dynamic array of int */
{
  int stamp;
  char name[NAME_L];
  int  max,                     /* max. array size */
       curr;                    /* current occupation */
  int* i;
};

struct double_array             /* dynamic array of double */
{
  int stamp;
  int max,                      /* max. array size */
      curr;                     /* current occupation */
  double* a;
};

struct char_array_list
{
  char name[NAME_L];
  int stamp;
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  struct char_array** ca;
};

// interface

struct char_array*      new_char_array(int length);
struct char_p_array*    new_char_p_array(int length);
struct int_array*       new_int_array(int length);
struct double_array*    new_double_array(int length);
struct char_array_list* new_char_array_list(int size);

int    addto_char_p_array(struct char_p_array* ch_p_arr, struct char_array* ch_arr);

struct char_p_array*    clone_char_p_array(struct char_p_array* p);
struct int_array*       clone_int_array(struct int_array* p);
struct double_array*    clone_double_array(struct double_array* p);
struct char_array*      delete_char_array(struct char_array* pa);
struct char_p_array*    delete_char_p_array(struct char_p_array* pa, int flag);
struct int_array*       delete_int_array(struct int_array* i);
struct double_array*    delete_double_array(struct double_array* a);
void    dump_char_array(struct char_array* a);
void    dump_char_p_array(struct char_p_array* p);
void    dump_int_array(struct int_array* ia);
void    grow_char_array(struct char_array* p);
void    grow_char_p_array(struct char_p_array* p);
void    grow_int_array(struct int_array* p);
void    grow_double_array(struct double_array* p);
void    grow_char_array_list(struct char_array_list* p);

void    copy_double(double* source, double* target, int n);
int     char_p_pos(char* name, struct char_p_array* p);
void    ftoi_array(struct double_array* da, struct int_array* ia);
int     int_in_array(int k, int n, int* array);

#endif // MAD_ARRAY_H

