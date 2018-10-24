#ifndef MAD_NAME_H
#define MAD_NAME_H

// types

struct name_list;
struct double_array;

struct name_list /* contains list of index sorted names plus int inform. */
{
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  int* index;                   /* index for alphabetic access */
  int* inform;                  /* array parallel to names with integer */
  int stamp;
  const char** names;           /* element names for sort */
};

struct vector_list              /* contains named vectors */
{
  int curr, max;
  struct name_list* names;
  struct double_array** vectors;
};

// interface

char*   get_new_name(void);
// double  find_value(char* name, int ntok, char** toks);
int     name_list_pos(const char* p, struct name_list* vlist);

struct name_list*  new_name_list(const char* list_name, int length);
struct name_list*  clone_name_list(struct name_list*);
struct name_list*  delete_name_list(struct name_list*);
void               dump_name_list(struct name_list*);
void               copy_name_list(struct name_list* out, struct name_list* in);
void               grow_name_list(struct name_list*);
int                add_to_name_list(const char* name, int inf, struct name_list*);
int                remove_from_name_list(const char* name, struct name_list* nl);

struct vector_list* new_vector_list(int length);
struct vector_list* delete_vector_list(struct vector_list*);
void                grow_vector_list(struct vector_list*);

#if 0 // kept for debugging
#define add_to_name_list(name, inf, list) \
  ((void)(fprintf(stderr, "**** calling add_to_name_list name=%s, %s:%d\n", (name), __FILE__, __LINE__)), \
    add_to_name_list(name, inf, list))
#endif

#if 0 // kept for debugging
#define name_list_pos(var, list) \
  ((void)(fprintf(stderr, "**** calling name_list_pos name=%s, list=%s from %s:%d\n", \
   (var), (list)?(list)->name:"null", __FILE__, __LINE__)), \
    name_list_pos(var, list))
#endif

#endif // MAD_NAME_H

