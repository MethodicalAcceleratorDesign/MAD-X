#ifndef MAD_MACRO_H
#define MAD_MACRO_H

// types

struct char_array;
struct char_p_array;
struct name_list;

struct macro     /* stores one line or macro definition */
{
  char name[NAME_L];
  int n_formal;                 /* no. of formal parameters */
  int dead;                     /* set to 1 to prevent line from expansion */
  struct char_p_array* formal;  /* list of formal parameters */
  struct char_p_array* tokens;  /* token pointers into body if split (line) */
  struct char_array* body;      /* contains all statements */
  int stamp;
};

struct macro_list
{
  int stamp;
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  struct name_list* list;
  struct macro** macros;
};

// interface

struct macro*       clone_macro(struct macro* org);
struct macro*       new_macro(int n_formal, int length, int p_length);
struct macro_list*  new_macro_list(int length);
struct macro*       delete_macro(struct macro* macro);
void  grow_macro_list(struct macro_list* p);
void  dump_macro(struct macro* m);
void  dump_macro_list(struct macro_list* ml);
int   make_macro(char* statement);
void  exec_macro(struct in_cmd* cmd, int pos);
void  add_to_macro_list(struct macro* macro, struct macro_list* nll);
void  disable_line(char* name, struct macro_list* nll);
int   remove_from_name_list(char* name, struct name_list* nl);
void  replace_lines(struct macro* org, int replace, char** reps);

#endif // MAD_MACRO_H

