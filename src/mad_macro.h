#ifndef MAD_MACRO_H
#define MAD_MACRO_H

// types

struct char_array;
struct char_p_array;
struct name_list;
struct in_cmd;

struct macro     /* stores one line or macro definition */
{
  char name[NAME_L];
  int n_formal;                 /* no. of formal parameters */
  int dead;                     /* set to 1 to prevent line from expansion */
  struct char_p_array* formal;  /* list of formal parameters */
  struct char_p_array* tokens;  /* token pointers into body if split (line) */
  struct char_array* body;      /* contains all statements */
  int stamp;
  struct char_array* original;
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

int           make_macro(char* statement);
struct macro* new_macro(int n_formal, int length, int p_length);

struct macro_list* new_macro_list(int length);
void               add_to_macro_list(struct macro*, struct macro_list*);

void  disable_line(char* name, struct macro_list*);
void  replace_lines(struct macro*, int replace, char** reps);
void  save_macros2file(const char *);
void  exec_macro(struct in_cmd* cmd, int pos);
#endif // MAD_MACRO_H

