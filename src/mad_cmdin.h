#ifndef MAD_CMDIN_H
#define MAD_CMDIN_H

// types

struct command;
struct name_list;
struct char_p_array;

struct in_cmd          /* contains information about classified command */
{
  char name[NAME_L];
  char* label;         /* pointer to label: if != NULL then buffer this */
  int type;            /*    0 command from list;
                             1 element definiton outside sequence;
                             2 variable definition;
                             3 start or end of sequence;
                             4 element definition; */
  int sub_type;        /* position in cmd_match_base */
  int stamp;
  int decl_start;      /* start of declarative part in tok_list */
  int clone_flag;      /* if zero, clone can be dropped after decoding */
  struct char_p_array* tok_list; /* contains pointers to tokens */
  struct command* cmd_def;       /* points to command definition */
  struct command* clone;         /* points to clone of command definition */
};

struct in_cmd_list /* contains list of in_cmd pointers sorted by label */
{
  int stamp;
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  struct name_list* labels;     /* index list of labels */
  struct in_cmd** in_cmds;      /* in_cmd pointer list */
};

// interface

struct in_cmd* new_in_cmd(int length);            // many uses
struct in_cmd* delete_in_cmd(struct in_cmd*);     // many uses
struct in_cmd* buffered_cmd (struct in_cmd*);     // used by mad_cmd.c
struct in_cmd* buffer_in_cmd(struct in_cmd*);     // used by mad_eval.c

struct in_cmd_list* new_in_cmd_list(int length);  // used by mad_core.c
void                grow_in_cmd_list(struct in_cmd_list*);
// missing delete?

void  scan_in_cmd (struct in_cmd*);
void  dump_in_cmd (struct in_cmd*); // for debugging

#endif // MAD_CMDIN_H
