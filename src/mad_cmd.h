#ifndef MAD_CMD_H
#define MAD_CMD_H

// types

struct name_list;
struct command_parameter_list;

struct command                     /* holds one command */
{
  char name[NAME_L];
  char module[NAME_L];                 /* name of module it belongs to */
  char group[NAME_L];                  /* command group it belongs to */
  int stamp;
  int link_type;                       /* 0 none, 1 start, 2 end of group */
  int mad8_type;                       /* 0 none, else mad-8 element code */
  int beam_def;                        /* beam commands: 1 if defined */
  struct name_list* par_names;         /* names + input flag of parameters */
  struct command_parameter_list* par;  /* parameter pointer list */
};

struct command_list /* contains list of command pointers sorted by name */
{
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  struct name_list* list;       /* index list of names */
  int stamp;
  struct command** commands;    /* command pointer list */
};

struct command_list_list /* contains list of command lists */
{
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  struct name_list* list;       /* index list of names */
  struct command_list** command_lists;    /* command_list pointer list */
  int stamp;
};

// interface

struct command*           new_command(char* name, int nl_length, int pl_length, char* module, char* group, int link, int mad_8);
struct command_list*      new_command_list(char* l_name, int length);
struct command_list_list* new_command_list_list(int length);
struct command*           clone_command(struct command* p);
struct command*           delete_command(struct command* cmd);
struct command_list*      delete_command_list(struct command_list* cl);
struct command_list_list* delete_command_list_list(struct command_list_list* ll);
void    exec_command(void);
int     decode_command(void);
void    control(struct in_cmd* cmd);
int     make_line(char* statement);
int     get_stmt(FILE* file, int supp_flag);
int     cmd_match(int cnt, char** toks, int* cmd_pos, int* decl_start);
struct command* find_command(char* name, struct command_list* cl);
struct command_list* find_command_list(char* name, struct command_list_list* sl);
void    remove_from_command_list(char* label, struct command_list* list);
void    store_command_def(char* cmd_string);  /* processes command definition */
void    get_defined_commands(void);
void    add_to_command_list(char* label, struct command* comm, struct command_list* cl, int flag);
void    add_to_command_list_list(char* label, struct command_list* cl, struct command_list_list* sl);
void    dump_command(struct command* cmd);
void    grow_command_list(struct command_list* p);
void    grow_command_list_list(struct command_list_list* p);
void    print_command(struct command* cmd);

#endif // MAD_CMD_H
