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

struct command* new_command(const char* name, int nl_length, int pl_length, const char* module, const char* group, int link, int mad_8);
struct command* delete_command(struct command*);
struct command* clone_command(struct command*);
struct command* clone_command_flat(struct command*);
struct command* find_command(const char* name, struct command_list*);

struct command_list* new_command_list(const char* l_name, int length);
struct command_list* delete_command_list(struct command_list* cl);
struct command_list* find_command_list(const char* name, struct command_list_list*);
void                 grow_command_list(struct command_list*);
void                 add_to_command_list(const char* label, struct command*, struct command_list*, int flag);

struct command_list_list* new_command_list_list(int length);
//struct command_list_list* delete_command_list_list(struct command_list_list*); // never used...
void                      add_to_command_list_list(char* label, struct command_list*, struct command_list_list*);

void    exec_command(void);
int     decode_command(void);
void    dump_command(struct command*); // for debugging
void    print_command(struct command*);
void    store_command_def(char* cmd_string);  /* processes command definition */
int     make_line(char* statement);
int     get_stmt(FILE* file, int supp_flag);
void    get_defined_commands(char *);
void    remove_from_command_list(char* label, struct command_list*);
void    exec_add_expression(struct in_cmd* cmd);

#endif // MAD_CMD_H
