#ifndef MAD_SELECT_H
#define MAD_SELECT_H

// types

struct node;
struct node_list;
struct table;
struct sequence;
struct command;
struct command_list;

// interface

int   pass_select(char* name, struct command* sc);
int   pass_select_list(char* name, struct command_list* cl);
void  get_select_t_ranges(struct command_list* select, struct command_list* deselect, struct table* t);
int   get_select_ranges(struct sequence* sequ, struct command_list* select, struct node_list* s_ranges);
int   get_select_ex_ranges(struct sequence* sequ, struct command_list* select, struct node_list* s_ranges);
int   get_ex_range(char* range, struct sequence* sequ, struct node** nodes);
int   get_sub_range(char* range, struct sequence* sequ, struct node** nodes);
int   get_range(char* range, struct sequence* sequ, struct node** nodes);
void  set_selected_errors(void);
void  set_range(char* range, struct sequence* sequ);
void  set_sector(void);
void  store_deselect(struct in_cmd* cmd);
void  store_select(struct in_cmd* cmd);

#endif // MAD_SELECT_H

