#ifndef MAD_SDDS_H
#define MAD_SDDS_H

// types

struct table;
struct in_cmd;
struct command;
struct command_list;
struct char_p_array;

// interface

void sel_table(char* tname, struct table* t);
void set_selected_rows_tab(struct table* t, struct command_list* select, struct command_list* deselect);
void pro_sdds(struct in_cmd* cmd);
int  sdds_ior(struct in_cmd* cmd);
int  sdds_iow(struct in_cmd* cmd);
int  sdds_readt(char *filename, char *tfsname);
int  sdds_writet_sel(char *filename, struct table *tfstab);
int  head_split(char* buf, struct char_p_array* list);
void sel_table(char* tname, struct table* t);
void set_selected_rows_tab(struct table* t, struct command_list* select, struct command_list* deselect);
int  pass_select_tab(char* name, struct command* sc);

#endif // MAD_SDDS_H
