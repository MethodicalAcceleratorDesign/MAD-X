#ifndef MAD_EXEC_H
#define MAD_EXEC_H

// types

struct in_cmd;

// interface

void exec_option(void);
void exec_help(struct in_cmd*);
void exec_assign(struct in_cmd*);
void exec_call(struct in_cmd*);
void exec_cmd_delete(struct in_cmd*);
void exec_print(struct in_cmd*);
void exec_printf(struct in_cmd*);
void exec_show(struct in_cmd*);
void exec_dump(struct in_cmd*);
void exec_shrink_table(struct in_cmd* cmd);
void exec_fill_table(struct in_cmd*);
void exec_create_table(struct in_cmd*);
void exec_removefile(struct in_cmd*);
void exec_renamefile(struct in_cmd*);
void exec_copyfile(struct in_cmd* cmd);
void exec_setvars_table(struct in_cmd*);
void exec_setvars_lin_table(struct in_cmd*);
void exec_addknob_table(struct in_cmd*);
void exec_store_coguess(struct in_cmd*);
void exec_delete_table(const char* name);

#endif // MAD_EXEC_H
