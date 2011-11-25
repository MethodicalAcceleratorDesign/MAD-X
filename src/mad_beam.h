#ifndef MAD_BEAM_H
#define MAD_BEAM_H

// types

struct in_cmd;
struct command;
struct sequence;

// interface

void    exec_beam(struct in_cmd* cmd, int flag);
double  get_beam_value(char* name, char* par);
void    save_beam(struct sequence* sequ, FILE* file);
void    show_beam(char* tok);
void    update_beam(struct command* comm);
void    adjust_beam(void);
int     attach_beam(struct sequence* sequ);
void    expand_line(struct char_p_array* l_buff);

#endif // MAD_BEAM_H

