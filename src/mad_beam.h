#ifndef MAD_BEAM_H
#define MAD_BEAM_H

// types

struct in_cmd;
struct command;
struct sequence;

// interface

void    exec_beam(struct in_cmd* cmd, int flag);        // used by mad_cmd.c
void    save_beam(struct sequence* sequ, FILE* file, int noexpr);   // used by mad_seq.c
void    show_beam(char* tok);                           // used by mad_exec.c
void    update_beam(struct command* comm);              // used by mad_option.c
void    adjust_beam(void);                              // many uses
int     attach_beam(struct sequence* sequ);             // many uses
void    adjust_probe(double delta_p);                   // many uses
void    adjust_rfc(void);                               // many uses
void    print_rfc(void);
void    print_probe(void);

#endif // MAD_BEAM_H

