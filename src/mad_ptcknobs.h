#ifndef MAD_MATCHPTCKNOBS_H
#define MAD_MATCHPTCKNOBS_H

// types

struct in_cmd;

// interface

void madx_mpk_run(struct in_cmd*);
void madx_mpk_prepare(void);
void madx_mpk_addvariable(struct in_cmd*);
void madx_mpk_addconstraint(const char* constr);
void madx_mpk_end(void);

void madx_mpk_setcreateuniverse(struct in_cmd*);
void madx_mpk_setcreatelayout(struct in_cmd*);
void madx_mpk_setsetswitch(struct in_cmd*);
void madx_mpk_setcalc(struct in_cmd*);

#endif // MAD_MATCHPTCKNOBS_H
