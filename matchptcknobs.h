#ifndef MATCHPTCKNOBS_H
#define MATCHPTCKNOBS_H

#include <stdio.h>

#define NAME_L 24           /* internal name length */
#define E_D_MAX 100         /* max. length of extra displacement tables */
#include "madx.h"

void madx_mpk_run(struct in_cmd* cmd);
void madx_mpk_prepare();
void madx_mpk_addvariable(struct in_cmd* cmd);
void madx_mpk_addconstraint(const char* constr);
void madx_mpk_end();



void madx_mpk_setcreateuniverse(struct in_cmd* cmd);
void madx_mpk_setcreatelayout(struct in_cmd* cmd);
void madx_mpk_setsetswitch(struct in_cmd* cmd);
void madx_mpk_setcalc(struct in_cmd* cmd);
#endif
