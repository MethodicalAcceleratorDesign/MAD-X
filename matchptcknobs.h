#ifndef MATCHPTCKNOBS_H
#define MATCHPTCKNOBS_H

#include <stdio.h>

#define NAME_L 24           /* internal name length */
#define E_D_MAX 500  /* IA: defined in madxl.h*/         /* max. length of extra displacement tables */ 
//#define E_DM_MAX 500         /* max. length of extra displacement table per element */
#define FIELD_MAX 42        /* field error array length */ /*defined in madxl.h*/

#include "madx.h"
#include "math.h"

void madx_mpk_run(struct in_cmd* cmd);
void madx_mpk_prepare();
int  madx_mpk_init();
void madx_mpk_addvariable(struct in_cmd* cmd);
void madx_mpk_addconstraint(const char* constr);
void madx_mpk_end();



void madx_mpk_setcreateuniverse(struct in_cmd* cmd);
void madx_mpk_setcreatelayout(struct in_cmd* cmd);
void madx_mpk_setsetswitch(struct in_cmd* cmd);
void madx_mpk_setcalc(struct in_cmd* cmd);
#endif
