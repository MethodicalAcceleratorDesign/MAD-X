#ifndef ARGS_H
#define ARGS_H

/*
 o---------------------------------------------------------------------o
 |
 | Numdiff
 |
 | Copyright (c) 2012+ laurent.deniau@cern.ch
 | Gnu General Public License
 |
 o---------------------------------------------------------------------o
  
   Purpose:
    manage arguments and options
    display the help
    run unit tests (immediate)
 
 o---------------------------------------------------------------------o
*/

#include <stdio.h>
#include <time.h>
#include "types.h"

struct option {
  int check, debug, serie, list, blank, utest, keep, reset, lgopt, trunc;
  const char *suite, *test;
  const char *fmt, *sfmt;
  const char *chr;
  const char *out_e, *ref_e, *cfg_e;
  const char *unzip[3];
  char lhs_file[FILENAME_MAX];
  char rhs_file[FILENAME_MAX];
  char cfg_file[FILENAME_MAX];
  int  lhs_zip, rhs_zip, cfg_zip;
  int  argi;

  const char *accum;
  time_t dat_t0;
  double clk_t0, clk_t1;
};

extern struct option option;

void usage(void);
void invalid_file(const char*);
void invalid_option(const char*);
void parse_args(int argc, const char *argv[]);

static inline bool
is_option(const char *arg)
{
  return arg[0] == '-' && (arg[1] == '-' || !arg[1] || !option.lgopt);
}

#endif
