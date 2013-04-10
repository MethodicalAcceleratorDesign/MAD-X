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

struct option {
  int check, debug, serie, list, blank, utest, keep, reset, lgopt;
  const char *suite, *test;
  const char *fmt, *sfmt;
  const char *chr;
  const char *out_e, *ref_e, *cfg_e;
  char   current_filename[FILENAME_MAX];
  char reference_filename[FILENAME_MAX];
  int argi;

  const char *accum;
  time_t dat_t0;
  double clk_t0, clk_t1;
};

extern struct option option;

void usage(void);
void invalid_file(const char*);
void invalid_option(const char*);
void parse_args(int argc, const char *argv[]);

#endif
