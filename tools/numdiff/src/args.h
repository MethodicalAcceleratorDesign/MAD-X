#ifndef ARGS_H
#define ARGS_H

/*
 o---------------------------------------------------------------------o
 |
 | Numdiff
 |
 | Copyright (c) 2012+ Laurent Deniau
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
  int check, debug, serie, list, blank, utest, keep, largerr;
  const char *suite, *test;
  const char *fmt;
  const char *chr;
  const char *out_e, *ref_e, *cfg_e;
  char indexed_filename[FILENAME_MAX+100];
  int argi;

  const char *acc;
  time_t dat_t0;
  double clk_t0, clk_t1;
};

extern struct option option;

void usage(void);
void invalid_file(const char*);
void invalid_option(const char*);
void parse_args(int argc, const char *argv[]);

#endif
