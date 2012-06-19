#ifndef ARGS_H
#define ARGS_H

#include <stdio.h>

struct option {
  int check, debug, serie, list, blank, utest;
  const char *suite, *test;
  const char *fmt;
  const char *out_e, *ref_e, *cfg_e;
  char indexed_filename[FILENAME_MAX+100];
  int argi;
};

extern struct option option;

void usage(void);
void invalid(void);
void parse_args(int argc, const char *argv[]);

#endif
