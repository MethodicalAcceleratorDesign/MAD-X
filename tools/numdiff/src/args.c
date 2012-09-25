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

#include <stdlib.h>
#include <string.h>

#include "args.h"
#include "utest.h"
#include "error.h"
#include "ndiff.h"
#include "context.h"

struct option option = {
  // names and series numbering
  .fmt = "%d",

  // punctuation part of identifiers
  .chr = "._",

  // number of diff displayed
  .keep = 25,

  // file extensions
  .out_e  = ".out", .ref_e =  ".ref", .cfg_e = ".cfg"
};

static void
run_utest(void)
{
  struct utest *ut = utest_alloc(0);

  inform("Running unit tests (incomplete)");

  // list of unit tests: TODO
  context_utest(ut);
  ndiff_utest(ut);

  // stat
  utest_stat(ut);
  utest_free(ut);
}

void
invalid_option(const char *str)
{
  warning("invalid program options or arguments '%s'", str);
  usage();
}

void
invalid_file(const char *str)
{
  warning("invalid filename argument '%s'", str);
  usage();
}

void
usage(void)
{
  logmsg_config.level = inform_level;

  inform("usage:");
  inform("\tnumdiff [options] lhs_file rhs_file [cfg_file]");
  inform("\tnumdiff [options] -list file1 file2 ...");
  inform("options:");
  inform("\t-a    -accum file   accumulate tests information in file");
  inform("\t-b    -blank        toggle ignore/no-ignore blank spaces (space and tabs)");
  inform("\t-c    -check        enable check mode");
  inform("\t-d    -debug        enable debug mode (include check mode)");
  inform("\t-e    -cfgext ext   specify the config file extension, default is .cfg");
  inform("\t-f    -format fmt   specify the (printf) format fmt for indexes, default is \"%%d\"");
  inform("\t-h    -help         display this help");
  inform("\t-k    -keep num     specify the number of diffs to display per file, default is 25");
  inform("\t-l    -list         enable list mode (list of filenames)");
  inform("\t-n    -serie        enable series mode (indexed filenames)");
  inform("\t-o    -outext ext   specify the output file extension, default is .out");
  inform("\t-p    -punct chrs   punctuation characters part of identifiers, default is \"._\"");
  inform("\t-q    -quiet        enable quiet mode (no output if no diff)");
  inform("\t-r    -refext ext   specify the reference file extension, default is .ref");
  inform("\t-s    -suite name   set testsuite name for output message (title)");
  inform("\t-t    -test name    set test name for output message (item)");
  inform("\t-u    -utest        run the test suite");

  exit(EXIT_FAILURE);
}

void
parse_args(int argc, const char *argv[])
{
  logmsg_config.level = inform_level;

  // parse command line arguments
  for (option.argi = 1; option.argi < argc; option.argi++) {

    // display help [immediate]
    if (!strcmp(argv[option.argi], "-h") || !strcmp(argv[option.argi], "-help")) {
      usage();
      continue;
    }

    // run utests [immediate]
    if (!strcmp(argv[option.argi], "-u") || !strcmp(argv[option.argi], "-utest")) {
      run_utest();
      option.utest += 1;
      continue;
    }

    // set debug mode [setup]
    if (!strcmp(argv[option.argi], "-trace")) {
      logmsg_config.level = trace_level;
      logmsg_config.locate = 1;
      debug("trace mode on");
      option.debug = 1;
      continue;
    }

    // set debug mode [setup]
    if (!strcmp(argv[option.argi], "-d") || !strcmp(argv[option.argi], "-debug")) {
      logmsg_config.level = debug_level;
      logmsg_config.locate = 1;
      debug("debug mode on");
      option.debug = 1;
      option.check = 1;
      continue;
    }

    // set check mode [setup]
    if (!strcmp(argv[option.argi], "-c") || !strcmp(argv[option.argi], "-check")) {
      debug("check mode on");
      option.check = 1;
      continue;
    }

    // set info mode [setup]
    if (!strcmp(argv[option.argi], "-i") || !strcmp(argv[option.argi], "-info")) {
      debug("info mode on");
      logmsg_config.level = inform_level;
      logmsg_config.locate = 0;
      continue;
    }

    // set blank mode [setup]
    if (!strcmp(argv[option.argi], "-b") || !strcmp(argv[option.argi], "-blank")) {
      debug("blank space ignored");
      option.blank = !option.blank;
      continue;
    }

    // set quiet mode [setup]
    if (!strcmp(argv[option.argi], "-q") || !strcmp(argv[option.argi], "-quiet")) {
      debug("quiet mode on");
      logmsg_config.level = warning_level;
      logmsg_config.locate = 0;
      continue;
    }

    // set serie mode [setup]
    if (!strcmp(argv[option.argi], "-n") || !strcmp(argv[option.argi], "-serie")) {
      debug("serie mode on");
      option.serie = 1;
      continue;
    }

    // set list mode [setup]
    if (!strcmp(argv[option.argi], "-l") || !strcmp(argv[option.argi], "-list")) {
      debug("list mode on");
      option.list = 1;
      continue;
    }

    // set suite name [setup]
    if (!strcmp(argv[option.argi], "-s") || !strcmp(argv[option.argi], "-suite")) {
      option.suite = argv[++option.argi];
      debug("suite name set to '%s'", option.suite);
      fprintf(stdout, " [ %s ]\n", option.suite);
      continue;
    }

    // set test name [setup]
    if (!strcmp(argv[option.argi], "-t") || !strcmp(argv[option.argi], "-test")) {
      option.test = argv[++option.argi];
      debug("test name set to '%s'", option.test);
      continue;
    }

    // set keep number [setup]
    if (!strcmp(argv[option.argi], "-k") || !strcmp(argv[option.argi], "-keep")) {
      option.keep = strtoul(argv[++option.argi],0,0);
      debug("keep set to %d", option.keep);
      continue;
    }

    // set serie format [setup]
    if (!strcmp(argv[option.argi], "-f") || !strcmp(argv[option.argi], "-format")) {
      if (option.serie) {
        option.fmt = argv[++option.argi];
        debug("format set to '%s'", option.fmt);
      } else {
        option.argi++;
        inform("serie mode is off, format ignored");
      }
      continue;
    }

    // set punctuation characters [setup]
    if (!strcmp(argv[option.argi], "-p") || !strcmp(argv[option.argi], "-punct")) {
      option.chr = argv[++option.argi]; 
      debug("punctuation characters set to '%s'", option.chr);
      continue;
    }

    // set output extension [setup]
    if (!strcmp(argv[option.argi], "-o") || !strcmp(argv[option.argi], "-outext")) {
      option.out_e = argv[++option.argi]; 
      debug("output extension set to '%s'", option.out_e);
      continue;
    }

    // set reference extension [setup]
    if (!strcmp(argv[option.argi], "-r") || !strcmp(argv[option.argi], "-refext")) {
      option.ref_e = argv[++option.argi]; 
      debug("reference extension set to '%s'", option.ref_e);
      continue;
    }

    // set config extension [setup]
    if (!strcmp(argv[option.argi], "-e") || !strcmp(argv[option.argi], "-cfgext")) {
      option.cfg_e = argv[++option.argi]; 
      debug("config extension set to '%s'", option.cfg_e);
      continue;
    }

    // set accumulation filename [setup]
    if (!strcmp(argv[option.argi], "-a") || !strcmp(argv[option.argi], "-accum")) {
      option.acc = argv[++option.argi]; 
      debug("accumulation filename set to '%s'", option.acc);
      continue;
    }

    // unknown option
    if (argv[option.argi][0] == '-')
      invalid_option(argv[option.argi]);

    break;
  }
}

