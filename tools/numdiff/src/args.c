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

#include <stdlib.h>
#include <string.h>

#include "args.h"
#include "utils.h"
#include "utest.h"
#include "error.h"
#include "ndiff.h"
#include "context.h"

#ifndef PUNCTCHRS
#define PUNCTCHRS "._$"
#endif

#ifndef SUITEFMT
#define SUITEFMT "[ %s ]"
#endif

#ifndef SERIEFMT
#define SERIEFMT "%d"
#endif

#ifndef MAXKEEP
#define MAXKEEP 25
#endif

#ifndef OUTFILEEXT
#define OUTFILEEXT ".out"
#endif

#ifndef REFFILEEXT
#define REFFILEEXT ".ref"
#endif

#ifndef CFGFILEEXT
#define CFGFILEEXT ".cfg"
#endif

struct option option = {
  // index of processed option
  .argi = 1,

  // names and series numbering
  .fmt = SERIEFMT,

  // names and series numbering
  .sfmt = SUITEFMT,

  // punctuation part of identifiers
  .chr = PUNCTCHRS,

  // number of diff displayed
  .keep = MAXKEEP,

  // file extensions
  .out_e  = OUTFILEEXT, .ref_e =  REFFILEEXT, .cfg_e = CFGFILEEXT
};

static void
run_utest(void)
{
  struct utest *ut = utest_alloc(0);

  inform("Running unit tests (incomplete)");

  // list of unit tests: TODO: more utests
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
  inform("\tnumdiff [options] --list file1 file2 ...");
  inform("\tnumdiff [options] --list --test '1st' file1 file2 --test '2nd' file3 ...");

  inform("");
  inform("options:");
  inform("\t-a   --accum file    accumulate tests information in file");
  inform("\t-b   --blank         ignore blank spaces (space and tabs)");
  inform("\t-c   --check         enable check mode (algorithms crosscheck)");
  inform("\t-d   --debug         enable debug mode (include check mode)");
  inform("\t-g   --cfgext ext    specify the config file extension, default is \"%s\"", option.cfg_e);
  inform("\t-h   --help          display this help");
  inform("\t-i   --info          enable info mode (default)");
  inform("\t-k   --keep num      specify the number of diffs to display per file, default is %d", option.keep);
  inform("\t-l   --list          enable list mode (list of filenames)");
  inform("\t-n   --serie         enable series mode (indexed filenames)");
  inform("\t-f   --seriefmt fmt  specify the (printf) format fmt for indexes, default is \"%s\"", option.fmt);
  inform("\t-o   --outext ext    specify the output file extension, default is \"%s\"", option.out_e);
  inform("\t-p   --punct chrs    punctuation characters part of identifiers, default is \"%s\"", option.chr);
  inform("\t-q   --quiet         enable quiet mode (no output if no diff)");
  inform("\t-r   --refext ext    specify the reference file extension, default is \"%s\"", option.ref_e);
  inform("\t-s   --suite name    set testsuite name for output message (title)");
  inform("\t     --suitefmt fmt  specify the (printf) format fmt for testsuite, default is \"%s\"", option.sfmt);
  inform("\t-t   --test name     set test name for output message (item)");
  inform("\t     --trace         enable trace mode (very verbose, include debug mode)");
  inform("\t-u   --utest         run the test suite (still incomplete)");
  inform("\t-z   --reset         reset accumulated information");

  inform("");
  inform("rules (%s):", option.cfg_e);
  inform("\t#row-range  column-range  commands");
  inform("\t 1-5        *             skip                            # banner");
  inform("\t *          2-$           any abs=1e-15 rel=1e-12 dig=1.5 # global");
  inform("\t 41         *             goto='penalty function'         # jump");
  inform("\t 109:20/5   2-8/3         abs=1e-8                        # specific");

  inform("ranges:");
  inform("\tnum                 row or column number, num >= 0");
  inform("\trange               start - end  [/ stride]");
  inform("\tslice               start : size [/ stride]");
  inform("\t$, *                last row or column, alias for 0-$");

  inform("commands:");
  inform("\tabs=num             absolute error (num < 1)");
  inform("\tall                 conjunctive constraints (default, qualifier)");
  inform("\tany                 disjunctive constraints (qualifier)");
  inform("\tdig=num             input-defined relative error (num > 1)");
  inform("\tequ                 strict numbers equality (same text)");
  inform("\tgoto='tag'          skip lines until 'tag' is found (action)");
  inform("\tign                 ignore numbers");
  inform("\tlarge               allow num >= 1 in abs and rel (qualifier)");
  inform("\tomit='tag'          ignore digits if preceded by 'tag'");
  inform("\trel=num             relative error (num < 1)");
  inform("\tskip                skip lines (action)");
  inform("\ttrace               trace rules (debug, qualifier)");

  inform("");
  inform("information:\thttp://cern.ch/mad");
  inform("contact    :\tmad@cern.ch");
  inform("version    :\t2013.03.22");

  exit(EXIT_FAILURE);
}

void
parse_args(int argc, const char *argv[])
{
  // parse command line arguments
  for (; option.argi < argc; option.argi++) {

    // not an option
    if (argv[option.argi][0] != '-') return;

// ---- [action]

    // display help [action]
    if (!strcmp(argv[option.argi], "-h") || !strcmp(argv[option.argi], "--help")) {
      usage();
      continue;
    }

    // run utests [action]
    if (!strcmp(argv[option.argi], "-u") || !strcmp(argv[option.argi], "--utest")) {
      run_utest();
      option.utest += 1;
      continue;
    }

    // reset accumulation information [action]
    if (!strcmp(argv[option.argi], "-z") || !strcmp(argv[option.argi], "--reset")) {
      ensure(option.accum, "no accumulation file specified");
      debug("reseting file '%s'", option.accum);
      option.reset = 1;
      accum_summary(0, 0, 0, 0);
      continue;
    }

// ---- [setup]

    // set trace mode [setup]
    if (!strcmp(argv[option.argi], "--trace")) {
      logmsg_config.level = trace_level;
      logmsg_config.locate = 1;
      debug("trace mode on");
      option.debug = 1;
      option.check = 1;
      continue;
    }

    // set debug mode [setup]
    if (!strcmp(argv[option.argi], "-d") || !strcmp(argv[option.argi], "--debug")) {
      logmsg_config.level = debug_level;
      logmsg_config.locate = 1;
      debug("debug mode on");
      option.debug = 1;
      option.check = 1;
      continue;
    }

    // set info mode [setup]
    if (!strcmp(argv[option.argi], "-i") || !strcmp(argv[option.argi], "--info")) {
      debug("info mode on");
      logmsg_config.level = inform_level;
      logmsg_config.locate = 0;
      continue;
    }

    // set quiet mode [setup]
    if (!strcmp(argv[option.argi], "-q") || !strcmp(argv[option.argi], "--quiet")) {
      debug("quiet mode on");
      logmsg_config.level = warning_level;
      logmsg_config.locate = 0;
      continue;
    }

    // set check mode [setup]
    if (!strcmp(argv[option.argi], "-c") || !strcmp(argv[option.argi], "--check")) {
      debug("check mode on");
      option.check = 1;
      continue;
    }

    // set blank mode [setup]
    if (!strcmp(argv[option.argi], "-b") || !strcmp(argv[option.argi], "--blank")) {
      debug("blank spaces ignored");
      option.blank = 1;
      continue;
    }

    // set accumulation filename [setup]
    if (!strcmp(argv[option.argi], "-a") || !strcmp(argv[option.argi], "--accum")) {
      option.accum = argv[++option.argi]; 
      debug("accumulation filename set to '%s'", option.accum);
      continue;
    }

    // set suite name [setup]
    if (!strcmp(argv[option.argi], "-s") || !strcmp(argv[option.argi], "--suite")) {
      option.suite = argv[++option.argi];
      debug("suite name set to '%s'", option.suite);
      continue;
    }

    // set suite format [setup]
    if (!strcmp(argv[option.argi], "--suitefmt")) {
      option.sfmt = argv[++option.argi];
      debug("suite format set to '%s'", option.sfmt);
      continue;
    }

    // set serie mode [setup]
    if (!strcmp(argv[option.argi], "-n") || !strcmp(argv[option.argi], "--serie")) {
      debug("serie mode on");
      option.serie = 1;
      continue;
    }

    // set serie format [setup]
    if (!strcmp(argv[option.argi], "--seriefmt")) {
      option.fmt = argv[++option.argi];
      debug("serie format set to '%s'", option.fmt);
      continue;
    }

    // set list mode [setup]
    if (!strcmp(argv[option.argi], "-l") || !strcmp(argv[option.argi], "--list")) {
      debug("list mode on");
      option.list = 1;
      continue;
    }

    // set test name [setup]
    if (!strcmp(argv[option.argi], "-t") || !strcmp(argv[option.argi], "--test")) {
      option.test = argv[++option.argi];
      debug("test name set to '%s'", option.test);
      continue;
    }

    // set keep number [setup]
    if (!strcmp(argv[option.argi], "-k") || !strcmp(argv[option.argi], "--keep")) {
      option.keep = strtoul(argv[++option.argi],0,0);
      debug("keep set to %d", option.keep);
      continue;
    }

    // set punctuation characters [setup]
    if (!strcmp(argv[option.argi], "-p") || !strcmp(argv[option.argi], "--punct")) {
      option.chr = argv[++option.argi]; 
      debug("punctuation characters set to '%s'", option.chr);
      continue;
    }

    // set output extension [setup]
    if (!strcmp(argv[option.argi], "-o") || !strcmp(argv[option.argi], "--outext")) {
      option.out_e = argv[++option.argi]; 
      debug("output extension set to '%s'", option.out_e);
      continue;
    }

    // set reference extension [setup]
    if (!strcmp(argv[option.argi], "-r") || !strcmp(argv[option.argi], "--refext")) {
      option.ref_e = argv[++option.argi]; 
      debug("reference extension set to '%s'", option.ref_e);
      continue;
    }

    // set config extension [setup]
    if (!strcmp(argv[option.argi], "-g") || !strcmp(argv[option.argi], "--cfgext")) {
      option.cfg_e = argv[++option.argi]; 
      debug("config extension set to '%s'", option.cfg_e);
      continue;
    }

    // unknown option
    invalid_option(argv[option.argi]);
  }
}

