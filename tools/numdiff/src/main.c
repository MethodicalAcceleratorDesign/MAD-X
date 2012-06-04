/* Compilation:
   gcc -std=c99 -W -Wall -pedantic -O3 *.c -o numdiff
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "context.h"
#include "constraint.h"
#include "ndiff.h"
#include "error.h"
#include "utils.h"
#include "utest.h"

static void
run_utest(void)
{
  struct utest *ut = utest_alloc(0);

  fprintf(stderr, "Running unit tests (incomplete)\n");

  // list of unit tests: TODO
  context_utest(ut);
  ndiff_utest(ut);

  // stat
  utest_stat(ut);
  utest_free(ut);
}

static void
usage(void)
{
  logmsg_config.level = inform_level;

  inform("usage:");
  inform("\tnumdiff [options] lhs_file rhs_file [cfg_file]");
  inform("options:");
  inform("\t-s    -serie        enable serie mode (indexed filenames)");
  inform("\t-f    -format fmt   specify the (printf) format fmt for indexes, default is \"%%d\"");
  inform("\t-q    -quiet        enable quiet mode (no output if no diff)");
  inform("\t-d    -debug        enable debug mode");
  inform("\t-u    -utest        run the test suite");
  inform("\t-h    -help         display this help");

  exit(EXIT_FAILURE);
}

static void
invalid(void)
{
  warning("invalid program options or arguments");
  usage();
}

static void
diff_summary(const struct ndiff *dif)
{
  int n, c;
  ndiff_getInfo(dif, &n, 0, &c);
  inform("% 6d lines have been diffed", n);
  inform("% 6d diffs have been detected", c);
}

int
main(int argc, const char* argv[])
{
  int debug = 0, serie = 0;
  const char *fmt = "%d";
  const char *lhs_s=0, *rhs_s=0, *cfg_s=0;
  logmsg_config.level = inform_level;

  // parse command line arguments
  for (int i = 1; i < argc; i++) {

    // display help [immediate]
    if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-help")) {
      usage();
      continue;
    }
    // run utests [immediate]
    if (!strcmp(argv[i], "-u") || !strcmp(argv[i], "-utest")) {
      run_utest();
      continue;
    }

    // set debug mode [setup]
    if (!strcmp(argv[i], "-d") || !strcmp(argv[i], "-debug")) {
      logmsg_config.level = debug_level;
      debug("debug mode on");
      debug = 1;
      continue;
    }

    // set info mode [setup]
    if (!strcmp(argv[i], "-i") || !strcmp(argv[i], "-info")) {
      debug("info mode on");
      logmsg_config.level = inform_level;
      continue;
    }

    // set quiet mode [setup]
    if (!strcmp(argv[i], "-q") || !strcmp(argv[i], "-quiet")) {
      debug("quiet mode on");
      logmsg_config.level = warning_level;
      continue;
    }

    // set serie mode [setup]
    if (!strcmp(argv[i], "-s") || !strcmp(argv[i], "-serie")) {
      debug("serie mode on");
      continue;
    }

    // set serie format [setup]
    if (!strcmp(argv[i], "-f") || !strcmp(argv[i], "-format")) {
      if (serie) {
        fmt = argv[++i];
        debug("format set to '%s'", fmt);
      } else {
        i += 1;
        inform("serie mode is off, format ignored");
      }
      continue;
    }

    // setup filenames [incremental]
         if (!lhs_s) lhs_s = argv[i];
    else if (!rhs_s) rhs_s = argv[i];
    else if (!cfg_s) cfg_s = argv[i];
  }

  // checks
  if (!lhs_s || !rhs_s) invalid();

  // serie loop
  for (int i = 0; i >= 0; i++) {
    FILE *lhs_fp=0, *rhs_fp=0, *cfg_fp=0;

    // open files
               lhs_fp = open_indexedFile(lhs_s, i, fmt, serie && i);
               rhs_fp = open_indexedFile(rhs_s, i, fmt, 0);
    if (cfg_s) cfg_fp = open_indexedFile(cfg_s, i, fmt, 0);

    // serie number
    if (debug && i)
      fprintf(stderr, "serie #%d\n", i);

    // create context of constraints
    struct context *cxt = context_alloc(0);

    // load constraints
    if (cfg_fp) cxt = context_scan(cxt, cfg_fp);

    // show constraints
    if (debug) {
      fprintf(stderr, "rules list:\n");
      context_print(cxt, stderr);
    }

    // numdiff loop
    struct ndiff *dif = ndiff_alloc(lhs_fp, rhs_fp, 0);
    ndiff_loop(dif, cxt, debug);

    // print summary
    diff_summary(dif);

    // destroy components
    ndiff_free(dif);
    context_free(cxt);

    // close files
    fclose(lhs_fp);
    fclose(rhs_fp);
    if (cfg_fp) fclose(cfg_fp);

    // not a serie, stop
    if (!serie) break;
  }

  return EXIT_SUCCESS;
}

