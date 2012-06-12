/* Compilation:
   gcc -std=c99 -W -Wall -pedantic -O3 *.c -o numdiff
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
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

  inform("Running unit tests (incomplete)");

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
  inform("\t-s    -suite name   set testsuite name for output message (title)");
  inform("\t-t    -test name    set test name for output message (item)");
  inform("\t-n    -serie        enable series mode (indexed filenames)");
  inform("\t-f    -format fmt   specify the (printf) format fmt for indexes, default is \"%%d\"");
  inform("\t-b    -blank        toggle ignore/no-ignore blank spaces (space and tabs)");
  inform("\t-q    -quiet        enable quiet mode (no output if no diff)");
  inform("\t-c    -check        enable check mode");
  inform("\t-d    -debug        enable debug mode (include check mode)");
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

static int
diff_summary(const struct ndiff *dif, int ns)
{
  static const char *msg[] = {
    "% 6d lines have been diffed",
    "% 6d lines have been diffed in serie #%d",
    "% 6d diffs have been detected",
    "% 6d diffs have been detected in serie #%d"
  };
  int n, c;
  ndiff_getInfo(dif, &n, 0, &c);

  if (c) {
    warning(msg[0 + (ns > 0)], n, ns);
    warning(msg[2 + (ns > 0)], c, ns);
  } else {
    inform (msg[0 + (ns > 0)], n, ns);
  }

  return c;
}

int
main(int argc, const char* argv[])
{
  int failed = 0, check = 0, debug = 0, serie = 0, blank = 0, utest = 0;
  const char *suite = 0, *test = 0;
  const char *fmt = "%d";
  const char *lhs_s=0, *rhs_s=0, *cfg_s=0;
  logmsg_config.level = inform_level;

  // parse command line arguments
  {
    int i;
    for (i = 1; i < argc; i++) {

      // display help [immediate]
      if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-help")) {
        usage();
        continue;
      }
      // run utests [immediate]
      if (!strcmp(argv[i], "-u") || !strcmp(argv[i], "-utest")) {
        run_utest();
        utest += 1;
        continue;
      }

      // set debug mode [setup]
      if (!strcmp(argv[i], "-trace")) {
        logmsg_config.level = trace_level;
        logmsg_config.locate = 1;
        debug("trace mode on");
        debug = 1;
        continue;
      }

      // set debug mode [setup]
      if (!strcmp(argv[i], "-d") || !strcmp(argv[i], "-debug")) {
        logmsg_config.level = debug_level;
        logmsg_config.locate = 1;
        debug("debug mode on");
        debug = 1;
        check = 1;
        continue;
      }

      // set check mode [setup]
      if (!strcmp(argv[i], "-c") || !strcmp(argv[i], "-check")) {
        debug("check mode on");
        check = 1;
        continue;
      }

      // set info mode [setup]
      if (!strcmp(argv[i], "-i") || !strcmp(argv[i], "-info")) {
        debug("info mode on");
        logmsg_config.level = inform_level;
        logmsg_config.locate = 0;
        continue;
      }

      // set blank mode [setup]
      if (!strcmp(argv[i], "-b") || !strcmp(argv[i], "-blank")) {
        debug("blank space ignored");
        blank = !blank;
        continue;
      }

      // set quiet mode [setup]
      if (!strcmp(argv[i], "-q") || !strcmp(argv[i], "-quiet")) {
        debug("quiet mode on");
        logmsg_config.level = warning_level;
        logmsg_config.locate = 0;
        continue;
      }

      // set serie mode [setup]
      if (!strcmp(argv[i], "-n") || !strcmp(argv[i], "-serie")) {
        debug("serie mode on");
        serie = 1;
        continue;
      }

      // set suite name [setup]
      if (!strcmp(argv[i], "-s") || !strcmp(argv[i], "-suite")) {
        suite = argv[++i];
        debug("suite name set to '%s'", suite);
        continue;
      }

      // set test name [setup]
      if (!strcmp(argv[i], "-t") || !strcmp(argv[i], "-test")) {
        test = argv[++i];
        debug("test name set to '%s'", test);
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
    if (!lhs_s || !rhs_s) {
      if (i == argc-utest) invalid();
      else exit(EXIT_SUCCESS);
    }
  }

  // testsuite
  if (suite)
    fprintf(stdout, " [ %s ]\n", suite);

  // serie loop
  {
    double t0 = clock();
    int n = 0;

    while (1) {
      FILE *lhs_fp=0, *rhs_fp=0, *cfg_fp=0;

      // open files
      lhs_fp = open_indexedFile(lhs_s, n, fmt, !(serie && n));
      if (!lhs_fp) break;
      rhs_fp = open_indexedFile(rhs_s, n, fmt, 1);
      if (cfg_s) cfg_fp = open_indexedFile(cfg_s, n, fmt, 1);

      // serie number
      if (n) inform("serie #%d", n);

      // create context of constraints
      struct context *cxt = context_alloc(0);

      // load constraints
      if (cfg_fp) cxt = context_scan(cxt, cfg_fp);

      // show constraints
      if (debug) {
        inform("rules list:");
        context_print(cxt, stderr);
      }

      // numdiff loop
      struct ndiff *dif = ndiff_alloc(lhs_fp, rhs_fp, 0);
      ndiff_loop(dif, cxt, blank, check);

      // print summary
      if (diff_summary(dif, n) > 0) ++failed;

      // destroy components
      ndiff_free(dif);
      context_free(cxt);

      // close files
      fclose(lhs_fp);
      fclose(rhs_fp);
      if (cfg_fp) fclose(cfg_fp);

      n += 1;

      // not a serie, stop
      if (!serie) break;
    }

    double t1 = clock();
    double t = (t1 - t0) / CLOCKS_PER_SEC;

    if (test)
      fprintf(stdout, " + %-50s (%.2f s) - %2d/%2d : %s\n", test, t, n-failed, n,
  #ifdef _WIN32
              failed ? "FAIL" : "PASS");
  #else
              failed ? "\033[31mFAIL\033[0m" : "\033[32mPASS\033[0m");
  #endif
  }

  return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}

