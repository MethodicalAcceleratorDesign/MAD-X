#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include "args.h"
#include "utils.h"
#include "error.h"
#include "ndiff.h"
#include "context.h"
#include "constraint.h"

static int
diff_summary(const struct ndiff *dif)
{
  int n, c;
  ndiff_getInfo(dif, &n, 0, &c);

  if (c) {
    warning("% 6d lines have been diffed   in %s", n, option.indexed_filename);
    warning("% 6d diffs have been detected in %s", c, option.indexed_filename);
  } else
    inform ("% 6d lines have been diffed in %s", n, option.indexed_filename);

  return c;
}

int
main(int argc, const char* argv[])
{
  // start timer
  double t0 = clock();

  // parse arguments
  parse_args(argc, argv);

  // test counter
  int total = 0, failed = 0;

  trace("arguments: total=%d, left=%d, right=%d, curr=%s",
        argc, option.argi, argc-option.argi, argv[option.argi]);

  // file list loop
  while (option.argi < argc) {
    const char *lhs_s = 0, *rhs_s = 0, *cfg_s = 0;
    int n = 0;

    // setup filenames [incremental]
    if (!option.list) {
      if (option.argi < argc) lhs_s = argv[option.argi++];
      if (option.argi < argc) rhs_s = argv[option.argi++];
      if (option.argi < argc) cfg_s = argv[option.argi++];
    } else
      if (option.argi < argc) lhs_s = rhs_s = cfg_s = argv[option.argi++];

    trace("arguments: total=%d, left=%d, right=%d, curr=%s",
          argc, option.argi, argc-option.argi, argv[option.argi]);

    // checks
    if (!lhs_s || !rhs_s) {
      if (option.argi == argc-option.utest) invalid();
      else exit(EXIT_SUCCESS);
    }

    // serie loop
    while (option.serie || !n) {
      FILE *lhs_fp=0, *rhs_fp=0, *cfg_fp=0;

      // open files
      lhs_fp = open_indexedFile(lhs_s, n, option.out_e, 1, 0);
      if (!lhs_fp && n) break;
      rhs_fp = open_indexedFile(rhs_s, n, option.ref_e, 0, 1);
      if (cfg_s) cfg_fp = open_indexedFile(cfg_s, n, option.cfg_e, 0, 0);

      if (!lhs_fp || !rhs_fp) invalid();

      // create context of constraints
      struct context *cxt = context_alloc(0);

      // load constraints
      if (cfg_fp) cxt = context_scan(cxt, cfg_fp);

      // show constraints
      if (option.debug) {
        debug("rules list:");
        context_print(cxt, stderr);
      }

      // numdiff loop
      struct ndiff *dif = ndiff_alloc(lhs_fp, rhs_fp, 0);
      ndiff_loop(dif, cxt, option.blank, option.check);

      // print summary
      if (diff_summary(dif) > 0) ++failed;

      // destroy components
      ndiff_free(dif);
      context_free(cxt);

      // close files
      fclose(lhs_fp);
      fclose(rhs_fp);
      if (cfg_fp) fclose(cfg_fp);

      n += 1;

      // not a serie, stop this loop
      if (!option.serie) break;
    }

    total += n;
  }

  double t1 = clock();
  double t = (t1 - t0) / CLOCKS_PER_SEC;

  if (option.test)
    printf(" + %-50s (%.2f s) - %2d/%2d : %s\n", option.test, t, total-failed, total,
#ifdef _WIN32
            failed ? "FAIL" : "PASS");
#else
            failed ? "\033[31mFAIL\033[0m" : "\033[32mPASS\033[0m");
#endif

  return EXIT_SUCCESS;
}

