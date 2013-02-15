/*
 o---------------------------------------------------------------------o
 |
 | Numdiff
 |
 | Copyright (c) 2012+ Laurent Deniau
 |
 o---------------------------------------------------------------------o
  
   Purpose:
    main file and loop
 
 o---------------------------------------------------------------------o
*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>

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

  if (!ndiff_feof(dif, 1)) {
    c += 1;
    warning("diff of files %s ended prematurely (truncated output?)", option.indexed_filename);
  }

  if (c) {
    if (option.test)
    warning("files %s from %s differ", option.indexed_filename, option.test);
    warning("% 6d lines have been diffed   in files %s", n, option.indexed_filename);
    warning("% 6d diffs have been detected in files %s", c, option.indexed_filename);
  } else {
    if (option.test)
    inform ("files %s from %s do not differ", option.indexed_filename, option.test);
    inform ("% 6d lines have been diffed in files %s", n, option.indexed_filename);
  }

  return c;
}

int
main(int argc, const char* argv[])
{
  // start timers
  option.dat_t0 = time(0);
  option.clk_t0 = clock();

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
      if (option.argi == argc-option.utest)
        invalid_option(argv[option.argi-1]);
      else exit(EXIT_SUCCESS);
    }

    // serie loop
    while (option.serie || !n) {
      FILE *lhs_fp=0, *rhs_fp=0, *cfg_fp=0;

      // open files
      lhs_fp = open_indexedFile(lhs_s, n, option.out_e, 1, 0);
      if (!lhs_fp && n) break;
      rhs_fp = open_indexedFile(rhs_s, n, option.ref_e, !option.list, 1);
      if (cfg_s) cfg_fp = open_indexedFile(cfg_s, n, option.cfg_e, !option.list, 0);

      if (!lhs_fp) {
        if (option.list) {
          warning("output file '%s[.out]' not found, skipping diff", lhs_s);
          if (rhs_fp) fclose(rhs_fp);
          if (cfg_fp) fclose(cfg_fp);
          ++failed;
          break;
        } else
          invalid_file(lhs_s);
      }

      if (!rhs_fp) {
        if (option.list) {
          warning("reference file '%s.ref' not found, skipping diff", rhs_s);
          if (lhs_fp) fclose(lhs_fp);
          if (cfg_fp) fclose(cfg_fp);
          ++failed;
          break;
        } else
          invalid_file(rhs_s);
      }

      // create context of constraints (using default size)
      struct context *cxt = context_alloc(0);

      // add rule #0: "* * abs=DBL_MIN"
      const struct constraint c =
        constraint_init(slice_initAll(), slice_initAll(), eps_init(eps_abs, DBL_MIN), 0);
      context_add(cxt, &c);

      // load constraints
      if (cfg_fp) cxt = context_scan(cxt, cfg_fp);

      // show constraints
      if (option.debug) {
        debug("rules list:");
        context_print(cxt, stderr);
      }

      // numdiff loop
      struct ndiff *dif = ndiff_alloc(lhs_fp, rhs_fp, cxt, 0);
      ndiff_option(dif, &option.keep, &option.blank, &option.check);
      ndiff_loop(dif);

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

  option.clk_t1 = clock();

  if (option.test) {
    double t = (option.clk_t1 - option.clk_t0) / CLOCKS_PER_SEC;
    printf(" + %-50s (%.2f s) - %2d/%2d : %s\n", option.test, t, total-failed, total,
            failed ? CSTR_RED("FAIL") : CSTR_GREEN("PASS"));
  }

  if (option.acc)
    accum_summary(total, failed);

  return EXIT_SUCCESS;
}

