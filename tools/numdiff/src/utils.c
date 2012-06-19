#include <stdlib.h>
#include <string.h>
#include "error.h"
#include "utils.h"
#include "args.h"

static const double pow10_tbl[2*99+1] = { 
  1e-99, 1e-98, 1e-97, 1e-96, 1e-95, 1e-94, 1e-93, 1e-92, 1e-91, 1e-90,
  1e-89, 1e-88, 1e-87, 1e-86, 1e-85, 1e-84, 1e-83, 1e-82, 1e-81, 1e-80,
  1e-79, 1e-78, 1e-77, 1e-76, 1e-75, 1e-74, 1e-73, 1e-72, 1e-71, 1e-70,
  1e-69, 1e-68, 1e-67, 1e-66, 1e-65, 1e-64, 1e-63, 1e-62, 1e-61, 1e-60,
  1e-59, 1e-58, 1e-57, 1e-56, 1e-55, 1e-54, 1e-53, 1e-52, 1e-51, 1e-50,
  1e-49, 1e-48, 1e-47, 1e-46, 1e-45, 1e-44, 1e-43, 1e-42, 1e-41, 1e-40,
  1e-39, 1e-38, 1e-37, 1e-36, 1e-35, 1e-34, 1e-33, 1e-32, 1e-31, 1e-30,
  1e-29, 1e-28, 1e-27, 1e-26, 1e-25, 1e-24, 1e-23, 1e-22, 1e-21, 1e-20,
  1e-19, 1e-18, 1e-17, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10,
  1e-09, 1e-08, 1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 1e-02, 1e-01,

  1e+00, 1e+01, 1e+02, 1e+03, 1e+04, 1e+05, 1e+06, 1e+07, 1e+08, 1e+09,
  1e+10, 1e+11, 1e+12, 1e+13, 1e+14, 1e+15, 1e+16, 1e+17, 1e+18, 1e+19,
  1e+20, 1e+21, 1e+22, 1e+23, 1e+24, 1e+25, 1e+26, 1e+27, 1e+28, 1e+29,
  1e+30, 1e+31, 1e+32, 1e+33, 1e+34, 1e+35, 1e+36, 1e+37, 1e+38, 1e+39,
  1e+40, 1e+41, 1e+42, 1e+43, 1e+44, 1e+45, 1e+46, 1e+47, 1e+48, 1e+49,
  1e+50, 1e+51, 1e+52, 1e+53, 1e+54, 1e+55, 1e+56, 1e+57, 1e+58, 1e+59,
  1e+60, 1e+61, 1e+62, 1e+63, 1e+64, 1e+65, 1e+66, 1e+67, 1e+68, 1e+69,
  1e+70, 1e+71, 1e+72, 1e+73, 1e+74, 1e+75, 1e+76, 1e+77, 1e+78, 1e+79,
  1e+80, 1e+81, 1e+82, 1e+83, 1e+84, 1e+85, 1e+86, 1e+87, 1e+88, 1e+89,
  1e+90, 1e+91, 1e+92, 1e+93, 1e+94, 1e+95, 1e+96, 1e+97, 1e+98, 1e+99,
};

const double *const pow10_table99 = &pow10_tbl[99];

FILE*
open_indexedFile(const char* str, int idx, const char *ext, int optext)
{
  char buf[FILENAME_MAX+100];

  strncpy(buf, str, sizeof buf);

  // remove extension, if any
  const char *dot = strrchr(buf, '.');
  int pos = dot ? dot-buf : (int)strlen(buf);
  buf[pos] = 0;

  // add formatted index
  if (idx > 0) pos += sprintf(buf+pos, option.fmt, idx);

  // copy filename into option
  strncpy(option.indexed_filename, buf, sizeof option.indexed_filename);

  // add extension
  strncat(buf+pos, ext, sizeof buf - pos);

  // try to open
  FILE *fp = fopen(buf, "r");
  if (!fp && !optext) ensure(fp, "failed to open %s", buf);

  // extension if optional, try again upon failure
  if (!fp && optext) {
    buf[pos] = 0;
    fp = fopen(buf, "r");
  }

  // resize buffer for faster read
  if (fp && BUFSIZ < 65536 && setvbuf(fp, 0, _IOFBF, 65536)) {
    fclose(fp);
    error("unable to resize the stream buffer size");
  }

  // debug information
  if (fp) {
    if (optext) inform("processing %s", option.indexed_filename);
    debug("file %s open for reading", buf);
  } else
    trace("<-open_indexedFile: unable to open file %s for reading", buf);

  return fp;
}

