#ifndef MAD_MAC_H
#define MAD_MAC_H

#include <assert.h>
#include <limits.h>
#include <string.h>

static inline char*
mad_strncpy(char *dst, const char *src, size_t siz)
{
  assert(dst && src && siz < INT_MAX);
  if (siz > 0) {
    dst[0] = '\0';
    strncat(dst, src, siz-1);
  }
  return dst;
}

static inline char*
mad_strncat(char *dst, const char *src, size_t siz)
{
  assert(dst && src && siz < INT_MAX);
  return strncat(dst, src, siz);
}

#undef  strncpy
#define strncpy(dst,src,siz) mad_strncpy(dst,src,siz)

#undef  strncat
#define strncat(dst,src,siz) mad_strncat(dst,src,siz)

#endif // MAD_MAC_H
