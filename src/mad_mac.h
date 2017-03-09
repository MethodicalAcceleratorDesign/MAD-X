#ifndef MAD_MAC_H
#define MAD_MAC_H

#include <string.h>

static inline char*
mad_strncpy(char *dst, const char *src, size_t siz)
{
  if (siz > 0) {
    strncpy(dst, src, siz-1);
    dst[siz-1] = '\0';
  }
  return dst;
}

static inline char*
mad_strncat(char *dst, const char *src, size_t siz)
{
  if (siz > 0) {
    size_t len = strlen(dst);
    strncat(dst, src, siz-1);
    dst[len+siz-1] = '\0';
  }
  return dst;
}

#undef  strncpy
#define strncpy(dst,src,siz) mad_strncpy(dst,src,siz)

#undef  strncat
#define strncat(dst,src,siz) mad_strncat(dst,src,siz)

#endif // MAD_MAC_H
