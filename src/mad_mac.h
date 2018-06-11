#ifndef MAD_MAC_H
#define MAD_MAC_H

#include <assert.h>
#include <limits.h>
#include <string.h>

static inline size_t
mad_strlen(const char *file, int line, const char *src)
{
  (void)file, (void)line;
//  printf("DEBUG:strlen:%s:%d:'%s'\n", file, line, src);
  assert(src);
  return strlen(src);
}

static inline char* // Safer copy (close with '\0')
mad_strncpy(const char *file, int line, char *dst, const char *src, size_t siz)
{
  (void)file, (void)line;
  // printf("DEBUG:strncpy:%s:%d:'%s'[%lu] -> %p\n", file, line, src, siz, dst);
  assert(dst && src && siz < INT_MAX);
  *dst = 0;
  if (siz > 0) strncat(dst, src, siz-1);
  return dst;
}

static inline char* // Safer concat (close with '\0')
mad_strncat(const char *file, int line, char *dst, const char *src, size_t siz)
{
  (void)file, (void)line;
//  printf("DEBUG:strncat:%s:%d:'%s'[%lu]\n", file, line, src, siz);
  assert(dst && src && siz < INT_MAX);
  return strncat(dst, src, siz);
}

static inline char* // Fortran copy (close with spaces)
mad_strfcpy(const char *file, int line, char *dst, const char *src, size_t siz)
{
  (void)file, (void)line;
//  printf("DEBUG:strfcpy:%s:%d:'%s'[%lu]\n", file, line, src, siz);
  assert(dst && src && siz < INT_MAX);
  *dst = 0;
  if (siz > 0) {
    strncat(dst, src, siz-1);
    size_t len = strlen(dst);
    memset(dst+len, ' ', siz-len);
  }
  return dst;
}

#undef  strlen
#define strlen(src)          mad_strlen (__FILE__,__LINE__,src)

#undef  strncpy
#define strncpy(dst,src,siz) mad_strncpy(__FILE__,__LINE__,dst,src,siz)

#undef  strncat
#define strncat(dst,src,siz) mad_strncat(__FILE__,__LINE__,dst,src,siz)

#undef  strfcpy
#define strfcpy(dst,src,siz) mad_strfcpy(__FILE__,__LINE__,dst,src,siz)

#endif // MAD_MAC_H
