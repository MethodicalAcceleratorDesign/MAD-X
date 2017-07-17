#ifndef MAD_UTILS_H
#define MAD_UTILS_H

// functions

int intrac(void);

///similar to strstr but rejects matchings inside words
const char * strword(const char * str, const char * word); 

// inliners

static inline int
imax(int a, int b) {
  return a > b ? a : b;
}

static inline int
imin(int a, int b) {
  return a < b ? a : b;
}

#endif // MAD_UTILS_H

