#ifndef MAD_UTILS_H
#define MAD_UTILS_H

// functions

int intrac(void);

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

