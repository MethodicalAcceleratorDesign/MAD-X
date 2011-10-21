#ifndef MAD_UTILS_H
#define MAD_UTILS_H

// functions

int intrac(void);

// inliners

static inline int
mymax(int a, int b) {
  return a > b ? a : b;
}

static inline int
mymin(int a, int b) {
  return a < b ? a : b;
}

#endif // MAD_UTILS_H

