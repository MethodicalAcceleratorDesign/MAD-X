#ifndef MAD_PORT_H
#define MAD_PORT_H

#if !defined(_WIN32) || defined(__MINGW32__)
// problem with C99 compliance on Windows
#include <stdint.h>
#endif

#ifdef __MINGW32__
// problem with unistd compliance on Cygwin
typedef long long off64_t;
// For Windows _mkdir
#include <direct.h>
#endif

#ifdef _WIN32
// problem with C99 compliance on Windows
#if !__STDC__ || __STDC_VERSION__ < 199901L
extern double rint(double);
extern double erf (double);
extern double erfc(double);
#endif
#endif

// problem with non-standard Intel names in math.h
#ifdef _ICC
#define compound(a,b) compound_intel(a,b)
#include <math.h>
#undef  compound
#endif

#endif // MAD_PORT_H
