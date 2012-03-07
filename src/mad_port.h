#ifndef MAD_PORT_H
#define MAD_PORT_H

#ifndef _WIN32
#include <stdint.h>
#endif

#ifdef _ICC
// problem with non-standard Intel names in math.h
#define compound(a,b) compound_intel(a,b)
#include <math.h>
#undef  compound
#else
#include <math.h>
#endif

#endif
