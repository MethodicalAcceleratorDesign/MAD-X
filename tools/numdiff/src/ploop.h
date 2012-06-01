#ifndef PLOOP_H
#define PLOOP_H

#include <stddef.h>

#define ploop(I,N,...) \
  do { \
    size_t ploop_n_   = N; \
    size_t ploop_rem_ = ploop_n_ % 8; \
    size_t ploop_cnt_ = ploop_n_ - ploop_rem_; \
\
    switch (ploop_rem_) \
      do {  ploop_cnt_ -= 8; \
            { size_t I = ploop_cnt_+7; __VA_ARGS__; } \
    case 7: { size_t I = ploop_cnt_+6; __VA_ARGS__; } \
    case 6: { size_t I = ploop_cnt_+5; __VA_ARGS__; } \
    case 5: { size_t I = ploop_cnt_+4; __VA_ARGS__; } \
    case 4: { size_t I = ploop_cnt_+3; __VA_ARGS__; } \
    case 3: { size_t I = ploop_cnt_+2; __VA_ARGS__; } \
    case 2: { size_t I = ploop_cnt_+1; __VA_ARGS__; } \
    case 1: { size_t I = ploop_cnt_+0; __VA_ARGS__; } \
    case 0: ; \
      } while(ploop_cnt_); \
  } while(0)

#endif
