#ifndef MASKED_H
#define MASKED_H
#include "params.h"
#include <stdint.h>

#define NSHARES 2
#define COEF_NBITS 12

#define BSSIZE 32

typedef uint32_t BsBBit[NSHARES];   // dense
typedef BsBBit BsBCoef[COEF_NBITS]; // dense
typedef int16_t Coef;
typedef Coef ACoef[NSHARES];             // dense
typedef Coef APoly[KYBER_N][NSHARES];    // dense
typedef APoly APolyVec[KYBER_K];         // dense
typedef Coef StrAPoly[NSHARES][KYBER_N]; // strided
typedef StrAPoly StrAPolyVec[KYBER_K];   // strided

#include "hal.h"
#include <stdio.h>
#define BAIL(...)                                                              \
  do {                                                                         \
    char bail_buf[100];                                                        \
    sprintf(bail_buf, __VA_ARGS__);                                            \
    hal_send_str(bail_buf);                                                    \
    hal_send_str("#");                                                         \
  } while (0)

#endif
