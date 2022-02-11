/* Copyright 2022 UCLouvain, Belgium and PQM4 contributors
 *
 * This file is part of pqm4_masked.
 *
 * pqm4_masked is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, version 3.
 *
 * pqm4_masked is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * pqm4_masked. If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef MASKED_H
#define MASKED_H
#include "params.h"
#include <stdint.h>

#ifndef NSHARES
#error "NSHARES missing"
#define NSHARES 4
#endif

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
