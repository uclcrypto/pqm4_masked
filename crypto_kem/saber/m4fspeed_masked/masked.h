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
#include "SABER_params.h"
#include <stdint.h>

#ifndef NSHARES
#error "NSHARES missing"
#define NSHARES 4
#endif

#define COEF_NBITS 13

#define BSSIZE 32

#define POLY_N SABER_N

typedef uint32_t BsBBit[NSHARES];   // dense
typedef BsBBit BsBCoef[COEF_NBITS]; // dense
typedef uint16_t Coef;
typedef Coef Poly[SABER_N] __attribute__ ((aligned(4)));
typedef Coef ACoef[NSHARES];             // dense
typedef Coef APoly[SABER_N][NSHARES] __attribute__ ((aligned(4)));    // dense
typedef APoly APolyVec[SABER_L];         // dense
typedef Coef StrAPoly[NSHARES][SABER_N] __attribute__ ((aligned(4))); // strided
typedef StrAPoly StrAPolyVec[SABER_L];   // strided

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
