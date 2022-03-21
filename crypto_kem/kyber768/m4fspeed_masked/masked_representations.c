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
#include "masked_representations.h"
#include "bench.h"
#include "masked.h"
#include <stdint.h>

#include <stdint.h>
#include <stdio.h>
#include <string.h>

/*************************************************
 * Name: transpose32
 *
 * Description: Transpose 32 x 32 bit matrix
 *
 * Arguments:
 * - uint32_t[32] a: the matrix, updated such that the new value of (a[i] >> j)
 *   & 0x1 is the initial value of (a[j] >> i) & 0x1.
 *
 * Adapted from "Hacker's Delight" by Henry S. Warren, Jr.
 * at http://www.icodeguru.com/Embedded/Hacker's-Delight/
 * **************************************************/
#if 1
void transpose32(uint32_t a[32]) {
  int j, k;
  uint32_t m, t;
  m = 0x0000FFFF;
  for (j = 16; j != 0; j = j >> 1, m = m ^ (m << j)) {
    for (k = 0; k < 32; k = (k + j + 1) & ~j) {
      t = (a[k + j] ^ (a[k] >> j)) & m;
      a[k + j] = a[k + j] ^ t;
      a[k] = a[k] ^ (t << j);
    }
  }
}
#else
#include <assert.h>
#define REP5(x) x x x x x
#define FOR5(init, cond, next, body)                                           \
  do {                                                                         \
    init;                                                                      \
    REP5({ {body} next; }) assert(!(cond));                                    \
  } while (0);
void transpose32_unrolled(uint32_t a[32]) {
  int j, k;
  uint32_t m, t;
  m = 0x0000FFFF;
  FOR5(j = 16, j != 0, (j = j >> 1, m = m ^ (m << j)), {
    for (k = 0; k < 32; k = (k + j + 1) & ~j) {
      t = (a[k + j] ^ (a[k] >> j)) & m;
      a[k + j] = a[k + j] ^ t;
      a[k] = a[k] ^ (t << j);
    }
  })
}
#endif

/*************************************************
 * Name:        StrAPoly2APoly
 *
 * Description: Maps strided polynomial into a dense representation
 *
 * Arguments:
 *           - APoly out : dense polynomial
 *           - StrAPoly in: strided polynomial
 * **************************************************/
void StrAPoly2APoly(APoly out, const StrAPoly in) {
  int i, d;
  for (i = 0; i < KYBER_N; i++) {
    for (d = 0; d < NSHARES; d++) {
      out[i][d] = in[d][i];
    }
  }
}

/*************************************************
 * Name:        APoly2StrAPoly
 *
 * Description: Maps dense polynomial into a strided representation
 *
 * Arguments:
 *           - StrAPoly out: strided polynomial
 *           - APoly in : dense polynomial
 * **************************************************/
void APoly2StrAPoly(StrAPoly out, const APoly in) {

  int i, d;
  for (i = 0; i < KYBER_N; i++) {
    for (d = 0; d < NSHARES; d++) {
      out[d][i] = in[i][d];
    }
  }
}

/*************************************************
 * Name:        masked_dense2bitslice
 *
 * Description: maps a dense reprensentation to a bitlisce one
 *
 * Arguments:
 *           - uint32_t *bitslice[]: output bitslice representation. Table of
 * coeffs_size x nshares.
 *           - int16_t *dense[]: input dense reprensetation. Table of n_coeffs x
 * nshares.
 *           - size_t coeffs_size: number of bits to represent the dense
 * coefficients
 *           - size_t n_coeffs: number of coefficients
 *           - size_t nshares: number of shares
 * **************************************************/
void masked_dense2bitslice(size_t nshares, size_t n_coeffs, size_t coeffs_size,
                           uint32_t *bitslice, size_t bitslice_msk_stride,
                           size_t bitslice_data_stride, const int16_t *dense,
                           size_t dense_msk_stride, size_t dense_data_stride) {

  start_bench(my_dense2bs);
  size_t d, c, b;
  for (b = 0; b < coeffs_size; b++) {
    for (d = 0; d < nshares; d++) {
      bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] = 0;
    }
  }

  if ((n_coeffs & 0x3) == 0) {
    for (d = 0; d < nshares; d++) {
      for (c = 0; c < n_coeffs; c += 4) {
        int16_t xd0 = dense[c * dense_data_stride + d * dense_msk_stride];
        int16_t xd1 = dense[(c + 1) * dense_data_stride + d * dense_msk_stride];
        int16_t xd2 = dense[(c + 2) * dense_data_stride + d * dense_msk_stride];
        int16_t xd3 = dense[(c + 3) * dense_data_stride + d * dense_msk_stride];

        for (b = 0; b < coeffs_size; b++) {

          bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] |=
              (xd0 & 0x1) << (c + 0);
          bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] |=
              (xd1 & 0x1) << (c + 1);
          bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] |=
              (xd2 & 0x1) << (c + 2);
          bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] |=
              (xd3 & 0x1) << (c + 3);
          xd0 = xd0 >> 1;
          xd1 = xd1 >> 1;
          xd2 = xd2 >> 1;
          xd3 = xd3 >> 1;
        }
      }
    }

  } else {
    for (d = 0; d < nshares; d++) {
      for (c = 0; c < n_coeffs; c++) {
        int16_t xd = dense[c * dense_data_stride + d * dense_msk_stride];
        for (b = 0; b < coeffs_size; b++) {
          bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] |=
              (xd & 0x1) << c;
          xd = xd >> 1;
        }
      }
    }
  }
  stop_bench(my_dense2bs);
}
void masked_dense2bitslice_u32(size_t nshares, size_t n_coeffs,
                               size_t coeffs_size, uint32_t *bitslice,
                               size_t bitslice_msk_stride,
                               size_t bitslice_data_stride,
                               const uint32_t *dense, size_t dense_msk_stride,
                               size_t dense_data_stride) {
  start_bench(my_dense2bs);
  size_t d, c, b;
  for (b = 0; b < coeffs_size; b++) {
    for (d = 0; d < nshares; d++) {
      bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] = 0;
    }
  }
  if ((n_coeffs & 0x3) == 0) {
    for (d = 0; d < nshares; d++) {
      for (c = 0; c < n_coeffs; c += 4) {
        uint32_t xd0 = dense[c * dense_data_stride + d * dense_msk_stride];
        uint32_t xd1 =
            dense[(c + 1) * dense_data_stride + d * dense_msk_stride];
        uint32_t xd2 =
            dense[(c + 2) * dense_data_stride + d * dense_msk_stride];
        uint32_t xd3 =
            dense[(c + 3) * dense_data_stride + d * dense_msk_stride];

        for (b = 0; b < coeffs_size; b++) {

          bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] |=
              (xd0 & 0x1) << (c + 0);
          bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] |=
              (xd1 & 0x1) << (c + 1);
          bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] |=
              (xd2 & 0x1) << (c + 2);
          bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] |=
              (xd3 & 0x1) << (c + 3);

          xd0 = xd0 >> 1;
          xd1 = xd1 >> 1;
          xd2 = xd2 >> 1;
          xd3 = xd3 >> 1;
        }
      }
    }

  } else {
    for (d = 0; d < nshares; d++) {
      for (c = 0; c < n_coeffs; c++) {
        uint32_t xd = dense[c * dense_data_stride + d * dense_msk_stride];
        for (b = 0; b < coeffs_size; b++) {
          bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] |=
              (xd & 0x1) << c;
          xd = xd >> 1;
        }
      }
    }
  }
  stop_bench(my_dense2bs);
}

/*************************************************
 * Name:        masked_bitslice2dense
 *
 * Description: maps a bitslice reprensentation to a dense one
 *
 * Arguments:
 *           - uint32_t *bitslice[]: input bitslice representation. Table of
 * coeffs_size x nshares.
 *           - int16_t *dense[]: output dense reprensetation. Table of n_coeffs
 * x nshares.
 *           - size_t coeffs_size: number of bits to represent the dense
 * coefficients
 *           - size_t n_coeffs: number of coefficients
 *           - size_t nshares: number of shares
 * **************************************************/
void masked_bitslice2dense(size_t nshares, size_t n_coeffs, size_t coeffs_size,
                           int16_t *dense, size_t dense_msk_stride,
                           size_t dense_data_stride, const uint32_t *bitslice,
                           size_t bitslice_msk_stride,
                           size_t bitslice_data_stride) {

  start_bench(my_bs2dense);
  size_t d, c, b;
  for (c = 0; c < n_coeffs; c += 2) {
    for (d = 0; d < nshares; d++) {
      int16_t xd0 = 0;
      int16_t xd1 = 0;
      for (b = 0; b < coeffs_size; b++) {
        uint32_t y =
            bitslice[b * bitslice_data_stride + d * bitslice_msk_stride];
        xd0 |= ((y >> c) & 0x1) << b;
        xd1 |= ((y >> (c + 1)) & 0x1) << b;
      }
      dense[c * dense_data_stride + d * dense_msk_stride] = xd0;
      dense[(c + 1) * dense_data_stride + d * dense_msk_stride] = xd1;
    }
  }
  stop_bench(my_bs2dense);
}

/*************************************************
 * Name: masked_dense2bitslice_opt
 *
 * Description: maps a dense reprensentation to a bitlisce one
 * REMARK:
 *   This is not equivalent to calling masked_dense2bitslice for two blocks!
 *   The slices of the two blocks are interleaved in the result. This should
 *   however have no impact is the corresponding bs2dense is used.
 *
 *   masked_dense2bitslice_opt(nshares, c, bitslice, bms, bds, dense, dms)
 *   is equivalent to
 *   masked_dense2bitslice(
 *       nshares, BSSIZE, c, bitslice, bms, bds, (int16_t *) dense, dms, 2
 *   );
 *   masked_dense2bitslice(
 *       nshares, BSSIZE, c, bitslice+c, bms, bds, ((int16_t *) dense)+1, dms, 2
 *   );
 *
 *
 * Arguments:
 * - uint32_t *bitslice: output bitslice representation. Array of length
 *   2 x coeffs_size x nshares.
 * - const uint32_t *dense: input dense reprensetation. Array of length BSSIZE
 *   x nshares.  Each uint32_t contains 2 coefficients (little-endian encoded).
 * - size_t coeffs_size: number of bits to represent the dense coefficients
 *   (<=16)
 * - size_t nshares: number of shares
 *
 *  number of coefficients: 2*BSSIZE
 *  dense data stride: 1
 * **************************************************/
void masked_dense2bitslice_opt(size_t nshares, size_t coeffs_size,
                               uint32_t *bitslice, size_t bitslice_msk_stride,
                               size_t bitslice_data_stride,
                               const int16_t *dense, size_t dense_msk_stride,
                               size_t dense_data_stride) {
  start_bench(my_dense2bs);
  uint32_t a[32];
  for (size_t d = 0; d < nshares; d++) {
    for (size_t i = 0; i < 32; i++) {
      a[i] = (dense[i * dense_data_stride + d * dense_msk_stride] << 0) |
             (dense[(i + 32) * dense_data_stride + d * dense_msk_stride] << 16);
    }
    transpose32(a);
    for (size_t i = 0; i < coeffs_size; i++) {
      bitslice[d * bitslice_msk_stride + i * bitslice_data_stride] = a[i];
      bitslice[d * bitslice_msk_stride +
               (i + coeffs_size) * bitslice_data_stride] = a[i + 16];
    }
  }
  stop_bench(my_dense2bs);
}
/*************************************************
 * Name: masked_dense2bitslice_opt_u32
 *
 * Description: maps a dense reprensentation to a bitlisce one
 *
 * Arguments:
 * - uint32_t *bitslice: output bitslice representation. Array of length
 *   coeffs_size x nshares.
 * - const uint32_t *dense: input dense reprensetation. Array of length BSSIZE
 *   x nshares.  Each uint32_t contains 1 coefficient.
 * - size_t coeffs_size: number of bits to represent the dense coefficients
 *   (<=32)
 * - size_t nshares: number of shares
 *
 *  number of coefficients: BSSIZE
 *  dense data stride: 1
 * **************************************************/
void masked_dense2bitslice_opt_u32(
    size_t nshares, size_t coeffs_size, uint32_t *bitslice,
    size_t bitslice_msk_stride, size_t bitslice_data_stride,
    const uint32_t *dense, size_t dense_msk_stride, size_t dense_data_stride) {
  start_bench(my_dense2bs_u32);
  uint32_t a[32];
  for (size_t d = 0; d < nshares; d++) {
    for (size_t i = 0; i < 32; i++) {
      a[i] = dense[i * dense_data_stride + d * dense_msk_stride];
    }
    transpose32(a);
    for (size_t i = 0; i < coeffs_size; i++) {
      bitslice[d * bitslice_msk_stride + i * bitslice_data_stride] = a[i];
    }
  }
  stop_bench(my_dense2bs_u32);
}

/*************************************************
 * Name:        masked_bitslice2dense_opt
 *
 * Description: maps a bitslice reprensentation to a dense one
 *
 * Arguments:
 * - const uint32_t *bitslice: input bitslice representation. Array of length
 *   2 x coeffs_size x nshares.
 * - uint32_t *dense: output dense reprensetation. Array of length BSSIZE x
 *   nshares.  Each uint32_t contains 2 coefficients (little-endian ordered).
 * - size_t coeffs_size: number of bits to represent the dense coefficients
 *   (<=16)
 *
 *  number of coefficients: BSSIZE
 *  dense data stride: 1
 * **************************************************/
void masked_bitslice2dense_opt(size_t nshares, size_t coeffs_size,
                               int16_t *dense, size_t dense_msk_stride,
                               size_t dense_data_stride,
                               const uint32_t *bitslice,
                               size_t bitslice_msk_stride,
                               size_t bitslice_data_stride) {
  start_bench(my_bs2dense);
  uint32_t a[32];
  for (size_t d = 0; d < nshares; d++) {
    for (size_t i = 0; i < coeffs_size; i++) {
      a[i] = bitslice[d * bitslice_msk_stride + i * bitslice_data_stride];
      a[i + 16] = bitslice[d * bitslice_msk_stride +
                           (i + coeffs_size) * bitslice_data_stride];
    }
    // Avoid uninitialized vars -> UB :(
    for (size_t i = coeffs_size; i < 16; i++) {
      a[i] = 0;
      a[i + 16] = 0;
    }
    transpose32(a);
    for (size_t i = 0; i < 32; i++) {
      dense[d * dense_msk_stride + i * dense_data_stride] =
          (a[i] >> 0) & ((1 << 16) - 1);
      dense[d * dense_msk_stride + (i + 32) * dense_data_stride] = a[i] >> 16;
    }
  }
  stop_bench(my_bs2dense);
}
