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
#include "masked_poly.h"
#include "bench.h"
#include "masked.h"
#include "masked_representations.h"
#include "poly.h"

#include "cbd.h"
#include "gadgets.h"
#include "masked_symmetric-fips202.h"
#include "ntt.h"
#include "params.h"
#include "symmetric.h"

#include <stdint.h>

/*************************************************
 * Name:        masked_poly_ntt
 *
 * Description: Computes negacyclic number-theoretic transform (NTT) of
 *              a polynomial in place;
 *              inputs assumed to be in normal order, output in bitreversed
 *order
 *
 * Arguments:   - uint16_t *r: pointer to in/output polynomial
 **************************************************/
void masked_poly_ntt(StrAPoly r) {
  start_bench(my_ntt);
  for (int d = 0; d < NSHARES; d++) {
    ntt(r[d]);
  }
  stop_bench(my_ntt);
}

/*************************************************
 * Name:        masked_poly_invntt
 *
 * Description: Computes inverse of negacyclic number-theoretic transform (NTT)
 *of a polynomial in place; inputs assumed to be in bitreversed order, output in
 *normal order
 *
 * Arguments:   - uint16_t *a: pointer to in/output polynomial
 **************************************************/
void masked_poly_invntt(StrAPoly r) {
  start_bench(my_ntt);
  for (int d = 0; d < NSHARES; d++) {
    invntt(r[d]);
  }
  stop_bench(my_ntt);
}

/*************************************************
 * Name:        masked_poly_getnoise
 *
 * Description: Sample a polynomial deterministically from a seed and a nonce,
 *              with output polynomial close to centered binomial distribution
 *              with parameter KYBER_ETA
 *
 * Arguments:   - poly *r:                   pointer to output polynomial
 *              - const unsigned char *seed: pointer to input seed (pointing to
 *array of length KYBER_SYMBYTES bytes)
 *              - unsigned char nonce:       one-byte input nonce
 *              - int add:                   boolean to indicate to accumulate
 *into r
 **************************************************/
void masked_poly_noise(StrAPoly r, const unsigned char *masked_seed,
                       size_t seed_msk_stride, size_t seed_data_stride,
                       unsigned char nonce, int add) {
  size_t kappa = KYBER_ETA;
  unsigned char buf_masked[(KYBER_ETA * KYBER_N / 4) * NSHARES];

  uint32_t a[kappa * NSHARES * 2];
  uint32_t b[kappa * NSHARES * 2];
  int16_t out[NSHARES * BSSIZE * 2];
  uint32_t bytes_off, by;

  masked_shake256_prf(buf_masked, KYBER_ETA * KYBER_N / 4,
                      KYBER_ETA * KYBER_N / 4, 1, masked_seed, seed_msk_stride,
                      seed_data_stride, nonce);
  // all the bitslice. 32*4 bits =
  for (uint32_t i = 0; i < KYBER_N / (2 * BSSIZE); i++) {
    // all the bits
    // 32*4 bits = 128 bits = 16 bytes
    for (uint32_t j = 0; j < kappa * NSHARES * 2; j++) {
      a[j] = 0;
      b[j] = 0;
    }

    for (int n = 0; n < 16; n += 1) {
      for (uint32_t d = 0; d < NSHARES; d++) {
        bytes_off = i * 32 + n;
        by = buf_masked[bytes_off + d * (kappa * KYBER_N / 4)];

        a[d] =
            (a[d] << 2) | (((by >> 0) & 0x1) << 1) | (((by >> 4) & 0x1) << 0);
        a[NSHARES + d] = (a[NSHARES + d] << 2) | (((by >> 1) & 0x1) << 1) |
                         (((by >> 5) & 0x1) << 0);

        b[d] =
            (b[d] << 2) | (((by >> 2) & 0x1) << 1) | (((by >> 6) & 0x1) << 0);
        b[NSHARES + d] = (b[NSHARES + d] << 2) | (((by >> 3) & 0x1) << 1) |
                         (((by >> 7) & 0x1) << 0);

        bytes_off = i * 32 + n + 16;
        by = buf_masked[bytes_off + d * (kappa * KYBER_N / 4)];

        a[d + kappa * NSHARES] = (a[d + kappa * NSHARES] << 2) |
                                 (((by >> 0) & 0x1) << 1) |
                                 (((by >> 4) & 0x1) << 0);
        a[NSHARES + d + kappa * NSHARES] =
            (a[NSHARES + d + kappa * NSHARES] << 2) | (((by >> 1) & 0x1) << 1) |
            (((by >> 5) & 0x1) << 0);

        b[d + kappa * NSHARES] = (b[d + kappa * NSHARES] << 2) |
                                 (((by >> 2) & 0x1) << 1) |
                                 (((by >> 6) & 0x1) << 0);
        b[NSHARES + d + kappa * NSHARES] =
            (b[NSHARES + d + kappa * NSHARES] << 2) | (((by >> 3) & 0x1) << 1) |
            (((by >> 7) & 0x1) << 0);
      }
    }
    masked_cbd(NSHARES, 2, KYBER_Q, COEF_NBITS, out, 1, NSHARES, a, 1, NSHARES,
               b, 1, NSHARES);

    for (uint32_t n = 0; n < BSSIZE; n++) {
      for (uint32_t j = 0; j < NSHARES; j++) {
        if (add) {
          r[j][(i * 64) + 31 - n] =
              (r[j][(i * 64) + 31 - n] + out[n * NSHARES + j]) % KYBER_Q;

          r[j][(i * 64) + 31 - n + 32] =
              (r[j][(i * 64) + 31 - n + 32] + out[(n + 32) * NSHARES + j]) %
              KYBER_Q;

        } else {
          r[j][(i * 64) + 31 - n] = out[n * NSHARES + j];
          r[j][(i * 64) + 31 - n + 32] = out[(n + 32) * NSHARES + j];
        }
      }
    }
  }
}

void masked_poly_tomsg(unsigned char *m, StrAPoly str_r) {
  start_bench(my_tomsg);
  APoly r;
  size_t i, j, d;
  uint32_t bits[NSHARES];

  StrAPoly2APoly(r, str_r);
  for (i = 0; i < KYBER_N; i += BSSIZE) {

    seccompress(NSHARES, KYBER_Q, 1, bits, 1, NSHARES, r[i], 1, NSHARES);

    for (d = 0; d < NSHARES; d++) {
      for (j = 0; j < BSSIZE / 8; j++) {
        m[d * KYBER_INDCPA_MSGBYTES + (i / 8) + j] =
            (bits[d] >> (j * 8)) & 0xFF;
      }
    }
  }
  stop_bench(my_tomsg);
}

/*************************************************
 * Name:       masked_poly_cmp
 *
 * Description: Compares masked polynomial with reference polynomial
 *
 * Arguments: - size_t c compression factor
 *            - uint32_t *rc: check bits array. Must be set to 0xFFFFFFFF if
 *              all the bits in two polynomials are equal.
 *            - StrAPoly mp: masked polynomial
 *            - poly *ref: reference polynomial
 **************************************************/

void masked_poly_cmp(size_t c, uint32_t *rc, const StrAPoly mp,
                     const poly *ref) {

  if(c == KYBER_DU){ start_bench(comp_du);}
  if(c == KYBER_DV){ start_bench(comp_dv);}
  start_bench(my_masked_poly_cmp);
  APoly r;
  size_t i, b;
  uint32_t bits[2 * NSHARES * c];
  uint32_t bits_ref[2 * c];

  StrAPoly2APoly(r, mp);

  for (i = 0; i < KYBER_N; i += BSSIZE * 2) {

    // compress masked polynomial
    seccompress(NSHARES, KYBER_Q, c, bits, 1, NSHARES, r[i], 1, NSHARES);
    seccompress(NSHARES, KYBER_Q, c, &bits[NSHARES * c], 1, NSHARES,
                r[i + BSSIZE], 1, NSHARES);

    // map public polynomial to bitslice
    masked_dense2bitslice_opt(1, c, bits_ref, 1, 1, &(ref->coeffs[i]), 1, 1);

    for (b = 0; b < 2 * c; b++) {
      // public polynomial and public one
      bits[b * NSHARES] ^= bits_ref[b] ^ 0xFFFFFFFF;

      masked_and(NSHARES, rc, 1, rc, 1, &bits[b * NSHARES], 1);
    }
  }
  stop_bench(my_masked_poly_cmp);
  if(c == KYBER_DU){ stop_bench(comp_du);}
  if(c == KYBER_DV){ stop_bench(comp_dv);}

}

/*************************************************
 * Name:       finalize
 *
 * Description: inplace ANDs all the 32-bits within bits
 *
 * Arguments: - uint32_t *bits: bits to ANDs
 * **************************************************/
void finalize_cmp(uint32_t *bits) {
  start_bench(my_cmp_finalize);
  uint32_t other[NSHARES];
  int d, shift;
  for (shift = 16; shift > 0; shift = shift >> 1) {
    for (d = 0; d < NSHARES; d++) {
      other[d] = bits[d] >> shift;
    }
    masked_and(NSHARES, bits, 1, bits, 1, other, 1);
  }
  stop_bench(my_cmp_finalize);
}

void masked_poly_frommsg(StrAPoly y,
                         const uint8_t m[KYBER_INDCPA_MSGBYTES * (NSHARES)],
                         size_t m_msk_stride, size_t m_data_stride) {

  start_bench(my_frommsg);
  uint32_t t1[NSHARES];
  int16_t t2[NSHARES];

  for (int i = 0; i < KYBER_N / 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      for (int k = 0; k < NSHARES; ++k)
        t1[k] = (m[i * m_data_stride + k * m_msk_stride] >> j) & 1;
      secb2a_1bit(NSHARES, t2, 1, t1, 1);

      for (int k = 0; k < NSHARES; ++k)
        y[k][i * 8 + j] = (t2[k] * ((KYBER_Q + 1) / 2)) % KYBER_Q;
    }
  }
  stop_bench(my_frommsg);
}
