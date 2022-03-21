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
/* Based on the implementation "libkeccak-tiny" by David Leon Gil.
 * available at https://github.com/coruus/keccak-tiny under CC0 License.
 * */

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "bench.h"
#include "gadgets.h"
#include "masked.h"
#include "masked_fips202.h"

#define SHAKE256_RATE 136
#define SHA3_512_RATE 72

#define NROUNDS 24
#define ROL(a, offset) ((a << offset) ^ (a >> (64 - offset)))
#define KECCAK_NWORDS 25
#define Plen 200

typedef union {
  uint64_t w[NSHARES][KECCAK_NWORDS];
  uint32_t h[NSHARES][2 * KECCAK_NWORDS];
} MaskedKeccakState;

/******** The Keccak-f[1600] permutation ********/

/*** Constants. ***/
static const uint8_t rho[24] = {1,  3,  6,  10, 15, 21, 28, 36, 45, 55, 2,  14,
                                27, 41, 56, 8,  25, 43, 62, 18, 39, 61, 20, 44};
static const uint8_t pi[24] = {10, 7,  11, 17, 18, 3, 5,  16, 8,  21, 24, 4,
                               15, 23, 19, 13, 12, 2, 20, 14, 22, 9,  6,  1};
static const uint64_t RC[24] = {1ULL,
                                0x8082ULL,
                                0x800000000000808aULL,
                                0x8000000080008000ULL,
                                0x808bULL,
                                0x80000001ULL,
                                0x8000000080008081ULL,
                                0x8000000000008009ULL,
                                0x8aULL,
                                0x88ULL,
                                0x80008009ULL,
                                0x8000000aULL,
                                0x8000808bULL,
                                0x800000000000008bULL,
                                0x8000000000008089ULL,
                                0x8000000000008003ULL,
                                0x8000000000008002ULL,
                                0x8000000000000080ULL,
                                0x800aULL,
                                0x800000008000000aULL,
                                0x8000000080008081ULL,
                                0x8000000000008080ULL,
                                0x80000001ULL,
                                0x8000000080008008ULL};

/*** Helper macros to unroll the permutation. ***/
#define rol(x, s) (((x) << s) | ((x) >> (64 - s)))
#define REPEAT6(e) e e e e e e
#define REPEAT24(e) REPEAT6(e e e e)
#define REPEAT5(e) e e e e e
#define FOR5(v, s, e)                                                          \
  v = 0;                                                                       \
  REPEAT5(e; v += s;)

void masked_keccak(MaskedKeccakState *state) {
  uint8_t x, y;
  for (int i = 0; i < NROUNDS; i++) {
    // Sharewise implementation for Theta, Rho and phi
    for (int j = 0; j < NSHARES; j++) {
      uint64_t *a = &state->w[j][0];
      uint64_t b[5];
      uint64_t t = 0;
      // Theta
      FOR5(x, 1, b[x] = 0; FOR5(y, 5, b[x] ^= a[x + y];))
      FOR5(x, 1,
           FOR5(y, 5, a[y + x] ^= b[(x + 4) % 5] ^ ROL(b[(x + 1) % 5], 1);))
      // Rho and pi
      t = a[1];
      x = 0;
      REPEAT24(b[0] = a[pi[x]]; a[pi[x]] = ROL(t, rho[x]); t = b[0]; x++;)
    }
    // Chi: non-linear -> not sharewise.
    // Masked gadgets are implemented on 32-bit words and Chi does not contain
    // rotations, so we can work on 32-bit words
    for (y = 0; y < 25; y += 5) {
      for (int off = 0; off < 2; off++) {
        uint32_t sb_state[5 * NSHARES];
        size_t sb_state_msk_stride = 1;        // in 32-bit words
        size_t sb_state_data_stride = NSHARES; // in 32-bit words
        uint32_t *sb_in = &state->h[0][2 * y + off];
        size_t sb_in_data_stride = 2;     // in 32-bit words
        size_t sb_in_msk_stride = 2 * 25; // in 32-bit words

        for (x = 0; x < 5; x++) {
          copy_sharing(
              NSHARES, sb_state + x * sb_state_data_stride, sb_state_msk_stride,
              sb_in + ((x + 1) % 5) * sb_in_data_stride, sb_in_msk_stride);
          sb_state[x * sb_state_data_stride] =
              ~sb_state[x * sb_state_data_stride]; // NOT: on a single share
          masked_and(
              NSHARES, sb_state + x * sb_state_data_stride, sb_state_msk_stride,
              sb_state + x * sb_state_data_stride, sb_state_msk_stride,
              sb_in + ((x + 2) % 5) * sb_in_data_stride, sb_in_msk_stride);
        }
        for (x = 0; x < 5; x++) {
          masked_xor(NSHARES, sb_in + x * sb_in_data_stride, sb_in_msk_stride,
                     sb_in + x * sb_in_data_stride, sb_in_msk_stride,
                     sb_state + x * sb_state_data_stride, sb_state_msk_stride);
        }
      }
    }
    // Iota
    // Add constant: on a single share
    state->w[0][0] ^= RC[i];
  }
}

#define XORU64(value, address, byte)                                           \
  do {                                                                         \
    (value)[(address) >> 3] ^= (((uint64_t)(byte)) << 8 * ((address)&0x7));    \
  } while (0)
#define ExtractU64(value, address)                                             \
  (((value)[(address) >> 3] >> 8 * ((address)&0x7)) & 0xFF)
void masked_hash_keccak(uint8_t *out, size_t outlen, size_t out_msk_stride,
                        size_t out_data_stride, const uint8_t *in, size_t inlen,
                        size_t in_msk_stride, size_t in_data_stride,
                        size_t rate, uint8_t delim) {
  start_bench(keccak);
  MaskedKeccakState state;
  memset(&state.w[0][0], 0, sizeof(state));
  uint64_t *msk_a = &state.w[0][0];
  // Absorb input.
  while (inlen >= rate) {
    for (size_t i = 0; i < rate; i++) {
      for (size_t j = 0; j < NSHARES; j++) {
        XORU64(msk_a, i + j * Plen, in[j * in_msk_stride + i * in_data_stride]);
      }
    }
    masked_keccak(&state);
    in += rate * in_data_stride;
    inlen -= rate;
  }
  // Xor in the last block.
  for (size_t i = 0; i < inlen; i++) {
    for (size_t j = 0; j < NSHARES; j++) {
      XORU64(msk_a, i + j * Plen, in[j * in_msk_stride + i * in_data_stride]);
    }
  }
  // Xor in the DS and pad frame.
  XORU64(msk_a, inlen, delim);
  XORU64(msk_a, rate - 1, 0x80);
  // Apply P
  masked_keccak(&state);
  // Squeeze output.
  while (outlen >= rate) {
    for (size_t i = 0; i < rate; i++) {
      for (size_t j = 0; j < NSHARES; j++) {
        out[i * out_data_stride + j * out_msk_stride] =
            ExtractU64(msk_a, i + j * Plen);
      }
    }
    masked_keccak(&state);
    out += rate * out_data_stride;
    outlen -= rate;
  }
  for (size_t i = 0; i < outlen; i++) {
    for (size_t j = 0; j < NSHARES; j++) {
      out[i * out_data_stride + j * out_msk_stride] =
          ExtractU64(msk_a, i + j * Plen);
    }
  }
  stop_bench(keccak);
}

/*************************************************
 * Name:        masked_shake256
 *
 * Description: SHAKE256 XOF with non-incremental API
 *
 * Arguments:   - uint8_t *output:      pointer to output
 *              - size_t outlen:        requested output length in bytes
 *              - size_t out_msk_stride: stride of output shares.
 *              - size_t out_data_stride: stride of output data bytes.
 *              - const uint8_t *input: pointer to input
 *              - size_t inlen:         length of input in bytes
 *              - size_t in_msk_stride: stride of input shares.
 *              - size_t in_data_stride: stride of input data bytes.
 **************************************************/
void masked_shake256(uint8_t *output, size_t outlen, size_t out_msk_stride,
                     size_t out_data_stride, const uint8_t *input, size_t inlen,
                     size_t in_msk_stride, size_t in_data_stride) {
  masked_hash_keccak(output, outlen, out_msk_stride, out_data_stride, input,
                     inlen, in_msk_stride, in_data_stride, 136, 0x1f);
}

/*************************************************
 * Name:        masked_sha3_512
 *
 * Description: SHA3-512 with non-incremental API
 *
 * Arguments:   - uint8_t *output:      pointer to output
 *              - size_t out_msk_stride: stride of output shares.
 *              - size_t out_data_stride: stride of output data bytes.
 *              - const uint8_t *input: pointer to input
 *              - size_t inlen:         length of input in bytes
 *              - size_t in_msk_stride: stride of input shares.
 *              - size_t in_data_stride: stride of input data bytes.
 **************************************************/
void masked_sha3_512(uint8_t *output, size_t out_msk_stride,
                     size_t out_data_stride, const uint8_t *input, size_t inlen,
                     size_t in_msk_stride, size_t in_data_stride) {
  masked_hash_keccak(output, 64, out_msk_stride, out_data_stride, input, inlen,
                     in_msk_stride, in_data_stride, 72, 0x06);
}
