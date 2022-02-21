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
#include "cbd.h"
#include "mNTT.h"
#include "masked_cbd.h"
#include <string.h>
#include "bench.h"
#include "fips202.h"
#include "gadgets.h"
#include "masked.h"
#include "masked_fips202.h"
#include "masked_representations.h"
#include "masked_utils.h"
#include "pack_unpack.h"

#define h1 (1 << (SABER_EQ - SABER_EP - 1))
#define h2                                                                     \
  ((1 << (SABER_EP - 2)) - (1 << (SABER_EP - SABER_ET - 1)) +                  \
   (1 << (SABER_EQ - SABER_EP - 1)))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

extern void __asm_poly_add_16(uint16_t *des, uint16_t *src1, uint16_t *src2);
extern void __asm_poly_add_32(uint32_t *des, uint32_t *src1, uint32_t *src2);

static inline shake128incctx
shake128_absorb_seed(const uint8_t seed[SABER_SEEDBYTES]) {
  shake128incctx ctx;
  shake128_inc_init(&ctx);
  shake128_inc_absorb(&ctx, seed, SABER_SEEDBYTES);
  shake128_inc_finalize(&ctx);
  return ctx;
}

/*************************************************
 * Name:        masked_InnerProdDecNTT
 *
 * Description: Performs Inner product between masked key and public 
 *            ciphertext.
 *
 * Arguments: - uint8_t *m: compression of inner product. Size of SABER_KEYBYTES
 *            - size_t m_msk_stide: m shares stride
 *            - size_t m_data_stride: m data stride
 *            - const uint8_t *: public ciphertext
 *            - StrAPolyVec: masked secret key
 **************************************************/
void masked_InnerProdDecNTT(uint8_t *m, size_t m_msk_stide,
                            size_t m_data_stride,
                            const uint8_t ciphertext[SABER_BYTES_CCA_DEC],
                            StrAPolyVec sk_masked) {

  // NTT of ciphertext
  uint32_t c_NTT_32[SABER_L][SABER_N];
  uint16_t c_NTT_16[SABER_L][SABER_N];
  StrAPoly m_poly;
  uint16_t cm[SABER_N];
  uint32_t masked_bs[NSHARES * SABER_EP * 2]; // TODO
  size_t i, j, b;

  // decompress and NTT public ciphertext.
  for (i = 0; i < SABER_L; i++) {
    uint16_t poly[SABER_N];
    BS2POLp(ciphertext + i * SABER_POLYCOMPRESSEDBYTES, poly);
    
    start_bench(my_ntt);
    NTT_forward_32(c_NTT_32[i], poly);
    NTT_forward_16(c_NTT_16[i], poly);
    stop_bench(my_ntt);
  }

  start_bench(my_ntt);
  // inner product: Ciphertex * sk
  for (j = 0; j < NSHARES; j++) {
    uint32_t acc_NTT_32[SABER_N];
    uint16_t acc_NTT_16[SABER_N];
    for (i = 0; i < SABER_L; i++) {
      uint32_t poly_NTT_32[SABER_N];
      uint16_t poly_NTT_16[SABER_N];

      NTT_forward_32(poly_NTT_32, sk_masked[i][j]);
      NTT_forward_16(poly_NTT_16, sk_masked[i][j]);
      if (i == 0) {
        NTT_mul_32(acc_NTT_32, poly_NTT_32, c_NTT_32[i]);
        NTT_mul_16(acc_NTT_16, poly_NTT_16, c_NTT_16[i]);
      } else {
        NTT_mul_acc_32(acc_NTT_32, poly_NTT_32, c_NTT_32[i]);
        NTT_mul_acc_16(acc_NTT_16, poly_NTT_16, c_NTT_16[i]);
      }
    }
    NTT_inv_32(acc_NTT_32);
    NTT_inv_16(acc_NTT_16);
    solv_CRT(m_poly[j], acc_NTT_32, acc_NTT_16);
  }
  stop_bench(my_ntt);

  start_bench(my_tomsg);
  BS2POLT(ciphertext + SABER_POLYVECCOMPRESSEDBYTES, cm);
  for (i = 0; i < SABER_N; i++) {
    m_poly[0][i] =
        (SABER_P + m_poly[0][i] + h2 - (cm[i] << (SABER_EP - SABER_ET))) & ((1<<SABER_EP)-1);
  }

  // compression of mpoly
  for (i = 0; i < SABER_N; i += 2 * BSSIZE) {
    masked_dense2bitslice_opt(NSHARES, SABER_EP, masked_bs, 1, NSHARES,
                              &m_poly[0][i], SABER_N, 1);
    seca2b(NSHARES, SABER_EP, masked_bs, 1, NSHARES);
    seca2b(NSHARES, SABER_EP, &masked_bs[NSHARES * SABER_EP], 1, NSHARES);
    masked_bitslice2dense_opt(NSHARES, 1, &m_poly[0][i], SABER_N, 1,
                              &masked_bs[(SABER_EP - 1) * NSHARES], 1,
                              NSHARES * SABER_EP);
  }

  // map compressed mpoly into m
  for (j = 0; j < NSHARES; j++) {
    for (i = 0; i < SABER_KEYBYTES; i++) {
      m[i * m_data_stride + j * m_msk_stide] = 0;
      for (b = 0; b < 8; b++) {
        m[i * m_data_stride + j * m_msk_stide] |= m_poly[j][i * 8 + b] << b;
      }
    }
  }
  stop_bench(my_tomsg);
}

/*************************************************
 * Name:        masked_MatrixVectorMulEncNTT_cmp
 *
 * Description: Performs matrix product and compare the result with ct1 and ct0
 *
 * Arguments: - uint8_t *ct0: reference ciphertext
 *            - uint8_t *ct1: reference cipehrtext
 *            - uint8_t *seed_s: buffer containing seed for CBD sampling
 *            - size_t seed_s_msk_stride: seed_s msk stride
 *            - size_t seed_s_data_stride: seed_s data stride
 *            - uint8_t *seed_A: buffer containing seed to derive public matrix
 *            - uint8_t *pk: public key
 *            - uint8_t *m: Message to decrypt. Size of SABER_KEYBYTES
 *            - size_t m_msk_stide: m shares stride
 *            - size_t m_data_stride: m data stride
 **************************************************/
uint32_t masked_MatrixVectorMulEncNTT_cmp(
    uint8_t ct0[SABER_POLYVECCOMPRESSEDBYTES],
    uint8_t ct1[SABER_SCALEBYTES_KEM], const uint8_t *seed_s,
    size_t seed_s_msk_stride, size_t seed_s_data_stride,
    const uint8_t seed_A[SABER_SEEDBYTES],
    const uint8_t pk[SABER_INDCPA_PUBLICKEYBYTES], const uint8_t *m,
    size_t m_msk_stide, size_t m_data_stride) {

  // s ac and A in NTT domain.
  uint32_t acc_NTT_32[NSHARES][SABER_N];
  uint32_t A_NTT_32[SABER_N];
  uint32_t s_NTT_32[NSHARES][SABER_L * SABER_N];

  uint16_t acc_NTT_16[NSHARES][SABER_N];
  uint16_t A_NTT_16[SABER_N];
  uint16_t s_NTT_16[NSHARES][SABER_L * SABER_N];

  uint16_t poly[SABER_N]; // public polynomial (from A)
  uint16_t m_poly[NSHARES][SABER_N];
  Poly myref; // reference polynomial to compare
  StrAPoly masked_acc;

  uint8_t shake_out[MAX(SABER_POLYBYTES, SABER_POLYCOINBYTES)];
  uint8_t masked_shake_out[MAX(SABER_POLYBYTES, SABER_POLYCOINBYTES) * NSHARES];

  size_t i, j, d;

  // comparison flags
  uint32_t correct = 0;
  uint32_t masked_correct[NSHARES];
  memset(masked_correct, 0, sizeof(masked_correct));
  masked_correct[0] = 0xFFFFFFFF;

  MaskedShakeCtx masked_shake_s_ctx;

  masked_shake128_inc_init(&masked_shake_s_ctx, seed_s, SABER_SEEDBYTES,
                           seed_s_msk_stride, seed_s_data_stride);

  for (i = 0; i < SABER_L; i++) {
    masked_shake128_squeeze(&masked_shake_s_ctx, masked_shake_out,
                            SABER_POLYCOINBYTES, SABER_POLYCOINBYTES, 1);

    masked_cbd_seed(NSHARES, &m_poly[0][0], SABER_N, 1, masked_shake_out,
                    SABER_POLYCOINBYTES, 1);

    start_bench(my_ntt);
    for (d = 0; d < NSHARES; d++) {
      NTT_forward_32(s_NTT_32[d] + i * SABER_N, m_poly[d]);
      NTT_forward_16(s_NTT_16[d] + i * SABER_N, m_poly[d]);
    }
    stop_bench(my_ntt);
  }

  // init Shake to generate A matrix
  start_bench(my_matacc);
  shake128incctx shake_A_ctx = shake128_absorb_seed(seed_A);
  stop_bench(my_matacc);


  for (i = 0; i < SABER_L; i++) {

    for (j = 0; j < SABER_L; j++) {

      start_bench(my_matacc);
      shake128_inc_squeeze(shake_out, SABER_POLYBYTES, &shake_A_ctx);
      stop_bench(my_matacc);
      BS2POLq(shake_out, poly);

      start_bench(my_ntt);
      NTT_forward_32(A_NTT_32, poly);
      NTT_forward_16(A_NTT_16, poly);

      for (d = 0; d < NSHARES; d++) {
        if (j == 0) {
          NTT_mul_32(acc_NTT_32[d], A_NTT_32, s_NTT_32[d] + j * SABER_N);
          NTT_mul_16(acc_NTT_16[d], A_NTT_16, s_NTT_16[d] + j * SABER_N);
        } else {
          NTT_mul_acc_32(acc_NTT_32[d], A_NTT_32, s_NTT_32[d] + j * SABER_N);
          NTT_mul_acc_16(acc_NTT_16[d], A_NTT_16, s_NTT_16[d] + j * SABER_N);
        }
      }
      stop_bench(my_ntt);
    }

    start_bench(my_ntt);
    for (d = 0; d < NSHARES; d++) {
      NTT_inv_32(acc_NTT_32[d]);
      NTT_inv_16(acc_NTT_16[d]);
      solv_CRT(masked_acc[d], acc_NTT_32[d], acc_NTT_16[d]);
    }
    stop_bench(my_ntt);

    for (j = 0; j < SABER_N; j++) {
      masked_acc[0][j] = (masked_acc[0][j] + h1) % SABER_Q;
    }

    // compare acc with ct0
    BS2POLp(ct0 + i * SABER_POLYCOMPRESSEDBYTES, myref);
    for (j = 0; j < SABER_N; j++) {
      myref[j] = myref[j] % SABER_P;
    }
    masked_poly_cmp(SABER_EQ - SABER_EP, SABER_EQ, SABER_EQ, masked_correct,
                    &masked_acc[0][0], SABER_N, 1, myref);
  }


  shake128_inc_ctx_release(&shake_A_ctx);

  start_bench(my_ntt);
  for (j = 0; j < SABER_L; j++) {

    BS2POLp(pk + j * SABER_POLYCOMPRESSEDBYTES, poly);

    NTT_forward_32(A_NTT_32, poly);
    NTT_forward_16(A_NTT_16, poly);

    for (d = 0; d < NSHARES; d++) {
      if (j == 0) {
        NTT_mul_32(acc_NTT_32[d], A_NTT_32, s_NTT_32[d] + j * SABER_N);
        NTT_mul_16(acc_NTT_16[d], A_NTT_16, s_NTT_16[d] + j * SABER_N);
      } else {
        NTT_mul_acc_32(acc_NTT_32[d], A_NTT_32, s_NTT_32[d] + j * SABER_N);
        NTT_mul_acc_16(acc_NTT_16[d], A_NTT_16, s_NTT_16[d] + j * SABER_N);
      }
    }
  }
  
  for (d = 0; d < NSHARES; d++) {
    NTT_inv_32(acc_NTT_32[d]);
    NTT_inv_16(acc_NTT_16[d]);
    solv_CRT(masked_acc[d], acc_NTT_32[d], acc_NTT_16[d]);
  }
  stop_bench(my_ntt);

  for (j = 0; j < SABER_N; j++) {
    // work in SABER_Q as for NTT. Could be done in SABER_P.
    masked_acc[0][j] =
        (masked_acc[0][j] -
         (((m[(j / 8) * m_data_stride] >> (j & 0x7)) & 0x1) << (SABER_EP - 1)) +
         h1) %
        SABER_Q;
    for (d = 1; d < NSHARES; d++) {
      masked_acc[d][j] =
          (masked_acc[d][j] -
           (((m[d * m_msk_stide + (j / 8) * m_data_stride] >> (j & 0x7)) & 0x1)
            << (SABER_EP - 1))) %
          SABER_Q;
    }
  }

  // compare acc with ct1
  BS2POLT(ct1, myref);
  for (j = 0; j < SABER_N; j++) {
    myref[j] = myref[j] % (1 << SABER_ET);
  }
  masked_poly_cmp(SABER_EP - SABER_ET, SABER_EP, SABER_EP, masked_correct,
                  &masked_acc[0][0], SABER_N, 1, myref);

  // finalize the comparison
  finalize_cmp(masked_correct);
  correct = 0;
  for (d = 0; d < NSHARES; d++) {
    correct ^= masked_correct[d];
  }
  return !correct;
}

/*************************************************
 * Name:       masked_poly_cmp 
 *
 * Description: Compares masked polynomial with reference polynomial 
 *
 * Arguments: - size_t b_start: first bit to compare
 *            - size_t b_end: last bit to compare
 *            - size_t coeffs_size: number of bits in polynomial modulus
 *            - uint32_t *rc: check bits array
 *            - const uint16_t *mp: masked polynomial 
 *            - size_t mp_msk_stide: m shares stride
 *            - size_t mp_data_stride: m data stride
 *            - Poly ref: reference polynomial
 **************************************************/
void masked_poly_cmp(size_t b_start, size_t b_end, size_t coeffs_size,
                     uint32_t *rc,
                     const uint16_t *mp,size_t mp_msk_stride, size_t mp_data_stride, 
                     Poly ref) {

  start_bench(my_masked_poly_cmp);
  size_t i, b;
  uint32_t bits[2 * NSHARES * coeffs_size];
  uint32_t bits_ref[2 * coeffs_size];

  for (i = 0; i < SABER_N; i += BSSIZE * 2) {

    // convert masked poly
    masked_dense2bitslice_opt(NSHARES, coeffs_size, bits, 1, NSHARES, mp,
                              mp_msk_stride, mp_data_stride);

    seca2b(NSHARES, coeffs_size, bits, 1, NSHARES);
    seca2b(NSHARES, coeffs_size, &bits[NSHARES * coeffs_size], 1, NSHARES);

    // map public polynomial to bitslice
    masked_dense2bitslice_opt(1, coeffs_size, bits_ref, 1, 1, ref, 1, 1);

    for (b = 0; b < b_end - b_start; b++) {

      // public polynomial and public one
      bits[(b + b_start) * NSHARES] ^= ~bits_ref[b];
      masked_and(NSHARES, rc, 1, rc, 1, &bits[(b + b_start) * NSHARES], 1);

      bits[(b + b_start + coeffs_size) * NSHARES] ^= ~bits_ref[b + coeffs_size];
      masked_and(NSHARES, rc, 1, rc, 1,
                 &bits[(b + b_start + coeffs_size) * NSHARES], 1);
    }
  }
  stop_bench(my_masked_poly_cmp);
}


void finalize_cmp(uint32_t *bits) {

  start_bench(my_cmp_finalize);
  uint32_t other[NSHARES];
  int d;
  for (d = 0; d < NSHARES; d++) {
    other[d] = bits[d] >> 16;
  }
  masked_and(NSHARES, bits, 1, bits, 1, other, 1);

  for (d = 0; d < NSHARES; d++) {
    other[d] = bits[d] >> 8;
  }
  masked_and(NSHARES, bits, 1, bits, 1, other, 1);

  for (d = 0; d < NSHARES; d++) {
    other[d] = bits[d] >> 4;
  }
  masked_and(NSHARES, bits, 1, bits, 1, other, 1);

  for (d = 0; d < NSHARES; d++) {
    other[d] = bits[d] >> 2;
  }
  masked_and(NSHARES, bits, 1, bits, 1, other, 1);
  for (d = 0; d < NSHARES; d++) {
    other[d] = bits[d] >> 1;
  }
  masked_and(NSHARES, bits, 1, bits, 1, other, 1);
  stop_bench(my_cmp_finalize);
}
