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
#include "bench.h"
#include "indcpa.h"
#include "masked.h"
#include "masked_poly.h"
#include "masked_utils.h"
#include "ntt.h"
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "randombytes.h"
#include "symmetric.h"
#include "gadgets_legacy.h"

#include <stdint.h>
#include <string.h>

static unsigned mswitch(unsigned x, unsigned q_start, unsigned q_end) {
  return (2 * q_end * x + q_start) / (2 * q_start);
}

static unsigned compress(unsigned x, unsigned q, unsigned d) {
  return mswitch(x, q, 1 << d) % (1 << d);
}

extern void doublebasemul_asm_acc(int16_t *r, const int16_t *a,
                                  const int16_t *b, int16_t zeta);
/*************************************************
 * Name:        masked matacc
 *
 * Description: Multiplies a row of A or A^T, generated on-the-fly,
 *              with a vector of polynomials and accumulates into the result.
 *
 * Arguments:   - StrAPoly *r_masked:                    pointer to output
 *polynomial to accumulate in
 *              - StrAPolyVec *b_masked:                 pointer to input vector
 *of polynomials to multiply with
 *              - unsigned char i:            byte to indicate the index <
 *KYBER_K of the row of A or A^T
 *              - const unsigned char *seed:  pointer to the public seed used to
 *generate A
 *              - int transposed:             boolean indicatin whether A or A^T
 *is generated
 **************************************************/
static void masked_matacc(StrAPoly r_masked, StrAPolyVec b_masked,
                          unsigned char i, const unsigned char *seed,
                          int transposed) {
  start_bench(my_matacc);
  unsigned char buf[XOF_BLOCKBYTES + 2];
  unsigned int buflen, off;
  xof_state state;
  unsigned int ctr, pos, k, l, d;
  uint16_t val0, val1;
  int16_t c[4];

  for (d = 0; d < NSHARES; d++) {
    for (l = 0; l < KYBER_N; l++) {
      r_masked[d][l] = 0;
    }
  }

  for (int j = 0; j < KYBER_K; j++) {
    ctr = pos = 0;
    if (transposed)
      xof_absorb(&state, seed, i, j);
    else
      xof_absorb(&state, seed, j, i);

    xof_squeezeblocks(buf, 1, &state);
    buflen = XOF_BLOCKBYTES;

    k = 0;
    while (ctr < KYBER_N / 4) {
      val0 = ((buf[pos + 0] >> 0) | ((uint16_t)buf[pos + 1] << 8)) & 0xFFF;
      val1 = ((buf[pos + 1] >> 4) | ((uint16_t)buf[pos + 2] << 4)) & 0xFFF;
      pos += 3;

      if (val0 < KYBER_Q) {
        c[k++] = (int16_t)val0;
        if (k == 4) {
          for (d = 0; d < NSHARES; d++) {
            doublebasemul_asm_acc(&r_masked[d][4 * ctr],
                                  &b_masked[j][d][4 * ctr], c, zetas[ctr]);
          }
          ctr++;
          k = 0;
        }
      }

      if (val1 < KYBER_Q && ctr < KYBER_N / 4) {
        c[k++] = (int16_t)val1;
        if (k == 4) {
          for (d = 0; d < NSHARES; d++) {
            doublebasemul_asm_acc(&r_masked[d][4 * ctr],
                                  &b_masked[j][d][4 * ctr], c, zetas[ctr]);
          }
          ctr++;
          k = 0;
        }
      }

      if (pos + 3 > buflen && ctr < KYBER_N / 4) {
        off = buflen % 3;
        for (l = 0; l < off; l++)
          buf[l] = buf[buflen - off + l];
        xof_squeezeblocks(buf + off, 1, &state);
        buflen = off + XOF_BLOCKBYTES;
        pos = 0;
      }
    }
  }
  stop_bench(my_matacc);
}

/*************************************************
 * Name:        masked_indcpa_enc_cmp
 *
 * Description: Re-encryption function.
 *              Compares the re-encypted ciphertext with the original ciphertext
 *byte per byte. The comparison is performed in a constant time manner.
 *
 *
 * Arguments:   - unsigned char *ct:         pointer to input ciphertext to
 *compare the new ciphertext with (of length KYBER_INDCPA_BYTES bytes)
 *              - const unsigned char *m:    pointer to input message (of length
 *KYBER_INDCPA_MSGBYTES bytes)
 *              - const unsigned char *pk:   pointer to input public key (of
 *length KYBER_INDCPA_PUBLICKEYBYTES bytes)
 *              - const unsigned char *coin: pointer to input random coins used
 *as seed (of length KYBER_SYMBYTES bytes) to deterministically generate all
 *randomness Returns:     - boolean byte indicating that re-encrypted ciphertext
 *is NOT equal to the original ciphertext
 **************************************************/
unsigned char masked_indcpa_enc_cmp(const unsigned char *c,
                                    const unsigned char *m_masked,
                                    size_t m_msk_stride, size_t m_data_stride,
                                    const unsigned char *pk,
                                    const unsigned char *masked_coins,
                                    size_t coins_msk_stride,
                                    size_t coins_data_stride) {

  Masked to_compare[KYBER_N * (KYBER_K+1)];
  uint16_t ref_to_compare[KYBER_N * (KYBER_K+1)];

  poly bp;
  polyvec c_ref;
  poly *pkp = &bp;
  const unsigned char *seed = pk + KYBER_POLYVECBYTES;
  int i, d;
  unsigned char nonce = 0;

  StrAPolyVec masked_sp;
  StrAPoly masked_bp;
  StrAPoly masked_v;
  StrAPoly masked_k;

  for (i = 0; i < KYBER_K; i++) {
    masked_poly_noise(masked_sp[i], masked_coins, coins_msk_stride,
                      coins_data_stride, nonce++, 0);
    masked_poly_ntt(masked_sp[i]);
  }

  // unpack reference ciphertext
  polyvec_decompress(&c_ref, c);
  for (i = 0; i < KYBER_K; i++) {
    for (int j = 0; j < KYBER_N; j++) {
      c_ref.vec[i].coeffs[j] =
          compress(c_ref.vec[i].coeffs[j], KYBER_Q, KYBER_DU);
    }
  }
  poly v_ref;
  poly_decompress(&v_ref, c + KYBER_POLYVECCOMPRESSEDBYTES);
  for (int j = 0; j < KYBER_N; j++) {
    v_ref.coeffs[j] = compress(v_ref.coeffs[j], KYBER_Q, KYBER_DV);
    ref_to_compare[ KYBER_K * KYBER_N  + j] = v_ref.coeffs[j];
  }

  for (i = 0; i < KYBER_K; i++) {

    masked_matacc(masked_bp, masked_sp, i, seed, 1);
    masked_poly_invntt(masked_bp);
    masked_poly_noise(masked_bp, masked_coins, coins_msk_stride,
                      coins_data_stride, nonce++, 1);
  

    for(int j = 0; j < KYBER_N; j++){
      for(int d = 0; d < NSHARES; d++){
        to_compare[(i * KYBER_N) + j].shares[d] =  (masked_bp[d][j] + KYBER_Q)%KYBER_Q;
      }
      ref_to_compare[(i * KYBER_N)  + j] = c_ref.vec[i].coeffs[j];
    }
  }

  // multiply sp vector with public key vector
  poly_frombytes(pkp, pk);
  for (d = 0; d < NSHARES; d++) {
    poly_basemul_i16(masked_v[d], pkp->coeffs, masked_sp[0][d]);
  }
  for (i = 1; i < KYBER_K; i++) {
    poly_frombytes(pkp, pk + i * KYBER_POLYBYTES);
    for (d = 0; d < NSHARES; d++) {
      poly_basemul_acc_i16(masked_v[d], pkp->coeffs, masked_sp[i][d]);
    }
  }

  masked_poly_invntt(masked_v);

  masked_poly_noise(masked_v, masked_coins, coins_msk_stride, coins_data_stride,
                    nonce++, 1);

  masked_poly_frommsg(masked_k, m_masked, m_msk_stride, m_data_stride);
  for (d = 0; d < NSHARES; d++) {
    poly_add_i16(masked_v[d], masked_v[d], masked_k[d]);
    poly_reduce_i16(masked_v[d]);
  }

  for(int j = 0; j < KYBER_N; j ++){
    for(int d = 0; d < NSHARES; d ++){
      to_compare[KYBER_K * KYBER_N + j].shares[d] = masked_v[d][j];
    }
  }

  int rc;
  rc = kyber_poly_comp_hybrid(to_compare,ref_to_compare);
  return (unsigned char)(!rc);
}

/*************************************************
 * Name:        masked_indcpa_dec
 *
 * Description: Decryption function of the CPA-secure
 *              public-key encryption scheme underlying Kyber.
 *
 * Arguments:   - unsigned char *m:        pointer to output decrypted message
 *(of length KYBER_INDCPA_MSGBYTES)
 *              - const unsigned char *c:  pointer to input ciphertext (of
 *length KYBER_INDCPA_BYTES)
 *              - const unsigned char *sk: pointer to input secret key (of
 *length KYBER_INDCPA_SECRETKEYBYTES)
 **************************************************/
void __attribute__((noinline))
masked_indcpa_dec(unsigned char *m, // secret
                  size_t o_msk_stride, size_t o_data_stride,
                  const unsigned char *c,    // public
                  const unsigned char *sk) { // secret
  poly mp, bp, psk;
  poly *v = &bp;
  StrAPoly masked_psk, masked_mp, masked_bp;
  int d;
  poly_unpackdecompress(&mp, c, 0);
  poly_ntt(&mp);

  // TODO store masked version of the private key
  // mask private key
  poly_frombytes(&psk, sk);
  masked_poly(masked_psk, &psk);

  for (d = 0; d < NSHARES; d++) {
    poly_basemul_i16(masked_mp[d], mp.coeffs, masked_psk[d]);
  }

  unmasked_poly(&mp, masked_mp);

  for (int i = 1; i < KYBER_K; i++) {
    poly_unpackdecompress(&bp, c, i);
    poly_ntt(&bp);

    poly_frombytes(&psk, sk + i * KYBER_POLYBYTES);
    masked_poly(masked_psk, &psk);

    for (d = 0; d < NSHARES; d++) {
      poly_basemul_i16(masked_bp[d], bp.coeffs, masked_psk[d]);
      poly_add_i16(masked_mp[d], masked_mp[d], masked_bp[d]);
    }
  }

  masked_poly_invntt(masked_mp);

  poly_decompress(v, c + KYBER_POLYVECCOMPRESSEDBYTES);

  poly_sub_i16(masked_mp[0], v->coeffs, masked_mp[0]);
  poly_reduce_i16(masked_mp[0]);
  poly_zeroize(v);
  for (d = 1; d < NSHARES; d++) {
    // compute -mp for all the other shares
    poly_sub_i16(masked_mp[d], v->coeffs, masked_mp[d]);
    poly_reduce_i16(masked_mp[d]);
  }

  unsigned char m_masked[KYBER_INDCPA_MSGBYTES * NSHARES];
  masked_poly_tomsg(m_masked, masked_mp);

  for (d = 0; d < NSHARES; d++) {
    for (int i = 0; i < KYBER_INDCPA_MSGBYTES; i++) {
      m[i * o_data_stride + d * o_msk_stride] =
          m_masked[d * KYBER_INDCPA_MSGBYTES + i];
    }
  }

  return;
}
