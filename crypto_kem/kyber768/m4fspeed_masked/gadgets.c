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
#include "gadgets.h"
#include "masked.h"
#include "masked_representations.h"
#include "masked_utils.h"
#include <stdint.h>

static inline uint32_t pini_and_core(uint32_t a, uint32_t b, uint32_t r) {
  uint32_t temp;
  uint32_t s;
  asm("eor %[temp], %[b], %[r]\n\t"
      "and %[temp], %[a], %[temp]\n\t"
      "bic %[s], %[r], %[a]\n\t"
      "eor %[s], %[s], %[temp]"
      : [ s ] "=r"(s),
        [ temp ] "=&r"(
            temp) /* outputs, use temp as an arbitrary-location clobber */
      : [ a ] "r"(a), [ b ] "r"(b), [ r ] "r"(r) /* inputs */
  );
  return s;
}

void masked_and(size_t nshares, uint32_t *z, size_t z_stride, const uint32_t *a,
                size_t a_stride, const uint32_t *b, size_t b_stride) {
  uint32_t ztmp[nshares];
  uint32_t r;
  uint32_t i, j;

  for (i = 0; i < nshares; i++) {
    ztmp[i] = a[i * a_stride] & b[i * b_stride];
  }

  for (i = 0; i < (nshares - 1); i++) {
    for (j = i + 1; j < nshares; j++) {
      r = rand32();
      // PINI
      ztmp[i] ^= pini_and_core(a[i * a_stride], b[j * b_stride], r);
      ztmp[j] ^= pini_and_core(a[j * a_stride], b[i * b_stride], r);
    }
  }
  for (i = 0; i < nshares; i++) {
    z[i * z_stride] = ztmp[i];
  }
}

void masked_xor(size_t nshares, uint32_t *out, size_t out_stride,
                const uint32_t *ina, size_t ina_stride, const uint32_t *inb,
                size_t inb_stride) {
  for (size_t i = 0; i < nshares; i++) {
    out[i * out_stride] = ina[i * ina_stride] ^ inb[i * inb_stride];
  }
}

uint32_t unmask_boolean(size_t nshares, const uint32_t *in, size_t in_stride) {
  uint32_t out, d;
  out = 0;
  for (d = 0; d < nshares; d++) {
    out ^= in[d * in_stride];
  }
  return out;
}

void copy_sharing(size_t nshares, uint32_t *out, size_t out_stride,
                  const uint32_t *in, size_t in_stride) {
  for (size_t i = 0; i < nshares; i++) {
    out[i * out_stride] = in[i * in_stride];
  }
}

void secadd(size_t nshares, size_t kbits, size_t kbits_out, uint32_t *out,
            size_t out_msk_stride, size_t out_data_stride, const uint32_t *in1,
            size_t in1_msk_stride, size_t in1_data_stride, const uint32_t *in2,
            size_t in2_msk_stride, size_t in2_data_stride) {
  
  if(nshares == NSHARES){
    start_bench(my_secadd);
  }

  size_t i, d;
  uint32_t carry[nshares];
  uint32_t xpy[nshares];
  uint32_t xpc[nshares];

  masked_and(nshares, carry, 1, &in1[0 * in1_data_stride], in1_msk_stride,
             &in2[0 * in2_data_stride], in2_msk_stride);

  masked_xor(nshares, &out[0 * out_data_stride], out_msk_stride,
             &in1[0 * in1_data_stride], in1_msk_stride,
             &in2[0 * in2_data_stride], in2_msk_stride);

  for (i = 1; i < kbits; i++) {

    // xpy = in2 ^ in1
    // xpc = in1 ^ carry
    // out = xpy ^ carry
    for (d = 0; d < nshares; d++) {
      xpy[d] = in1[i * in1_data_stride + d * in1_msk_stride] ^
               in2[i * in2_data_stride + d * in2_msk_stride];
      xpc[d] = in1[i * in1_data_stride + d * in1_msk_stride] ^ carry[d];
      out[i * out_data_stride + d * out_msk_stride] = xpy[d] ^ carry[d];
    }

    if ((i == (kbits - 1)) && (i == (kbits_out - 1))) {
      break;
    } else if (i == (kbits - 1)) {
      masked_and(nshares, carry, 1, xpy, 1, xpc, 1);
      masked_xor(nshares, &out[(kbits)*out_data_stride], out_msk_stride, carry,
                 1, &in1[i * in1_data_stride], in1_msk_stride);
      break;
    }

    masked_and(nshares, carry, 1, xpy, 1, xpc, 1);
    masked_xor(nshares, carry, 1, carry, 1, &in1[i * in1_data_stride],
               in1_msk_stride);
  }

  if(nshares == NSHARES){
    stop_bench(my_secadd);
  }
}

void secadd_constant_bmsk(size_t nshares, size_t kbits, size_t kbits_out,
                          uint32_t *out, size_t out_msk_stride,
                          size_t out_data_stride, const uint32_t *in1,
                          size_t in1_msk_stride, size_t in1_data_stride,
                          uint32_t constant, const uint32_t *bmsk,
                          size_t bmsk_msk_stride) {

  size_t i, d;
  uint32_t carry[nshares];
  uint32_t xpy[nshares];
  uint32_t xpc[nshares];

  if (constant & 0x1) {
    masked_and(nshares, carry, 1, &in1[0 * in1_data_stride], in1_msk_stride,
               bmsk, bmsk_msk_stride);
    masked_xor(nshares, &out[0 * out_data_stride], out_msk_stride,
               &in1[0 * in1_data_stride], in1_msk_stride, bmsk,
               bmsk_msk_stride);
  } else {
    for (d = 0; d < nshares; d++) {
      carry[d] = 0;
    }
    copy_sharing(nshares, out, out_msk_stride, in1, in1_msk_stride);
  }
  for (i = 1; i < kbits; i++) {
    // xpy = in2 ^ in1
    // xpc = in1 ^ carry
    // out = xpy ^ carry
    if ((constant >> i) & 0x1) {
      for (d = 0; d < nshares; d++) {
        xpy[d] = in1[i * in1_data_stride + d * in1_msk_stride] ^
                 bmsk[d * bmsk_msk_stride];
        xpc[d] = in1[i * in1_data_stride + d * in1_msk_stride] ^ carry[d];
        out[i * out_data_stride + d * out_msk_stride] = xpy[d] ^ carry[d];
      }

      if ((i == (kbits - 1)) && (i == (kbits_out - 1))) {
        return;
      } else if (i == (kbits - 1)) {
        masked_and(nshares, carry, 1, xpy, 1, xpc, 1);
        masked_xor(nshares, &out[(kbits)*out_data_stride], out_msk_stride,
                   carry, 1, &in1[i * in1_data_stride], in1_msk_stride);
        return;
      }

      masked_and(nshares, carry, 1, xpy, 1, xpc, 1);
      masked_xor(nshares, carry, 1, carry, 1, &in1[i * in1_data_stride],
                 in1_msk_stride);
    } else {
      // compute the carry
      masked_xor(nshares, &out[i * out_data_stride], out_msk_stride, carry, 1,
                 &in1[i * in1_data_stride], 1);

      if ((i == (kbits - 1)) && (i == (kbits_out - 1))) {
        return;
      } else if (i == (kbits - 1)) {
        masked_and(nshares, &out[(kbits)*out_data_stride], out_msk_stride,
                   carry, 1, &in1[i * in1_data_stride], in1_msk_stride);
        return;
      }
      masked_and(nshares, carry, 1, carry, 1, &in1[i * in1_data_stride],
                 in1_msk_stride);
    }
  }
}

void secadd_constant(size_t nshares, size_t kbits, size_t kbits_out,
                     uint32_t *out, size_t out_msk_stride,
                     size_t out_data_stride, const uint32_t *in1,
                     size_t in1_msk_stride, size_t in1_data_stride,
                     uint32_t constant) {

  size_t i, d;
  uint32_t carry[nshares];
  uint32_t xpy[nshares];
  uint32_t xpc[nshares];

  if (constant & 0x1) {
    for (d = 0; d < nshares; d++) {
      carry[d] = in1[d * in1_msk_stride];
    }
    copy_sharing(nshares, out, out_msk_stride, in1, in1_msk_stride);
    out[0] ^= 0xFFFFFFFF;
  } else {
    for (d = 0; d < nshares; d++) {
      carry[d] = 0;
    }
    copy_sharing(nshares, out, out_msk_stride, in1, in1_msk_stride);
  }

  for (i = 1; i < kbits; i++) {
    if ((constant >> i) & 0x1) {
      for (d = 0; d < nshares; d++) {
        xpy[d] = in1[i * in1_data_stride + d * in1_msk_stride];
        xpc[d] = in1[i * in1_data_stride + d * in1_msk_stride] ^ carry[d];
        out[i * out_data_stride + d * out_msk_stride] = xpy[d] ^ carry[d];
      }
      xpy[0] ^= 0xFFFFFFFF;
      out[i * out_data_stride] ^= 0xFFFFFFFF;

      if ((i == (kbits - 1)) && (i == (kbits_out - 1))) {
        return;
      } else if (i == (kbits - 1)) {
        masked_and(nshares, carry, 1, xpy, 1, xpc, 1);
        masked_xor(nshares, &out[(kbits)*out_data_stride], out_msk_stride,
                   carry, 1, &in1[i * in1_data_stride], in1_msk_stride);

        // add the kbits_out of the constant
        out[(kbits)*out_data_stride] ^=
            0xFFFFFFFF * ((constant >> kbits) & 0x1);

        return;
      }

      masked_and(nshares, carry, 1, xpy, 1, xpc, 1);
      masked_xor(nshares, carry, 1, carry, 1, &in1[i * in1_data_stride],
                 in1_msk_stride);
    } else {
      // compute the carry
      masked_xor(nshares, &out[i * out_data_stride], out_msk_stride, carry, 1,
                 &in1[i * in1_data_stride], 1);

      if ((i == (kbits - 1)) && (i == (kbits_out - 1))) {
        return;
      } else if (i == (kbits - 1)) {
        masked_and(nshares, &out[(kbits)*out_data_stride], out_msk_stride,
                   carry, 1, &in1[i * in1_data_stride], in1_msk_stride);

        // add the kbits_out of the constant
        out[(kbits)*out_data_stride] ^=
            0xFFFFFFFF * ((constant >> kbits) & 0x1);
        return;
      }
      masked_and(nshares, carry, 1, carry, 1, &in1[i * in1_data_stride],
                 in1_msk_stride);
    }
  }
}

void secadd_modp(size_t nshares, size_t kbits, uint32_t q, uint32_t *out,
                 size_t out_msk_stride, size_t out_data_stride,
                 const uint32_t *in1, size_t in1_msk_stride,
                 size_t in1_data_stride, const uint32_t *in2,
                 size_t in2_msk_stride, size_t in2_data_stride) {

  uint32_t s[(kbits + 1) * nshares];
  uint32_t sp[(kbits + 1) * nshares];

  secadd(nshares, kbits, kbits + 1, s, 1, nshares, in1, in1_msk_stride,
         in1_data_stride, in2, in2_msk_stride, in2_data_stride);

  secadd_constant(nshares, kbits + 1, kbits + 1, sp, 1, nshares, s, 1, nshares,
                  (1 << (kbits + 1)) - q);

  secadd_constant_bmsk(nshares, kbits, kbits, out, out_msk_stride,
                       out_data_stride, sp, 1, nshares, q, &sp[kbits * nshares],
                       1);
}

void seca2b(size_t nshares, size_t kbits, uint32_t *in, size_t in_msk_stride,
            size_t in_data_stride) {

  // TODO optimize inplace ?
  if(nshares == NSHARES)
    start_bench(my_seca2b);

  size_t i, d;

  if (nshares == 1) {
    return;
  }

  size_t nshares_low = nshares / 2;
  size_t nshares_high = nshares - nshares_low;

  seca2b(nshares_low, kbits, in, in_msk_stride, in_data_stride);
  seca2b(nshares_high, kbits, &in[nshares_low], in_msk_stride, in_data_stride);

  uint32_t expanded_low[kbits * nshares];
  uint32_t expanded_high[kbits * nshares];

  for (i = 0; i < kbits; i++) {
    for (d = 0; d < nshares_low; d++) {
      expanded_low[i * nshares + d] =
          in[i * in_data_stride + d * in_msk_stride];
      expanded_high[i * nshares + d] = 0;
    }
    for (d = nshares_low; d < nshares; d++) {
      expanded_high[i * nshares + d] =
          in[i * in_data_stride + d * in_msk_stride];
      expanded_low[i * nshares + d] = 0;
    }
  }

  secadd(nshares, kbits, kbits, in, in_msk_stride, in_data_stride, expanded_low,
         1, nshares, expanded_high, 1, nshares);
 
  if(nshares == NSHARES)
    stop_bench(my_seca2b);
}

void seca2b_modp(size_t nshares, size_t kbits, uint32_t p, uint32_t *in,
                 size_t in_msk_stride, size_t in_data_stride) {

  size_t i, d;
  if(nshares == NSHARES)
    start_bench(my_seca2b_modp);

  if (nshares == 1) {
    return;
  }

  size_t nshares_low = nshares / 2;
  size_t nshares_high = nshares - nshares_low;

  seca2b_modp(nshares_low, kbits, p, in, in_msk_stride, in_data_stride);
  seca2b_modp(nshares_high, kbits, p, &in[nshares_low], in_msk_stride,
              in_data_stride);

  uint32_t expanded_low[(kbits + 1) * nshares];
  uint32_t expanded_high[(kbits + 1) * nshares];
  uint32_t u[(kbits + 1) * nshares];

  secadd_constant(nshares_low, kbits, kbits + 1, expanded_low, 1, nshares, in,
                  in_msk_stride, in_data_stride, (1 << (kbits + 1)) - p);

  for (i = 0; i < (kbits + 1); i++) {
    for (d = 0; d < nshares_low; d++) {
      // has already been written by secadd_constant_bmsk
      // expanded_low[i*nshares + d] = in[i*in_data_stride + d*in_msk_stride];
      expanded_high[i * nshares + d] = 0;
    }
    for (d = nshares_low; d < nshares; d++) {
      // kbits + 1 within in is unset
      expanded_high[i * nshares + d] =
          (i < (kbits)) ? in[i * in_data_stride + d * in_msk_stride] : 0;
      expanded_low[i * nshares + d] = 0;
    }
  }

  secadd(nshares, kbits + 1, kbits + 1, u, 1, nshares, expanded_high, 1,
         nshares, expanded_low, 1, nshares);

  secadd_constant_bmsk(nshares, kbits, kbits, in, in_msk_stride, in_data_stride,
                       u, 1, nshares, p, &u[kbits * nshares], 1);

  if(nshares == NSHARES)
    stop_bench(my_seca2b_modp);

}

// operates on 64 words the time
void secb2a_modp(size_t nshares,
                 //   size_t kbits, // MUST BE EQUAL TO COEF_NBITS
                 uint32_t p, uint32_t *in, size_t in_msk_stride,
                 size_t in_data_stride) {

  int16_t z_dense[2*BSSIZE * nshares];
  int16_t zp_dense[2*BSSIZE * nshares];
  uint32_t zp_str[2*COEF_NBITS * nshares];
  uint32_t b_str[2*COEF_NBITS * nshares];
  int16_t r[2];
  size_t d, i;

  // generate uniform sharing for z
  // zp = p - z;
  for (d = 0; d < nshares - 1; d++) {
    for (i = 0; i < BSSIZE; i += 2) {
      rand_q(r);
      z_dense[i * nshares + d] = r[0];
      zp_dense[i * nshares + d] = p - r[0];

      z_dense[(i + 1) * nshares + d] = r[1];
      zp_dense[(i + 1) * nshares + d] = p - r[1];

      rand_q(r);
      z_dense[i * nshares + d + (BSSIZE*nshares)] = r[0];
      zp_dense[i * nshares + d+ (BSSIZE*nshares)] = p - r[0];

      z_dense[(i + 1) * nshares + d+ (BSSIZE*nshares)] = r[1];
      zp_dense[(i + 1) * nshares + d+ (BSSIZE*nshares)] = p - r[1];

    }
#if ((BSSIZE & 0x1) == 0x1)
    rand_q(r);
    z_dense[(BSSIZE - 1) * nshares + d] = r[0];
    zp_dense[(BSSIZE - 1) * nshares + d] = p - r[0];

    rand_q(r);
    z_dense[(BSSIZE - 1) * nshares + d+ (BSSIZE*nshares)] = r[0];
    zp_dense[(BSSIZE - 1) * nshares + d+ (BSSIZE*nshares)] = p - r[0];
#endif
  }

  // map zp to bitslice representation
  masked_dense2bitslice_opt(nshares - 1, COEF_NBITS, zp_str, 1, nshares,
                        zp_dense, 1, nshares);

  // last shares of zp set to zero
  for (i = 0; i < COEF_NBITS; i++) {
    zp_str[i * nshares + (nshares - 1)] = 0;
    zp_str[i * nshares + (nshares - 1) + COEF_NBITS*nshares] = 0;
  }

  // last shares of zp_str to zero
  seca2b_modp(nshares, COEF_NBITS, p, zp_str, 1, nshares);
  seca2b_modp(nshares, COEF_NBITS, p, &zp_str[COEF_NBITS*nshares], 1, nshares);

  secadd_modp(nshares, COEF_NBITS, p, b_str, 1, nshares, in, in_msk_stride,
              in_data_stride, zp_str, 1, nshares);
  secadd_modp(nshares, COEF_NBITS, p, &b_str[COEF_NBITS*nshares], 1, nshares, &in[COEF_NBITS*nshares], in_msk_stride,
              in_data_stride, &zp_str[COEF_NBITS*nshares], 1, nshares);


  // map z to bistlice in output buffer
  masked_dense2bitslice_opt(nshares - 1, COEF_NBITS, in, in_msk_stride,
                        in_data_stride, z_dense, 1, nshares);

  // unmask b_str and set to the last share of the output
  for (i = 0; i < COEF_NBITS; i++) {
    RefreshIOS_rec(&b_str[i * nshares], nshares);
    RefreshIOS_rec(&b_str[i * nshares + COEF_NBITS*nshares], nshares);

    in[i * in_data_stride + (nshares - 1) * in_msk_stride] = 0;
    in[i * in_data_stride + (nshares - 1) * in_msk_stride + COEF_NBITS*nshares] = 0;
    for (d = 0; d < nshares; d++) {
      in[i * in_data_stride + (nshares - 1) * in_msk_stride] ^=
          b_str[i * nshares + d];
      in[i * in_data_stride + (nshares - 1) * in_msk_stride + COEF_NBITS*nshares] ^=
          b_str[i * nshares + d + COEF_NBITS*nshares];
    }
  }
}

void RefreshIOS_rec(uint32_t *x, uint32_t d) {
  uint32_t r;
  if (d == 1) {
  } else if (d == 2) {
    r = get_random();
    x[0] ^= r;
    x[1] ^= r;
  } else {
    RefreshIOS_rec(x, d / 2);
    RefreshIOS_rec(x + d / 2, d - d / 2);
    for (unsigned int i = 0; i < d / 2; i++) {
      r = rand32();
      x[i] ^= r;
      x[i + d / 2] ^= r;
    }
  }
}

void seccompress(size_t nshares, size_t ncoeffs, uint32_t q, uint32_t c,
                 uint32_t *out, size_t out_msk_stride, size_t out_data_stride,
                 const int16_t *in, size_t in_msk_stride,
                 size_t in_data_stride) {

  size_t i, d;
  uint32_t ell = 0;
  uint32_t prod = q * nshares;
  while (prod > 0) {
    prod = prod >> 1;
    ell++;
  }

  uint32_t in_expanded[ncoeffs * nshares];
  uint32_t bs_expanded[(ell + c) * nshares];
  uint32_t bs_expanded_ref[(ell + c) * nshares];

  // map mod q to mod 2^ell.
  uint32_t tmp32;
  uint64_t tmp64;
  for (i = 0; i < ncoeffs; i++) {
    tmp64 = (in[i * in_data_stride] + q) % q;
    tmp32 = ((((tmp64 << (c + ell + 1)) + q) >> 1) / (q));
    tmp32 += (1 << (ell - 1));
    in_expanded[i * nshares] = tmp32 & ((1 << (ell + c)) - 1);
    for (d = 1; d < nshares; d++) {
      tmp64 = (q + in[i * in_data_stride + d * in_msk_stride]) % q;
      tmp32 = ((((tmp64 << (c + ell + 1)) + q) >> 1) / (q));
      in_expanded[i * nshares + d] = tmp32 & ((1 << (ell + c)) - 1);
    }
  }

  // map to bitslice
  masked_dense2bitslice_opt_u32(NSHARES, ell + c, bs_expanded, 1, NSHARES,
                            in_expanded, 1, NSHARES);
  
  // convert A2B
  seca2b(NSHARES, ell + c, bs_expanded, 1, NSHARES);

  // map to the output
  for (i = 0; i < c; i++) {
    for (d = 0; d < nshares; d++) {
      out[i * out_data_stride + d * out_msk_stride] =
          bs_expanded[(ell + i) * nshares + d];
    }
  }
}

void secfulladd(size_t nshares, uint32_t *co, size_t co_msk_stride, uint32_t *w,
                size_t w_msk_stide, uint32_t *ci, size_t ci_msk_stride,
                uint32_t *x, size_t x_msk_stride, uint32_t *y,
                size_t y_msk_stride) {
  uint32_t a[nshares];
  uint32_t b[nshares];

  masked_xor(nshares, a, 1, x, x_msk_stride, y, y_msk_stride);

  masked_xor(nshares, b, 1, x, x_msk_stride, ci, ci_msk_stride);

  // compute carry out
  masked_and(nshares, a, 1, b, 1, a, 1);

  masked_xor(nshares, a, 1, a, 1, x, 1);

  // compute w
  masked_xor(nshares, w, w_msk_stide, b, 1, y, y_msk_stride);

  // write carry out
  for (size_t d = 0; d < nshares; d++) {
    co[d * co_msk_stride] = a[d];
  }
}

void masked_cbd(size_t nshares, size_t eta, size_t n_coeffs, size_t p,
                size_t kbits, int16_t *z, size_t z_msk_stride,
                size_t z_data_stride, uint32_t *a, size_t a_msk_stride,
                size_t a_data_stride, uint32_t *b, size_t b_msk_stride,
                size_t b_data_stride) {

  start_bench(my_cbd);
  size_t np = 2 * eta;
  size_t i, d, k, j, s;
  uint32_t sp[nshares * np], z_str_full[nshares * COEF_NBITS *2];
  uint32_t *a_in,*b_in,*z_str;

  //  k = ceil(log2(2*eta +1))
  k = 0;
  i = 2 * eta + 1;
  while (i != 0) {
    k++;
    i >>= 1;
  }

  // compte HW(a)-HW(b) for all 64 input coefficients 
  // levaraging 32 bus size
  for(s=0;s<2;s++){
    // copy input puts
    a_in = &a[s*(eta)*a_data_stride];
    b_in = &b[s*(eta)*b_data_stride];
    z_str = &z_str_full[s*nshares * COEF_NBITS];

    // compute HW(a) - HW(b) for current 32-bit slices
    for (i = 0; i < eta; i++) {
      for (d = 0; d < nshares; d++) {
        sp[(i * 2) * nshares + d] = a_in[i * a_data_stride + d * a_msk_stride];
        sp[(i * 2) * nshares + nshares + d] =
            b_in[i * b_data_stride + d * b_msk_stride];
        sp[(i * 2) * nshares + nshares + d] ^= (d == 0) ? 0xFFFFFFFF : 0;
      }
    }

    uint32_t *c;
    np = 2*eta;
    for (i = 0; i < k; i++) { // iterate on output bits

      c = &z_str[i * nshares];

      // init the carry
      for (d = 0; d < nshares; d++) {
        c[d] = (np & 0x1) ? sp[(np - 1) * nshares + d] : 0;
      }

      // sum np bits
      for (j = 0; j < np / 2; j++) {
        secfulladd(nshares, &sp[j * nshares], 1, c, 1, c, 1,
                   &sp[(1 + 2 * j) * nshares], 1, &sp[(2 * j) * nshares], 1);
      }
      np /= 2;
    }

    for (i = k; i < COEF_NBITS; i++) {
      for (d = 0; d < nshares; d++) {
        z_str[i * nshares + d] = 0;
      }
    }
  }

  secb2a_modp(nshares, p, z_str_full, 1, nshares);

  masked_bitslice2dense_opt(nshares, kbits, z, z_msk_stride,
                        z_data_stride, z_str_full, 1, nshares);

  for (i = 0; i < 2*n_coeffs; i++) {
    z[i * nshares] = (z[i * nshares] + p - eta) % p;
  }

  stop_bench(my_cbd);
}

// algo 7 in SOPG
static void secb2a_qbit_n(size_t nshares, int16_t *c, size_t c_msk_stride,
                          int16_t *a, size_t a_msk_stride, uint32_t x) {

  int16_t b[nshares];
  int16_t r[2];
  size_t j;

  b[nshares - 1] = 0;

  for (j = 0; j < nshares - 2; j += 2) {
    rand_q(r);

    b[j] = (a[j * a_msk_stride] - r[0] + KYBER_Q) % KYBER_Q;
    b[nshares - 1] = (b[nshares - 1] + r[0]) % KYBER_Q;

    b[j + 1] = (a[(j + 1) * a_msk_stride] - r[1] + KYBER_Q) % KYBER_Q;
    b[nshares - 1] = (b[nshares - 1] + r[1]) % KYBER_Q;
  }

  if (((nshares - 1) & 0x1) == 0x1) {
    rand_q(r);
    j = nshares - 2;
    b[j] = (a[j * a_msk_stride] - r[0] + KYBER_Q) % KYBER_Q;
    b[nshares - 1] = (b[nshares - 1] + r[0]) % KYBER_Q;
  }

  for (j = 0; j < nshares; j++) {
    c[j * c_msk_stride] = b[j] + 2 * KYBER_Q;
    c[j * c_msk_stride] -= (2 * b[j] * x);
    c[j * c_msk_stride] = c[j * c_msk_stride] % KYBER_Q;
  }
  c[0] = (c[0] + x) % KYBER_Q;
}

// algo 6
static void b2a_qbit(size_t nshares, int16_t *a, size_t a_msk_stride,
                     uint32_t *x, size_t x_msk_stride) {

  a[0] = x[0];
  for (size_t i = 1; i < nshares; i++) {
    secb2a_qbit_n(i + 1, a, a_msk_stride, a, a_msk_stride, x[i * x_msk_stride]);
  }
}

// algo 8
static void refresh_add(size_t nshares, int16_t *a, size_t a_msk_stride) {

  size_t i, j;
  int16_t rnd;
  size_t REFRESH_ADD_POOL_SIZE = (nshares * (nshares + 1) / 2);
  int16_t pool[REFRESH_ADD_POOL_SIZE];

  for (i = 0; i < REFRESH_ADD_POOL_SIZE; i += 2) {
    rand_q(&pool[i]);
  }

  if (((REFRESH_ADD_POOL_SIZE)&0x1) == 1) {
    int16_t r[2];
    rand_q(r);
    pool[REFRESH_ADD_POOL_SIZE - 1] = r[0];
  }

  int cpt = 0;
  for (i = 0; i < NSHARES - 1; i++) {
    for (j = i + 1; j < NSHARES; j++) {
      rnd = pool[cpt];
      cpt++;
      a[i * a_msk_stride] = (a[i * a_msk_stride] + rnd) % KYBER_Q;
      a[j * a_msk_stride] = (a[j * a_msk_stride] + KYBER_Q - rnd) % KYBER_Q;
    }
  }
}
// algo 5
void secb2a_1bit(size_t nshares, int16_t *a, size_t a_msk_stride, uint32_t *x,
                 size_t x_msk_stride) {

  b2a_qbit(nshares, a, a_msk_stride, x, x_msk_stride);
  refresh_add(nshares, a, a_msk_stride);
}
