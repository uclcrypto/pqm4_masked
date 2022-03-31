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
#include "gadgets.h"
#include "bench.h"
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

/*************************************************
 * Name:        masked_and
 *
 * Description: Performs masked AND (z = a & b ) gate with nshares.
 *
 * Arguments:   - size_t nshares: number of shares
 *            - uint32_t *z: output buffer
 *            - size_t z_stride: output buffer stride
 *            - uint32_t *a: first input buffer
 *            - size_t a_stride: a buffer stride
 *            - uint32_t *b: second input buffer
 *            - size_t b_stride: b buffer stride
 **************************************************/
void masked_and_c(size_t nshares, uint32_t *z, size_t z_stride, const uint32_t *a,
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

/*************************************************
 * Name:        masked_xor
 *
 * Description: Performs masked XOR (z = a ^ b ) gate with nshares.
 *
 * Arguments:   - size_t nshares: number of shares
 *            - uint32_t *z: output buffer
 *            - size_t z_stride: output buffer stride
 *            - uint32_t *a: first input buffer
 *            - size_t a_stride: a buffer stride
 *            - uint32_t *b: second input buffer
 *            - size_t b_stride: b buffer stride
 **************************************************/
void masked_xor_c(size_t nshares, uint32_t *out, size_t out_stride,
                const uint32_t *ina, size_t ina_stride, const uint32_t *inb,
                size_t inb_stride) {
  for (size_t i = 0; i < nshares; i++) {
    out[i * out_stride] = ina[i * ina_stride] ^ inb[i * inb_stride];
  }
}

/*************************************************
 * Name:        unmask_boolean
 *
 * Description: Unmask sharing by XORing all the shares
 *
 * Arguments: - size_t nshares: number of shares
 *            - uint32_t *in: input buffer
 *            - size_t in_stride: in buffer stride
 * Returns:   - uint32_t out: masked value
 **************************************************/
uint32_t unmask_boolean(size_t nshares, const uint32_t *in, size_t in_stride) {
  uint32_t out, d;
  out = 0;
  for (d = 0; d < nshares; d++) {
    out ^= in[d * in_stride];
  }
  return out;
}

/*************************************************
 * Name:        copy_sharing
 *
 * Description: Copy input sharing to output sharing 
 *
 * Arguments: - size_t nshares: number of shares
 *            - uint32_t *out: output buffer
 *            - size_t out_stride: out buffer stride
 *            - uint32_t *in: input buffer
 *            - size_t in_stride: in buffer stride
 **************************************************/
void copy_sharing(size_t nshares, uint32_t *out, size_t out_stride,
                  const uint32_t *in, size_t in_stride) {
  for (size_t i = 0; i < nshares; i++) {
    out[i * out_stride] = in[i * in_stride];
  }
}

/*************************************************
 * Name:        secadd
 *
 * Description: Performs masked addition of two bitslice words. 
 *              out = (in1 + in2)%(2**kbits_out) 
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t kbits: number of bits in input words
 *            - size_t kbits_out: number of bits in the output word.
 *                kbits = kbits_out or kbits = kbits_out - 1
 *            - uint32_t *out: output buffer
 *            - size_t out_msk_stride: stride between shares
 *            - size_t out_data_stride: stride between masked bits
 *            - uint32_t *in1: first input buffer
 *            - size_t in1_msk_stride: stride between shares
 *            - size_t in1_data_stride: stride between masked bits
 *            - uint32_t *in2: second input buffer
 *            - size_t in2_msk_stride: stride between shares
 *            - size_t in2_data_stride: stride between masked bits
 **************************************************/
void secadd(size_t nshares, size_t kbits, size_t kbits_out, uint32_t *out,
            size_t out_msk_stride, size_t out_data_stride, const uint32_t *in1,
            size_t in1_msk_stride, size_t in1_data_stride, const uint32_t *in2,
            size_t in2_msk_stride, size_t in2_data_stride) {

  if (nshares == NSHARES) {
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

  if (nshares == NSHARES) {
    stop_bench(my_secadd);
  }
}

/*************************************************
 * Name:        seca2b
 *
 * Description: Inplace arithmetic to boolean masking conversion:
 *            sum(in_i)%(2**kbits) = XOR(in_i) 
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t kbits: number of bits in input words. Arithmetic
 *            masking is on 2**kbits.
 *            - uint32_t *in: input buffer
 *            - size_t in_msk_stride: stride between shares
 *            - size_t in_data_stride: stride between masked bits
 **************************************************/
void seca2b(size_t nshares, size_t kbits, uint32_t *in, size_t in_msk_stride,
            size_t in_data_stride) {

  if (nshares == NSHARES)
    start_bench(my_seca2b);

  size_t i, d;

  if (nshares == 1) {
    return;
  }

  size_t nshares_low = nshares / 2;
  size_t nshares_high = nshares - nshares_low;

  seca2b(nshares_low, kbits, in, in_msk_stride, in_data_stride);
  seca2b(nshares_high, kbits, &in[nshares_low * in_msk_stride], in_msk_stride,
         in_data_stride);

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

  if (nshares == NSHARES)
    stop_bench(my_seca2b);
}


/*************************************************
 * Name:        secb2a
 *
 * Description: Inplace boolean to arithmetic masking conversion:
 *            sum(in_i)%((1<<kbits)) = XOR(in_i) 
 *
 * /!\ This function is operating on two slices such that 64 conversions are
 * performed in parallel.
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t kbits: number of bits in input words. kbits < 16
 *            - uint32_t p: modulus of the arithmetic masking
 *            - uint32_t *in: input buffer
 *            - size_t in_msk_stride: stride between shares
 *            - size_t in_data_stride: stride between masked bits
 **************************************************/
void secb2a(size_t nshares,
            size_t kbits,
            uint32_t *in, size_t in_msk_stride, size_t in_data_stride) {

  uint16_t z_dense[2 * BSSIZE * nshares];
  uint16_t zp_dense[2 * BSSIZE * nshares];
  uint32_t zp_str[2 * kbits * nshares];
  uint32_t b_str[2 * kbits * nshares];
  size_t d, i;
  uint32_t r;
  uint16_t r0, r1;
  uint32_t p = 1 << kbits;
  // generate uniform sharing for z
  // zp = p - z;
  for (d = 0; d < nshares - 1; d++) {
    for (i = 0; i < BSSIZE; i += 2) {
      r = rand32();
      r0 = r & ((1 << kbits) - 1);
      r1 = (r >> 16) & ((1 << kbits) - 1);
      z_dense[i * nshares + d] = r0;
      zp_dense[i * nshares + d] = p - r0;

      z_dense[(i + 1) * nshares + d] = r1;
      zp_dense[(i + 1) * nshares + d] = p - r1;

      r = rand32();
      r0 = r & ((1 << kbits) - 1);
      r1 = (r >> 16) & ((1 << kbits) - 1);
      z_dense[i * nshares + d + (BSSIZE * nshares)] = r0;
      zp_dense[i * nshares + d + (BSSIZE * nshares)] = p - r0;

      z_dense[(i + 1) * nshares + d + (BSSIZE * nshares)] = r1;
      zp_dense[(i + 1) * nshares + d + (BSSIZE * nshares)] = p - r1;
    }
#if ((BSSIZE & 0x1))
#error "Can only handle even number of slices."
#endif
  }

  // map zp to bitslice representation
  masked_dense2bitslice_opt(nshares - 1, kbits, zp_str, 1, nshares, zp_dense, 1,
                            nshares);

  // last shares of zp set to zero
  for (i = 0; i < kbits; i++) {
    zp_str[i * nshares + (nshares - 1)] = 0;
    zp_str[i * nshares + (nshares - 1) + kbits * nshares] = 0;
  }

  // last shares of zp_str to zero
  seca2b(nshares, kbits, zp_str, 1, nshares);
  seca2b(nshares, kbits, &zp_str[kbits * nshares], 1, nshares);

  secadd(nshares, kbits, kbits, b_str, 1, nshares, in, in_msk_stride,
         in_data_stride, zp_str, 1, nshares);
  secadd(nshares, kbits, kbits, &b_str[kbits * nshares], 1, nshares,
         &in[kbits * nshares], in_msk_stride, in_data_stride,
         &zp_str[kbits * nshares], 1, nshares);

  // map z to bistlice in output buffer
  masked_dense2bitslice_opt(nshares - 1, kbits, in, in_msk_stride,
                            in_data_stride, z_dense, 1, nshares);

  // unmask b_str and set to the last share of the output
  for (i = 0; i < kbits; i++) {
    RefreshIOS_rec(nshares,nshares,&b_str[i * nshares], 1);
    RefreshIOS_rec(nshares,nshares,&b_str[i * nshares + kbits * nshares],1);

    in[i * in_data_stride + (nshares - 1) * in_msk_stride] = 0;
    in[i * in_data_stride + (nshares - 1) * in_msk_stride + kbits * nshares] =
        0;
    for (d = 0; d < nshares; d++) {
      in[i * in_data_stride + (nshares - 1) * in_msk_stride] ^=
          b_str[i * nshares + d];
      in[i * in_data_stride + (nshares - 1) * in_msk_stride +
         kbits * nshares] ^= b_str[i * nshares + d + kbits * nshares];
    }
  }
}

/*************************************************
 * Name:        RefreshIOS_rec
 *
 * Description: IOS refresh on boolean sharing 
 * 
 * Arguments: - size_t nshares: number of shares
 *            - size_t d: current recursion shars:
 *            - uint32_t *x: input buffer
 *            - size_t x_msk_stride: stride between shares
 **************************************************/
void RefreshIOS_rec(size_t nshares, size_t d,
    uint32_t *x, size_t x_msk_stride) {
  uint32_t r;
  if (d == 1) {
  } else if (d == 2) {
    r = get_random();
    x[0*x_msk_stride] ^= r;
    x[1*x_msk_stride] ^= r;
  } else {
    RefreshIOS_rec(nshares, d / 2, x, x_msk_stride);
    RefreshIOS_rec(nshares, d - d / 2, &x[(d / 2)*x_msk_stride], x_msk_stride);
    for (unsigned int i = 0; i < d / 2; i++) {
      r = rand32();
      x[i * x_msk_stride] ^= r;
      x[(i + d / 2) * x_msk_stride] ^= r;
    }
  }
}

/*************************************************
 * Name:        secfulladd
 *
 * Description: Performs a 3-bit full adder. See for details: 
 *              https://eprint.iacr.org/2022/158.pdf, Algo 5.
 *
 * Arguments: - size_t nshares: number of shares
 *            - uint32_t *co: carry out buffer
 *            - size_t co_msk_stride: carry out shares stride
 *            - uint32_t *w: output bit
 *            - size_t w_msk_stride: w shares stride
 *            - uint32_t *ci: carry in buffer
 *            - size_t ci_msk_stride: carry in shares stride
 *            - uint32_t *x: first input bit
 *            - size_t x_msk_stride: x shares stride
 *            - uint32_t *y: second input bit
 *            - size_t y_msk_stride: y shares stride
 **************************************************/
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

/*************************************************
 * Name:        masked_cbd
 *
 * Description: Computes central binomial distribution.
 *              z = (HW(a) - HW(b))%(1<<kbits). Output is an arithmetic masking
 *              https://eprint.iacr.org/2022/158.pdf, Algo 6.
 *
 * /!\ operates on 64 coefficients at the time.
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t eta: number of bits in a and b
 *            - size_t kbits: number of bits in modulus
 *            - int16_t *z: output central binomial distribution
 *            - size_t z_msk_stride: z shares stride
 *            - size_t z_data_stride: z data stride
 *            - uint32_t *a: input uniform bits 
 *            - size_t a_msk_stride: a shares stride
 *            - size_t a_data_stride: a data stride
 *            - uint32_t *b: input uniform bits 
 *            - size_t b_msk_stride: b shares stride
 *            - size_t b_data_stride: b data stride
 **************************************************/
void masked_cbd(size_t nshares, size_t eta, size_t kbits, uint16_t *z,
                size_t z_msk_stride, size_t z_data_stride, uint32_t *a,
                size_t a_msk_stride, size_t a_data_stride, uint32_t *b,
                size_t b_msk_stride, size_t b_data_stride) {

  start_bench(my_cbd);
  size_t np = 2 * eta;
  size_t i, d, k, j, s;
  uint32_t sp[nshares * np], z_str_full[nshares * kbits * 2];
  uint32_t *a_in, *b_in, *z_str;

  //  k = ceil(log2(2*eta +1))
  k = 0;
  i = 2 * eta + 1;
  while (i != 0) {
    k++;
    i >>= 1;
  }
  uint32_t z_32[NSHARES*BSSIZE];
  for (s = 0; s < 2; s++) {
    sec_sampler2(z_32,&a[(s*eta)*a_data_stride],&b[(s*eta)*b_data_stride],4,SABER_EQ);
    for(i = s*32; i < ((s*32) + 32); i++){
      for( d = 0; d < NSHARES; d ++){
        z[ i * z_data_stride + d * z_msk_stride] = z_32[ (i%32) * NSHARES + d]; 
      }
    }
  }
  stop_bench(my_cbd);
}
