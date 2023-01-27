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
void masked_and_c(size_t nshares, uint32_t *z, size_t z_stride,
                  const uint32_t *a, size_t a_stride, const uint32_t *b,
                  size_t b_stride) {
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
void copy_sharing_c(size_t nshares, uint32_t *out, size_t out_stride,
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

    masked_xor(nshares, xpy, 1, &in1[i * in1_data_stride], in1_msk_stride,
               &in2[i * in2_data_stride], in2_msk_stride);
    masked_xor(nshares, xpc, 1, &in1[i * in1_data_stride], in1_msk_stride,
               carry, 1);
    masked_xor(nshares, &out[i * out_data_stride], out_msk_stride, xpy, 1,
               carry, 1);

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
 * Name:        secadd_constant_bmsk
 *
 * Description: Performs masked addition of an input with a masked
 *              constant sucht that:
 *              out = (in + bmsk*constant)%(2**kbits_out)
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t kbits: number of bits in input words
 *            - size_t kbits_out: number of bits in the output word.
 *                kbits = kbits_out or kbits = kbits_out - 1
 *            - uint32_t *out: output buffer
 *            - size_t out_msk_stride: stride between shares
 *            - size_t out_data_stride: stride between masked bits
 *            - uint32_t *in: first input buffer
 *            - size_t in_msk_stride: stride between shares
 *            - size_t in_data_stride: stride between masked bits
 *            - uint32_t constant: public constant
 *            - uint32_t *bmsk: masked bit buffer
 *            - size_t bmsk_msk_stride: shares stride of masked bit
 **************************************************/
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
      masked_xor(nshares, xpy, 1, &in1[i * in1_data_stride], in1_msk_stride,
                 bmsk, bmsk_msk_stride);
      masked_xor(nshares, xpc, 1, &in1[i * in1_data_stride], in1_msk_stride,
                 carry, 1);
      masked_xor(nshares, &out[i * out_data_stride], out_msk_stride, xpy, 1,
                 carry, 1);

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

/*************************************************
 * Name:        secadd_constant
 *
 * Description: Performs masked addition of an input with a
 *              public constant such that
 *              out = (in + constant)%(2**kbits_out)
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t kbits: number of bits in input words
 *            - size_t kbits_out: number of bits in the output word.
 *                kbits = kbits_out or kbits = kbits_out - 1
 *            - uint32_t *out: output buffer
 *            - size_t out_msk_stride: stride between shares
 *            - size_t out_data_stride: stride between masked bits
 *            - uint32_t *in: first input buffer
 *            - size_t in_msk_stride: stride between shares
 *            - size_t in_data_stride: stride between masked bits
 *            - uint32_t constant: public constant
 **************************************************/
void secadd_constant(size_t nshares, size_t kbits, size_t kbits_out,
                     uint32_t *out, size_t out_msk_stride,
                     size_t out_data_stride, const uint32_t *in1,
                     size_t in1_msk_stride, size_t in1_data_stride,
                     uint32_t constant) {

  size_t i, d;
  uint32_t carry[nshares];
  uint32_t xpy[nshares];
  uint32_t xpc[nshares];
  uint32_t dummy = 0;

  if (constant & 0x1) {
    copy_sharing(nshares, carry, 1, in1, in1_msk_stride);
    copy_sharing(nshares, out, out_msk_stride, in1, in1_msk_stride);
    secxor_cst(out, 0xFFFFFFFF, &dummy);
  } else {
    for (d = 0; d < nshares; d++) {
      carry[d] = 0;
    }
    copy_sharing(nshares, out, out_msk_stride, in1, in1_msk_stride);
  }

  for (i = 1; i < kbits; i++) {
    if ((constant >> i) & 0x1) {
      copy_sharing(nshares, xpy, 1, &in1[i * in1_data_stride], in1_msk_stride);
      masked_xor(nshares, xpc, 1, &in1[i * in1_data_stride], in1_msk_stride,
                 carry, 1);
      masked_xor(nshares, &out[i * out_data_stride], out_msk_stride, xpy, 1,
                 carry, 1);

      secxor_cst(xpy, 0xFFFFFFFF, &dummy);
      secxor_cst(out + i * out_data_stride, 0xFFFFFFFF, &dummy);

      if ((i == (kbits - 1)) && (i == (kbits_out - 1))) {
        return;
      } else if (i == (kbits - 1)) {
        masked_and(nshares, carry, 1, xpy, 1, xpc, 1);
        masked_xor(nshares, &out[(kbits)*out_data_stride], out_msk_stride,
                   carry, 1, &in1[i * in1_data_stride], in1_msk_stride);

        // add the kbits_out of the constant
        secxor_cst(out + kbits*out_data_stride, 0xFFFFFFFF * ((constant >> kbits) & 0x1), &dummy);

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
        secxor_cst(out + kbits*out_data_stride, 0xFFFFFFFF * ((constant >> kbits) & 0x1), &dummy);
        return;
      }
      masked_and(nshares, carry, 1, carry, 1, &in1[i * in1_data_stride],
                 in1_msk_stride);
    }
  }
}

/*************************************************
 * Name:        secadd_modp
 *
 * Description: Performs masked addition of two bitslice words with
 *              arbitrary modulo such that
 *              out = (in1 + in2)%(p)
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t kbits: number of bits in input words:
 *                requires kbits = ceil(log(p))
 *            - uint32_t: modulo for the reduction
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
void secadd_modp(size_t nshares, size_t kbits, uint32_t q, uint32_t *out,
                 size_t out_msk_stride, size_t out_data_stride,
                 const uint32_t *in1, size_t in1_msk_stride,
                 size_t in1_data_stride, const uint32_t *in2,
                 size_t in2_msk_stride, size_t in2_data_stride) {

  start_bench(my_secaddmodp_12);
  uint32_t s[(kbits + 1) * nshares];
  uint32_t sp[(kbits + 1) * nshares];

  start_bench(my_secadd_12);
  secadd(nshares, kbits, kbits + 1, s, 1, nshares, in1, in1_msk_stride,
         in1_data_stride, in2, in2_msk_stride, in2_data_stride);
  stop_bench(my_secadd_12);

  secadd_constant(nshares, kbits + 1, kbits + 1, sp, 1, nshares, s, 1, nshares,
                  (1 << (kbits + 1)) - q);

  secadd_constant_bmsk(nshares, kbits, kbits, out, out_msk_stride,
                       out_data_stride, sp, 1, nshares, q, &sp[kbits * nshares],
                       1);

  stop_bench(my_secaddmodp_12);
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
    copy_sharing(nshares_low, &expanded_low[i * nshares], 1,
                 &in[i * in_data_stride], in_msk_stride);
    copy_sharing(nshares_high, &expanded_high[i * nshares + nshares_low], 1,
                 &in[i * in_data_stride + nshares_low * in_msk_stride],
                 in_msk_stride);

    for (d = 0; d < nshares_low; d++) {
      expanded_high[i * nshares + d] = 0;
    }
    for (d = nshares_low; d < nshares; d++) {
      expanded_low[i * nshares + d] = 0;
    }
  }

  secadd(nshares, kbits, kbits, in, in_msk_stride, in_data_stride, expanded_low,
         1, nshares, expanded_high, 1, nshares);

  if (nshares == NSHARES)
    stop_bench(my_seca2b);
}

/*************************************************
 * Name:        seca2b_modp
 *
 * Description: Inplace arithmetic to boolean masking conversion:
 *            sum(in_i)%(p) = XOR(in_i)
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t kbits: number of bits in input words. kbits =
 *ceil(log(p))
 *            - uint32_t p: modulus of the arithmetic masking
 *            - uint32_t *in: input buffer
 *            - size_t in_msk_stride: stride between shares
 *            - size_t in_data_stride: stride between masked bits
 **************************************************/
void seca2b_modp(size_t nshares, size_t kbits, uint32_t p, uint32_t *in,
                 size_t in_msk_stride, size_t in_data_stride) {

  size_t i, d;
  if (nshares == NSHARES)
    start_bench(my_seca2b_modp);

  if (nshares == 1) {
    return;
  }

  size_t nshares_low = nshares / 2;
  size_t nshares_high = nshares - nshares_low;

  seca2b_modp(nshares_low, kbits, p, in, in_msk_stride, in_data_stride);
  seca2b_modp(nshares_high, kbits, p, &in[nshares_low * in_msk_stride],
              in_msk_stride, in_data_stride);

  uint32_t expanded_low[(kbits + 1) * nshares];
  uint32_t expanded_high[(kbits + 1) * nshares];
  uint32_t u[(kbits + 1) * nshares];

  secadd_constant(nshares_low, kbits, kbits + 1, expanded_low, 1, nshares, in,
                  in_msk_stride, in_data_stride, (1 << (kbits + 1)) - p);

  for (i = 0; i < (kbits + 1); i++) {
    if (i < kbits) {
      copy_sharing(nshares_high, &expanded_high[i * nshares + nshares_low], 1,
                   &in[i * in_data_stride + nshares_low * in_msk_stride],
                   in_msk_stride);
    }
    for (d = 0; d < nshares_low; d++) {
      // has already been written by secadd_constant_bmsk
      // expanded_low[i*nshares + d] = in[i*in_data_stride + d*in_msk_stride];
      expanded_high[i * nshares + d] = 0;
    }
    for (d = nshares_low; d < nshares; d++) {
      // kbits + 1 within in is unset
      if (i >= kbits) {
        expanded_high[i * nshares + d] = 0;
      }
      expanded_low[i * nshares + d] = 0;
    }
  }

  secadd(nshares, kbits + 1, kbits + 1, u, 1, nshares, expanded_high, 1,
         nshares, expanded_low, 1, nshares);

  secadd_constant_bmsk(nshares, kbits, kbits, in, in_msk_stride, in_data_stride,
                       u, 1, nshares, p, &u[kbits * nshares], 1);

  if (nshares == NSHARES)
    stop_bench(my_seca2b_modp);
}

/*************************************************
 * Name:        secb2a_modp
 *
 * Description: Inplace boolean to arithmetic masking conversion:
 *            sum(in_i)%(p) = XOR(in_i)
 *
 * /!\ This function is operating on two slices such that 64 conversions are
 * performed in parallel.
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t kbits: number of bits in input words. kbits =
 *ceil(log(p))
 *            - uint32_t p: modulus of the arithmetic masking
 *            - uint32_t *in: input buffer
 *            - size_t in_msk_stride: stride between shares
 *            - size_t in_data_stride: stride between masked bits
 **************************************************/
void secb2a_modp(size_t nshares,
                 //   size_t kbits, // MUST BE EQUAL TO COEF_NBITS
                 uint32_t p, uint32_t *in, size_t in_msk_stride,
                 size_t in_data_stride) {

  int16_t z_dense[2 * BSSIZE * nshares];
  int16_t zp_dense[2 * BSSIZE * nshares];
  uint32_t zp_str[2 * COEF_NBITS * nshares];
  uint32_t b_str[2 * COEF_NBITS * nshares];
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
      z_dense[i * nshares + d + (BSSIZE * nshares)] = r[0];
      zp_dense[i * nshares + d + (BSSIZE * nshares)] = p - r[0];

      z_dense[(i + 1) * nshares + d + (BSSIZE * nshares)] = r[1];
      zp_dense[(i + 1) * nshares + d + (BSSIZE * nshares)] = p - r[1];
    }
#if ((BSSIZE & 0x1) == 0x1)
    rand_q(r);
    z_dense[(BSSIZE - 1) * nshares + d] = r[0];
    zp_dense[(BSSIZE - 1) * nshares + d] = p - r[0];

    rand_q(r);
    z_dense[(BSSIZE - 1) * nshares + d + (BSSIZE * nshares)] = r[0];
    zp_dense[(BSSIZE - 1) * nshares + d + (BSSIZE * nshares)] = p - r[0];
#endif
  }

  // map zp to bitslice representation
  masked_dense2bitslice_opt(nshares - 1, COEF_NBITS, zp_str, 1, nshares,
                            zp_dense, 1, nshares);

  // last shares of zp set to zero
  for (i = 0; i < COEF_NBITS; i++) {
    zp_str[i * nshares + (nshares - 1)] = 0;
    zp_str[i * nshares + (nshares - 1) + COEF_NBITS * nshares] = 0;
  }

  // last shares of zp_str to zero
  seca2b_modp(nshares, COEF_NBITS, p, zp_str, 1, nshares);
  seca2b_modp(nshares, COEF_NBITS, p, &zp_str[COEF_NBITS * nshares], 1,
              nshares);

  secadd_modp(nshares, COEF_NBITS, p, b_str, 1, nshares, in, in_msk_stride,
              in_data_stride, zp_str, 1, nshares);
  secadd_modp(nshares, COEF_NBITS, p, &b_str[COEF_NBITS * nshares], 1, nshares,
              &in[COEF_NBITS * nshares], in_msk_stride, in_data_stride,
              &zp_str[COEF_NBITS * nshares], 1, nshares);

  // map z to bistlice in output buffer
  masked_dense2bitslice_opt(nshares - 1, COEF_NBITS, in, in_msk_stride,
                            in_data_stride, z_dense, 1, nshares);

  // unmask b_str and set to the last share of the output
  for (i = 0; i < COEF_NBITS; i++) {
    RefreshIOS_rec(nshares, nshares, &b_str[i * nshares], 1);
    RefreshIOS_rec(nshares, nshares, &b_str[i * nshares + COEF_NBITS * nshares],
                   1);

    in[i * in_data_stride + (nshares - 1) * in_msk_stride] = 0;
    in[i * in_data_stride + (nshares - 1) * in_msk_stride +
       COEF_NBITS * nshares] = 0;
    for (d = 0; d < nshares; d++) {
      in[i * in_data_stride + (nshares - 1) * in_msk_stride] ^=
          b_str[i * nshares + d];
      in[i * in_data_stride + (nshares - 1) * in_msk_stride +
         COEF_NBITS * nshares] ^= b_str[i * nshares + d + COEF_NBITS * nshares];
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
void RefreshIOS_rec(size_t nshares, size_t d, uint32_t *x,
                    size_t x_msk_stride) {
  uint32_t r;
  if (d == 1) {
  } else if (d == 2) {
    r = get_random();
    x[0 * x_msk_stride] ^= r;
    x[1 * x_msk_stride] ^= r;
  } else {
    RefreshIOS_rec(nshares, d / 2, x, x_msk_stride);
    RefreshIOS_rec(nshares, d - d / 2, &x[(d / 2) * x_msk_stride],
                   x_msk_stride);
    for (unsigned int i = 0; i < d / 2; i++) {
      r = rand32();
      x[i * x_msk_stride] ^= r;
      x[(i + d / 2) * x_msk_stride] ^= r;
    }
  }
}

/*************************************************
 * Name:        seccompress
 *
 * Description: Performs polynomial coefficients. See for details:
 *              https://eprint.iacr.org/2021/1615.pdf, Algo 6.
 *
 * /!\ performs compression on 32 coefficients at once.
 *
 * Arguments: - size_t nshares: number of shares
 *            - uint32_t q: prime
 *            - uint32_t c: compression factor
 *            - uint32_t *out: output bitslice buffer. c-bits buffer.
 *            - size_t out_msk_stride: stride between shares
 *            - size_t out_data_stride: stride between masked bits
 *            - int16_t *in: input buffer. Contains 32 coefficients
 *            - size_t in_msk_stride: stride between shares
 *            - size_t in_data_stride: stride between masked bits
 **************************************************/
void seccompress(size_t nshares, uint32_t q, uint32_t c, uint32_t *out,
                 size_t out_msk_stride, size_t out_data_stride,
                 const int16_t *in, size_t in_msk_stride,
                 size_t in_data_stride) {

  size_t i, d;
  uint32_t ell = 0;
  uint32_t prod = q * nshares;
  while (prod > 0) {
    prod = prod >> 1;
    ell++;
  }

  uint32_t in_expanded[BSSIZE * nshares];
  uint32_t bs_expanded[(ell + c) * nshares];

  // map mod q to mod 2^ell.
  uint32_t tmp32;
  uint64_t tmp64;
  for (i = 0; i < BSSIZE; i++) {
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
  masked_dense2bitslice_opt_u32(nshares, ell + c, bs_expanded, 1, nshares,
                                in_expanded, 1, nshares);

  // convert A2B
  seca2b(nshares, ell + c, bs_expanded, 1, nshares);

  // map to the output
  for (i = 0; i < c; i++) {
    for (d = 0; d < nshares; d++) {
      out[i * out_data_stride + d * out_msk_stride] =
          bs_expanded[(ell + i) * nshares + d];
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
 *              z = (HW(a) - HW(b))%p. Output is an arithmetic masking
 *              https://eprint.iacr.org/2022/158.pdf, Algo 6.
 *
 * /!\ operates on 64 coefficients at the time.
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t eta: number of bits in a and b
 *            - size_t p: modulus of the
 *            - size_t kbits: number of bits in p ceil(log2(p))
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
void masked_cbd(size_t nshares, size_t eta, size_t p, size_t kbits, int16_t *z,
                size_t z_msk_stride, size_t z_data_stride, uint32_t *a,
                size_t a_msk_stride, size_t a_data_stride, uint32_t *b,
                size_t b_msk_stride, size_t b_data_stride) {

  start_bench(my_cbd);
  size_t np = 2 * eta;
  size_t i, d, k, j, s;
  uint32_t sp[nshares * np], z_str_full[nshares * COEF_NBITS * 2];
  uint32_t *a_in, *b_in, *z_str;

  //  k = ceil(log2(2*eta +1))
  k = 0;
  i = 2 * eta + 1;
  while (i != 0) {
    k++;
    i >>= 1;
  }

  start_bench(cbd_bool);
  // compte HW(a)-HW(b) for all 64 input coefficients
  // levaraging 32 bus size
  for (s = 0; s < 2; s++) {
    // copy input puts
    a_in = &a[s * (eta)*a_data_stride];
    b_in = &b[s * (eta)*b_data_stride];
    z_str = &z_str_full[s * nshares * COEF_NBITS];

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
    np = 2 * eta;
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

  stop_bench(cbd_bool);

  start_bench(cbd_b2a);
  secb2a_modp(nshares, p, z_str_full, 1, nshares);

  masked_bitslice2dense_opt(nshares, kbits, z, z_msk_stride, z_data_stride,
                            z_str_full, 1, nshares);

  stop_bench(cbd_b2a);
  for (i = 0; i < 2 * BSSIZE; i++) {
    z[i * nshares] = (z[i * nshares] + p - eta) % p;
  }

  stop_bench(my_cbd);
}

// algo 7 in https://eprint.iacr.org/2019/910
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

// algo 6 in https://eprint.iacr.org/2019/910
static void b2a_qbit(size_t nshares, int16_t *a, size_t a_msk_stride,
                     uint32_t *x, size_t x_msk_stride) {

  a[0] = x[0];
  for (size_t i = 1; i < nshares; i++) {
    secb2a_qbit_n(i + 1, a, a_msk_stride, a, a_msk_stride, x[i * x_msk_stride]);
  }
}

// algo 8 in https://eprint.iacr.org/2019/910
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

/*************************************************
 * Name:        secb2a_1bit
 *
 * Description: Performs boolean to arithmetic conversion
 *            of a single masked bit.
 *            See https://eprint.iacr.org/2019/910, algo 5
 *
 * /!\ operates on 64 coefficients at the time.
 *
 * Arguments: - size_t nshares: number of shares
 *            - int16_t *a: output masked bit (with arith. masking)
 *            - size_t a_msk_stride: a shares stride
 *            - uint32_t *b: input masked bit (bool. masking)
 *            - size_t b_msk_stride: b masked stride
 **************************************************/
void secb2a_1bit(size_t nshares, int16_t *a, size_t a_msk_stride, uint32_t *x,
                 size_t x_msk_stride) {
  b2a_qbit(nshares, a, a_msk_stride, x, x_msk_stride);
  refresh_add(nshares, a, a_msk_stride);
}
