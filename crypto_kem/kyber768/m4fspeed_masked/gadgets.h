#ifndef GADGETS_H
#define GADGETS_H

#include <stdint.h>
#include <stddef.h>
#include "masked.h"

// Atomic gadgets
void masked_xor(
        size_t nshares,
        uint32_t *out, size_t out_stride,
        const uint32_t *ina, size_t ina_stride,
        const uint32_t *inb, size_t inb_stride
        );
void masked_and(
        size_t nshares,
        uint32_t *z, size_t z_stride,
        const uint32_t *a, size_t a_stride,
        const uint32_t *b, size_t b_stride
        );
void copy_sharing(
        size_t nshares,
        uint32_t *out, size_t out_stride,
        const uint32_t *in, size_t in_stride
        );
void RefreshIOS_rec(uint32_t *x, uint32_t d);

// Adders
void secfulladd(size_t nshares,
        uint32_t *co, size_t co_msk_stride,
        uint32_t *w, size_t w_msk_stide,
        uint32_t *ci, size_t ci_msk_stride,
        uint32_t *x, size_t x_msk_stride,
        uint32_t *y, size_t y_msk_stride);
void secadd(size_t nshares,
        size_t kbits,size_t kbits_out,
        uint32_t *out, size_t out_msk_stride, size_t out_data_stride,
        const uint32_t *in1, size_t in1_msk_stride, size_t in1_data_stride,
        const uint32_t *in2, size_t in2_msk_stride, size_t in2_data_stride);
void secadd_modp(size_t nshares,
        size_t kbits,
        uint32_t p,
        uint32_t *out, size_t out_msk_stride, size_t out_data_stride,
        const uint32_t *in1, size_t in1_msk_stride, size_t in1_data_stride,
        const uint32_t *in2, size_t in2_msk_stride, size_t in2_data_stride);
void secadd_constant_bmsk(size_t nshares,
        size_t kbits,
        size_t kbits_out,
        uint32_t *out, size_t out_msk_stride, size_t out_data_stride,
        const uint32_t *in1, size_t in1_msk_stride, size_t in1_data_stride,
        uint32_t constant, const uint32_t *bmsk, size_t bmsk_msk_stride);
void secadd_constant(size_t nshares,
        size_t kbits,
        size_t kbits_out,
        uint32_t *out, size_t out_msk_stride, size_t out_data_stride,
        const uint32_t *in1, size_t in1_msk_stride, size_t in1_data_stride,
        uint32_t constant);

// Conversions
void seca2b(size_t nshares,
        size_t kbits,
        uint32_t *in, size_t in_msk_stride, size_t in_data_stride);
void seca2b_modp(size_t nshares,
        size_t kbits,
        uint32_t p,
        uint32_t *in, size_t in_msk_stride, size_t in_data_stride);
void secb2a_1bit(
        size_t nshares,
        int16_t *a, size_t a_msk_stride,
        uint32_t *x, size_t x_msk_stride);
void secb2a_modp(size_t nshares,
        uint32_t p,
        uint32_t *in, size_t in_msk_stride, size_t in_data_stride);

// Lattice-based KEM specific
void seccompress(size_t nshares,
        size_t ncoeffs,
        uint32_t q,
        uint32_t c,
        uint32_t *out, size_t out_msk_stride, size_t out_data_stride,
        const int16_t *in, size_t in_msk_stride, size_t in_data_stride);

void masked_cbd(size_t nshares,
        size_t eta,
        size_t n_coeffs,
        size_t p, size_t kbits,
        int16_t *z, size_t z_msk_stride, size_t z_data_stride, 
        uint32_t *a, size_t a_msk_stride, size_t a_data_stride,
        uint32_t *b, size_t b_msk_stride, size_t b_data_stride);

#endif // GADGETS_H
