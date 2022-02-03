#ifndef MASKED_REPRESENTATION_H
#define MASKED_REPRESENTATION_H

#include "masked.h"
#include <stddef.h>
void APoly2StrAPoly(StrAPoly out, const APoly in);
void StrAPoly2APoly(APoly out, const StrAPoly in);

void masked_dense2bitslice(
        size_t nshares,
        size_t n_coeffs,
        size_t coeffs_size,
        uint32_t *bitslice,size_t bitslice_data_stride,size_t bitlice_msk_stride,
        int16_t *dense,size_t dense_data_stride,size_t dense_msk_stide);

void masked_bitslice2dense(
        size_t nshares,
        size_t n_coeffs,
        size_t coeffs_size,
        int16_t *dense,size_t dense_data_stride,size_t dense_msk_stide,
        uint32_t *bitslice,size_t bitslice_data_stride,size_t bitlice_msk_stride);

#endif
