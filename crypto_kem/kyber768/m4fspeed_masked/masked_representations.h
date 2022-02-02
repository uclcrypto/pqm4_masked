#ifndef MASKED_REPRESENTATION_H
#define MASKED_REPRESENTATION_H

#include "masked.h"
#include <stddef.h>
void APoly2StrAPoly(StrAPoly out, const APoly in);
void StrAPoly2APoly(APoly out, const StrAPoly in);

void masked_bitslice2dense(
        int16_t *dense[],
        uint32_t *bitslice[],
        size_t coeffs_size,
        size_t n_coeffs,
        size_t nshares);

void masked_dense2bitslice(
        uint32_t *bitslice[],
        int16_t *dense[],
        size_t coeffs_size,
        size_t n_coeffs,
        size_t nshares);

#endif
