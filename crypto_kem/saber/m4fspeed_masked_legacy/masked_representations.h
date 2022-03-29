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
#ifndef MASKED_REPRESENTATION_H
#define MASKED_REPRESENTATION_H

#include "masked.h"
#include <stddef.h>
void APoly2StrAPoly(StrAPoly out, const APoly in);
void StrAPoly2APoly(APoly out, const StrAPoly in);

void masked_dense2bitslice(size_t nshares, size_t n_coeffs, size_t coeffs_size,
                           uint32_t *bitslice, size_t bitslice_msk_stride,
                           size_t bitlice_data_stride, const int16_t *dense,
                           size_t dense_msk_stride, size_t dense_data_stide);

void masked_dense2bitslice_u32(size_t nshares, size_t n_coeffs,
                               size_t coeffs_size, uint32_t *bitslice,
                               size_t bitslice_msk_stride,
                               size_t bitlice_data_stride,
                               const uint32_t *dense, size_t dense_msk_stride,
                               size_t dense_data_stide);

void masked_bitslice2dense(size_t nshares, size_t n_coeffs, size_t coeffs_size,
                           int16_t *dense, size_t dense_msk_stride,
                           size_t dense_data_stide, const uint32_t *bitslice,
                           size_t bitslice_msk_stride,
                           size_t bitlice_data_stride);

void masked_dense2bitslice_opt(size_t nshares, size_t coeffs_size,
                               uint32_t *bitslice, size_t bitslice_msk_stride,
                               size_t bitslice_data_stride,
                               const uint16_t *dense, size_t dense_msk_stride,
                               size_t dense_data_stide);
void masked_dense2bitslice_opt_u32(
    size_t nshares, size_t coeffs_size, uint32_t *bitslice,
    size_t bitslice_msk_stride, size_t bitslice_data_stride,
    const uint32_t *dense, size_t dense_msk_stride, size_t dense_data_stide);
void masked_bitslice2dense_opt(size_t nshares, size_t coeffs_size,
                               uint16_t *dense, size_t dense_msk_stride,
                               size_t dense_data_stide,
                               const uint32_t *bitslice,
                               size_t bitslice_msk_stride,
                               size_t bitslice_data_stride);
#endif
