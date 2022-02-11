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
#ifndef MASKED_FIPS202_H
#define MASKED_FIPS202_H

#include <stddef.h>
#include <stdint.h>

void masked_shake256(uint8_t *output, size_t outlen, size_t out_msk_stride,
                     size_t out_data_stride, const uint8_t *input, size_t inlen,
                     size_t in_msk_stride, size_t in_data_stride);

void masked_sha3_512(uint8_t *output, size_t out_msk_stride,
                     size_t out_data_stride, const uint8_t *input, size_t inlen,
                     size_t in_msk_stride, size_t in_data_stride);

#endif
