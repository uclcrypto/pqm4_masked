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
#ifndef MASKED_SYM_FIPS202
#define MASKED_SYM_FIPS202

void masked_shake256_prf(unsigned char *output, size_t outlen,
                         size_t o_msk_stride, size_t o_data_stride,
                         const unsigned char *key, size_t k_msk_stride,
                         size_t k_data_stride, unsigned char nonce);

#endif
