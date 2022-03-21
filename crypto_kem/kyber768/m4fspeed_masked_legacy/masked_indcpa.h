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
#ifndef MASKED_INDCPA_H
#define MASKED_INDCPA_H

#include <stddef.h>

unsigned char
masked_indcpa_enc_cmp(const unsigned char *c, const unsigned char *m,
                      size_t m_msk_stride, size_t m_data_stride,
                      const unsigned char *pk, const unsigned char *coins,
                      size_t coins_msk_stride, size_t coins_data_stride);

void masked_indcpa_dec(unsigned char *m, size_t o_msk_stride,
                       size_t o_data_stride, const unsigned char *c,
                       const unsigned char *sk);
#endif
