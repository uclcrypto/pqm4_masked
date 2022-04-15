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

#include "SABER_params.h"
#include "masked.h"
#include <stdint.h>

uint8_t
masked_indcpa_kem_enc_cmp(const uint8_t *m, size_t m_msk_stride,
                          size_t m_data_stride, const uint8_t *seed_s,
                          size_t seed_s_msk_stride, size_t seed_s_data_stride,
                          const uint8_t pk[SABER_INDCPA_PUBLICKEYBYTES],
                          const uint8_t ciphertext[SABER_BYTES_CCA_DEC]);

void masked_indcpa_kem_dec(StrAPolyVec masked_sk,
                           const uint8_t ciphertext[SABER_BYTES_CCA_DEC],
                           uint8_t *m, size_t m_msk_stide,
                           size_t m_data_stride);
#endif
