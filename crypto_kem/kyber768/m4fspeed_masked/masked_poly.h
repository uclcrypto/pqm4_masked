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
#ifndef MASKEDPOLY_H
#define MASKEDPOLY_H
#include "masked.h"
#include "params.h"
#include "poly.h"

void masked_poly_ntt(StrAPoly r);
void masked_poly_invntt(StrAPoly r);
void masked_poly_tomsg(unsigned char *m, StrAPoly str_r);
void masked_poly_cmp(size_t c, uint32_t *rc, const StrAPoly mp,
                     const poly *ref);

void finalize_cmp(uint32_t *bits);

void masked_poly_noise(StrAPoly r, const unsigned char *seed,
                       size_t seed_msk_stride, size_t seed_data_stride,
                       unsigned char nonce, int add);
void masked_poly_frommsg(StrAPoly y,
                         const uint8_t m[KYBER_INDCPA_MSGBYTES * (NSHARES)],
                         size_t m_msk_stride, size_t m_data_stride);
#endif
