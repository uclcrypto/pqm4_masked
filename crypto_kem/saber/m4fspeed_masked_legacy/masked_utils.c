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
#include "masked_utils.h"
#include "SABER_params.h"
#include "masked.h"
#include <libopencm3/stm32/rng.h>

extern inline uint32_t get_random();
extern inline uint32_t rand32();
extern inline void rand_q(uint16_t v[2], uint16_t modulus);
extern inline uint16_t rand16();
uint16_t sampleq_id=0;
uint16_t sampleq_rs[2];
uint16_t r16_id = 0;
uint32_t r16_rs;

/*************************************************
 * Name:        masked_poly
 *
 * Description: masked a polynomial into a strided polynomial
 *
 * Arguments:   - uint16_t *p: pointer to in/output polynomial
 *              - StrAPoly mp: strided masked polynomial
 **************************************************/
void masked_poly(StrAPoly mp, const Poly p, uint16_t modulus) {
  for (int i = 0; i < SABER_N; i++) {
    mp[0][i] = (p[i] + modulus) % modulus;
  }
  mask_poly_inplace(mp, modulus);
}

// assumes first share contains correct value
void mask_poly_inplace(StrAPoly mp, uint16_t modulus) {
  uint16_t v[2];
  for (int d = 1; d < NSHARES; d++) {
    for (int i = 0; i < SABER_N; i += 2) {
      rand_q(v, modulus);
      mp[0][i] = (mp[0][i] + v[0]) % modulus;
      mp[0][i + 1] = (mp[0][i + 1] + v[1]) % modulus;

      mp[d][i] = (modulus - v[0]) % modulus;
      mp[d][i + 1] = (modulus - v[1]) % modulus;
    }
  }
}

/*************************************************
 * Name:        unmasked_poly
 *
 * Description: unmasked a strided masked polynomial
 *
 * Arguments:   - uint16_t *p: pointer to in/output polynomial
 *              - StrAPoly mp: strided masked polynomial
 **************************************************/
void unmasked_poly(Poly p, const StrAPoly mp, uint16_t modulus) {

  for (int i = 0; i < SABER_N; i++) {
    p[i] = mp[0][i];
  }

  for (int d = 1; d < NSHARES; d++) {
    for (int i = 0; i < SABER_N; i += 1) {
      p[i] = (p[i] + mp[d][i]) % modulus;
    }
  }
}
