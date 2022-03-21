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
#include "masked.h"
#include "poly.h"
#include <libopencm3/stm32/rng.h>

extern inline uint32_t get_random();
extern inline uint32_t rand32();
extern inline void rand_q(int16_t v[2]);

/*************************************************
 * Name:        masked_poly
 *
 * Description: masked a polynomial into a strided polynomial
 *
 * Arguments:   - uint16_t *p: pointer to in/output polynomial
 *              - StrAPoly mp: strided masked polynomial
 **************************************************/
void masked_poly(StrAPoly mp, const poly *p) {
  int16_t v[2];

  for (int i = 0; i < KYBER_N; i++) {
    mp[0][i] = p->coeffs[i];
  }

  for (int d = 1; d < NSHARES; d++) {
    for (int i = 0; i < KYBER_N; i += 2) {
      rand_q(v);
      mp[0][i] = (mp[0][i] + v[0]) % KYBER_Q;
      mp[0][i + 1] = (mp[0][i + 1] + v[1]) % KYBER_Q;

      mp[d][i] = (KYBER_Q - v[0]) % KYBER_Q;
      mp[d][i + 1] = (KYBER_Q - v[1]) % KYBER_Q;
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
void unmasked_poly(poly *p, const StrAPoly mp) {

  for (int i = 0; i < KYBER_N; i++) {
    p->coeffs[i] = mp[0][i];
  }

  for (int d = 1; d < NSHARES; d++) {
    for (int i = 0; i < KYBER_N; i += 1) {
      p->coeffs[i] = (p->coeffs[i] + mp[d][i] + KYBER_Q) % KYBER_Q;
    }
  }
}
