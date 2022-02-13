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
#ifndef MASKED_UTILS_H
#define MASKED_UTILS_H
#include "masked.h"
#include "SABER_params.h"
#include <libopencm3/stm32/rng.h>
#include <stdint.h>

#if BENCH_RND==1
extern uint64_t rng_cnt;
#endif

inline uint32_t get_random() {
#if BENCH_RND==1
  rng_cnt+=1;
#endif 
  while (1) {
    if ((RNG_SR & RNG_SR_DRDY) == 1) { // check if data is ready
      return RNG_DR;
    }
  }
  return 0;
}

inline uint32_t rand32() { return get_random(); }

// No rejection sampling -> requires modulus to be power of 2
inline void rand_q(uint16_t v[2], uint16_t modulus) {
  uint32_t r = rand32();
  v[0] = r % modulus;
  v[1] = (r / modulus) % modulus;
}

void masked_poly(StrAPoly mp, const Poly p, uint16_t modulus);
void unmasked_poly(Poly p, const StrAPoly mp, uint16_t modulus);
void mask_poly_inplace(StrAPoly mp, uint16_t modulus);
#endif
