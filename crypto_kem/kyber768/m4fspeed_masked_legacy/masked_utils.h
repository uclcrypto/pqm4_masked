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
#include "poly.h"
#include <libopencm3/stm32/rng.h>
#include <stdint.h>

#if BENCH_RND == 1
extern uint64_t rng_cnt;
#endif

inline uint32_t get_random() {
#if BENCH_RND == 1
  rng_cnt += 1;
#endif
  while (1) {
    if ((RNG_SR & RNG_SR_DRDY) == 1) { // check if data is ready
      return RNG_DR;
    }
  }
  return 0;
}

inline uint32_t rand32() { return get_random(); }

inline void rand_q(int16_t v[2]) {
  uint32_t r;
  do {
    r = rand32();
  } while (r > (387U * KYBER_Q * KYBER_Q));
  r = r % (KYBER_Q * KYBER_Q);
  v[0] = r % (KYBER_Q);
  v[1] = r / KYBER_Q;
}

extern uint16_t sampleq_id;
extern int16_t sampleq_rs[2];
inline int16_t sampleq(){
    sampleq_id ^= 0x1;
    if((sampleq_id&0x1) == 0x1){
        rand_q(sampleq_rs);
    }
    return sampleq_rs[sampleq_id];
}

extern uint16_t r16_id;
extern uint32_t r16_rs;
inline uint16_t rand16(){
    r16_id ^= 0x1;
    if((r16_id&0x1) == 0x1){
        r16_rs = rand32();
    }
    return (r16_rs >> (r16_id * 16)) & 0xFFFF; 
}
inline uint64_t rand64(){
  return ((uint64_t)rand32() << 32) + (uint64_t)rand32();
}
void masked_poly(StrAPoly mp, const poly *p);
void unmasked_poly(poly *p, const StrAPoly mp);
#endif
