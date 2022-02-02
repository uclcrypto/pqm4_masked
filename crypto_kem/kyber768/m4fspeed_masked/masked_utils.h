#ifndef MASKED_UTILS_H
#define MASKED_UTILS_H
#include <stdint.h>
#include <libopencm3/stm32/rng.h>
#include "masked.h"
#include "poly.h"

inline uint32_t get_random(){
    while (1) {
        if ((RNG_SR & RNG_SR_DRDY) == 1){  // check if data is ready
            //return RNG_DR;
            return 0;
        }
    }
    return 0;
}

inline uint32_t rand32(){
    return get_random();
}

inline void rand_q(uint16_t v[2]){
  uint32_t r;
  do{
    r = rand32();
  } while(r > (387U * KYBER_Q*KYBER_Q));
  r = r%(KYBER_Q*KYBER_Q);
  v[0] = r%(KYBER_Q);
  v[1] = r/KYBER_Q;
}

void masked_poly(StrAPoly mp, const poly *p);
void unmasked_poly(poly *p, const StrAPoly mp);
#endif
