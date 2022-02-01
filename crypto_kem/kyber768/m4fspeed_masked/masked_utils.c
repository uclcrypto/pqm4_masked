#include "masked_utils.h"
#include "masked.h"
#include <libopencm3/stm32/rng.h>
#include "poly.h"

inline uint32_t get_random(){
    while (1) {
        if ((RNG_SR & RNG_SR_DRDY) == 1){  // check if data is ready
            return RNG_DR;
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

void masked_poly(StrAPoly mp, const poly *p){
    uint16_t v[2];

    for(int i=0;i<KYBER_N;i++){
        mp[0][i] = p->coeffs[i];
    }
    
    for(int d=1;d<NSHARES;d++){
        for(int i=0;i<KYBER_N;i+=2){
            rand_q(v);
            mp[0][i] = (mp[0][i] + v[0])%KYBER_Q;
            mp[0][i+1] = (mp[0][i+1] + v[1])%KYBER_Q;

            mp[d][i] = (mp[0][i] + KYBER_Q - v[0])%KYBER_Q;
            mp[d][i+1] = (mp[0][i+1] + KYBER_Q - v[1])%KYBER_Q;
        }
    }
}


