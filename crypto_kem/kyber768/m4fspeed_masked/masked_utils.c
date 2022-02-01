#include "masked_utils.h"
#include "masked.h"
#include <libopencm3/stm32/rng.h>
#include "poly.h"

extern inline uint32_t get_random();
extern inline uint32_t rand32();
extern inline void rand_q(uint16_t v[2]);

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


