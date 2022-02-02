#include <stdint.h>
#include "masked.h"
#include "masked_representations.h"

void StrAPoly2APoly(APoly out, const StrAPoly in){
    int i,d;
    for(i=0;i<KYBER_N;i++){
        for(d=0;d<NSHARES;d++){
            out[i][d] = in[d][i];
        }
    }
}

void APoly2StrAPoly(StrAPoly out, const APoly in){

    int i,d;
    for(i=0;i<KYBER_N;i++){
        for(d=0;d<NSHARES;d++){
            out[d][i] = in[i][d];
        }
    }
}
