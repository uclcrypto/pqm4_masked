#include <stdint.h>
#include "masked.h"
#include "masked_representations.h"

/*************************************************
* Name:        StrAPoly2APoly
*
* Description: Maps strided polynomial into a dense representation
*
* Arguments:  
*           - APoly out : dense polynomial
*           - StrAPoly in: strided polynomial
* **************************************************/
void StrAPoly2APoly(APoly out, const StrAPoly in){
    int i,d;
    for(i=0;i<KYBER_N;i++){
        for(d=0;d<NSHARES;d++){
            out[i][d] = in[d][i];
        }
    }
}

/*************************************************
* Name:        APoly2StrAPoly
*
* Description: Maps dense polynomial into a strided representation
*
* Arguments:  
*           - StrAPoly out: strided polynomial
*           - APoly in : dense polynomial
* **************************************************/
void APoly2StrAPoly(StrAPoly out, const APoly in){

    int i,d;
    for(i=0;i<KYBER_N;i++){
        for(d=0;d<NSHARES;d++){
            out[d][i] = in[i][d];
        }
    }
}

/*************************************************
* Name:        masked_dense2bitslice
*
* Description: maps a dense reprensentation to a bitlisce one 
*
* Arguments:  
*           - uint32_t *bitslice[]: output bitslice representation. Table of coeffs_size x nshares.
*           - int16_t *dense[]: input dense reprensetation. Table of n_coeffs x nshares.
*           - size_t coeffs_size: number of bits to represent the dense coefficients
*           - size_t n_coeffs: number of coefficients
*           - size_t nshares: number of shares
* **************************************************/
void masked_dense2bitslice(
        uint32_t *bitslice[],
        int16_t *dense[],
        size_t coeffs_size,
        size_t n_coeffs,
        size_t nshares){

    size_t d,c,b;
    for(b=0;b<coeffs_size;b++){
        for(d=0;d<nshares;d++){
            bitslice[b][d] = 0;
        }
    }
    
    
    for(d=0;d<nshares;d++){
        for(c=0; c<n_coeffs;c++){
            int16_t xd = dense[c][d];
            for(b=0; b<coeffs_size;b++){
                bitslice[b][d] |= (xd&0x1)<<c;
                xd = xd >> 1;
            }
        }
    }
}

/*************************************************
* Name:        masked_bitslice2dense
*
* Description: maps a bitslice reprensentation to a dense one 
*
* Arguments:  
*           - uint32_t *bitslice[]: input bitslice representation. Table of coeffs_size x nshares.
*           - int16_t *dense[]: output dense reprensetation. Table of n_coeffs x nshares.
*           - size_t coeffs_size: number of bits to represent the dense coefficients
*           - size_t n_coeffs: number of coefficients
*           - size_t nshares: number of shares
* **************************************************/
void masked_bitslice2dense(
        int16_t *dense[],
        uint32_t *bitslice[],
        size_t coeffs_size,
        size_t n_coeffs,
        size_t nshares){

    size_t d,c,b;
    for(d=0;d<nshares;d++){
        for(c=0; c<n_coeffs;c++){
            int16_t xd = 0;
            for(b=0; b<coeffs_size;b++){
                xd |= ((((bitslice[b][d])>>c)&0x1)<<b);
            }
            dense[c][d] = xd;
        }
    }
}
