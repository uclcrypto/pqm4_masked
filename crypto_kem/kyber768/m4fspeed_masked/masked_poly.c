#include "poly.h"
#include "masked_poly.h"
#include "masked.h"
#include "masked_representations.h"

#include "cbd.h"
#include "ntt.h"
#include "params.h"
#include "symmetric.h"
#include "gadgets.h"

#include <stdint.h>

/*************************************************
* Name:        masked_poly_ntt
*
* Description: Computes negacyclic number-theoretic transform (NTT) of
*              a polynomial in place;
*              inputs assumed to be in normal order, output in bitreversed order
*
* Arguments:   - uint16_t *r: pointer to in/output polynomial
**************************************************/
void masked_poly_ntt(StrAPoly r) {
    for(int d=0;d<NSHARES;d++){
        ntt(r[d]);
    }
}

/*************************************************
* Name:        masked_poly_invntt
*
* Description: Computes inverse of negacyclic number-theoretic transform (NTT) of
*              a polynomial in place;
*              inputs assumed to be in bitreversed order, output in normal order
*
* Arguments:   - uint16_t *a: pointer to in/output polynomial
**************************************************/
void masked_poly_invntt(StrAPoly r) {
    for(int d=0; d<NSHARES;d++){
        invntt(r[d]);
    }
}

void masked_poly_tomsg(unsigned char *m, StrAPoly str_r){
    APoly r;
    size_t i,j,d;
    uint32_t bits[NSHARES];

    StrAPoly2APoly(r,str_r);
    for(i=0;i<KYBER_N;i+=BSSIZE){

        seccompress(NSHARES,BSSIZE,KYBER_Q,1,
                bits,1,NSHARES,
                r[i],1,NSHARES);

        for(d=0;d<NSHARES;d++){
            for(j=0;j<BSSIZE/8;j++){
                m[d*KYBER_N + (i/8)+j] = (bits[d]>>(j*8))&0xFF;
            }
        }
    }
}

void masked_poly_cmp(
        size_t c,
        uint32_t *rc,
        const StrAPoly mp,
        const poly *ref){

    APoly r;
    size_t i,b;
    uint32_t bits[NSHARES*c];
    uint32_t bits_ref[c];

    StrAPoly2APoly(r, mp);
 
    for(i=0;i<KYBER_N;i+=BSSIZE){

        // compress masked polynomial
        seccompress(NSHARES,BSSIZE,KYBER_Q,c,
                bits,1,NSHARES,
                r[i],1,NSHARES);

        // map public polynomial to bitslice
        masked_dense2bitslice(1,
                BSSIZE,
                c,
                bits_ref,1,1,
                &(ref->coeffs[i]),1,1);

        for(b=0;b<c;b++){
            // public polynomial and public one
            bits[b*NSHARES] ^= bits_ref[b] ^ 0xFFFFFFFF;

            masked_and(NSHARES,
                    rc,1,
                    rc,1,
                    &bits[b*NSHARES],1);
        }
    }
}
void finalize_cmp(uint32_t *bits){
    uint32_t other[NSHARES];
    int d;
    for(d=0;d<NSHARES;d++){
        other[d] = bits[d] >> 16;
    }
    masked_and(NSHARES,
            bits,1,
            bits,1,
            other,1);

    for(d=0;d<NSHARES;d++){
        other[d] = bits[d] >> 8;
    }
    masked_and(NSHARES,
            bits,1,
            bits,1,
            other,1);
    for(d=0;d<NSHARES;d++){
        other[d] = bits[d] >> 8;
    }
    masked_and(NSHARES,
            bits,1,
            bits,1,
            other,1);

    for(d=0;d<NSHARES;d++){
        other[d] = bits[d] >> 4;
    }
    masked_and(NSHARES,
            bits,1,
            bits,1,
            other,1);

    for(d=0;d<NSHARES;d++){
        other[d] = bits[d] >> 2;
    }
    masked_and(NSHARES,
            bits,1,
            bits,1,
            other,1);

    for(d=0;d<NSHARES;d++){
        other[d] = bits[d] >> 1;
    }
    masked_and(NSHARES,
            bits,1,
            bits,1,
            other,1);
}
