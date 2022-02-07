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

/*************************************************
* Name:        masked_poly_getnoise
*
* Description: Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter KYBER_ETA
*
* Arguments:   - poly *r:                   pointer to output polynomial
*              - const unsigned char *seed: pointer to input seed (pointing to array of length KYBER_SYMBYTES bytes)
*              - unsigned char nonce:       one-byte input nonce
*              - int add:                   boolean to indicate to accumulate into r
**************************************************/
void masked_poly_noise(StrAPoly r, const unsigned char *seed, unsigned char nonce, int add) {
    size_t kappa = KYBER_ETA;
    unsigned char buf[KYBER_ETA * KYBER_N / 4];
    unsigned char buf_masked[(KYBER_ETA * KYBER_N / 4)*NSHARES];

    uint32_t a[kappa*NSHARES];
    uint32_t b[kappa*NSHARES];
    int16_t out[NSHARES*BSSIZE];

    prf(buf, KYBER_ETA * KYBER_N / 4, seed, nonce);
    for(int i=0;i<(KYBER_ETA * KYBER_N/4);i++){
        buf_masked[(0*KYBER_ETA * KYBER_N / 4) + i] = buf[i];
        for(int d=1;d<NSHARES;d++){
            buf_masked[(d*KYBER_ETA * KYBER_N / 4) + i] = 0;
        }
    }

    // all the bitslice. 32*4 bits =  
    for(uint32_t i=0;i<KYBER_N/BSSIZE;i++){
        // all the bits
        //32*4 bits = 128 bits = 16 bytes
        for(uint32_t j=0;j<kappa*NSHARES;j++){
            a[j] = 0; b[j] = 0;
        }
        for(int n=0;n<16;n++){
            for(uint32_t d=0;d<NSHARES;d++){
                uint32_t bytes_off = i*16 + n;
                uint32_t by = buf_masked[bytes_off + d*(kappa*KYBER_N/4)];
                
                a[d] = (a[d] << 2) | (((by>>0)&0x1)<<1) | (((by>>4)&0x1)<<0);
                a[NSHARES+d] = (a[NSHARES+d] << 2) | (((by>>1)&0x1)<<1) | (((by>>5)&0x1)<<0);

                b[d] = (b[d] << 2) | (((by>>2)&0x1)<<1) | (((by>>6)&0x1)<<0);
                b[NSHARES+d] = (b[NSHARES+d] << 2) | (((by>>3)&0x1)<<1) | (((by>>7)&0x1)<<0);
            }
        }
        masked_cbd(NSHARES,
                2,
                BSSIZE,
                KYBER_Q,COEF_NBITS,
                out,1,NSHARES,
                a,1,NSHARES,
                b,1,NSHARES);


        for(uint32_t n=0;n<BSSIZE;n++){
            for(uint32_t j=0;j<NSHARES;j++){
                if(add){
                    r[j][(i*32)+31-n] = (r[j][(i*32)+31-n] + out[n*NSHARES+j])%KYBER_Q;
                }else{
                    r[j][(i*32)+31-n] = out[n*NSHARES+j];
                }
            }
        }
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
