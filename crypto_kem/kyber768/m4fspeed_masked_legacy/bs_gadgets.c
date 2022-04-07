#include <stdint.h>
#include <stdlib.h>
#include "gadgets.h"
#include "bs_gadgets.h"
#include "masked_utils.h"

#define KYBER_MASKING_ORDER (D-1)
#define D NSHARES
#define Q 3329
#define GET_BIT_SHARING(i,j) ((i)*(j)) 

//#define SNI
void RefreshXOR(uint32_t *z,
        uint32_t d){

    uint32_t r;
    for(uint32_t i=0;i<(d-1);i++){
        for(uint32_t j=i+1; j<d;j++){
            r = get_random();
            z[i] ^= (r);
            z[j] ^= (r); 
        }
    }
}

void RefreshIOS(uint32_t *z,
        uint32_t d){
    
    uint32_t r;
    for(uint32_t i=0;i<(d-1);i++){
        for(uint32_t j=i+1; j<d;j++){
            r = get_random();
            z[i] ^= (r);
            z[j] ^= (r); 
        }
    }
}


void RefreshXOR16(uint32_t *z,
        uint32_t d){

    uint32_t r;
    for(uint32_t i=0;i<(d-1);i++){
        for(uint32_t j=i+1; j<d;j++){
            r = rand16();
            z[i] ^= (r);
            z[j] ^= (r); 
        }
    }
}


void FullRefreshXOR(uint32_t *x,
        uint32_t d){

    uint32_t r;
    
    uint32_t x0 = x[0];
    for(uint32_t i=0;i<(d);i++){
        for(uint32_t j=1; j<d;j++){
            r = get_random();
            x0 ^= (r);
            x[j] ^= (r); 
        }
    }
    x[0] = x0;
}

static inline uint32_t pini_and_core(uint32_t a, uint32_t b, uint32_t r) {
    uint32_t temp;
    uint32_t s;
    asm(
    "eor %[temp], %[b], %[r]\n\t"
    "and %[temp], %[a], %[temp]\n\t"
    "bic %[s], %[r], %[a]\n\t"
    "eor %[s], %[s], %[temp]"
    :[s]"=r"(s), [temp]"=&r"(temp) /* outputs, use temp as an arbitrary-location clobber */
    :[a]"r"(a),  [b]"r"(b), [r]"r"(r)  /* inputs */
    );
    return s;
}

void SecAnd16(uint32_t *x,
        uint32_t *y,
        uint32_t *z,
        uint32_t d){
    masked_and(d,z,1,x,1,y,1);
}

void SecAnd(uint32_t *x,
        uint32_t *y,
        uint32_t *z,
        uint32_t d){
    masked_and(d,z,1,x,1,y,1);
}
void SecANDbs(uint32_t *z,
        const uint32_t *a,
        const uint32_t *b){
    masked_and(NSHARES,z,1,a,1,b,1);
}

void SecXORbs(uint32_t *z,
        const uint32_t *a,
        const uint32_t *b){
    masked_xor(NSHARES,z,1,a,1,b,1);
}
