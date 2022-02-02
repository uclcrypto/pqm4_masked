#include <stdint.h>
#include "masked.h"
#include "gadgets.h"
#include "masked_utils.h"

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

void masked_and(
        size_t nshares,
        uint32_t *z, size_t z_stride,
        const uint32_t *a, size_t a_stride,
        const uint32_t *b, size_t b_stride
        ) {
    uint32_t ztmp[nshares];
    uint32_t r;
    uint32_t i,j;
    
    for(i=0;i<nshares;i++){
        ztmp[i] = a[i*a_stride] & b[i*b_stride];
    }

    for(i=0;i<(nshares-1);i++){
        for(j=i+1; j<nshares;j++){
            r = rand32();
            // PINI
            ztmp[i] ^= pini_and_core(a[i*a_stride],b[j*b_stride],r);
            ztmp[j] ^= pini_and_core(a[j*a_stride],b[i*b_stride],r);
        }
    }
    for(i=0;i<nshares;i++){
        z[i*z_stride] = ztmp[i];
    }
}

void masked_xor(
        size_t nshares,
        uint32_t *out, size_t out_stride,
        const uint32_t *ina, size_t ina_stride,
        const uint32_t *inb, size_t inb_stride
        ) {
    for (size_t i=0; i<nshares; i++){
        out[i*out_stride] = ina[i*ina_stride] ^ inb[i*inb_stride];
    }
}

void copy_sharing(
        size_t nshares,
        uint32_t *out, size_t out_stride,
        const uint32_t *in, size_t in_stride
        ) {
    for (size_t i=0; i<nshares; i++){
        out[i*out_stride] = in[i*in_stride];
    }
}
