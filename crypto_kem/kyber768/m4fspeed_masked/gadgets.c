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

uint32_t unmask_boolean(
        size_t nshares,
        const uint32_t *in, size_t in_stride){
    uint32_t out,d;
    out =0;
    for(d=0;d<nshares;d++){
        out ^= in[d*in_stride];
    }
    return out;
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


void secadd(size_t nshares,
        size_t kbits,
        uint32_t *out, size_t out_msk_stride, size_t out_data_stride,
        const uint32_t *in1, size_t in1_msk_stride, size_t in1_data_stride,
        const uint32_t *in2, size_t in2_msk_stride, size_t in2_data_stride){
    
    size_t i,d;
    if(nshares==1){
        uint32_t c = in1[0] & in2[0];
        out[0] = in1[0] ^ in2[0];
        for(i=1;i<kbits-1;i++){
            uint32_t tmp = in1[i*in1_data_stride] ^ in2[i*in2_data_stride] ^ c;
            c = (in1[i*in1_data_stride]&in2[i*in2_data_stride]) ^ (c & (in1[i*in1_data_stride] ^ in2[i*in2_data_stride]));
            out[i*out_data_stride] = tmp;
        }
        out[(kbits-1)*out_data_stride] = c ^ in1[(kbits-1)*in1_data_stride] ^ in2[(kbits-1)*in2_data_stride];
    }else{

        uint32_t carry[nshares];
        uint32_t xpy[nshares];
        uint32_t xpc[nshares];
        
        masked_and(nshares,
                    carry,1,
                    &in1[0*in1_data_stride],in1_msk_stride,
                    &in2[0*in2_data_stride],in2_msk_stride);
        
        masked_xor(nshares,
                &out[0*out_data_stride],out_msk_stride,
                &in1[0*in1_data_stride],in1_msk_stride,
                &in2[0*in2_data_stride],in2_msk_stride);

        for(i=1;i<kbits;i++){

            // xpy = in2 ^ in1
            // xpc = in1 ^ carry
            // out = xpy ^ carry
            for(d= 0;d<nshares;d++){
                xpy[d] = in1[i*in1_data_stride + d*in1_msk_stride] ^ in2[i*in2_data_stride + d*in2_msk_stride];
                xpc[d] = in1[i*in1_data_stride + d*in1_msk_stride] ^ carry[d];
                out[i*out_data_stride + d*out_msk_stride] = xpy[d] ^ carry[d];
            }
            
            if(i == (kbits-1)){ 
                return;
            }else{
                masked_and(nshares,
                        carry,1,
                        xpy,1,
                        xpc,1);
            }
            masked_xor(nshares,
                    carry,1,
                    carry,1,
                    &in1[i*in1_data_stride],in1_msk_stride);
        }
    }
}

void seca2b(size_t nshares,
                size_t kbits,
                uint32_t *in, size_t in_msk_stride, size_t in_data_stride){

    // TODO optimize inplace ? 
    size_t i,d;

    if(nshares==1){
        return;
    }

    size_t nshares_low = nshares/2;
    size_t nshares_high = nshares - nshares_low;

    seca2b(nshares_low,kbits,in,in_msk_stride,in_data_stride);
    seca2b(nshares_high,kbits,&in[nshares_low],in_msk_stride,in_data_stride);

    uint32_t expanded_low[kbits*nshares];
    uint32_t expanded_high[kbits*nshares];

    for(i=0;i<kbits;i++){
        for(d=0;d<nshares_low;d++){
            expanded_low[i*nshares + d] = in[i*in_data_stride + d*in_msk_stride];
            expanded_high[i*nshares + d] = 0;
        }
        for(d=nshares_low;d<nshares;d++){
            expanded_high[i*nshares + d] = in[i*in_data_stride + d*in_msk_stride];
            expanded_low[i*nshares + d] = 0;
        }
    }

    secadd(nshares,kbits,
            in,in_msk_stride,in_data_stride,
            expanded_low,1,nshares,
            expanded_high,1,nshares);
}

