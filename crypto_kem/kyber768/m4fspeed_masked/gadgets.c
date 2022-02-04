#include <stdint.h>
#include "masked.h"
#include "gadgets.h"
#include "masked_utils.h"
#include "masked_representations.h"

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
        size_t kbits_out,
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
            
            if((i == (kbits-1)) && (i == (kbits_out-1))){ 
                return;
            }else if(i == (kbits-1)){
                masked_and(nshares,
                        carry,1,
                        xpy,1,
                        xpc,1);
                masked_xor(nshares,
                        &out[(kbits)*out_data_stride],out_msk_stride,
                        carry,1,
                        &in1[i*in1_data_stride],in1_msk_stride);
                return;
            }

            masked_and(nshares,
                    carry,1,
                    xpy,1,
                    xpc,1);
            masked_xor(nshares,
                    carry,1,
                    carry,1,
                    &in1[i*in1_data_stride],in1_msk_stride);
        }
    }
}

void secadd_constant_bmsk(size_t nshares,
        size_t kbits,
        size_t kbits_out,
        uint32_t *out, size_t out_msk_stride, size_t out_data_stride,
        const uint32_t *in1, size_t in1_msk_stride, size_t in1_data_stride,
        uint32_t constant, const uint32_t *bmsk, size_t bmsk_msk_stride){
    
    size_t i,d;
    uint32_t carry[nshares];
    uint32_t xpy[nshares];
    uint32_t xpc[nshares];
   
    if(constant & 0x1){
        masked_and(nshares,
                    carry,1,
                    &in1[0*in1_data_stride],in1_msk_stride,
                    bmsk,bmsk_msk_stride);
        masked_xor(nshares,
                &out[0*out_data_stride],out_msk_stride,
                &in1[0*in1_data_stride],in1_msk_stride,
                bmsk,bmsk_msk_stride);
    }else{
        for(d=0;d<nshares;d++){
            carry[d] = 0;
        }
        copy_sharing(nshares,
                out,out_msk_stride,
                in1,in1_msk_stride);
    }
    for(i=1;i<kbits;i++){
        // xpy = in2 ^ in1
        // xpc = in1 ^ carry
        // out = xpy ^ carry
        if((constant >> i)&0x1){
            for(d= 0;d<nshares;d++){
                xpy[d] = in1[i*in1_data_stride + d*in1_msk_stride] ^ bmsk[d*bmsk_msk_stride];
                xpc[d] = in1[i*in1_data_stride + d*in1_msk_stride] ^ carry[d];
                out[i*out_data_stride + d*out_msk_stride] = xpy[d] ^ carry[d];
            }
            
            if((i == (kbits-1)) && (i == (kbits_out-1))){ 
                return;
            }else if(i == (kbits-1)){
                masked_and(nshares,
                        carry,1,
                        xpy,1,
                        xpc,1);
                masked_xor(nshares,
                        &out[(kbits)*out_data_stride],out_msk_stride,
                        carry,1,
                        &in1[i*in1_data_stride],in1_msk_stride);
                return;
            }

            masked_and(nshares,
                    carry,1,
                    xpy,1,
                    xpc,1);
            masked_xor(nshares,
                    carry,1,
                    carry,1,
                    &in1[i*in1_data_stride],in1_msk_stride);
        }else{
            // compute the carry
            masked_xor(nshares,
                    &out[i*out_data_stride],out_msk_stride,
                    carry,1,
                    &in1[i*in1_data_stride],1);
            
            if((i == (kbits-1)) && (i == (kbits_out-1))){ 
                return;
            }else if(i == (kbits-1)){
                masked_and(nshares,
                        &out[(kbits)*out_data_stride],out_msk_stride,
                        carry,1,
                        &in1[i*in1_data_stride],in1_msk_stride);
                return;
            }
            masked_and(nshares,
                    carry,1,
                    carry,1,
                    &in1[i*in1_data_stride],in1_msk_stride);
        }
    }
}

void secadd_constant(size_t nshares,
        size_t kbits,
        size_t kbits_out,
        uint32_t *out, size_t out_msk_stride, size_t out_data_stride,
        const uint32_t *in1, size_t in1_msk_stride, size_t in1_data_stride,
        uint32_t constant){
    
    size_t i,d;
    uint32_t carry[nshares];
    uint32_t xpy[nshares];
    uint32_t xpc[nshares];
   
    if(constant & 0x1){
        for(d=0; d<nshares; d++){
            carry[d] = in1[d*in1_msk_stride];
        }
        copy_sharing(nshares,
                out,out_msk_stride,
                in1,in1_msk_stride);
        out[0] ^= 0xFFFFFFFF;
    }else{
        for(d=0;d<nshares;d++){
            carry[d] = 0;
        }
        copy_sharing(nshares,
                out,out_msk_stride,
                in1,in1_msk_stride);
    }

    for(i=1;i<kbits;i++){
        if((constant >> i)&0x1){
            for(d= 0;d<nshares;d++){
                xpy[d] = in1[i*in1_data_stride + d*in1_msk_stride];
                xpc[d] = in1[i*in1_data_stride + d*in1_msk_stride] ^ carry[d];
                out[i*out_data_stride + d*out_msk_stride] = xpy[d] ^ carry[d];
            }
            xpy[0] ^= 0xFFFFFFFF;
            out[i*out_data_stride] ^= 0xFFFFFFFF;

            if((i == (kbits-1)) && (i == (kbits_out-1))){ 
                return;
            }else if(i == (kbits-1)){
                masked_and(nshares,
                        carry,1,
                        xpy,1,
                        xpc,1);
                masked_xor(nshares,
                        &out[(kbits)*out_data_stride],out_msk_stride,
                        carry,1,
                        &in1[i*in1_data_stride],in1_msk_stride);
                return;
            }

            masked_and(nshares,
                    carry,1,
                    xpy,1,
                    xpc,1);
            masked_xor(nshares,
                    carry,1,
                    carry,1,
                    &in1[i*in1_data_stride],in1_msk_stride);
        }else{
            // compute the carry
            masked_xor(nshares,
                    &out[i*out_data_stride],out_msk_stride,
                    carry,1,
                    &in1[i*in1_data_stride],1);
            
            if((i == (kbits-1)) && (i == (kbits_out-1))){ 
                return;
            }else if(i == (kbits-1)){
                masked_and(nshares,
                        &out[(kbits)*out_data_stride],out_msk_stride,
                        carry,1,
                        &in1[i*in1_data_stride],in1_msk_stride);
                return;
            }
            masked_and(nshares,
                    carry,1,
                    carry,1,
                    &in1[i*in1_data_stride],in1_msk_stride);
        }
    }
}



void secadd_modq(size_t nshares,
        size_t kbits,
        uint32_t q,
        uint32_t *out, size_t out_msk_stride, size_t out_data_stride,
        const uint32_t *in1, size_t in1_msk_stride, size_t in1_data_stride,
        const uint32_t *in2, size_t in2_msk_stride, size_t in2_data_stride){

    uint32_t s[(kbits+1)*NSHARES];
    uint32_t sp[(kbits+1)*NSHARES];
    uint32_t bmsk_ones[NSHARES];

    secadd(nshares,
            kbits,kbits+1,
            s,1,NSHARES,
            in1,in1_msk_stride,in1_data_stride,
            in2,in2_msk_stride,in2_data_stride);
    bmsk_ones[0] = 0xFFFFFFFF;
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

    secadd(nshares,kbits,kbits,
            in,in_msk_stride,in_data_stride,
            expanded_low,1,nshares,
            expanded_high,1,nshares);
}

void seccompress(size_t nshares,
                    size_t ncoeffs,
                    uint32_t q,
                    uint32_t c,
                    uint32_t *out, size_t out_msk_stride, size_t out_data_stride,
                    const int16_t *in, size_t in_msk_stride, size_t in_data_stride){

    size_t i,d;
    uint32_t ell=0;
    uint32_t prod = q * nshares;
    while(prod > 0){
        prod = prod >> 1;
        ell++;
    }

    uint32_t in_expanded[ncoeffs*nshares];
    uint32_t bs_expanded[(ell+c)*nshares];

    // map mod q to mod 2^ell.
    uint32_t tmp32;
    uint64_t tmp64;
    for(i=0;i<ncoeffs;i++){
        tmp64 = (in[i*in_data_stride]+q)%q;
        tmp32 = ((((tmp64<<(c+ell+1)) +q)>>1)/(q));
        tmp32 += (1<<(ell-1));
        in_expanded[i*nshares] = tmp32 & ((1<<(ell+c))-1);
        for(d=1;d<nshares;d++){
            tmp64 = (q+in[i*in_data_stride + d*in_msk_stride])%q;
            tmp32 = ((((tmp64<<(c+ell+1)) + q)>>1)/(q));
            in_expanded[i*nshares+d] = tmp32 & ((1<<(ell+c))-1);
        }
    }
    
    // map to bitslice
    masked_dense2bitslice_u32(
            nshares,
            ncoeffs,
            ell+c,
            bs_expanded,1,nshares,
            in_expanded,1,nshares);

    // convert A2B
    seca2b(nshares,ell+c,
            bs_expanded,1,nshares);

    // map to the output
    for(i=0;i<c;i++){
        for(d=0;d<nshares;d++){
            out[i*out_data_stride + d * out_msk_stride] = bs_expanded[(ell+i)*nshares + d];
        }
    }
}

/*static void init_tables(void);
// init a tables with q and -q in it, bitslice form
static uint32_t is_init = 0;
static uint32_t bs_q[NSHARES*(COEF_NBITS+1)];
static uint32_t bs_negq[NSHARES*(COEF_NBITS+1)];

static void init_tables(){

    int16_t tmp[1];

    tmp[0] = (1<<COEF_NBITS) - KYBER_Q;
    map_dense2bitslice(
            1,
            1,
            COEF_NBITS+1,
            bs_negq,1,NSHARES,
            tmp,1,1);
    
    tmp[0] = KYBER_Q;
    map_dense2bitslice(
            1,
            1,
            COEF_NBITS+1,
            bs_q,1,NSHARES,
            tmp,1,1);  

    for(int i=0;i<COEF_NBITS+1;i++){
        bs_q[i*NSHARES] *= 0xFFFFFFFF;
        bs_negq[i*NSHARES] *= 0xFFFFFFFF;
        for(int d=1;d<NSHARES;d++){
            bs_q[i*NSHARES + d] = 0;
            bs_negq[i*NSHARES + d] = 0;
        }
    }
    is_init = 1;
    return;
}
*/

