#include <string.h>
#include "masked_poly.h"
#include "masked_cbd.h"
#include "cbd.h"
//#include "mNTT.h"
#include "NTT.h"
#include "fips202.h"
#include "pack_unpack.h"
#include "masked_utils.h"
#include "masked_representations.h"
#include "gadgets.h"
#include "masked.h"

#define h1 (1 << (SABER_EQ - SABER_EP - 1))
#define h2 ((1 << (SABER_EP - 2)) - (1 << (SABER_EP - SABER_ET - 1)) + (1 << (SABER_EQ - SABER_EP - 1)))
#define MAX(a,b) (((a)>(b))?(a):(b))

//#define MY_DEBUG
extern void __asm_poly_add_16(uint16_t *des, uint16_t *src1, uint16_t *src2);
extern void __asm_poly_add_32(uint32_t *des, uint32_t *src1, uint32_t *src2);

static inline shake128incctx shake128_absorb_seed(const uint8_t seed[SABER_SEEDBYTES]){

    shake128incctx ctx;
    shake128_inc_init(&ctx);
    shake128_inc_absorb(&ctx, seed, SABER_SEEDBYTES);
    shake128_inc_finalize(&ctx);

    return ctx;

}

uint32_t masked_MatrixVectorMulEncNTT(uint8_t ct0[SABER_POLYVECCOMPRESSEDBYTES], 
                uint8_t ct1[SABER_SCALEBYTES_KEM], 
                const uint8_t seed_s[SABER_NOISE_SEEDBYTES], 
                const uint8_t seed_A[SABER_SEEDBYTES], 
                const uint8_t pk[SABER_INDCPA_PUBLICKEYBYTES], 
                const uint8_t m[SABER_KEYBYTES], int compare){

    uint32_t acc_NTT[SABER_N];
    uint32_t A_NTT[SABER_N];
    uint32_t s_NTT[SABER_L * SABER_N];

    uint16_t poly[SABER_N];
    uint16_t poly_ref[SABER_N];
    uint16_t m_poly[NSHARES*SABER_N];
    uint16_t acc[SABER_N];
    uint16_t myref[SABER_N];
    uint16_t masked_acc[NSHARES*SABER_N];

    uint8_t shake_out[MAX(SABER_POLYBYTES, SABER_POLYCOINBYTES)];
    uint8_t masked_shake_out[MAX(SABER_POLYBYTES, SABER_POLYCOINBYTES)*NSHARES];

    uint16_t *mp = poly;

    size_t i, j;
    uint32_t fail = 0;

    shake128incctx shake_s_ctx = shake128_absorb_seed(seed_s);

    for(i = 0; i < SABER_L; i++){
        shake128_inc_squeeze(shake_out, SABER_POLYCOINBYTES, &shake_s_ctx);
        memset(masked_shake_out,0,sizeof(masked_shake_out));
        for(j=0;j<SABER_POLYCOINBYTES;j++)
          masked_shake_out[0*SABER_POLYCOINBYTES + j] = shake_out[j]; 
        
        masked_cbd_seed(NSHARES,
                    m_poly,SABER_N,1,
                    masked_shake_out,SABER_POLYCOINBYTES,1);
        
        cbd(poly_ref,shake_out);

        unmasked_poly(poly,m_poly,SABER_Q);

#ifdef MY_DEBUG
        char buf_x [128];
        hal_send_str("----------");
        for(int n = 0; n<SABER_N;n++){
          sprintf(buf_x,"n%d- > %d %d",n,(uint16_t) poly[n],(uint16_t) poly_ref[n]%SABER_Q);
          hal_send_str(buf_x);
        }
#endif
        NTT_forward_32(s_NTT + i * SABER_N, poly_ref);
    }

    shake128_inc_ctx_release(&shake_s_ctx);

    shake128incctx shake_A_ctx = shake128_absorb_seed(seed_A);

    uint32_t rc[NSHARES];
    memset(rc,0,sizeof(rc));
    rc[0] = 0xFFFFFFFF;

    for (i = 0; i < SABER_L; i++) {

        for (j = 0; j < SABER_L; j++) {

            shake128_inc_squeeze(shake_out, SABER_POLYBYTES, &shake_A_ctx);
            BS2POLq(shake_out, poly);
            
            NTT_forward_32(A_NTT, poly);

            if (j == 0) {
                NTT_mul_32(acc_NTT, A_NTT, s_NTT + j * SABER_N);
            } else {
                NTT_mul_32(A_NTT, A_NTT, s_NTT + j * SABER_N);
                __asm_poly_add_32(acc_NTT, acc_NTT, A_NTT);
            }
        }

        NTT_inv_32(acc, acc_NTT);

        for (j = 0; j < SABER_N; j++) {
            acc[j] = (acc[j] + h1);
        }
        masked_poly(masked_acc,acc,SABER_Q);
        for (j = 0; j < SABER_N; j++) {
            acc[j] = ((acc[j]) >> (SABER_EQ - SABER_EP));
        }

        
        if (compare) {
            BS2POLp(ct0 + i*SABER_POLYCOMPRESSEDBYTES,myref);
            for(j=0;j<SABER_N;j++){
              myref[j] = myref[j] % SABER_P;
            }
            masked_poly_cmp(SABER_EQ-SABER_EP,SABER_EQ,SABER_EQ,rc,masked_acc,myref);
        } else {
            POLp2BS(ct0 + i * SABER_POLYCOMPRESSEDBYTES, acc);
        }
    }

    shake128_inc_ctx_release(&shake_A_ctx);

    for(j = 0; j < SABER_L; j++){

        BS2POLp(pk + j * SABER_POLYCOMPRESSEDBYTES, poly);

        NTT_forward_32(A_NTT, poly);

        if(j == 0){
            NTT_mul_32(acc_NTT, A_NTT, s_NTT + j * SABER_N);
        }else{
            NTT_mul_32(A_NTT, A_NTT, s_NTT + j * SABER_N);
            __asm_poly_add_32(acc_NTT, acc_NTT, A_NTT);
        }
    }

    NTT_inv_32(acc, acc_NTT);

    BS2POLmsg(m, mp);

    for(j = 0; j < SABER_N; j++){
        acc[j] = (acc[j] - (mp[j] << (SABER_EP - 1)) + h1);
    }

    masked_poly(masked_acc,acc,SABER_P);
    for(j = 0; j < SABER_N; j++){
        acc[j] = (acc[j]) >> (SABER_EP - SABER_ET);
    }
 
    if(compare){
        fail |= POLT2BS_cmp(ct1, acc);
        BS2POLT(ct1,myref);
        for(j=0;j<SABER_N;j++){
            myref[j] = myref[j] % (1<<SABER_ET);
        }           
        masked_poly_cmp(SABER_EP-SABER_ET,SABER_EP,SABER_EP,rc,masked_acc,myref);
    }else{
        POLT2BS(ct1, acc);
    }
    
    finalize_cmp(rc);
    fail = 0;
    for(int d=0;d<NSHARES;d++){
      fail ^= rc[d];
    }
    return ~fail;

}
