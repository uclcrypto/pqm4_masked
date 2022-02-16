#include <string.h>
#include "masked_poly.h"
#include "masked_cbd.h"
#include "cbd.h"
#include "mNTT.h"
//#include "NTT.h"
#include "fips202.h"
#include "masked_fips202.h"
#include "pack_unpack.h"
#include "masked_utils.h"
#include "masked_representations.h"
#include "gadgets.h"
#include "masked.h"

#define h1 (1 << (SABER_EQ - SABER_EP - 1))
#define h2 ((1 << (SABER_EP - 2)) - (1 << (SABER_EP - SABER_ET - 1)) + (1 << (SABER_EQ - SABER_EP - 1)))
#define MAX(a,b) (((a)>(b))?(a):(b))

extern void __asm_poly_add_16(uint16_t *des, uint16_t *src1, uint16_t *src2);
extern void __asm_poly_add_32(uint32_t *des, uint32_t *src1, uint32_t *src2);

static inline shake128incctx shake128_absorb_seed(const uint8_t seed[SABER_SEEDBYTES]){

    shake128incctx ctx;
    shake128_inc_init(&ctx);
    shake128_inc_absorb(&ctx, seed, SABER_SEEDBYTES);
    shake128_inc_finalize(&ctx);
    return ctx;
}


void masked_InnerProdDecNTT(uint8_t m[SABER_KEYBYTES],
    const uint8_t ciphertext[SABER_BYTES_CCA_DEC], const StrAPolyVec sk_masked){

    // NTT of ciphertext
    uint32_t c_NTT_32[SABER_L][SABER_N];
    uint16_t c_NTT_16[SABER_L][SABER_N];
    for (size_t i = 0; i < SABER_L; i++) {
        uint16_t poly[SABER_N];
        BS2POLp(ciphertext + i * SABER_POLYCOMPRESSEDBYTES, poly);
        NTT_forward_32(c_NTT_32[i], poly);
        NTT_forward_16(c_NTT_16[i], poly);
    }

    // Ciphertex * sk
    StrAPoly m_poly;
    for (size_t j=0; j<NSHARES; j++) {
        uint32_t acc_NTT_32[SABER_N];
        uint16_t acc_NTT_16[SABER_N];
        for (size_t i = 0; i < SABER_L; i++) {
            uint32_t poly_NTT_32[SABER_N];
            uint16_t poly_NTT_16[SABER_N];
            NTT_forward_32(poly_NTT_32, sk_masked[i][j]);
            NTT_forward_16(poly_NTT_16, sk_masked[i][j]);
            if (i == 0) {
                NTT_mul_32(acc_NTT_32, poly_NTT_32, c_NTT_32[i]);
                NTT_mul_16(acc_NTT_16, poly_NTT_16, c_NTT_16[i]);
            } else {
                NTT_mul_acc_32(acc_NTT_32, poly_NTT_32, c_NTT_32[i]);
                NTT_mul_acc_16(acc_NTT_16, poly_NTT_16, c_NTT_16[i]);
            }
        }
        NTT_inv_32(acc_NTT_32);
        NTT_inv_16(acc_NTT_16);
        solv_CRT(m_poly[j], acc_NTT_32, acc_NTT_16);
    }

    uint16_t cm[SABER_N];
    BS2POLT(ciphertext + SABER_POLYVECCOMPRESSEDBYTES, cm);
    for (size_t i = 0; i < SABER_N; i++) {
        m_poly[0][i] = (SABER_P + m_poly[0][i] + h2 - (cm[i] << (SABER_EP - SABER_ET))) % SABER_P;
    }

    for (size_t i = 0; i < SABER_N; i+=2*BSSIZE) {
        uint32_t masked_bs[NSHARES*SABER_EP*2];
        masked_dense2bitslice_opt(
                NSHARES, SABER_EP,
                masked_bs, 1, NSHARES,
                &m_poly[0][i], SABER_N, 1
                );
        seca2b(NSHARES, SABER_EP, masked_bs, 1, NSHARES);
        seca2b(NSHARES, SABER_EP, &masked_bs[NSHARES*SABER_EP], 1, NSHARES);
        masked_bitslice2dense_opt(
                NSHARES, 1,
                &m_poly[0][i], SABER_N, 1,
                &masked_bs[(SABER_EP-1)*NSHARES], 1, NSHARES*SABER_EP);    
    }
    
    Poly poly;
    unmasked_poly(poly,m_poly,2);
   
    POLmsg2BS(m, poly);
}

uint32_t masked_MatrixVectorMulEncNTT_cmp(uint8_t ct0[SABER_POLYVECCOMPRESSEDBYTES], 
                uint8_t ct1[SABER_SCALEBYTES_KEM], 
                const uint8_t seed_s[SABER_NOISE_SEEDBYTES*NSHARES], 
                const uint8_t seed_A[SABER_SEEDBYTES], 
                const uint8_t pk[SABER_INDCPA_PUBLICKEYBYTES], 
                const uint8_t m[SABER_KEYBYTES]){

    uint32_t acc_NTT_32[NSHARES][SABER_N];
    uint32_t A_NTT_32[SABER_N];
    uint32_t s_NTT_32[NSHARES][SABER_L * SABER_N];

    uint16_t acc_NTT_16[NSHARES][SABER_N];
    uint16_t A_NTT_16[SABER_N];
    uint16_t s_NTT_16[NSHARES][SABER_L * SABER_N];

    uint16_t poly[SABER_N];
    uint16_t poly_ref[SABER_N];
    uint16_t m_poly[NSHARES][SABER_N];
    uint16_t acc[SABER_N];
    uint16_t myref[SABER_N];
    uint16_t masked_acc[NSHARES][SABER_N];

    uint8_t shake_out[MAX(SABER_POLYBYTES, SABER_POLYCOINBYTES)];
    uint8_t masked_shake_out[MAX(SABER_POLYBYTES, SABER_POLYCOINBYTES)*NSHARES];

    uint16_t *mp = poly;

    size_t i, j;
    uint32_t fail = 0;

    uint8_t masked_seed_s[NSHARES*SABER_SEEDBYTES];
    memset(masked_seed_s,0,sizeof(masked_seed_s));
    memcpy(masked_seed_s,seed_s,SABER_SEEDBYTES);

    MaskedShakeCtx masked_shake_s_ctx;
    masked_shake128_inc_init(&masked_shake_s_ctx,
        masked_seed_s,SABER_SEEDBYTES,SABER_SEEDBYTES,1);

    for(i = 0; i < SABER_L; i++){
        masked_shake128_squeeze(&masked_shake_s_ctx,
            masked_shake_out,SABER_POLYCOINBYTES,SABER_POLYCOINBYTES,1);
      
        masked_cbd_seed(NSHARES,
                    m_poly,SABER_N,1,
                    masked_shake_out,SABER_POLYCOINBYTES,1);
        
#ifdef MY_DEBUG
        char buf_x [128];
        hal_send_str("----------");
        for(int n = 0; n<SABER_N;n++){
          sprintf(buf_x,"n%d- > %d %d",n,(uint16_t) poly[n],(uint16_t) poly_ref[n]%SABER_Q);
          hal_send_str(buf_x);
        }
#endif
        for(int d = 0; d<NSHARES; d++){
          NTT_forward_32(s_NTT_32[d] + i * SABER_N, m_poly[d]);
          NTT_forward_16(s_NTT_16[d] + i * SABER_N, m_poly[d]);
        }
    }

    shake128incctx shake_A_ctx = shake128_absorb_seed(seed_A);

    uint32_t rc[NSHARES];
    memset(rc,0,sizeof(rc));
    rc[0] = 0xFFFFFFFF;

    for (i = 0; i < SABER_L; i++) {

        for (j = 0; j < SABER_L; j++) {

            shake128_inc_squeeze(shake_out, SABER_POLYBYTES, &shake_A_ctx);
            BS2POLq(shake_out, poly);
            
            NTT_forward_32(A_NTT_32, poly);
            NTT_forward_16(A_NTT_16, poly);

            // TODO
            for(int d = 0; d<NSHARES; d++){
              if (j == 0) {
                  NTT_mul_32(acc_NTT_32[d], A_NTT_32, s_NTT_32[d] + j * SABER_N);
                  NTT_mul_16(acc_NTT_16[d], A_NTT_16, s_NTT_16[d] + j * SABER_N);
              } else {
                  NTT_mul_acc_32(acc_NTT_32[d], A_NTT_32, s_NTT_32[d] + j * SABER_N);
                  NTT_mul_acc_16(acc_NTT_16[d], A_NTT_16, s_NTT_16[d] + j * SABER_N);
              }
            }
        }

        for(int d = 0; d<NSHARES; d ++){
          NTT_inv_32(acc_NTT_32[d]);
          NTT_inv_16(acc_NTT_16[d]);
          solv_CRT(masked_acc[d], acc_NTT_32[d], acc_NTT_16[d]);
        }

        for (j = 0; j < SABER_N; j++) {
            masked_acc[0][j] = (masked_acc[0][j] + h1)%SABER_Q;
        }
       
        // compare acc with ct0
        BS2POLp(ct0 + i*SABER_POLYCOMPRESSEDBYTES,myref);
        for(j=0;j<SABER_N;j++){
          myref[j] = myref[j] % SABER_P;
        }
        masked_poly_cmp(SABER_EQ-SABER_EP,SABER_EQ,SABER_EQ,rc,masked_acc,myref);
    }

    shake128_inc_ctx_release(&shake_A_ctx);

    for(j = 0; j < SABER_L; j++){

        BS2POLp(pk + j * SABER_POLYCOMPRESSEDBYTES, poly);

        NTT_forward_32(A_NTT_32, poly);
        NTT_forward_16(A_NTT_16, poly);

        for(int d =0; d<NSHARES; d++){
          if(j == 0){
              NTT_mul_32(acc_NTT_32[d], A_NTT_32, s_NTT_32[d] + j * SABER_N);
              NTT_mul_16(acc_NTT_16[d], A_NTT_16, s_NTT_16[d] + j * SABER_N);
          }else{
              NTT_mul_acc_32(acc_NTT_32[d], A_NTT_32, s_NTT_32[d] + j * SABER_N);
              NTT_mul_acc_16(acc_NTT_16[d], A_NTT_16, s_NTT_16[d] + j * SABER_N);
          }
        }
    }

    for(int d = 0; d<NSHARES; d ++){
      NTT_inv_32(acc_NTT_32[d]);
      NTT_inv_16(acc_NTT_16[d]);
      solv_CRT(masked_acc[d], acc_NTT_32[d], acc_NTT_16[d]);
    }

    BS2POLmsg(m, mp);
    for (j = 0; j < SABER_N; j++) {
        // work in SABER_Q as for NTT. Could be done in SABER_P.
        masked_acc[0][j] = (masked_acc[0][j] - (mp[j] << (SABER_EP-1)) + h1)%SABER_Q;
    }

    // compare acc with ct1
    BS2POLT(ct1,myref);
    for(j=0;j<SABER_N;j++){
        myref[j] = myref[j] % (1<<SABER_ET);
    }           
    masked_poly_cmp(SABER_EP-SABER_ET,SABER_EP,SABER_EP,rc,masked_acc,myref);
 
    // finalize the comparison
    finalize_cmp(rc);
    fail = 0;
    for(int d=0;d<NSHARES;d++){
      fail ^= rc[d];
    }
    return !fail;
}


void masked_poly_cmp(size_t b_start, size_t b_end, size_t coeffs_size, uint32_t *rc, const uint16_t *mp,
                     int16_t *ref) {

  size_t i, b;
  uint32_t bits[2* NSHARES * coeffs_size];
  uint32_t bits_ref[2*coeffs_size];

  for (i = 0; i < SABER_N; i += BSSIZE*2) {

    // convert masked poly
    masked_dense2bitslice_opt(NSHARES,coeffs_size,
        bits,1,NSHARES,
        mp,SABER_N,1);

    seca2b(NSHARES,  coeffs_size, bits, 1, NSHARES);
    seca2b(NSHARES,  coeffs_size, &bits[NSHARES*coeffs_size], 1, NSHARES);

    // map public polynomial to bitslice
    masked_dense2bitslice_opt(1, coeffs_size, bits_ref, 1, 1, ref, 1,
                          1);

    for (b = 0; b < b_end-b_start; b++){

      // public polynomial and public one
      bits[(b+b_start) * NSHARES] ^= bits_ref[b] ^ 0xFFFFFFFF;
      masked_and(NSHARES, rc, 1, rc, 1, &bits[(b+b_start) * NSHARES], 1);

      bits[(b+b_start+coeffs_size) * NSHARES] ^= bits_ref[b+coeffs_size] ^ 0xFFFFFFFF;
      masked_and(NSHARES, rc, 1, rc, 1, &bits[(b+b_start+coeffs_size) * NSHARES], 1);
    }
  }
}
void finalize_cmp(uint32_t *bits) {

  uint32_t other[NSHARES];
  int d;
  for (d = 0; d < NSHARES; d++) {
    other[d] = bits[d] >> 16;
  }
  masked_and(NSHARES, bits, 1, bits, 1, other, 1);

  for (d = 0; d < NSHARES; d++) {
    other[d] = bits[d] >> 8;
  }
  masked_and(NSHARES, bits, 1, bits, 1, other, 1);

  for (d = 0; d < NSHARES; d++) {
    other[d] = bits[d] >> 4;
  }
  masked_and(NSHARES, bits, 1, bits, 1, other, 1);

  for (d = 0; d < NSHARES; d++) {
    other[d] = bits[d] >> 2;
  }
  masked_and(NSHARES, bits, 1, bits, 1, other, 1);
  for (d = 0; d < NSHARES; d++) {
    other[d] = bits[d] >> 1;
  }
  masked_and(NSHARES, bits, 1, bits, 1, other, 1);
}
