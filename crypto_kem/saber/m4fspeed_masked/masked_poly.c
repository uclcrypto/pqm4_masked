#include <string.h>
#include "masked_poly.h"
#include "cbd.h"
#include "mNTT.h"
#include "fips202.h"
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

#if 0
static inline shake128incctx shake128_absorb_seed(const uint8_t seed[SABER_SEEDBYTES]){

    shake128incctx ctx;
    shake128_inc_init(&ctx);
    shake128_inc_absorb(&ctx, seed, SABER_SEEDBYTES);
    shake128_inc_finalize(&ctx);

    return ctx;

}

uint32_t masked_MatrixVectorMulEncNTT(uint8_t ct0[SABER_POLYVECCOMPRESSEDBYTES], uint8_t ct1[SABER_SCALEBYTES_KEM], const uint8_t seed_s[SABER_NOISE_SEEDBYTES], const uint8_t seed_A[SABER_SEEDBYTES], const uint8_t pk[SABER_INDCPA_PUBLICKEYBYTES], const uint8_t m[SABER_KEYBYTES], int compare){

    uint32_t acc_NTT[SABER_N];
    uint32_t A_NTT[SABER_N];
    uint32_t s_NTT[SABER_L * SABER_N];

    uint16_t poly[SABER_N];
    uint16_t acc[SABER_N];

    uint8_t shake_out[MAX(SABER_POLYBYTES, SABER_POLYCOINBYTES)];

    uint16_t *mp = poly;

    size_t i, j;
    uint32_t fail = 0;

    shake128incctx shake_s_ctx = shake128_absorb_seed(seed_s);

    for(i = 0; i < SABER_L; i++){
        shake128_inc_squeeze(shake_out, SABER_POLYCOINBYTES, &shake_s_ctx);
        cbd(poly, shake_out);
        NTT_forward_32(s_NTT + i * SABER_N, poly);
    }

    shake128_inc_ctx_release(&shake_s_ctx);

    shake128incctx shake_A_ctx = shake128_absorb_seed(seed_A);

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
            acc[j] = ((acc[j] + h1) >> (SABER_EQ - SABER_EP));
        }

        if (compare) {
            fail |= POLp2BS_cmp(ct0 + i * SABER_POLYCOMPRESSEDBYTES, acc);
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
        acc[j] = (acc[j] - (mp[j] << (SABER_EP - 1)) + h1) >> (SABER_EP - SABER_ET);
    }

    if(compare){
        fail |= POLT2BS_cmp(ct1, acc);
    }else{
        POLT2BS(ct1, acc);
    }

    return fail;

}
#endif

void masked_InnerProdDecNTT(uint8_t m[SABER_KEYBYTES], const uint8_t ciphertext[SABER_BYTES_CCA_DEC], const uint8_t sk[SABER_INDCPA_SECRETKEYBYTES]){

    // TODO Store the secret key masked.
    StrAPolyVec sk_masked;
    for (size_t i = 0; i < SABER_L; i++) {
#ifdef SABER_COMPRESS_SECRETKEY
        BS2POLmu(sk + i * SABER_POLYSECRETBYTES, sk_masked[i][0]);
#else
        BS2POLq(sk + i * SABER_POLYSECRETBYTES, sk_masked[i][0]);
#endif
        mask_poly_inplace(sk_masked[i], SABER_P);
    }

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
        m_poly[0][i] = (m_poly[0][i] + h2 - (cm[i] << (SABER_EP - SABER_ET))) % SABER_P;
    }

    Poly poly;
    for (size_t i = 0; i < SABER_N; i+=2*BSSIZE) {
        uint32_t masked_bs[NSHARES][2*SABER_EP];
        masked_dense2bitslice_opt(
                NSHARES, SABER_EP,
                masked_bs[0], 2*SABER_EP, 1,
                &m_poly[0][i], SABER_N, 1
                );
        seca2b(NSHARES, SABER_EP, &masked_bs[0][0], 2*SABER_EP, 1);
        seca2b(NSHARES, SABER_EP, &masked_bs[0][SABER_EP], 2*SABER_EP, 1);
        masked_bitslice2dense_opt(
                NSHARES, 1,
                &m_poly[0][i], SABER_N, 1,
                &masked_bs[0][SABER_EP-1], 2*SABER_EP, SABER_EP
                );
    }
    unmasked_poly(poly, m_poly, SABER_P);
    POLmsg2BS(m, poly);
}




