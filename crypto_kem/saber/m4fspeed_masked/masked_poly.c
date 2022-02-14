#include <string.h>
#include "masked_poly.h"
#include "masked_cbd.h"
#include "cbd.h"
#include "mNTT.h"
//#include "NTT.h"
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