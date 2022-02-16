#include "SABER_indcpa.h"
#include "randombytes.h"
#include "fips202.h"
#include "poly.h"
#include "masked_poly.h"
#include <string.h>
#include <stdlib.h>

uint8_t masked_indcpa_kem_enc_cmp(const uint8_t m[SABER_KEYBYTES], const uint8_t seed_s[SABER_NOISE_SEEDBYTES], const uint8_t pk[SABER_INDCPA_PUBLICKEYBYTES], const uint8_t ciphertext[SABER_BYTES_CCA_DEC]){

    uint32_t fail = 0;
    const uint8_t *seed_A = pk + SABER_POLYVECCOMPRESSEDBYTES;
    const uint8_t *ct0 = ciphertext;
    const uint8_t *ct1 = ciphertext + SABER_POLYVECCOMPRESSEDBYTES;

    fail = masked_MatrixVectorMulEncNTT_cmp((uint8_t*)ct0, (uint8_t*)ct1, seed_s, seed_A, pk, m);
    return (uint8_t)fail;
}

void masked_indcpa_kem_dec(const uint8_t sk[SABER_INDCPA_SECRETKEYBYTES], const uint8_t ciphertext[SABER_BYTES_CCA_DEC], uint8_t m[SABER_KEYBYTES]){

    masked_InnerProdDecNTT(m, ciphertext, sk); // m <- Pack(Round(b'*s - cm))

}
