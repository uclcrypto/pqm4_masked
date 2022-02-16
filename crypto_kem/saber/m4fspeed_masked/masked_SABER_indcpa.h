#ifndef MASKED_INDCPA_H
#define MASKED_INDCPA_H

#include "SABER_params.h"
#include "masked.h"
#include <stdint.h>

uint8_t masked_indcpa_kem_enc_cmp(const uint8_t m[SABER_KEYBYTES], const uint8_t seed_sp[SABER_NOISE_SEEDBYTES], const uint8_t pk[SABER_INDCPA_PUBLICKEYBYTES], const uint8_t ciphertext[SABER_BYTES_CCA_DEC]);

void masked_indcpa_kem_dec(const StrAPolyVec masked_sk, 
    const uint8_t ciphertext[SABER_BYTES_CCA_DEC], 
    uint8_t m[SABER_KEYBYTES]);
#endif
