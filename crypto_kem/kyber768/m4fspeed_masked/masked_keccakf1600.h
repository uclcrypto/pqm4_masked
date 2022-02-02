#ifndef MASKED_KECCAKF1600_H
#define MASKED_KECCAKF1600_H

#include <stdint.h>
#include <stddef.h>
#include "masked.h"

#define KECCAK_NWORDS 25

typedef uint64_t KeccakState[KECCAK_NWORDS];
typedef KeccakState MaskedKeccakState[NSHARES];

void MaskedKeccakF1600_StateExtractBytes(
        const MaskedKeccakState state,
        unsigned char *data,
        unsigned int offset,
        unsigned int length,
        size_t data_msk_stride,
        size_t data_data_stride
        );
void MaskedKeccakF1600_StateXORPublicBytes(
        MaskedKeccakState state,
        const unsigned char *data,
        unsigned int offset,
        unsigned int length
        );
void MaskedKeccakF1600_StateXORBytes(
        MaskedKeccakState state,
        const unsigned char *data,
        unsigned int offset,
        unsigned int length,
        size_t data_msk_stride,
        size_t data_data_stride
        );
void MaskedKeccakF1600_StatePermute(
        MaskedKeccakState state
        );

#endif
