#ifndef MASKED_FIPS202_H
#define MASKED_FIPS202_H

#include <stddef.h>
#include <stdint.h>

// Copied from mupq/common/fips202.{h,c}
// Kept only functions that have to be masked in kyber768.


/* One-stop SHAKE256 call */
void masked_shake256(
        uint8_t *output, size_t outlen, size_t out_msk_stride, size_t out_data_stride,
        const uint8_t *input, size_t inlen, size_t in_msk_stride, size_t in_data_stride
        );


/* One-stop SHA3-512 shop */
void masked_sha3_512(
        uint8_t *output, size_t out_msk_stride, size_t out_data_stride,
        const uint8_t *input, size_t inlen, size_t in_msk_stride, size_t in_data_stride
        );

#endif
