#include "fips202.h"
#include "symmetric.h"
#include <stdlib.h>

/*************************************************
* Name:        shake256_prf
*
* Description: Usage of SHAKE256 as a PRF, concatenates secret and public input
*              and then generates outlen bytes of SHAKE256 output
*
* Arguments:   - unsigned char *output:      pointer to output
*              - size_t outlen:              number of requested output bytes
*              - const unsigned char * key:  pointer to the key (of length KYBER_SYMBYTES)
*              - const unsigned char nonce:  single-byte nonce (public PRF input)
**************************************************/
#include "masked.h"
#include "masked_utils.h"
#include "masked_fips202.h"
#include "masked_symmetric-fips202.h"

void masked_shake256_prf(unsigned char *output, size_t outlen, size_t o_msk_stride, size_t o_data_stride,
        const unsigned char *key, size_t k_msk_stride, size_t k_data_stride, unsigned char nonce) {

    unsigned char extkey[(KYBER_SYMBYTES + 1)*NSHARES];
    size_t i,d;

    for (d = 0; d < NSHARES; d++){
        for (i = 0; i < KYBER_SYMBYTES; i++) {
            extkey[d*(KYBER_SYMBYTES+1) + i] = key[d*k_msk_stride + i * k_data_stride];
        }
        extkey[d*(KYBER_SYMBYTES+1) + i] = (d == 0) ? nonce:0;
    }
    
    masked_shake256(output, outlen, o_msk_stride, o_data_stride, 
            extkey, KYBER_SYMBYTES+1, KYBER_SYMBYTES+1, 1);
}
