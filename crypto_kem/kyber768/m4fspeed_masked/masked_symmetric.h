#ifndef MASKED_SYMMETRIC_H
#define MASKED_SYMMETRIC_H

#include "masked_fips202.h"
#include "masked_symmetric-fips202.h"
#include "params.h"
#include <stddef.h>

//void kyber_shake128_absorb(shake128ctx *s, const unsigned char *input, unsigned char x, unsigned char y);
//void kyber_shake128_squeezeblocks(unsigned char *output, size_t nblocks, shake128ctx *s);
//void masked_shake256_prf(unsigned char *output, size_t outlen, const unsigned char *key, unsigned char nonce);
//void sha3_512_masked_check(uint8_t *output, const uint8_t *input, size_t inlen);

//#define hash_h(OUT, IN, INBYTES) sha3_256(OUT, IN, INBYTES)
//#define hash_g(OUT, IN, INBYTES) sha3_512_masked_check(OUT, IN, INBYTES)
//#define xof_absorb(STATE, IN, X, Y) kyber_shake128_absorb(STATE, IN, X, Y)
//#define xof_squeezeblocks(OUT, OUTBLOCKS, STATE) kyber_shake128_squeezeblocks(OUT, OUTBLOCKS, STATE)
//#define masked_prf(OUT, OUTBYTES, KEY, NONCE) masked_shake256_prf(OUT, OUTBYTES, KEY, NONCE)
//#define kdf(OUT, IN, INBYTES) shake256(OUT, KYBER_SSBYTES, IN, INBYTES)


#endif /* SYMMETRIC_H */
