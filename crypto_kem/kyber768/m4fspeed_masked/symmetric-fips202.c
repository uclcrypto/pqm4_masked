#include "fips202.h"
#include "symmetric.h"

#include <stdlib.h>

/*************************************************
 * Name:        kyber_shake128_absorb
 *
 * Description: Absorb step of the SHAKE128 specialized for the Kyber context.
 *
 * Arguments:   - shake128ctx *s:                  pointer to (uninitialized)
 *output Keccak state
 *              - const unsigned char *input:      pointer to KYBER_SYMBYTES
 *input to be absorbed into s
 *              - unsigned char i                  additional byte of input
 *              - unsigned char j                  additional byte of input
 **************************************************/
void kyber_shake128_absorb(shake128ctx *s, const unsigned char *input,
                           unsigned char x, unsigned char y) {
  unsigned char extseed[KYBER_SYMBYTES + 2];
  int i;

  for (i = 0; i < KYBER_SYMBYTES; i++) {
    extseed[i] = input[i];
  }
  extseed[i++] = x;
  extseed[i] = y;
  shake128_absorb(s, extseed, KYBER_SYMBYTES + 2);
}

/*************************************************
 * Name:        kyber_shake128_squeezeblocks
 *
 * Description: Squeeze step of SHAKE128 XOF. Squeezes full blocks of
 *SHAKE128_RATE bytes each. Modifies the state. Can be called multiple times to
 *keep squeezing, i.e., is incremental.
 *
 * Arguments:   - unsigned char *output:      pointer to output blocks
 *              - size_t nblocks:             number of blocks to be squeezed
 *(written to output)
 *              - shake128ctx *s:            pointer to in/output Keccak state
 **************************************************/
void kyber_shake128_squeezeblocks(unsigned char *output, size_t nblocks,
                                  shake128ctx *s) {
  shake128_squeezeblocks(output, nblocks, s);
}

/*************************************************
 * Name:        shake256_prf
 *
 * Description: Usage of SHAKE256 as a PRF, concatenates secret and public input
 *              and then generates outlen bytes of SHAKE256 output
 *
 * Arguments:   - unsigned char *output:      pointer to output
 *              - size_t outlen:              number of requested output bytes
 *              - const unsigned char * key:  pointer to the key (of length
 *KYBER_SYMBYTES)
 *              - const unsigned char nonce:  single-byte nonce (public PRF
 *input)
 **************************************************/
#include "masked.h"
#include "masked_fips202.h"
#include "masked_utils.h"

void shake256_prf(unsigned char *output, size_t outlen,
                  const unsigned char *key, unsigned char nonce) {
  unsigned char extkey[KYBER_SYMBYTES + 1];
  size_t i;

  for (i = 0; i < KYBER_SYMBYTES; i++) {
    extkey[i] = key[i];
  }
  extkey[i] = nonce;

  shake256(output, outlen, extkey, KYBER_SYMBYTES + 1);
}

void sha3_512_masked_check(uint8_t *output, const uint8_t *input,
                           size_t inlen) {
  uint8_t *masked_input = malloc(inlen * NSHARES);
  for (size_t i = 0; i < inlen; i++) {
    masked_input[i] = input[i];
    for (size_t j = 1; j < NSHARES; j++) {
      unsigned char rnd = get_random() & 0xFF;
      masked_input[i] ^= rnd;
      masked_input[i + j * inlen] = rnd;
    }
  }
  sha3_512(output, input, inlen);
  uint8_t masked_output[NSHARES * 64];
  masked_sha3_512(masked_output, 64, 1, masked_input, inlen, inlen, 1);
  for (size_t i = 0; i < 64; i++) {
    unsigned char o = output[i];
    output[i] = 0;
    for (size_t j = 0; j < NSHARES; j++) {
      output[i] ^= masked_output[i + j * 64];
    }
    if (o != output[i]) {
      BAIL("ERROR sha3_512 out %i %x %x", i, o, output[i]);
    }
  }
  free(masked_input);
}
