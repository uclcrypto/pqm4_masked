#include "SABER_indcpa.h"
#include "api.h"
#include "fips202.h"
#include "masked.h"
#include "masked_SABER_indcpa.h"
#include "masked_fips202.h"
#include "masked_poly.h"
#include "masked_utils.h"
#include "pack_unpack.h"
#include "randombytes.h"
#include "verify.h"
#include <string.h>

int crypto_kem_keypair(uint8_t *pk, uint8_t *sk) {
  indcpa_kem_keypair(pk, sk); // sk[0:SABER_INDCPA_SECRETKEYBYTES-1] <-- sk

  memcpy(sk + SABER_INDCPA_SECRETKEYBYTES, pk,
         SABER_INDCPA_PUBLICKEYBYTES); // sk[SABER_INDCPA_SECRETKEYBYTES:SABER_INDCPA_SECRETKEYBYTES+SABER_INDCPA_SECRETKEYBYTES-1]
                                       // <-- pk

  sha3_256(sk + SABER_SECRETKEYBYTES - 64, pk,
           SABER_INDCPA_PUBLICKEYBYTES); // Then hash(pk) is appended.

  randombytes(
      sk + SABER_SECRETKEYBYTES - SABER_KEYBYTES,
      SABER_KEYBYTES); // Remaining part of sk contains a pseudo-random number,
                       // this is output when check in crypto_kem_dec() fails.

  return (0);
}

int crypto_kem_enc(uint8_t *c, uint8_t *k, const uint8_t *pk) {
  uint8_t kr[64]; // Will contain key, coins
  uint8_t buf[64];

  randombytes(buf, 32);

  sha3_256(buf, buf,
           32); // BUF[0:31] <-- random message (will be used as the key for
                // client) Note: hash doesnot release system RNG output

  sha3_256(buf + 32, pk,
           SABER_INDCPA_PUBLICKEYBYTES); // BUF[32:63] <-- Hash(public key);
                                         // Multitarget countermeasure for coins
                                         // + contributory KEM

  sha3_512(kr, buf, 64); // kr[0:63] <-- Hash(buf[0:63]), K^ <-- kr[0:31],
                         // noiseseed (r) <-- kr[32:63]

  indcpa_kem_enc(
      buf, kr + 32, pk,
      c); // buf[0:31] contains message; kr[32:63] contains randomness r;

  sha3_256(kr + 32, c, SABER_BYTES_CCA_DEC);

  sha3_256(k, kr, 64); // hash concatenation of pre-k and h(c) to k

  return (0);
}

int crypto_kem_dec(uint8_t *k, const uint8_t *c, const uint8_t *sk) {
  uint8_t fail;
  uint8_t masked_buf[64 * NSHARES];

  uint8_t kr[64];                  // Will contain key, coins
  uint8_t masked_kr[64 * NSHARES]; // Will contain key, coins
  const uint8_t *pk = sk + SABER_INDCPA_SECRETKEYBYTES;
  const uint8_t *hpk =
      sk + SABER_SECRETKEYBYTES - 64; // Save hash by storing h(pk) in sk
  StrAPolyVec masked_sk;

  // TODO Store the secret key masked.
  for (size_t i = 0; i < SABER_L; i++) {
#ifdef SABER_COMPRESS_SECRETKEY
    BS2POLmu(sk + i * SABER_POLYSECRETBYTES, masked_sk[i][0]);
#else
    BS2POLq(sk + i * SABER_POLYSECRETBYTES, masked_sk[i][0]);
#endif
    mask_poly_inplace(masked_sk[i], SABER_P);
  }

  // masked_sk has no const flag.
  masked_indcpa_kem_dec(masked_sk, c, masked_buf, 64,
                        1); // buf[0:31] <-- message

  // masking the hpk
  memcpy(masked_buf + 32, hpk, 32);
  for (size_t d = 1; d < NSHARES; d++) {
    memset(masked_buf + 32 + d * 64, 0, 32);
  }
  // Hash decrypted to get coins
  masked_sha3_512(masked_kr, 64, 1, masked_buf, 64, 64, 1);

  // Compare with re-encrypted
  fail = masked_indcpa_kem_enc_cmp(
      masked_buf, 64, 1, &masked_kr[32], 64, 1, pk,
      c); // in-place verification of the re-encryption

  // Not masked: H(c)
  sha3_256(kr + 32, c, SABER_BYTES_CCA_DEC); // overwrite coins in kr with h(c)

  // No need to mask this, like for Kyber
  // -> unmask K in kr
  memset(kr, 0, 32);
  for (size_t d = 0; d < NSHARES; d++) {
    for (size_t i = 0; i < 32; i++) {
      kr[i] ^= masked_kr[d * 64 + i];
    }
  }
  cmov(kr, sk + SABER_SECRETKEYBYTES - SABER_KEYBYTES, SABER_KEYBYTES, fail);

  // No need to mask this, like for Kyber
  sha3_256(k, kr, 64); // hash concatenation of pre-k and h(c) to k

  return (0);
}
