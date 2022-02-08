#ifndef MASKED_SYM_FIPS202
#define MASKED_SYM_FIPS202

void masked_shake256_prf(unsigned char *output, size_t outlen, const unsigned char *key, unsigned char nonce);

#endif
