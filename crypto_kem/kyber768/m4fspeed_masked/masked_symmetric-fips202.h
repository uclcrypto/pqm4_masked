#ifndef MASKED_SYM_FIPS202
#define MASKED_SYM_FIPS202

void masked_shake256_prf(unsigned char *output, size_t outlen,
                         size_t o_msk_stride, size_t o_data_stride,
                         const unsigned char *key, size_t k_msk_stride,
                         size_t k_data_stride, unsigned char nonce);

#endif
