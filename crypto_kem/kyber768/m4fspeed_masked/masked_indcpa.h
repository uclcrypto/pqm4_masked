#ifndef MASKED_INDCPA_H
#define MASKED_INDCPA_H

#include <stddef.h>

unsigned char masked_indcpa_enc_cmp(const unsigned char *c,
                             const unsigned char *m, size_t m_msk_stride, size_t m_data_stride,
                             const unsigned char *pk,
                             const unsigned char *coins,size_t coins_msk_stride, size_t coins_data_stride);

void masked_indcpa_dec(unsigned char *m,
                size_t o_msk_stride,
                size_t o_data_stride,
                const unsigned char *c,
                const unsigned char *sk);
#endif
