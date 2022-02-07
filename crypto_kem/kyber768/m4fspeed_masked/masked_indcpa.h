#ifndef MASKED_INDCPA_H
#define MASKED_INDCPA_H

unsigned char masked_indcpa_enc_cmp(const unsigned char *ct,
                             const unsigned char *m,
                             const unsigned char *pk,
                             const unsigned char *coins);

void masked_indcpa_dec(unsigned char *m,
                const unsigned char *c,
                const unsigned char *sk);
#endif
