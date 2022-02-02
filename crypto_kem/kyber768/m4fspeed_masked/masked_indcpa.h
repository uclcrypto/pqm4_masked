#ifndef MASKED_INDCPA_H
#define MASKED_INDCPA_H

unsigned char masked_indcpa_enc_cmp(const unsigned char *ct,
                             const unsigned char *m,
                             const unsigned char *pk,
                             const unsigned char *coins);

#endif
