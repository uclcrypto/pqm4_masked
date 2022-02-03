#ifndef MASKEDPOLY_H
#define MASKEDPOLY_H
#include "params.h"
#include "masked.h"

void masked_poly_ntt(StrAPoly r);
void masked_poly_invntt(StrAPoly r);
void masked_poly_tomsg(unsigned char *m, StrAPoly str_r);
void masked_poly_cmp(
        size_t c,
        uint32_t *rc,
        const StrAPoly mp,
        const poly *ref);

void finalize_cmp(uint32_t *bits);
#endif
