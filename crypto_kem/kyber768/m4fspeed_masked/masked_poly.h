#ifndef MASKEDPOLY_H
#define MASKEDPOLY_H
#include "params.h"
#include "masked.h"

void masked_poly_ntt(StrAPoly r);
void masked_poly_invntt(StrAPoly r);
void masked_poly_tomsg(unsigned char *m, StrAPoly str_r);
#endif
