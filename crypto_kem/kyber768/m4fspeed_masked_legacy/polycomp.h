#ifndef __POLYCOMP__
#define __POLYCOMP__
#include "gadgets_legacy.h"
#include "params.h"
#include "polycomp.h"
#include <stdint.h>
#define D NSHARES
typedef struct u32Masked {uint32_t shares[KYBER_MASKING_ORDER+1];} u32Masked;

uint32_t GoubinAB(uint32_t A,uint32_t r,int k);
void convert_A2B_CGV14_32bits(u32Masked* x, u32Masked* y, unsigned k);
int zero_test_mod_mult(Masked* a, int q);
int zero_testing_prime_multi(Masked* ppoly, int q, const int SIZE);
int zero_test_poly_mul(Masked* ppoly, int q, int lambda, const int SIZE);
int zero_test_poly_mul_with_reduction(Masked* ppoly, int q, int kappa, const int SIZE);
void boolean_zero_test(Masked* x, Masked* y, int k, int logk);
void bool_poly_zero_test(Masked* ppoly, Masked* b, int k, int logk, const int SIZE);
void high_order_compress(Masked* x, Masked* y, unsigned q, unsigned k, unsigned ell);
int kyber_poly_comp_hybrid(Masked* mmasked_poly, uint16_t* ppoly);

void ConvertAB(uint32_t *A,uint32_t *z,int k,int n);
void range_compare(Masked* x, Masked* z, unsigned c, unsigned k, unsigned q);

#endif
