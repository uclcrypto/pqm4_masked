#ifndef __GADGETS_LEGACY_H__
#define __GADGETS_LEGACY_H__

#define KYBER_MASKING_ORDER (NSHARES-1)
#define KYBER_Q 3329
typedef struct Masked {int32_t shares[KYBER_MASKING_ORDER+1];} Masked;


void linear_arithmetic_refresh(Masked* x, unsigned q);
void convert_2_l_to_1bit_bool(Masked* x, Masked* b, unsigned l);
void convert_B2A(Masked* x, Masked* y, unsigned k, unsigned q);
void modulus_switch(Masked* x, unsigned q, unsigned shift);
void kyber_decryption(Masked* x, Masked* b);
#endif
