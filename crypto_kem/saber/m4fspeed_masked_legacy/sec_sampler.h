#include <stdint.h>
#include <stddef.h>

#ifndef __SECSAMPLES__
#define __SECSAMPLES__
void secbitadd(uint32_t *z, const uint32_t *x, const size_t lambda, const size_t kappa);
void sec_sampler2(uint32_t *a, const uint32_t *x, 
                        const uint32_t *y,
                        const uint32_t kappa,const uint32_t kmodulus);

void sec_b2a_q(uint32_t *a, uint32_t *x,size_t k, uint32_t kb);
void refresh_add(uint32_t *a,uint32_t kb);
void sec_b2a_qbit(uint32_t *a, uint32_t *x, uint32_t kb);
void b2a_qbit(uint32_t *a,uint32_t *x, const uint32_t kb);
void sec_b2a_qbit_n(uint32_t *c,const uint32_t *a,const uint32_t xn,const uint32_t n, const uint32_t kb);
void secconstadd2(uint32_t *z);
void secconstadd4(uint32_t *z);
void secbitsub(uint32_t *z, 
        const uint32_t *x, 
        const size_t lambda, const size_t kappa);
#endif
