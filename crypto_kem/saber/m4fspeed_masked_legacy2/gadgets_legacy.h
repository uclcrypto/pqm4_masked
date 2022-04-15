#ifndef __GADGETS_LEGACY__H
#define __GADGETS_LEGACY__H
void my_SecAdd(uint32_t *x,
        uint32_t *y,
        uint32_t *z,
        uint32_t d,
        uint32_t k);
void SecA2BModpow2(uint32_t *a,
        uint32_t *out,
        uint32_t d,
        uint32_t k);
#endif
