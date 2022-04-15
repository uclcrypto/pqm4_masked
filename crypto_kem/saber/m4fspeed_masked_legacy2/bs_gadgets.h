#include <stdint.h>
#ifndef __BSGADGETS__
#define __BSGADGETS__
void SecAnd(uint32_t *x,
        uint32_t *y,
        uint32_t *z,
        uint32_t d);
void RefreshXOR(uint32_t *x,
        uint32_t d);
void FullRefreshXOR(uint32_t *x,
        uint32_t d);
void SecANDbs(uint32_t *z,
        const uint32_t *a,
        const uint32_t *b);
void SecXORbs(uint32_t *z, 
        const uint32_t *a, 
        const uint32_t *b);
#endif
