#ifndef GADGETS_H
#define GADGETS_H

#include <stdint.h>
#include "masked.h"

void masked_xor(
        size_t nshares,
        uint32_t *ina, size_t ina_stride,
        uint32_t *inb, size_t inb_stride,
        uint32_t *out, size_t out_stride
        );
void masked_and(
        size_t nshares,
        uint32_t *z, size_t z_stride,
        uint32_t *a, size_t a_stride,
        uint32_t *b, size_t b_stride
        );
 

#endif // GADGETS_H
