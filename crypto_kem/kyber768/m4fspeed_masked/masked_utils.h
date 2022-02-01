#ifndef MASKED_UTILS_H
#define MASKED_UTILS_H
#include <stdint.h>
#include "masked.h"
#include "poly.h"

inline uint32_t get_random();
inline uint32_t rand32();
inline void rand_q(uint16_t v[2]);
void masked_poly(StrAPoly mp, const poly *p);
void unmasked_poly(poly *p, const StrAPoly mp);
#endif
