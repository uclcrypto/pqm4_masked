#ifndef MASKED_H
#define MASKED_H
#include <stdint.h>
#include "params.h"

#define NSHARES 3

typedef poly maskedpoly[NSHARES];

typedef struct {
    maskedpoly vec[KYBER_K];
} maskedpolyvec;

#endif
