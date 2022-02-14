#ifndef MCBD_H
#define MCBD_H

#include "SABER_params.h"
#include "masked.h"
#include <stdint.h>

void masked_cbd_seed(size_t nshares,
                uint16_t *s, 
                size_t s_msk_stride,
                size_t s_data_stride,
                const uint8_t *buf,
                size_t buf_msk_stride,
                size_t buf_data_stride);

#endif
