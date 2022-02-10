#include "bench.h"
#include "masked_representations.h"
#include "masked.h"
#include <stdint.h>

/*************************************************
 * Name:        StrAPoly2APoly
 *
 * Description: Maps strided polynomial into a dense representation
 *
 * Arguments:
 *           - APoly out : dense polynomial
 *           - StrAPoly in: strided polynomial
 * **************************************************/
void StrAPoly2APoly(APoly out, const StrAPoly in) {
  int i, d;
  for (i = 0; i < KYBER_N; i++) {
    for (d = 0; d < NSHARES; d++) {
      out[i][d] = in[d][i];
    }
  }
}

/*************************************************
 * Name:        APoly2StrAPoly
 *
 * Description: Maps dense polynomial into a strided representation
 *
 * Arguments:
 *           - StrAPoly out: strided polynomial
 *           - APoly in : dense polynomial
 * **************************************************/
void APoly2StrAPoly(StrAPoly out, const APoly in) {

  int i, d;
  for (i = 0; i < KYBER_N; i++) {
    for (d = 0; d < NSHARES; d++) {
      out[d][i] = in[i][d];
    }
  }
}

/*************************************************
 * Name:        masked_dense2bitslice
 *
 * Description: maps a dense reprensentation to a bitlisce one
 *
 * Arguments:
 *           - uint32_t *bitslice[]: output bitslice representation. Table of
 * coeffs_size x nshares.
 *           - int16_t *dense[]: input dense reprensetation. Table of n_coeffs x
 * nshares.
 *           - size_t coeffs_size: number of bits to represent the dense
 * coefficients
 *           - size_t n_coeffs: number of coefficients
 *           - size_t nshares: number of shares
 * **************************************************/
void masked_dense2bitslice(size_t nshares, size_t n_coeffs, size_t coeffs_size,
                           uint32_t *bitslice, size_t bitslice_msk_stride,
                           size_t bitslice_data_stride, 
                           int16_t *dense,
                           size_t dense_msk_stride, size_t dense_data_stride) {

  start_bench(my_dense2bs);
  size_t d, c, b;
  for (b = 0; b < coeffs_size; b++) {
    for (d = 0; d < nshares; d++) {
      bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] = 0;
    }
  }
  
  if((n_coeffs & 0x3) == 0){
   for (d = 0; d < nshares; d++) {
      for (c = 0; c < n_coeffs; c+=4) {
        int16_t xd0 = dense[c * dense_data_stride + d * dense_msk_stride];
        int16_t xd1 = dense[(c+1) * dense_data_stride + d * dense_msk_stride];
        int16_t xd2 = dense[(c+2) * dense_data_stride + d * dense_msk_stride];
        int16_t xd3 = dense[(c+3) * dense_data_stride + d * dense_msk_stride];
        
        for (b = 0; b < coeffs_size; b++) {

          bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] |=
              (xd0 & 0x1) << (c+0);
          bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] |=
              (xd1 & 0x1) << (c+1);
          bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] |=
              (xd2 & 0x1) << (c+2);
          bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] |=
              (xd3 & 0x1) << (c+3);
          xd0 = xd0 >> 1;
          xd1 = xd1 >> 1;
          xd2 = xd2 >> 1;
          xd3 = xd3 >> 1;
        }
      }
    }

  }
  else{
    for (d = 0; d < nshares; d++) {
      for (c = 0; c < n_coeffs; c++) {
        int16_t xd = dense[c * dense_data_stride + d * dense_msk_stride];
        for (b = 0; b < coeffs_size; b++) {
          bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] |=
              (xd & 0x1) << c;
          xd = xd >> 1;
        }
      }
    }
  }
  stop_bench(my_dense2bs);
}
void masked_dense2bitslice_u32(size_t nshares, size_t n_coeffs, size_t coeffs_size,
                           uint32_t *bitslice, size_t bitslice_msk_stride,
                           size_t bitslice_data_stride, uint32_t *dense,
                           size_t dense_msk_stride, size_t dense_data_stride) {
  start_bench(my_dense2bs);
  size_t d, c, b;
  for (b = 0; b < coeffs_size; b++) {
    for (d = 0; d < nshares; d++) {
      bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] = 0;
    }
  }
  if((n_coeffs & 0x3) == 0){
   for (d = 0; d < nshares; d++) {
      for (c = 0; c < n_coeffs; c+=4) {
        uint32_t xd0 = dense[c * dense_data_stride + d * dense_msk_stride];
        uint32_t xd1 = dense[(c+1) * dense_data_stride + d * dense_msk_stride];
        uint32_t xd2 = dense[(c+2) * dense_data_stride + d * dense_msk_stride];
        uint32_t xd3 = dense[(c+3) * dense_data_stride + d * dense_msk_stride];
        
        for (b = 0; b < coeffs_size; b++) {

          bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] |=
              (xd0 & 0x1) << (c+0);
          bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] |=
              (xd1 & 0x1) << (c+1);
          bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] |=
              (xd2 & 0x1) << (c+2);
          bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] |=
              (xd3 & 0x1) << (c+3);

          xd0 = xd0 >> 1;
          xd1 = xd1 >> 1;
          xd2 = xd2 >> 1;
          xd3 = xd3 >> 1;
        }
      }
    }

  }
  else{
    for (d = 0; d < nshares; d++) {
      for (c = 0; c < n_coeffs; c++) {
        uint32_t xd = dense[c * dense_data_stride + d * dense_msk_stride];
        for (b = 0; b < coeffs_size; b++) {
          bitslice[b * bitslice_data_stride + d * bitslice_msk_stride] |=
              (xd & 0x1) << c;
          xd = xd >> 1;
        }
      }
    }
  }
  stop_bench(my_dense2bs);
}

/*************************************************
 * Name:        masked_bitslice2dense
 *
 * Description: maps a bitslice reprensentation to a dense one
 *
 * Arguments:
 *           - uint32_t *bitslice[]: input bitslice representation. Table of
 * coeffs_size x nshares.
 *           - int16_t *dense[]: output dense reprensetation. Table of n_coeffs
 * x nshares.
 *           - size_t coeffs_size: number of bits to represent the dense
 * coefficients
 *           - size_t n_coeffs: number of coefficients
 *           - size_t nshares: number of shares
 * **************************************************/
void masked_bitslice2dense(size_t nshares, size_t n_coeffs, size_t coeffs_size,
                           int16_t *dense, size_t dense_msk_stride,
                           size_t dense_data_stride, uint32_t *bitslice,
                           size_t bitslice_msk_stride,
                           size_t bitslice_data_stride) {

  start_bench(my_bs2dense);
  size_t d, c, b;
  for (c = 0; c < n_coeffs; c+=2) {
    for (d = 0; d < nshares; d++) {
      int16_t xd0 = 0;
      int16_t xd1 = 0;
      for (b = 0; b < coeffs_size; b++) {
        uint32_t y = bitslice[b*bitslice_data_stride + d* bitslice_msk_stride];
          xd0 |= ((y >> c) & 0x1)<<b;
          xd1 |= ((y >> (c+1)) & 0x1)<<b;
      }
      dense[c * dense_data_stride + d * dense_msk_stride] = xd0;
      dense[(c+1) * dense_data_stride + d * dense_msk_stride] = xd1;
    }
  }
  stop_bench(my_bs2dense);
}
