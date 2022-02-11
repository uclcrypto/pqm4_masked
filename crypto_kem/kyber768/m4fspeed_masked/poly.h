/* Copyright 2022 UCLouvain, Belgium and PQM4 contributors
 *
 * This file is part of pqm4_masked.
 *
 * pqm4_masked is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, version 3.
 *
 * pqm4_masked is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * pqm4_masked. If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef POLY_H
#define POLY_H

#include "params.h"

#include <stdint.h>

#define poly_getnoise(p, seed, nonce) poly_noise(p, seed, nonce, 0)
#define poly_addnoise(p, seed, nonce) poly_noise(p, seed, nonce, 1)

/*
 * Elements of R_q = Z_q[X]/(X^n + 1). Represents polynomial
 * coeffs[0] + X*coeffs[1] + X^2*xoeffs[2] + ... + X^{n-1}*coeffs[n-1]
 */
typedef struct {
  int16_t coeffs[KYBER_N];
} poly;

void poly_compress(unsigned char *r, poly *a);
void poly_decompress(poly *r, const unsigned char *a);

void poly_packcompress(unsigned char *r, poly *a, int i);
void poly_unpackdecompress(poly *r, const unsigned char *a, int i);

int cmp_poly_compress(const unsigned char *r, poly *a);
int cmp_poly_packcompress(const unsigned char *r, poly *a, int i);

void poly_tobytes(unsigned char *r, poly *a);
void poly_frombytes(poly *r, const unsigned char *a);
void poly_frombytes_mul(poly *r, const unsigned char *a);

void poly_frommsg(poly *r, const unsigned char msg[KYBER_SYMBYTES]);
void poly_tomsg(unsigned char msg[KYBER_SYMBYTES], poly *a);

void poly_noise(poly *r, const unsigned char *seed, unsigned char nonce,
                int add);

void poly_ntt(poly *r);
void poly_invntt(poly *r);
void poly_basemul(poly *r, const poly *a, const poly *b);
void poly_basemul_i16(int16_t *r, const int16_t *a, const int16_t *b);
void poly_basemul_acc(poly *r, const poly *a, const poly *b);
void poly_basemul_acc_i16(int16_t *r, const int16_t *a, const int16_t *b);
void poly_frommont(poly *r);

void poly_reduce(poly *r);
void poly_reduce_i16(int16_t *r);

void poly_add(poly *r, const poly *a, const poly *b);
void poly_add_i16(int16_t *r, const int16_t *a, const int16_t *b);
void poly_sub(poly *r, const poly *a, const poly *b);
void poly_sub_i16(int16_t *r, const int16_t *a, const int16_t *b);

void poly_zeroize(poly *p);

#endif
