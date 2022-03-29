/* Written in 2021-2022 by Amin Abdulrahman and Jiun-Peng Chen and Yu-Jia Chen
 * Vincent Hwang and Matthias J. Kannwischer and Bo-Yin Yang and  UCLouvain 
 * 
 *  To the extent possible under law, the author(s) have dedicated all
 *  copyright and related and neighboring rights to this software to the
 *  public domain worldwide. This software is distributed without any
 *  warranty.
 *  You should have received a copy of the CC0 Public Domain Dedication along
 *  with this software. If not, see
 *  <http://creativecommons.org/publicdomain/zero/1.0/>. 
 *  see https://github.com/multi-moduli-ntt-saber/multi-moduli-ntt-saber
 */
#ifndef NTT_H
#define NTT_H

#include "mNTT_params.h"

static const int32_t CRT_constants[8] = {Q1half,  Q1, Q2invRmod,
                                         Q1prime, Q2, Q2bar};

static const int32_t streamlined_CT_negacyclic_table_Q1[NTT_N] = {
    -17815616, 20179042, 5173450,   -10771126, 16264029,  7771222,   22209164,
    20168288,  12843351, -21065524, 9549694,   -18485124, 13730545,  -3408104,
    -12756821, -2111760, 9083870,   12322149,  17357617,  15860553,  22071560,
    10680947,  528447,   -18818674, -7231958,  14373826,  -3110586,  -13073414,
    12881717,  17638719, -15509723, -10953473, -16918429, -19802953, 16381402,
    3052338,   -9581747, 8857767,   17005820,  -16701106, -20029869, 21170404,
    16873658,  18074717, -11294276, 18712276,  -5233444,  -14669727, -14263517,
    18655939,  19085761, 1633657,   13490703,  -4090736,  7989938,   -16270819,
    3219765,   742202,   2078156,   12706786,  19999052,  21125353,  406480,
    0};

static const int16_t
    streamlined_CT_negacyclic_table_Q2[(NTT_N - 1) + (1 << 0) + (1 << 3)] = {
        0,    -164, -81,  361,  186,  -3,   -250, -120, 0,    -308, -76,  -98,
        147,  -114, -272, 54,   0,    129,  36,   -75,  -2,   -124, -80,  -346,
        0,    -16,  -339, -255, 86,   -51,  364,  267,  0,    -223, 282,  -203,
        161,  -15,  288,  169,  0,    -362, -34,  199,  191,  307,  -50,  -24,
        0,    -143, 178,  270,  -170, 226,  121,  -188, 0,    131,  -10,  149,
        -380, 279,  180,  -375, 0,    -337, 369,  -192, -157, 263,  -128, -246};

static const int32_t mul_Rmod_table_Q1[NTT_N >> 1] = {
    9549694,   -18485124, 13730545,  -3408104,  12322149,  17357617,  15860553,
    22071560,  -7231958,  14373826,  -3110586,  -13073414, -10953473, -16918429,
    -19802953, 16381402,  17005820,  -16701106, -20029869, 21170404,  18712276,
    -5233444,  -14669727, -14263517, 13490703,  -4090736,  7989938,   -16270819,
    12706786,  19999052,  21125353,  406480};

static const int16_t mul_Rmod_table_Q2[NTT_N >> 1] = {
    147, -114, -272, 54,  -2,  -124, -80,  -346, 86,   -51,  364,
    267, 161,  -15,  288, 169, 191,  307,  -50,  -24,  -170, 226,
    121, -188, -380, 279, 180, -375, -157, 263,  -128, -246};

static const int32_t streamlined_inv_CT_negacyclic_table_Q1[NTT_N << 1] = {
    5361568,   5361568,   17815616,  5361568,   -5173450,  17815616,  -20179042,
    5361568,   5361568,   17815616,  5361568,   -5173450,  17815616,  -20179042,
    5211980,   10191609,  -18958318, -2536408,  -20370394, 6798571,   17653340,
    9256845,   -22209164, -3219765,  -18655939, -2078156,  -1633657,  -742202,
    -19085761, -3913294,  -9518846,  -6761239,  -22261278, -13080044, 4166159,
    -12172318, -22181021, -5173450,  -22209164, -7771222,  -3219765,  -16873658,
    -18655939, -3052338,  -9980694,  1477987,   -15907284, -20235687, 21681607,
    8100971,   -16574872, -8749497,  -16264029, -12881717, -10680947, 15509723,
    18818674,  -17638719, -528447,   17287360,  -15924680, -1954003,  12671720,
    -187509,   -5761119,  -18776506, -8460216,  17815616,  -5173450,  -20179042,
    -22209164, -16264029, -7771222,  10771126,  -2901148,  11565223,  21558827,
    -18751213, 21416902,  16292070,  14968358,  -19786940, -7771222,  -16873658,
    -3052338,  11294276,  -8857767,  -18074717, 9581747,   18051007,  -11458335,
    -4297738,  402981,    -654814,   19144908,  15910857,  -13448571, -20179042,
    -16264029, 10771126,  -12881717, 12756821,  -10680947, -20168288, -6753887,
    6198868,   5533187,   -16866197, -17400052, 10085233,  -11482435, 8453394,
    10771126,  12756821,  -20168288, -9083870,  21065524,  2111760,   -12843351,
    9582513,   -14672722, -16368664, -10330663, 17818646,  -7297800,  16215382,
    -9371731,  0};

static const int16_t
    streamlined_inv_CT_negacyclic_table_Q2[(NTT_N - 1) + (1 << 0) + (1 << 3) +
                                           NTT_N] = {
        0,    171,  171,  164,  171,  -361, 164,  81,   0,    171,  171,  164,
        171,  -361, 164,  81,   -228, -160, 225,  -4,   294,  -77,  -108, 248,
        0,    120,  337,  -131, 192,  -149, -369, 10,   -328, 269,  -162, 372,
        342,  240,  47,   6,    0,    -361, 120,  250,  337,  143,  -131, 362,
        -256, 360,  -314, 9,    -277, -19,  243,  211,  0,    3,    223,  16,
        203,  255,  -282, 339,  -31,  356,  -20,  -68,  384,  229,  -298, 371,
        0,    164,  -361, 81,   120,  3,    250,  -186, -193, -41,  322,  172,
        -338, 235,  30,   102,  0,    250,  143,  362,  -270, -199, -178, 34,
        262,  -32,  45,   153,  -95,  -323, 286,  -258, 0,    81,   3,    -186,
        223,  -129, 16,   308,  242,  -100, -340, 382,  376,  48,   317,  155,
        0,    -186, -129, 308,  75,   98,   -36,  76,   -205, 72,   91,   -152,
        -363, 150,  -259, 196};

extern void __asm_negacyclic_ntt_32_m(uint32_t *des, const int32_t *table,
                                      int32_t Qprime, int32_t Q, uint16_t *src,
                                      int32_t RmodQ);
extern void __asm_base_mul_32(uint32_t *des, const int32_t *table,
                              int32_t Qprime, int32_t Q, uint32_t *src1,
                              uint32_t *src2);
extern void __asm_base_mul_acc_32(uint32_t *des, const int32_t *table,
                                  int32_t Qprime, int32_t Q, uint32_t *src1,
                                  uint32_t *src2);
extern void __asm_negacyclic_intt_32_m(uint32_t *src, const int32_t *table,
                                       int32_t Qprime, int32_t Q);
extern void __asm_negacyclic_intt_central_32(uint32_t *src,
                                             const int32_t *table,
                                             int32_t Qprime, int32_t Q);

extern void __asm_negacyclic_ntt_16_light_m(uint16_t *des, const int16_t *table,
                                            int32_t QQprime, uint16_t *src);
extern void __asm_base_mul_16(uint16_t *des, const int16_t *table,
                              int32_t QQprime, uint16_t *src1, uint16_t *src2);
extern void __asm_base_mul_acc_16(uint16_t *des, const int16_t *table,
                                  int32_t QQprime, uint16_t *src1,
                                  uint16_t *src2);
extern void __asm_negacyclic_intt_16_light(uint16_t *des, const int16_t *table,
                                           int32_t QQprime, int32_t RmodQ);

extern void __asm_CRT(uint16_t *des, uint32_t *src_32, uint16_t *src_16,
                      const int32_t *CRT_const);

#define NTT_forward_32(out, in)                                                \
  __asm_negacyclic_ntt_32_m(out, streamlined_CT_negacyclic_table_Q1, Q1prime,  \
                            Q1, in, RmodQ1)
#define NTT_forward_16(out, in)                                                \
  __asm_negacyclic_ntt_16_light_m(out, streamlined_CT_negacyclic_table_Q2,     \
                                  Q2Q2prime, in)

#define NTT_mul_32(des, src1, src2)                                            \
  __asm_base_mul_32(des, mul_Rmod_table_Q1, Q1prime, Q1, src1, src2)
#define NTT_mul_16(des, src1, src2)                                            \
  __asm_base_mul_16(des, mul_Rmod_table_Q2, Q2Q2prime, src1, src2)

#define NTT_mul_acc_32(des, src1, src2)                                        \
  __asm_base_mul_acc_32(des, mul_Rmod_table_Q1, Q1prime, Q1, src1, src2)
#define NTT_mul_acc_16(des, src1, src2)                                        \
  __asm_base_mul_acc_16(des, mul_Rmod_table_Q2, Q2Q2prime, src1, src2)

#define NTT_inv_32(des)                                                        \
  __asm_negacyclic_intt_32_m(des, streamlined_inv_CT_negacyclic_table_Q1,      \
                             Q1prime, Q1)
#define NTT_inv_central_32(des)                                                \
  __asm_negacyclic_intt_central_32(                                            \
      des, streamlined_inv_CT_negacyclic_table_Q1, Q1prime, Q1)
#define NTT_inv_16(des)                                                        \
  __asm_negacyclic_intt_16_light(des, streamlined_inv_CT_negacyclic_table_Q2,  \
                                 Q2Q2prime, RmodQ2)
#define solv_CRT(des, src_32, src_16)                                          \
  __asm_CRT(des, src_32, src_16, CRT_constants)

#endif
