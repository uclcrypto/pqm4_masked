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

#include "macros.i"
#include "CT_butterflies.i"

#ifndef LOOP
#define LOOP
#endif

.macro sbfx_dual a, tmp
    sbfx.w \tmp, \a, #16, #13
    sbfx.w \a, \a, #0, #13
    pkhbt.w \a, \a, \tmp, lsl #16
.endm

.syntax unified
.cpu cortex-m4

.align 2
.global __asm_negacyclic_ntt_16_light_m
.type __asm_negacyclic_ntt_16_light_m, %function
__asm_negacyclic_ntt_16_light_m:
    push {r4-r12, lr}

    .equ width, 2
    .equ logq, 13

    mov.w r14, r3

    vldm.w r1!, {s4-s7}
    vmov.w s1, r1

#ifdef LOOP
    add.w r12, r0, #32*width
    vmov.w s2, r12
    _0_1_2_16_light:
#else
.rept 16
#endif

    ldrstrvecjump ldr.w, r14, r4, r5, r6, r7, r8, r9, r10, r11, #32*width, #64*width, #96*width, #128*width, #160*width, #192*width, #224*width, #2*width

    sbfx_dual r4, r12
    sbfx_dual r5, r12
    sbfx_dual r6, r12
    sbfx_dual r7, r12
    sbfx_dual r8, r12
    sbfx_dual r9, r12
    sbfx_dual r10, r12
    sbfx_dual r11, r12

    _3_layer_double_CT_16 r4, r5, r6, r7, r8, r9, r10, r11, s4, s5, s6, s7, r1, r2, r3, r12
    ldrstrvecjump str.w, r0, r4, r5, r6, r7, r8, r9, r10, r11, #32*width, #64*width, #96*width, #128*width, #160*width, #192*width, #224*width, #2*width

#ifdef LOOP
    vmov.w r12, s2
    cmp.w r0, r12
    bne.w _0_1_2_16_light
#else
.endr
#endif

    sub.w r0, r0, #32*width

#ifdef LOOP
    add.w r12, r0, #256*width
    vmov.w s2, r12
    _3_4_5_16_light:
#else
.rept 8
#endif

    vmov.w r1, s1
    vldm.w r1!, {s4-s7}
    vmov.w s1, r1

#ifdef LOOP
    add.w r14, r0, #4*width
    vmov.w s3, r14
    _3_4_5_16_inner_light:
#else
.rept 2
#endif

    ldrstrvec ldr.w, r0, r4, r5, r6, r7, r8, r9, r10, r11, #0*width, #4*width, #8*width, #12*width, #16*width, #20*width, #24*width, #28*width
    _3_layer_double_CT_16 r4, r5, r6, r7, r8, r9, r10, r11, s4, s5, s6, s7, r1, r2, r3, r12
    ldrstrvecjump str.w, r0, r4, r5, r6, r7, r8, r9, r10, r11, #4*width, #8*width, #12*width, #16*width, #20*width, #24*width, #28*width, #2*width

#ifdef LOOP
    vmov.w r14, s3
    cmp.w r0, r14
    bne.w _3_4_5_16_inner_light
#else
.endr
#endif

    add.w r0, r0, #28*width

#ifdef LOOP
    vmov.w r12, s2
    cmp.w r0, r12
    bne.w _3_4_5_16_light
#else
.endr
#endif

    pop {r4-r12, pc}

.align 2
.global __asm_negacyclic_ntt_32_m
.type __asm_negacyclic_ntt_32_m, %function
__asm_negacyclic_ntt_32_m:
    push {r4-r12, lr}

    .equ ldrwidth, 2
    .equ strwidth, 4

    vmov.w s0, r0
    ldr.w r0, [sp, #40]

    vldm.w r1!, {s4-s10}
    vmov.w s1, r1

#ifdef LOOP
    add.w r12, r0, #32*ldrwidth
    vmov.w s2, r12
    _0_1_2_32:
#else
.rept 16
#endif

.rept 2

    ldrstrvecjump ldrsh.w, r0, r4, r5, r6, r7, r8, r9, r10, r11, #32*ldrwidth, #64*ldrwidth, #96*ldrwidth, #128*ldrwidth, #160*ldrwidth, #192*ldrwidth, #224*ldrwidth, #ldrwidth

    sbfx.w r4, r4, #0, #13
    sbfx.w r5, r5, #0, #13
    sbfx.w r6, r6, #0, #13
    sbfx.w r7, r7, #0, #13
    sbfx.w r8, r8, #0, #13
    sbfx.w r9, r9, #0, #13
    sbfx.w r10, r10, #0, #13
    sbfx.w r11, r11, #0, #13

    _3_layer_CT_32 r4, r5, r6, r7, r8, r9, r10, r11, s4, s5, s6, s7, s8, s9, s10, r1, r2, r3, r12, r14

    vmov.w r14, s0
    ldrstrvecjump str.w, r14, r4, r5, r6, r7, r8, r9, r10, r11, #32*strwidth, #64*strwidth, #96*strwidth, #128*strwidth, #160*strwidth, #192*strwidth, #224*strwidth, #strwidth
    vmov.w s0, r14

.endr

#ifdef LOOP
    vmov.w r12, s2
    cmp.w r0, r12
    bne.w _0_1_2_32
#else
.endr
#endif

    sub.w r0, r14, #32*strwidth

#ifdef LOOP
    add.w r12, r0, #256*strwidth
    vmov.w s2, r12
    _3_4_5_32:
#else
.rept 8
#endif

    vmov.w r1, s1
    vldm.w r1!, {s4-s10}
    vmov.w s1, r1

#ifdef LOOP
    add.w r14, r0, #4*strwidth
    vmov.w s3, r14
    _3_4_5_inner_32:
#else
.rept 2
#endif

.rept 2

    ldrstrvec ldr.w, r0, r4, r5, r6, r7, r8, r9, r10, r11, #0*strwidth, #4*strwidth, #8*strwidth, #12*strwidth, #16*strwidth, #20*strwidth, #24*strwidth, #28*strwidth
    _3_layer_CT_32 r4, r5, r6, r7, r8, r9, r10, r11, s4, s5, s6, s7, s8, s9, s10, r1, r2, r3, r12, r14
    ldrstrvecjump str.w, r0, r4, r5, r6, r7, r8, r9, r10, r11, #4*strwidth, #8*strwidth, #12*strwidth, #16*strwidth, #20*strwidth, #24*strwidth, #28*strwidth, #strwidth

.endr

#ifdef LOOP
    vmov.w r14, s3
    cmp.w r0, r14
    bne.w _3_4_5_inner_32
#else
.endr
#endif

    add.w r0, r0, #28*strwidth

#ifdef LOOP
    vmov.w r12, s2
    cmp.w r0, r12
    bne.w _3_4_5_32
#else
.endr
#endif

    pop {r4-r12, pc}


















