
#include "macros.i"

#ifndef LOOP
#define LOOP
#endif

.syntax unified
.cpu cortex-m4

.align 2
.global __asm_CRT
.type __asm_CRT, %function
__asm_CRT:
    push {r4-r12, lr}

    .equ width, 2

    vmov.w s0, r0
    mov.w r0, r1
    mov.w r1, r2

    // r2 = Q1half; r3 = Q1;
    // r8 = Q2invRmod; r9 = Q1prime;
    // r10 = Q2; r11 = Q2bar;
    ldm.w r3, {r2-r3, r8-r11}

#ifdef LOOP
    add.w r12, r1, #256*width
    vmov.w s2, r12
    _CRT:
#else
.rept 16
#endif

.rept 8

    ldr.w r5, [r0, #2*width]
    ldr.w r4, [r0], #4*width
    ldrsh.w r7, [r1, #1*width]
    ldrsh.w r6, [r1], #2*width

    barrett_32 r6, r11, r10, r12
    barrett_32 r7, r11, r10, r12

    sub.w r4, r4, r6
    montgomery_mul_32 r4, r8, r9, r3, r12, r14
    central_reduce r4, r2, r3
    mla.w r4, r4, r10, r6
    sbfx.w r4, r4, #0, #13

    sub.w r5, r5, r7
    montgomery_mul_32 r5, r8, r9, r3, r12, r14
    central_reduce r5, r2, r3
    mla.w r5, r5, r10, r7
    sbfx.w r5, r5, #0, #13

    vmov.w r14, s0
    strh.w r5, [r14, #1*width]
    strh.w r4, [r14], #2*width
    vmov.w s0, r14

.endr

#ifdef LOOP
    vmov.w r12, s2
    cmp.w r1, r12
    bne.w _CRT
#else
.endr
#endif

    pop {r4-r12, pc}

