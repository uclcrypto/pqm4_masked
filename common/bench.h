#ifndef BENCH_H
#define BENCH_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <hal.h>

#define BENCH 1 // run bench
#define BENCH_RND 0 // count randomness instead of time

#define BENCH_CASES X(keypair) X(encaps) X(decaps) X(keccak)

#define X(x) #x,
static const char* bench_cases_names[] = { BENCH_CASES };
#undef X
#define N_BENCH_CASES (sizeof(bench_cases_names)/sizeof(const char *))

typedef enum {
#define X(x) x,
BENCH_CASES
#undef X
} bench_case;

void reset_bench(bench_case bench);
void reset_benches();
void start_bench(bench_case bench);
void stop_bench(bench_case bench);
uint32_t bench_tot_time(bench_case bench);
uint32_t bench_tot_calls(bench_case bench);
void cyccnt_start();
void print_bench(bench_case bench, const char *s);
void print_all_benches(const char *s);

#endif // BENCH_H
