#ifndef BENCH_H
#define BENCH_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <hal.h>

#ifndef BENCH
#error "BENCH unset"
#define BENCH 1 // run bench
#endif

#ifndef BENCH_RND
#define BENCH_RND 0 // count randomness instead of time
#endif

#define BENCH_CASES X(keypair) X(encaps) X(decaps) X(keccak) X(my_secadd) X(my_masked_poly_cmp) X(my_cbd) X(my_tomsg) X(my_frommsg) X(my_cmp_finalize) X(my_matacc) X(my_ntt) X(my_seca2b)   X(my_dense2bs) X(my_bs2dense) X(my_seca2b_modp)

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
