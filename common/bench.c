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

#include "bench.h"
#include "sendfn.h"
#include <stdio.h>
#include <string.h>

#define printcycles(S, U) send_unsignedll((S), (U))

#if (BENCH_RND)
#define GET_TIME rngcnt_get // TODO
#else
#define GET_TIME cyccnt_get
#endif

#define ERROR(...) do { \
    char bail_buf[100]; sprintf(bail_buf, __VA_ARGS__); hal_send_str(bail_buf); hal_send_str("#"); \
} while (0)


/* Don't use opencm3 here, since not all platforms might use opencm3, but all
   have these DWT registers */
#define DWT_CTRL           (*(volatile uint32_t*)(0xE0001000u + 0x00))
#define DWT_CYCCNT         (*(volatile uint32_t*)(0xE0001000u + 0x04))
#define DWT_CTRL_CYCCNTENA (1 << 0)
#define SCS_DEMCR          (*(volatile uint32_t*)(0xE000E000u + 0xDFC))
#define SCS_DEMCR_TRCENA   (1 << 24)
void cyccnt_start() {
    SCS_DEMCR |= SCS_DEMCR_TRCENA;
    DWT_CYCCNT = 0;
    DWT_CTRL |= DWT_CTRL_CYCCNTENA;
    DWT_CYCCNT = 0;
}
static inline uint32_t cyccnt_get() {
  return DWT_CYCCNT;
}
uint64_t rng_cnt = 0;
static inline uint32_t rngcnt_get() {
  return rng_cnt;
}
typedef struct {
    uint32_t start_time;
    uint32_t tot_time;
    uint32_t n_calls;
    bool running;
} bench_struct;

#define X(x) #x,
static const char* bench_cases_names[] = { BENCH_CASES };
#undef X
#define N_BENCH_CASES (sizeof(bench_cases_names)/sizeof(const char *))

bench_struct benches[N_BENCH_CASES] = { 0 };


void reset_bench(bench_case bench) {
    memset(&benches[bench], 0, sizeof(bench_struct));
}

void reset_benches() {
    memset(benches, 0, sizeof(benches));
}

void start_bench(bench_case bench) {
#if (BENCH)
    bench_struct *b = &benches[bench];
    if (b->running) {
        ERROR("bench running: %s", bench_cases_names[bench]);
    }
    b->running = true;
    b->start_time = GET_TIME();
#else
    (void) bench; // avoid unused variable warning
#endif
}

void stop_bench(bench_case bench) {
#if (BENCH)
    uint32_t t = GET_TIME();
    bench_struct *b = &benches[bench];
    if (!b->running) {
        ERROR("bench not running: %s", bench_cases_names[bench]);
    }
    b->running = false;
    b->tot_time += t- b->start_time;
    b->n_calls += 1;
#else
    (void) bench; // avoid unused variable warning
#endif
}

uint32_t bench_tot_time(bench_case bench) {
    return benches[bench].tot_time;
}

uint32_t bench_tot_calls(bench_case bench) {
    return benches[bench].n_calls;
}

void print_bench(bench_case bench, const char *s) {
    char buf[200];
    sprintf(buf, "%s,%s,%d,%ld,%ld",s,bench_cases_names[bench],NSHARES,bench_tot_calls(bench),bench_tot_time(bench));
    hal_send_str(buf);
}

void print_all_benches(const char *s) {
    for (size_t i=0; i<N_BENCH_CASES; i++) {
        print_bench(i, s);
    }
}
