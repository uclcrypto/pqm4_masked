#include "bench.h"
#include "sendfn.h"
#include <stdio.h>
#include <string.h>

#define printcycles(S, U) send_unsignedll((S), (U))

#if (BENCH_RND)
#define GET_TIME 0 // TODO
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

typedef struct {
    uint32_t start_time;
    uint32_t tot_time;
    uint32_t n_calls;
    bool running;
} bench_struct;

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
//    sprintf(buf, "%s %s calls:", s, bench_cases_names[bench]);
    sprintf(buf, "%s,%s,%d,%ld,%ld",s,bench_cases_names[bench],NSHARES,bench_tot_calls(bench),bench_tot_time(bench));
    hal_send_str(buf);
/*    printcycles(buf, bench_tot_calls(bench));
#if BENCH_RND
    sprintf(buf, "%s %s rnd:", s, bench_cases_names[bench]);
#else
    sprintf(buf, "%s %s cycles:", s, bench_cases_names[bench]);
#endif
    printcycles(buf, bench_tot_time(bench));
*/
}

void print_all_benches(const char *s) {
    for (size_t i=0; i<N_BENCH_CASES; i++) {
        print_bench(i, s);
    }
}
