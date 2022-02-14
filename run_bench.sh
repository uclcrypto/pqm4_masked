#!/bin/bash

rm benchmarks/* -rf 
for D in {2..16}
do
    rm -rf obj/ bin/
    echo "------------------------------"
    echo "BENCHMARK CYCLES $D SHARES"
    echo "------------------------------"
    CFLAGS="-DNSHARES=$D -DBENCH=1 -DBENCH_RND=0" python3 benchmarks.py -p nucleo-l4r5zi --uart /dev/ttyACM0 kyber768/m4fspeed_masked --subspeed -o speed 
done

echo "case,bench,shares,calls,perf" > bench_masked_cycles.csv
cat benchmarks/speed_sub/crypto_kem/kyber768/m4fspeed_masked/* >> bench_masked_cycles.csv

rm benchmarks/* -rf 
for D in {2..16}
do
    rm -rf obj/ bin/
    echo "------------------------------"
    echo "BENCHMARK RANDOMNES $D SHARES"
    echo "------------------------------"
    CFLAGS="-DNSHARES=$D -DBENCH=1 -DBENCH_RND=1" python3 benchmarks.py -p nucleo-l4r5zi --uart /dev/ttyACM0 kyber768/m4fspeed_masked --subspeed -o speed 
done

echo "case,bench,shares,calls,perf" > bench_masked_rnd.csv
cat benchmarks/speed_sub/crypto_kem/kyber768/m4fspeed_masked/* >> bench_masked_rnd.csv
