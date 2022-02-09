#!/bin/bash

rm benchmarks/* -rf 
for D in {2..8}
do
    rm -rf obj/ bin/
    CFLAGS="-DNSHARES=$D" python3 benchmarks.py -p nucleo-l4r5zi --uart /dev/ttyACM0 kyber768/m4fspeed_masked --nostack --nospeed --nohashing --nosize -o speed 
done

echo "case,bench,shares,calls,perf" > bench_masked.csv
cat benchmarks/speed_sub/crypto_kem/kyber768/m4fspeed_masked/* >> bench_masked.csv
