#!/bin/bash

SCHEME=kyber768
TARGET=$SCHEME/m4fspeed_masked
CYCLES_NAME=$SCHEME\_cycles.csv
RND_NAME=$SCHEME\_rnd.csv
echo $RND_NAME
rm benchmarks/* -rf 
for D in {2..16}
do
    rm -rf obj/ bin/
    echo "------------------------------"
    echo "BENCHMARK CYCLES $D SHARES"
    echo "------------------------------"
    CFLAGS="-DNSHARES=$D -DBENCH=1 -DBENCH_RND=0" python3 benchmarks.py -p nucleo-l4r5zi --uart /dev/ttyACM0 $TARGET --subspeed -o speed 
done

echo "case,bench,shares,calls,perf" > $CYCLES_NAME
cat benchmarks/speed_sub/crypto_kem/$TARGET/* >> $CYCLES_NAME 

rm benchmarks/* -rf 
for D in {2..16}
do
    rm -rf obj/ bin/
    echo "------------------------------"
    echo "BENCHMARK RANDOMNES $D SHARES"
    echo "------------------------------"
    CFLAGS="-DNSHARES=$D -DBENCH=1 -DBENCH_RND=1" python3 benchmarks.py -p nucleo-l4r5zi --uart /dev/ttyACM0 $TARGET --subspeed -o speed 
done

echo "case,bench,shares,calls,perf" > $RND_NAME
cat benchmarks/speed_sub/crypto_kem/$TARGET/* >> $RND_NAME
