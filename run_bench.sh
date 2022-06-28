#!/bin/bash
# Copyright 2022 UCLouvain, Belgium
#
# This file is part of pqm4_masked.
#
# pqm4_masked is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation, version 3.
#
# pqm4_masked is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along with
# pqm4_masked. If not, see <https://www.gnu.org/licenses/>.
#

SCHEME=$1
MAXD=8

function run_bench() {
    rm benchmarks/* -rf 
    for D in {2..16}
    do
        rm -rf obj/ bin/
        echo "------------------------------"
        echo "BENCHMARK CYCLES $D SHARES"
        echo "------------------------------"
        CFLAGS="-DNSHARES=$D -DBENCH=1 -DBENCH_RND=0 $CUSE"
        CFLAGS=$CFLAGS python3 benchmarks.py -p nucleo-l4r5zi --uart /dev/ttyACM0 $TARGET --subspeed -o speed 
    done
    echo "case,bench,shares,calls,perf" > $CYCLES_NAME
    cat benchmarks/speed_sub/crypto_kem/$TARGET/* >> $CYCLES_NAME 
}

TARGET=$SCHEME/m4fspeed_masked

# ASM
echo "Benchmarking ASM implementation"
CYCLES_NAME=$SCHEME\_asm_cycles.csv
CUSE=""
run_bench

# C
echo "Benchmarking C implementation"
CYCLES_NAME=$SCHEME\_c_cycles.csv
CUSE="-DUSEC"
run_bench
