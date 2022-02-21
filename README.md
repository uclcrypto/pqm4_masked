# pqm4_masked
Post-quantum masked crypto library for the ARM Cortex-M4

## Introduction
The **pqm4\_masked** library is based on [pqm4](https://github.com/mupq/pqm4).

The main differences are
* the addition of masked implementations of the following schemes:
    - Kyber768 (in `crypto_kem/kyber768/m4fspeed_masked`, tested on the `nucleo-l4r5zi` target);
    **/!\\This implementation is provided for demonstration and benchmarking purposes.
    It has NOT been subjected to any leakage assessment. Use it at your own risk./!\\**
    - SABER (in `crypto_kem/saber/m4fspeed_masked`, tested on the `nucleo-l4r5zi` target);
    **/!\\This implementation is provided for demonstration and benchmarking purposes.
    It has NOT been subjected to any leakage assessment. Use it at your own risk./!\\**

* some changes in the tooling to allow for finer-grained benchmarking and testing
  (for this, we use a slightly modified version of the `mupq` submodule located
  at <https//github.com/uclcrypto/mupq>).

## Setup/Installation

Same as [pqm4](https://github.com/mupq/pqm4).


## API documentation
For the purpose of simple testing and benchmarking, the **pqm4\_masked**
library uses the NIST/SUPERCOP/[PQClean
API](https://github.com/PQClean/PQClean).
(See [pqm4](https://github.com/mupq/pqm4).)

**/!\\As a result, the secret key is stored unmasked, which is insecure./!\\**

**/!\\The epheremal (encapsulated) is also unprotected similarly as [this paper](https://tches.iacr.org/index.php/TCHES/article/view/9064)./!\\**


## Running tests and benchmarks

In addition to pqm4 common testing (`testvectors.py` and benchmarking `benchmarking.py`), `pqm4\_masked` supports:
* testing a single implementation (e.g. `kyber768/m4fspeed_masked` as command-line parameter)
* benchmarking components of the masked implementation `subspeed` in `benchmarks.py`

### Example usages

Before running any code, run `cd libopencm3; make; cd ..`, otherwise the scripts fail.

```shell
# Testvectors
CFLAGS="-DNSHARES=2" python3 testvectors.py -p nucleo-l4r5zi --uart /dev/ttyACM0 kyber768/m4fspeed_masked
# Benchmark whole
CFLAGS="-DNSHARES=2 -DBENCH=0" python3 benchmarks.py -p nucleo-l4r5zi --uart /dev/ttyACM0 kyber768/m4fspeed_masked --speed
# Benchmark components
CFLAGS="-DNSHARES=2 -DBENCH=1" python3 benchmarks.py -p nucleo-l4r5zi --uart /dev/ttyACM0 kyber768/m4fspeed_masked --subspeed
```

### Benchmarks
In order to benchmark a specific function `my_function`, you should first edit `BENCH_CASES` in `common/bench.h` such that:

```c
#define BENCH_CASES X(my_function) X(my_other_function)
```

Then you can increase the performance counter of that specific function within the source files with:
```c 
start_bench(my_function);
my_function(...);
stop_bench(my_function);
```
See also `common/bench.c` and `common/speed_sub.c` for additional details about the benchmarking procedure. 

Note that the benchmark configurations is done through `CFLAGS`. To run them, the flag `BENCH=1` must be defined. In order to benchmark randomness usage instead of cycle count, `BENCH_RND=1` must be defined.
In order to run all benchmarks for `saber` or `kyber768`, run:

```shell
./run_bench.sh kyber768 
python3 parse_bench.py kyber768
```


## Bibliography

When referring to this framework in academic literature, please consider using the following bibTeX excerpt:

```
@misc{cryptoeprint:2022:158,
    author       = {Olivier Bronchain and
		    Gaëtan Cassiers},
    title        = {Bitslicing Arithmetic/Boolean Masking Conversions for Fun and Profit with Application to Lattice-Based KEMs},
    howpublished = {Cryptology ePrint Archive, Report 2022/158},
    year         = {2022},
    note         = {\url{https://ia.cr/2022/158}},
}
```

## License
Different parts of **pqm4\_masked** have different licenses. 
Each subdirectory containing implementations contains a LICENSE or COPYING file stating 
under what license that specific implementation is released. 
The files in common contain licensing information at the top of the file (and 
are currently either public domain or MIT). 
All other code in this repository is released under the conditions of [CC0](https://creativecommons.org/publicdomain/zero/1.0/).

Note: masked implementations are released under the conditions of the
[GPLv3](https://www.gnu.org/licenses/gpl-3.0.txt).
