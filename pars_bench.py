"""
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
"""
import csv
import pandas
import numpy as np
import sys
import matplotlib.pyplot as plt

ds = np.arange(2,16)
target = sys.argv[1]
rnd_name = target+"_rnd.csv"
cycles_name = target+"_cycles.csv"

if target == "kyber768":
    benchs = ['my_cmp_finalize',"decaps","keccak",
            'my_masked_poly_cmp',"my_cbd","my_frommsg","my_tomsg","my_matacc","my_ntt"]
    
    #benchs +=["my_dense2bs","my_bs2dense"]
    #benchs += ["my_secadd","my_seca2b"]
elif target == "saber":
    benchs = ["decaps","keccak",
            "my_cbd","my_cmp_finalize",
            "my_masked_poly_cmp","my_matacc","my_ntt","my_tomsg"] 

def get_perf(df,l,k):
    for x in l:
        df = df[df[x[0]].isin([x[1]])]
    return int(df[k])

def extract_all(df,ds,names,k):
    dic = {}
    for n in names:
        cycles_k = np.array([get_perf(df,
                        [("case","decaps"),
                        ("shares",d),
                        ("bench",n)],k) for d in ds])

        dic[n] = cycles_k

    return dic

if __name__ == "__main__":

    
    df = pandas.read_csv(cycles_name)
    cycles = extract_all(df,ds,benchs,"perf")

    df = pandas.read_csv(rnd_name)
    rnds = extract_all(df,ds,benchs,"perf")

    # print raw numbers
    plt.figure()
    for k,cycle in cycles.items():
        plt.semilogy(ds,cycle,label=k)
    plt.legend()
    plt.grid(True,which="both",ls="--")

    # print percentage
    plt.figure()
    full = 0
    for k,cycle in cycles.items():
        if k == "decaps":
            continue
        plt.plot(ds,cycle/cycles["decaps"],label=k)
        full += cycle
    plt.plot(ds,(cycles["decaps"]  - full) / cycles["decaps"],label="others")
    plt.legend()
    plt.grid(True,which="both",ls="--")

    # print rnd ratio numbers
    plt.figure()
    plt.title("rnd proportion")
    for (k,cycle),(_,rnd) in zip(cycles.items(),rnds.items()):
        plt.plot(ds,rnd*53/cycle,label=k)
    plt.legend()
    plt.grid(True,which="both",ls="--")

    plt.show()
