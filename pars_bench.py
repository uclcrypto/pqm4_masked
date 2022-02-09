import csv
import pandas
import numpy as np
import matplotlib.pyplot as plt

def get_perf(df,l,k):
    for x in l:
        df = df[df[x[0]].isin([x[1]])]

    print(l)
    print(k)
    print(df)
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

    ds = np.arange(2,9)
    k = "perf"

    df = pandas.read_csv("bench_masked.csv")
    cycles = extract_all(df,ds,["decaps","keccak",'my_cmp_finalize',
        'my_masked_poly_cmp',"my_cbd","my_frommsg","my_tomsg"],"perf")

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

    plt.show()
