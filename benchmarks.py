from my_graphs import *
import cProfile, pstats, StringIO
from collections import deque
import my_bliss
import timeit


G = graphs.CycleGraph(500)

def my_fun():
    my_bliss.canonical_form(G, certificate=True, return_graph=False)


if __name__ == "__main__":
    print timeit.timeit(stmt="my_fun()", setup="from __main__ import my_fun", number=100)
