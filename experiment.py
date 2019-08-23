### experiments from the paper
### not all of these made it into the paper, but they give a good idea of how
### the code works

from sage.all import *
import numpy.random
from collections import deque
import random
import numpy as np
import my_bliss
from my_graphs import *
import cProfile, pstats, StringIO
from test import *
import itertools
import time

# holes pigeons, m holes
def mk_pigeonhole_fg(n, m, order=True):
    w1 = 10000000
    w2 = 100000
    (g, (variables, factors)) = gen_pigeonhole_fg(n, m)
    def potential(state):
        total = 0.0
        # to see every pigeon in exactly one hole
        for p in range(0, n):
            # check the holes for the pigeons
            in_hole = False
            for h in range(0, m):
                if state[(p, h)]:
                    if in_hole:
                        return 0.000000001
                    else:
                        in_hole = True
            if in_hole:
                total += w1


        # check to see no no hole has 2 pigeons
        for h in range(0, m):
            for (p1, p2) in findsubsets(range(0, n), 2):
                if not state[(p1, h)] or not state[(p2, h)]:
                    total += w2

        return total
    return FactorGraph(g, variables, factors, potential)


### computes the total variation distance comparing different sampling methods
def complete_pairwise_dtv():
    model = gen_complete_pairwise_factorgraph(6)
    gibbs = model.gibbs_transition()
    # print(gibbs)
    # print(sum(gibbs))
    within_orbit = model.orbit_transition()
    orbitalmcmc = np.matmul(within_orbit, gibbs)
    # unif = model.uniform_transition()
    burnside = model.burnside_mh_transition(4)
    M = orbitalmcmc

    pv = model.brute_force_prob_vector()
    start = np.zeros([2**len(model.variables)])
    start[10] = 1
    # print(np.linalg.matrix_power(M, 5))
    print("-------------------")
    print("pure gibbs")
    model.total_variation(gibbs, start, 100)
    print("-------------------")
    print("lifted MCMC")
    model.total_variation(orbitalmcmc, start, 100)
    print("------------------")
    print("orbit jump MCMC")
    model.total_variation(burnside, start, 100)

def complete_pairwise_exact():
    model = gen_complete_pairwise_factorgraph(6)
    print("partition: %f" % model.partition())

def pigeonhole_dtv():
    model = mk_pigeonhole_fg(2,5)
    gibbs = model.gibbs_transition()
    # print(gibbs)
    # print(sum(gibbs))
    within_orbit = model.orbit_transition()
    orbitalmcmc = np.matmul(within_orbit, gibbs)
    # unif = model.uniform_transition()
    burnside = model.burnside_mh_transition(4)
    M = orbitalmcmc

    pv = model.brute_force_prob_vector()
    start = np.zeros([2**len(model.variables)])
    start[10] = 1
    # print(np.linalg.matrix_power(M, 5))
    print("-------------------")
    print("pure gibbs")
    model.total_variation(gibbs, start, 100)
    print("-------------------")
    print("lifted MCMC")
    model.total_variation(orbitalmcmc, start, 100)
    print("------------------")
    print("orbit jump MCMC")
    model.total_variation(burnside, start, 100)


if __name__ == "__main__":
    # run your test here if you want
    pigeonhole_dtv()
