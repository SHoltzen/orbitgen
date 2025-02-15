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


def pigeonhole_exact_lifted():
    model = mk_pigeonhole_fg(2,5)
    model.partition()


if __name__ == "__main__":
    # run your test here if you want
    pigeonhole_dtv()
