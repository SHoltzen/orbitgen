from sage.all import *
import numpy.random
import random
from my_graphs import *

def flip():
    return numpy.random.binomial(1, 0.5)

class MarkovModel:
    ### graph: a sage graph
    ### variables: a list of graph vertices which correspond to variables in the factor graph
    ### potential: state -> real: a function which evaluates the potential on a particular state
    ###    a state is a dictionary assigning variables to Boolean values
    def __init__(self, graph, variables, potential):
        self.graph = graph
        self.variables = variables
        self.potential = potential

    ### converts a variable state into a variable partition
    def state_to_partition(self, state):
        var_part_1 = []
        var_part_2 = []
        for v,val in state.iteritems():
            if val == True:
                var_part_1.append(v)
            else:
                var_part_2.append(v)
        return [var_part_1, var_part_2]

    ### computes a burnside process transition beginning from a state particular state
    ### n: number of moves
    def burnside(self, state, n):
        for j in range(0,n):
            var_part = self.state_to_partition(state)
            partition = var_part
            stab = self.graph.automorphism_group(partition=partition)
            p = stab.random_element().cycle_tuples(singletons=True)
            # stab = libgap(stab)
            # p = (libgap.PseudoRandom(stab).sage()).to_cycles(singletons=True)
            # if len(p) == 0:
            #     # this is the identity; for some reason, Sage doesn't do this one correctly
            #     for v in self.variables:
            #         p += [(v,)]

            # now convert p into a state
            # print(p)
            state = dict()
            for cyc in p:
                # sample a value for cyc and fill in the new current value
                v = flip()
                for idx in cyc:
                    if (idx) in self.variables:
                        #-1, zero indexed
                        state[idx] = v
        return state

    ### perform a single step of orbit jumping
    ### returns: a pair, (the ratio of transition probabilities, new state)
    ### n: number of burnside steps to take
    def orbitjump(self, state, n):
        hatx = self.burnside(state, n)
        probx = self.potential(state)
        probhatx = self.potential(hatx)
        orbx = self.graph.automorphism_group(partition=self.state_to_partition(state)).order()
        orbhatx = self.graph.automorphism_group(partition=self.state_to_partition(hatx)).order()
        transitionprob = (probhatx * orbx) / (probx * orbhatx)
        return (transitionprob, hatx)

    ### draw n samples according to orbit jump MCMC with no orbital MCMC
    ### burnsidesize: number of burnside steps to take
    ### returns a list of samples
    def orbitjumpmcmc(self, burnsidesize, n):
        samples = []
        # set up initial random state
        cur_state = dict()
        for v in self.variables:
            cur_state[v] = flip()

        for i in range(0, n):
            (ratio, new_state) = self.orbitjump(cur_state, burnsidesize)
            # accept with transition probability
            acceptprob = numpy.minimum(1, ratio)
            if numpy.random.binomial(1, acceptprob):
                samples.append(new_state)
                cur_state = new_state
            else:
                samples.append(cur_state)
        return samples



def gen_complete_pairwise_factorgraph(n):
    (g, (v, factors)) = gen_complete_pairwise_factor(n)
    def potential(state):
        p = 0.0
        for v in state.itervalues():
            if v:
                p += 1
        return p
    return FactorGraph(g, v, [factors], potential)

def run_burnside():
    cur_state = dict()
    counts = dict()
    g = graphs.CompleteGraph(2)
    f = MarkovModel(g, g.vertices(), None)
    for v in f.variables:
        cur_state[v] = flip()
        counts[v] = 0
    counts[len(f.variables)] = 0

    for i in range(0, 5000):
        cur_state = f.burnside(cur_state, 12)
        # print(cur_state)
        num_true = 0
        for k,v in cur_state.iteritems():
            if v:
                num_true += 1
        counts[num_true] = counts[num_true] + 1

    print(counts)


def main():
    # print some burnside samples
    # print(graph.orbitjumpmcmc(4, 10))
    None

if __name__ == "__main__":
    main()
