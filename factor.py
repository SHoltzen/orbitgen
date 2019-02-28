from sage.all import *
import numpy.random
from collections import deque
import random
import my_bliss
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
        self.graph_aut = graph.automorphism_group()
        self.graph_aut_order = self.graph_aut.order()

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

    def partition_to_state(self, part):
        state1 = dict()
        state2 = dict()
        for v in part[0]:
            state1[v] = True
            state2[v] = False
        for v in part[1]:
            state1[v] = False
            state2[v] = True
        if len(part[0]) < len(self.variables) / 2.0:
            return [state1, state2]
        else:
            return [state1]


    def query_enumerate(self, query):
        # queue of colorings yet to be considered
        queue = deque()
        queue.append([[], self.variables])
        g_order = self.graph_aut_order
        # set of representative graph colorings
        prob = 0.0
        Z = 0.0
        reps = set()
        while len(queue) != 0:
            # print("queue: %s" % queue)
            c = queue.popleft()
            # TODO: turn this into a single GI call by having it return both the
            # automorphism group and the canonical fora
            gcanon, cert = my_bliss.canonical_form(self.graph, partition=c, certificate=True,
                                                   return_graph=False)
            gcanon = tuple(gcanon)
            # convert the colors to their coloring in the canonical graph
            c_canon = [[], []]
            for c1 in c[0]:
                c_canon[0].append(cert[c1])
            for c1 in c[1]:
                c_canon[1].append(cert[c1])

            c_canon[0].sort()
            c_canon[1].sort()

            if (gcanon, (tuple(c_canon[0]), tuple(c_canon[1]))) in reps:
                continue

            A, orbits = self.graph.automorphism_group(partition=c, orbits=True,
                                                 algorithm="bliss")
            orbitsz = g_order / A.order()
            states = self.partition_to_state(c)
            for s in states:
                weight = orbitsz * self.potential(s)
                Z += weight
                if query(s):
                    prob += weight
            reps.add((gcanon, (tuple(c_canon[0]), tuple(c_canon[1]))))
            if len(c_canon[0]) + 1 <= len(self.variables) / 2.0:
                # if we need the full automorphism group, we can compute it here.
                # expand this node and add it to the queue
                for o in orbits:
                    e = o.pop() # pick an element in the orbit arbitrarily
                    # check if this orbit is already true, or if it is any of the fixed colors
                    if e in c[0]:
                        continue
                    # we know this is an uncolored vertex
                    newcolor = [c[0] + [e], [x for x in c[1] if x != e]]
                    queue.append(newcolor)
        return prob / Z


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
                    state[idx] = v
        return state

    ### perform a single step of orbit jumping
    ### returns: a pair, (the ratio of transition probabilities, new state)
    ### n: number of burnside steps to take
    def orbitjump(self, state, n):
        hatx = self.burnside(state, n)
        probx = self.potential(state)
        probhatx = self.potential(hatx)
        orbx = self.graph_aut_order / self.graph.automorphism_group(partition=self.state_to_partition(state)).order()
        orbhatx = self.graph_aut_order / self.graph.automorphism_group(partition=self.state_to_partition(hatx)).order()
        # print("x: %s, hatx: %s, probx: %s, probhatx: %s, orbx: %s, orbhatx: %s" %
        #       (state, hatx, probx, probhatx, orbx, orbhatx))
        try:
            transitionprob = (probhatx * orbhatx) / (probx * orbx)
            return (transitionprob, hatx)
        except:
            # divided by 0
            return (1.0, hatx)

    ### draw n samples according to orbit jump MCMC with no orbital MCMC
    ### burnsidesize: number of burnside steps to take
    ### query: state -> bool, evaluates a query on a state
    ### gamma: how often to orbit-jump
    ### burn: the burn in of the chain
    ### returns the probability of the query
    def orbitjumpmcmc(self, burnsidesize, n, query, burn):
        samples = []
        # set up initial random state
        cur_state = dict()
        for v in self.variables:
            cur_state[v] = flip()

        query_count = 0.0
        cur_step = 0
        counts = dict()
        for v in self.variables:
            counts[v] = 0
        counts[len(self.variables)] = 0

        for i in range(0, n * burn):
            (ratio, new_state) = self.orbitjump(cur_state, burnsidesize)
            num_true = 0
            for k,v in new_state.iteritems():
                if v:
                    num_true += 1
            counts[num_true] = counts[num_true] + 1

            # accept with transition probability
            acceptprob = numpy.minimum(1, ratio)
            if numpy.random.binomial(1, acceptprob):
                cur_state = new_state
            if cur_step % burn == 0:
                if query(cur_state):
                    query_count += 1
            cur_step += 1
        print("counts: %s " % counts)
        return query_count / n



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
    n = 100
    g = graphs.CompleteGraph(2)
    def potential(state):
        p = 0.0
        for v in state.itervalues():
            if v:
                p += 1
        return p
    def query(state):
        num_true = 0
        # print(r)
        for k,v in state.iteritems():
            if v:
                num_true += 1
        return num_true == 1


    graph = MarkovModel(g, g.vertices(), potential)
    print("exact: %s" % graph.query_enumerate(query))
    print("query: %s" % graph.orbitjumpmcmc(4, n, query, 10))


if __name__ == "__main__":
    main()
