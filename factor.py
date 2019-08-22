from sage.all import *
import numpy.random
from collections import deque
import random
import numpy as np
import my_bliss
from my_graphs import *
import cProfile, pstats, StringIO
import itertools
import time

def findsubsets(S,m):
    return set(itertools.combinations(S, m))

def flip():
    return numpy.random.binomial(1, 0.5)

def div_or_0(a, b):
    if b == 0.0:
        return 0
    else:
        return a / b

### returns a random element from a group g using product replacement
def fast_random_element(g):
    stab = libgap(g)
    return (libgap.PseudoRandom(stab).sage(parent=g)).cycle_tuples(singletons=True)


class FactorGraph:
    ### graph: a sage graph
    ### variables: a list of graph vertices which correspond to variables in the factor graph
    ### factors: a list of list of factors; each list is a single color
    ### potential: state -> real: a function which evaluates the potential on a particular state
    ###    a state is a dictionary assigning variables to Boolean values
    def __init__(self, graph, variables, factors, potential):
        self.graph = graph
        self.variables = variables
        self.factors = factors
        self.potential = potential
        g = graph.automorphism_group(partition=[variables] + factors)
        self.graph_aut = g
        # print self.graph_aut
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

    ### returns: a list of all possible states (as sorted tuples)
    def gen_all_states(self):
        def add_vertex(cur_lst, rst):
            if len(rst) == 0:
                return cur_lst
            cur_add = rst.pop()
            new_lst = []
            for itm in cur_lst:
                itm1 = itm.copy()
                itm2 = itm.copy()
                itm1[cur_add] = False
                itm2[cur_add] = True
                new_lst.append(itm1)
                new_lst.append(itm2)
            return add_vertex(new_lst, rst)
        tbl = add_vertex([dict()], self.variables)
        res = []
        for itm in tbl:
            l = itm.items()
            l.sort()
            res.append(tuple(l))
        return res


    ### compute partition via exhaustive enumeration using orbit generation
    def partition(self):
        # queue of colorings yet to be considered
        queue = deque()
        queue.append([[], self.variables])
        g_order = self.graph_aut_order
        # set of representative graph colorings
        Z = 0.0
        reps = set()
        while len(queue) != 0:
            c = queue.popleft()
            # TODO: turn this into a single GI call by having it return both the
            # automorphism group and the canonical form
            gcanon, cert = my_bliss.canonical_form(self.graph,
                                                   partition=c + self.factors,
                                                   certificate=True,
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

            A, orbits = self.graph.automorphism_group(partition=c+self.factors,
                                                      orbits=True,
                                                      algorithm="bliss")
            orbitsz = g_order / A.order()
            A = None

            states = self.partition_to_state(c)
            for s in states:
                weight = orbitsz * self.potential(s)
                Z += weight
            reps.add((gcanon, (tuple(c_canon[0]), tuple(c_canon[1]))))
            if len(c_canon[0]) + 1 <= len(self.variables) / 2.0:
                # if we need the full automorphism group, we can compute it here.
                # expand this node and add it to the queue
                for o in orbits:
                    e = o.pop() # pick an element in the orbit arbitrarily
                    # check if this orbit is already true, or if it is any of the fixed colors
                    if e in c[1]:
                        # we know this is an uncolored vertex
                        newcolor = [c[0] + [e], [x for x in c[1] if x != e]]
                        queue.append(newcolor)
        return Z


    def brute_force_partition(self):
        states = self.gen_all_states()
        # print "states: " + str(states)
        state_to_idx = dict()
        idx_to_state = dict()
        for (idx, st) in enumerate(states):
            state_to_idx[st] = idx
            idx_to_state[idx] = st
        Z = 0.0
        for st in states:
            Z += self.potential(dict(st))
        return Z

def gen_complete_pairwise_factorgraph(n):
    (g, (v, factors)) = gen_complete_pairwise_factor(n)
    def potential(state):
        p = 0.0
        for v in state.itervalues():
            if v:
                p += 1
        return p
    # make half the factors different colors
    # return FactorGraph(g, v, [factors[:len(factors)/2],factors[len(factors)/2:]], potential)
    return FactorGraph(g, v, [factors], potential)


if __name__ == "__main__":
    fg = gen_complete_pairwise_factorgraph(10)
    print fg.partition()
    print fg.brute_force_partition()
