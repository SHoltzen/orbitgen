### The main file that was used to run experiments in the paper.
### Defines a simple factor class, a brute force exact inference
### algorithm, lifted orbit generation and orbit-jump MCMC.

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
        g = graph.automorphism_group(partition=[variables] + self.factors)
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
        tbl = add_vertex([dict()], self.variables[:])
        res = []
        for itm in tbl:
            l = itm.items()
            l.sort()
            res.append(tuple(l))
        return res

    ### Given a group g, strips the factor generators from that group
    def strip_factors(self, g):
        stripped_gens = []
        for gen in g.gens():
            sgen = []
            for cycle in gen.cycle_tuples(singletons=True):
                # print "checking cycle " + str(cycle)
                # check this cycle contains a factor
                contains_factor = False
                for factor in self.factors:
                    for e in factor:
                        # print "checking if " + str(e) + " in " + str(cycle)
                        if e in cycle:
                            contains_factor = True
                if not contains_factor:
                    # print "no factor in " + str(cycle)
                    sgen.append(cycle)
            stripped_gens.append(tuple(sgen))
        return PermutationGroup(stripped_gens)



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


    ### brute forces the partition function by enumerating all possible states
    ### for testing and benchmarking purposes
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

    ### computes the transition matrix of a single step of the burnside process
    def burnside_transition(self):
        states = self.gen_all_states()
        state_to_idx = dict()
        idx_to_state = dict()
        for (idx, st) in enumerate(states):
            state_to_idx[st] = idx
            idx_to_state[idx] = st

        transition = np.zeros([len(states),len(states)])
        for (idx, s) in enumerate(states):
            var_part = self.state_to_partition(dict(s))
            partition = var_part
            stab = self.strip_factors(self.graph.automorphism_group(partition=partition + self.factors))
            order = stab.order()
            gelems = stab.list()
            for e in gelems:
                row_states = [dict()]
                # build of list of states which are transitioned to
                cyc_tuples = e.cycle_tuples(singletons=True)
                for cyc in cyc_tuples:
                    new_states = []
                    for state in row_states:
                        st1 = state.copy()
                        st2 = state.copy()
                        for i in cyc:
                            st1[i] = True
                            st2[i] = False
                        new_states.append(st1)
                        new_states.append(st2)
                    row_states = new_states
                for st in row_states:
                    l = st.items()
                    l.sort()
                    cur_idx = state_to_idx[tuple(l)]
                    transition[idx,cur_idx] += (1.0 / order) * (1.0 / (2**len(cyc_tuples)))
        return transition

    ### computes the transition matrix for a Gibbs step
    def gibbs_transition(self):
        states = self.gen_all_states()
        state_to_idx = dict()
        idx_to_state = dict()
        for (idx, st) in enumerate(states):
            state_to_idx[st] = idx
            idx_to_state[idx] = st

        transition = np.zeros([len(states),len(states)])
        for (idx, s) in enumerate(states):
            for v in self.variables:
                # compute two resulting new states
                st1 = dict(s)
                st2 = dict(s)
                st1[v] = True
                st2[v] = False

                st1prob = self.potential(st1)
                st1 = st1.items()
                st1.sort()
                st1 = tuple(st1)
                st1idx = state_to_idx[st1]

                st2prob = self.potential(st2)
                st2 = st2.items()
                st2.sort()
                st2 = tuple(st2)
                st2idx = state_to_idx[st2]
                if st1prob + st2prob > 0:
                    transition[idx, st1idx] += ((1.0/ len(self.variables)) * float(st1prob)/ (st1prob + st2prob))
                    transition[idx, st2idx] += ((1.0 / len(self.variables)) * float(st2prob)/ (st1prob + st2prob))
                else:
                    # stay
                    s_idx = state_to_idx[s]
                    transition[s_idx, st1idx] += 1.0 / len(self.variables)
        return transition

    ### compute transition matrix of a within-orbit step, for implementing
    ### lifted MCMC
    def orbit_transition(self):
        states = self.gen_all_states()
        state_to_idx = dict()
        idx_to_state = dict()
        for (idx, st) in enumerate(states):
            state_to_idx[st] = idx
            idx_to_state[idx] = st

        transition = np.zeros([len(states),len(states)])
        for (idx, s) in enumerate(states):
            st_dict = dict(s)
            stripped = self.strip_factors(self.graph_aut)
            for g in stripped.list():
                cur_s = dict()
                for cyc in g.cycle_tuples(singletons=True):
                    for i, var in enumerate(cyc):
                        cur_s[var] = st_dict[cyc[(i + 1) % len(cyc)]]
                st = cur_s.items()
                st.sort()
                st = tuple(st)

                new_idx = state_to_idx[st]
                transition[idx, new_idx] += 1.0 / self.graph_aut_order
        return transition

    ### computes the metropolis hastings transition matrix for a burnside
    ### proposal of `burnsidesteps` steps
    def burnside_mh_transition(self, burnsidesteps):
        # now apply the metropolis correction for the transition x -> y, where
        # each entry in `transition` contains the burnside transition
        # probability
        states = self.gen_all_states()
        state_to_idx = dict()
        idx_to_state = dict()
        for (idx, st) in enumerate(states):
            state_to_idx[st] = idx
            idx_to_state[idx] = st

        aut_cache = dict()
        B = np.linalg.matrix_power(self.burnside_transition(), burnsidesteps)
        for x in range(0, len(states)):
            for y in range(0, len(states)):
                st_x = idx_to_state[x]
                part_x = self.state_to_partition(dict(st_x))
                orbx = None
                if st_x in aut_cache:
                    orbx = aut_cache[st_x]
                else:
                    orbx = self.graph_aut_order / \
                           self.graph.automorphism_group(partition=self.state_to_partition(dict(st_x)) + self.factors).order()
                    aut_cache[st_x] = orbx

                prx = self.potential(dict(st_x))
                st_y = idx_to_state[y]
                orby = None
                if st_y in aut_cache:
                    orby = aut_cache[st_y]
                else:
                    orby = self.graph_aut_order / \
                           self.graph.automorphism_group(partition=self.state_to_partition(dict(st_y)) + self.factors).order()
                    aut_cache[st_y] = orby

                # orby = self.graph_aut_order / self.graph.automorphism_group(partition=self.state_to_partition(dict(st_y))).order()
                pry = self.potential(dict(st_y))
                if prx != 0.0:
                    B[x,y] = B[x,y] * np.minimum(1, (pry * orby) / (prx * orbx))
            # compute probability of staying in the current position
            prob_stay = 0.0
            for yp in range(0, len(states)):
                prob_stay += B[x,yp]
            prob_stay = 1.0 - prob_stay
            B[x,x] = B[x,x] + prob_stay
        return B

    def brute_force_prob_vector(self):
        states = self.gen_all_states()
        state_to_idx = dict()
        idx_to_state = dict()
        for (idx, st) in enumerate(states):
            state_to_idx[st] = idx
            idx_to_state[idx] = st
        Z = 0.0
        for st in states:
            Z += self.potential(dict(st))
        # return the vector
        vec = np.zeros(len(states))
        for st in states:
            vec[state_to_idx[st]] = self.potential(dict(st)) / Z

        return vec

    ### compute the total variation distance of the factor graph from its
    ### stationary distribution
    ### for each time step t, print d_tv(stationary_dist, transition_matrix^t * starting_point)
    ### transition_matrix: describes transition probabilities between states
    ### starting_point: initial probability vector on states
    ### num_steps: total number of steps for which to run
    def total_variation(self, transition_matrix, starting_point, num_steps):
        v = self.brute_force_prob_vector()
        T = transition_matrix
        for i in range(1, num_steps):
            T = np.matmul(T, transition_matrix)
            p = starting_point.dot(T)
            d = 0.0
            for idx in range(0, len(v)):
                d += abs(v[idx] - p[idx])
            print("%s\t%s" % (i, d))

    ### perform a single step of orbit jumping
    ### returns: a pair, (the ratio of transition probabilities, new state)
    ### n: number of burnside steps to take
    def orbitjump(self, state, n):
        hatx = self.burnside(state, n)
        probx = self.potential(state)
        probhatx = self.potential(hatx)
        orbx = self.graph_aut_order / \
               self.graph.automorphism_group(partition=self.state_to_partition(state) + self.factors).order()
        orbhatx = self.graph_aut_order / \
                  self.graph.automorphism_group(partition=self.state_to_partition(hatx) + self.factors).order()
        try:
            transitionprob = (probhatx * orbhatx) / (probx * orbx)
            return (transitionprob, hatx)
        except:
            # divided by 0
            return (1.0, hatx)

    ### perform a standard orbital MCMC gibbs update, a la Niepert 2012
    ### state: a state
    ### returns: the gibbs update state
    def orbitgibbs(self, state):
        v = sage.misc.prandom.choice(self.variables)
        # resample the variable
        v_true = state.copy()
        v_true[v] = True
        v_false = state.copy()
        v_false[v] = False
        prob = None
        try:
            prob = self.potential(v_true) / (self.potential(v_true) + self.potential(v_false))
        except:
            return state # evidence was not satisfied
        new_v = numpy.random.binomial(1, prob)
        state[v] = new_v

        # now walk along the orbit
        g = fast_random_element(self.graph_aut)
        # apply g to the state
        for cyc in g:
            fst = state[cyc[0]]
            for idx, var in enumerate(cyc):
                state[var] = state[cyc[(idx + 1) % len(cyc)]]
            state[cyc[-1]] = fst
        return state

    ### draw n samples according to orbit jump MCMC with no orbital MCMC
    ### burnsidesize: number of burnside steps to take
    ### n: int, total number of samples to take
    ### gamma: int, number of steps between taking orbit jumps
    ### burn: the burn in of the chain
    ### returns the probability of the query
    def orbitjumpmcmc(self, n, query, burnsidesize=10, gamma=10, burn=100):
        samples = []
        # set up initial random state
        cur_state = dict()
        for v in self.variables:
            cur_state[v] = flip()

        query_count = 0.0
        cur_step = 0

        for i in range(0, n * burn):
            if (cur_step % gamma) == 0:
                # do a jump update
                (ratio, new_state) = self.orbitjump(cur_state, burnsidesize)
                # accept with transition probability
                acceptprob = numpy.minimum(1, ratio)
                if numpy.random.binomial(1, acceptprob):
                    cur_state = new_state
            else:
                # do a gibbs update
                cur_state = self.orbitgibbs(cur_state)

            if cur_step % burn == 0:
                if query(cur_state):
                    query_count += 1
            cur_step += 1
        return query_count / n

if __name__ == "__main__":
    model = gen_complete_pairwise_factorgraph_half(6)
    print(model)
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
    model.total_variation(orbitalmcmc, start, 100)
    print("------------------")
    print("pure jump")
    model.total_variation(burnside, start, 100)


    # print fg.partition()
    # print fg.brute_force_partition()
