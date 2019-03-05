from sage.all import *
import numpy.random
from collections import deque
import random
import numpy as np
import my_bliss
from my_graphs import *
import cProfile, pstats, StringIO
import itertools

def findsubsets(S,m):
    return set(itertools.combinations(S, m))

def flip():
    return numpy.random.binomial(1, 0.5)

### returns a random element from a group g using product replacement
def fast_random_element(g):
    stab = libgap(g)
    return (libgap.PseudoRandom(stab).sage(parent=g)).cycle_tuples(singletons=True)


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


    ### perform a query via exhaustive enumeration using orbit generation
    ### query: a function state -> bool
    ### Z: if true, return a tuple (prob, Z)
    def query_enumerate(self, query, Z=False):
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
        if Z:
            return (prob / Z, Z)
        else:
            return prob / Z


    ### computes a burnside process transition beginning from a state particular state
    ### n: number of moves
    def burnside(self, state, n):
        for j in range(0,n):
            var_part = self.state_to_partition(state)
            partition = var_part
            stab = self.graph.automorphism_group(partition=partition)
            p = fast_random_element(stab)
            #p = stab.random_element()cycle_tuples(singletons=True)

            # now convert p into a state
            # print(p)
            state = dict()
            for cyc in p:
                # sample a value for cyc and fill in the new current value
                v = flip()
                for idx in cyc:
                    state[idx] = v
        return state

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
        tbl = add_vertex([dict()], self.graph.vertices())
        res = []
        for itm in tbl:
            l = itm.items()
            l.sort()
            res.append(tuple(l))
        return res


    ### computes the transition matrix of a single step of the burnside process
    ### Z: the normalizing constant for the distribution
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
            stab = self.graph.automorphism_group(partition=partition)
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

        B = np.linalg.matrix_power(self.burnside_transition(), burnsidesteps)
        for x in range(0, len(states)):
            for y in range(0, len(states)):
                st_x = idx_to_state[x]
                orbx = self.graph.automorphism_group(partition=self.state_to_partition(dict(st_x))).order()
                prx = self.potential(dict(st_x))
                st_y = idx_to_state[y]
                orby = self.graph.automorphism_group(partition=self.state_to_partition(dict(st_y))).order()
                pry = self.potential(dict(st_y))
                if prx != 0.0:
                    B[x,y] = B[x,y] * np.minimum(1, (pry * orby) / (prx * orbx))
        return B

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

    # def total_variation(self, burnsidesize):

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



    ### proportion of support of the distribution explored after n steps,
    ### starting in a random initial state.
    ###
    ### n: number of steps
    ### lag: number of steps between each coverage report
    ### burnsidesize: number of burnside steps to take
    ### gamma: int, number of steps between taking orbit jumps
    ### Z: optionally provide the normalizing constant, to avoid recomputing it
    ### returns: a list of proportions, one for each lag step
    def support_explored(self, n, burnsidesize=10, gamma=10, Z=None):
        def trivial_query(state):
            return True
        if Z is None:
            Z = self.query_enumerate(trivial_query)

        states = set()

        # set up initial random state
        cur_state = dict()
        for v in self.variables:
            cur_state[v] = flip()

        cur_step = 0
        supp_explored = 0.0

        for i in range(0, n):
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

            state_tup = cur_state.items()
            state_tup.sort()
            state_tup = tuple(state_tup)

            if state_tup not in states:
                states.add(state_tup)
                supp_explored += self.potential(cur_state) / Z
            cur_step += 1
        return supp_explored


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


def experiment_motivating():
    # print some burnside samples
    nlist = [25, 50, 75, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    # nlist = [25]

    g = graphs.CompleteGraph(6)
    def potential(state):
        num_t = 0.0
        for v in state.itervalues():
            if v:
                num_t += 1

        if num_t == 0 or num_t == 6:
            return 100
        if num_t == 1 or num_t == 5:
            return 6
        if num_t == 2 or num_t == 4:
            return 7
        if num_t == 3:
            return 20

    def trivial(state):
        return True

    graph = MarkovModel(g, g.vertices(), potential)
    prob, Z = graph.query_enumerate(trivial, Z=True)
    for n in nlist:
        res = []
        for i in range(0, 15):
            v1 = graph.support_explored(n, burnsidesize=4, gamma=20, Z=Z)
            res.append(v1)
        print("%s\t%s\t%s" % (n, numpy.average(res), numpy.std(res)))

def experiment_friends_smokers():
    # make a friends/smokers markov model
    g = Graph(sparse=true)
    num_smokers = 5
    # make n smoker vertices
    smokers = [x for x in range(0,num_smokers)]
    # connect all the smokers
    smokeredges = findsubsets(smokers, 2)
    # make friends
    friends = []
    friendedges = []
    count = num_smokers
    for (s1,s2) in findsubsets(smokers, 2):
        friends += [count]
        friendedges += [(s1, count), (s2, count)]
        count += 1

    g.add_vertices(smokers)
    g.add_vertices(friends)
    g.add_edges(friendedges)
    g.add_edges(smokeredges)


    # print some burnside samples
    nlist = [25, 50, 75, 100, 200, 300, 400, 500, 600, 700, 800, 900]

    def potential(state):
        count = num_smokers
        total = 0.0
        for (s1,s2) in findsubsets(smokers, 2):
            # friend vertex == counta
            # add 3 if S(s1) /\ F(s1,s2) => S(s2), 1 otherwise
            if state[s1] == state[s2]:
                total += 1000000
            else:
                total += 1
            count += 1
        return total

    def trivial(state):
        return True

    graph = MarkovModel(g, g.vertices(), potential)
    prob, Z = graph.query_enumerate(trivial, Z=True)
    for n in nlist:
        res = []
        for i in range(0, 15):
            v1 = graph.support_explored(n, burnsidesize=8, gamma=n+1, Z=Z)
            res.append(v1)
        print("%s\t%s\t%s" % (n, numpy.average(res), numpy.std(res)))



def experiment_pigeonhole():
    # make a friends/smokers markov model
    nholes = 2
    npigeons = 5
    g = gen_pigeonhole(nholes,npigeons)
    # print some burnside samples
    nlist = [50,
             75,
             100,
             200,
             300,
             400,
             500,
             600,
             700,
             800,
             900]

    def potential(state):
        total = 0.0
        distinct_holes = set()
        for p in range(0, npigeons):
            num_holes = 0
            for h in range(0, nholes):
                distinct_holes.add(h)
                if state[(h, p)]:
                    num_holes += 1
            if num_holes == 1:
                total += 100
            else:
                total -= 100

            # penalize if many distinct holes are occupied
            total -= 100 * len(distinct_holes)
        return math.exp(total)

    def trivial(state):
        return True

    graph = MarkovModel(g, g.vertices(), potential)
    prob, Z = graph.query_enumerate(trivial, Z=True)
    for n in nlist:
        res = []
        for i in range(0, 5):
            v1 = graph.support_explored(n, burnsidesize=8, gamma=5, Z=Z)
            res.append(v1)
        print("%s\t%s\t%s" % (n, numpy.average(res), numpy.std(res)))


def motivating_example():
    # print some burnside samples
    total_people = 3
    g = graphs.CompleteGraph(total_people)
    def potential(state):
        num_t = 0.0
        for v in state.itervalues():
            if v:
                num_t += 1

        if num_t == 0 or num_t == 6:
            return 100
        if num_t == 1 or num_t == 5:
            return 4
        if num_t == 2 or num_t == 4:
            return 5
        if num_t == 3:
            return 20

    # pr = cProfile.Profile()
    # pr.enable()
    graph = MarkovModel(g, g.vertices(), potential)
    for n in range(0, total_people+1):
        def query(state):
            num_true = 0
            # print(r)
            for k,v in state.iteritems():
                if v:
                    num_true += 1
            return (num_true == n)
        print("%s: %s" % (n, graph.query_enumerate(query)))
    # v1 = graph.orbitjumpmcmc(20, query, gamma=1, burn=2, burnsidesize=5)
    # print("exact: %s, approx: %s" % (exact, v1))

def test():
    def potential(state):
        num_t = 0.0
        for v in state.itervalues():
            if v:
                num_t += 1
        return num_t
    g = graphs.CompleteGraph(3)
    graph = MarkovModel(g, g.vertices(), potential)
    print(graph.burnside_mh_transition(4))

if __name__ == "__main__":
    test()
