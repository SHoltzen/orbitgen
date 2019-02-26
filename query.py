from sage.all import *
from my_graphs import *
import cProfile, pstats, StringIO
from collections import deque
import my_bliss
from orbitgen import *

### partition_function
### computes the partition of a fully symmetric MLN the specified parameter
###
### potential: a function which maps variable truth assignments to
### probabilities. This function *must* be invariant wrt. the automorphism
### group of the graph.
def partition(G, variables, factors, potential):
    # folding function
    print("G: %s, vars: %s, factors: %s" % (G, variables, factors))
    def sum_prob(acc, g, color, A):
        # TODO: avoid recomputing these two automorphism groups
        g_order = g.automorphism_group().order()
        aut_order = A.order()
        orbit_sz = g_order / aut_order # yay orbit stabilizer theorem!
        if len(color[0]) < len(g.vertices()) / 2.0:
            return acc + (orbit_sz * potential(color)) + (orbit_sz * potential(color[::-1]))

        return acc + (orbit_sz * potential(color))
    return bfs_foldrep(G, variables, factors, 0.0, sum_prob)


def compute_partition():
    G = gen_complete_pairwise_factor(50)
    def potential(assgn):
        # count potential: potential is # incoming nodes which are true
        p = 0.0
        for factor in G[1][1]:
            connected = G[0].neighbors(factor)
            for v in connected:
                if v in assgn[0]:
                    p += 1
        return p

    pr = cProfile.Profile()
    pr.enable()

    part = partition(G[0], G[1][0], [G[1][1]], potential)
    print("partition: %d" % part)
    # print("brute forced: %d" % (bruteforce_partition(G, partfun)))

    pr.disable()
    s = StringIO.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print s.getvalue()



def mpe(G, variables, factors, potential):
    # folding function
    print("G: %s, vars: %s, factors: %s" % (G, variables, factors))
    def sum_prob(acc, g, color):
        (cur_max_state, cur_max_value) = acc
        if len(color[0]) < len(g.vertices()) / 2.0:
            # check both possibilities
            pot = potential(color[::-1])
            if cur_max_value < pot:
                cur_max_value = pot
                cur_max_state = color[::-1]
        pot = potential(color)
        if cur_max_value < pot:
            cur_max_value = pot
            cur_max_state = color
        return (cur_max_state, cur_max_value)
    return bfs_foldrep(G, variables, factors, (None, 0.0), sum_prob)


def compute_mpe_complete_2factor():
    G = gen_complete_pairwise_factor(10)
    # add evidence: the first variable is false
    G[0].add_vertex(name="e")
    G[0].add_edge(("e", G[1][0][0]))
    def potential(assgn):
        p = 0.0
        # check for evidence
        if G[1][0][0] in assgn[0]:
            return 0.0
        for factor in G[1][1]:
            connected = G[0].neighbors(factor)
            if connected[0] != connected[1]:
                p += 0.3
            else:
                p += 0.2
            for v in connected:
                if v in assgn[0]:
                    p += 1
        return p

    pr = cProfile.Profile()
    pr.enable()

    res = mpe(G[0], G[1][0], [G[1][1] + ["e"]], potential)
    print("max state: %s, max value: %s" % (res[0], res[1]))
    # print("brute forced: %d" % (bruteforce_partition(G, partfun)))

    pr.disable()
    s = StringIO.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print s.getvalue()





if __name__ == "__main__":
    compute_mpe_complete_2factor()
