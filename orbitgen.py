from sage.all import *
from my_graphs import *
import cProfile, pstats, StringIO
from sage.groups.perm_gps.partn_ref.refinement_graphs import isomorphic
from collections import deque

### foldrep: applies f to each canonical representative, maintains an
### accumulator
### graph: a sage graph
### colors: a list of lists which partitions the vertices in G into colors
###   colors[0]: True vertices
###   colors[1]: False vertices
###   colors[3] is reserved for the algorithm
### fold: apply a fold method to each orbit
def foldrep(graph, colors, acc, f):
    # TODO: reduce number of isomorphism calls if possible
    acc = f(acc, graph, colors)
    if len(colors[0]) >= (len(graph.vertices()) / 2):
        return acc
    A, orbits = graph.automorphism_group(partition=colors, orbits=True)
    for o in orbits:
        e = o[0] # pick the first element in the orbit arbitrarily
        # check if this orbit is already true
        if e in colors[0]:
            continue
        # we know this is an uncolored vertex
        Dcolor = [colors[0], [c for c in colors[1] if c != e], [e]]
        D = graph.canonical_label(partition=Dcolor, algorithm="bliss")

        Gpcolor = [colors[0] + [e], [c for c in colors[1] if c != e]]
        Gprime, c = graph.canonical_label(partition=Gpcolor,
                                           certificate=True)

        # c is a map from old vertices to new vertices; we need to invert it
        inv_c = {v: k for k, v in c.iteritems()}

        # compute the canonical deletion of Gprime
        # find lexically first vertex whose value is true in the canonical
        # Gprime; default to vertex 0.
        first_t = None
        for v in Gprime.vertices(sort=True):
            if inv_c[v] in Gpcolor[0]:
                first_t = inv_c[v]
                break
        if first_t is None:
            first_t = Gprime.vertices(sort=True)[0]

        canoncolor = [[x for x in Gpcolor[0] if x != first_t],
                      Gpcolor[1],
                      [first_t]]
        Dcanon = graph.canonical_label(partition=canoncolor)
        if Dcanon == D:
            acc = foldrep(graph, Gpcolor, acc, f)
    return acc

def bfs_foldrep(graph, colors, acc, f):
    # queue of colorings yet to be considered
    queue = deque([colors])
    # set of representative graph colorings
    reps = set()
    while len(queue) != 0:
        c = queue.popleft()
        # print("----")
        # print("Current color: %s" % c)
        gcanon, cert = graph.canonical_label(partition=c, certificate=True)
        G_immut = Graph(gcanon, immutable=True)

        # convert the colors to their coloring in the canonical graph
        c_canon = [[], []]
        for c1 in c[0]:
            c_canon[0].append(cert[c1])
        for c1 in c[1]:
            c_canon[1].append(cert[c1])

        c_canon[0].sort()
        c_canon[1].sort()

        if (G_immut, (tuple(c_canon[0]), tuple(c_canon[1]))) in reps:
            continue
        acc = f(acc, graph, c)
        reps.add((G_immut, (tuple(c_canon[0]), tuple(c_canon[1]))))
        # print("added rep: %s" % reps)

        # print(c)
        if len(c_canon[0]) + 1 <= len(graph.vertices()) / 2.0:
            # if we can add more colors, try to
            A, orbits = graph.automorphism_group(partition=c, orbits=True)
            # print("orbits: %s" % orbits)
            # expand this node and add it to the queue
            for o in orbits:
                # print("Considering orbit %s" % o)
                e = o[0] # pick the first element in the orbit arbitrarily
                # check if this orbit is already true
                if e in c[0]:
                    continue
                # we know this is an uncolored vertex
                newcolor = [c[0] + [e], [x for x in c[1] if x != e]]
                queue.append(newcolor)
    return acc

### genrep: generates a representative of each orbit class of a graph
def genrep(G):
    colors = [[], [i for i in G.vertices()]]
    # folding function
    def add_vertex(acc, g, color):
        if len(color[0]) < len(g.vertices()) / 2.0:
            return acc + [color] + [color[::-1]]

        return acc + [color]
    return foldrep(G, colors, [], add_vertex)

def genrep_bfs(G):
    colors = [[], [i for i in G.vertices()]]
    # folding function
    def add_vertex(acc, g, color):
        if len(color[0]) < len(g.vertices()) / 2.0:
            return acc + [color] + [color[::-1]]

        return acc + [color]
    return bfs_foldrep(G, colors, [], add_vertex)


### partition_function
### computes the partition of a fully symmetric MLN the specified parameter
###
### potential: a function which maps variable truth assignments to
### probabilities. This function *must* be invariant wrt. the automorphism
### group of the graph.
def partition(G, potential):
    colors = [[], [i for i in G.vertices()]]
    # folding function
    def sum_prob(acc, g, color):
        # TODO: avoid recomputing these two automorphism groups
        g_order = g.automorphism_group().order()
        aut_order = g.automorphism_group(partition=color).order()
        orbit_sz = g_order / aut_order # yay orbit stabilizer theorem!
        if len(color[0]) < len(g.vertices()) / 2.0:
            return acc + (orbit_sz * potential(color)) + (orbit_sz * potential(color[::-1]))

        return acc + (orbit_sz * potential(color))
    return foldrep(G, colors, 0.0, sum_prob)

def bruteforce_partition(G, potential):
    # evaluate the potential on every assignment to variables
    def recurse(c1, c2):
        print(c1)
        if len(c1) == 0:
            return 0
        tot = potential([c1, c2])
        for x in c1:
            new_c1 = c1[:]
            new_c1.remove(x)
            tot += recurse(new_c1, c2 + [x])
        return tot
    return recurse(G.vertices(), [])

# given a graph g, count the number of distinct 2-colorings using Polya's
# theorem.
# TODO: This can in be done more efficiently without expanding the cycle index
# polynomial on 2 variables, but I can't figure out how.
def count_num_distinct(g):
    return g.automorphism_group().cycle_index().expand(2)((1,1))


def compute_partition():
    def partfun(assgn):
        return len(assgn[0])
    G = gen_complete_extra(3)
    print("partition: %d" % (partition(G, partfun)))
    print("brute forced: %d" % (bruteforce_partition(G, partfun)))

def find_representatives():
    # n=3
    # pr = cProfile.Profile()
    # pr.enable()

    # G = graphs.CompleteGraph(4)
    G = gen_friends_smokers(10)
    # G = graphs.CycleGraph(20)
    # G = gen_complete_extra(3)
    num_vert = len(G.vertices())
    r = genrep_bfs(G)
    for x in r:
        print(x)
    print("Found configurations: %d" % len(r))
    print("Expected configurations: %d" % count_num_distinct(G))

    # pr.disable()
    # s = StringIO.StringIO()
    # sortby = 'cumulative'
    # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    # ps.print_stats()
    # print s.getvalue()


if __name__ == "__main__":
    find_representatives()
