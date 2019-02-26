from sage.all import *
from my_graphs import *
import cProfile, pstats, StringIO
from collections import deque
import my_bliss

### bfs_foldrep
### visits all possible isomorphic graph colorings in `colors`. The colors in
### `fixcolors` are fixed throughout.
def bfs_foldrep(graph, colors, fixcolors, acc, f):
    # queue of colorings yet to be considered
    queue = deque()
    queue.append([[], colors])
    # set of representative graph colorings
    reps = set()
    while len(queue) != 0:
        c = queue.popleft()
        # print("----")
        # print("Current color: %s" % (c + fixcolors))
        # TODO: turn this into a single GI call by having it return both the
        # automorphism group and the canonical form
        gcanon, cert = my_bliss.canonical_form(graph, partition=c + fixcolors, certificate=True,
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
        acc = f(acc, graph, c)
        reps.add((gcanon, (tuple(c_canon[0]), tuple(c_canon[1]))))

        if len(c_canon[0]) + 1 <= len(colors) / 2.0:
            # if we can add more colors, try to

            gens = my_bliss.raw_automorphism_generators(graph, partition=c + fixcolors)
            # big speed hack: compute orbits from generators. avoid calling GAP for this.
            vset = set(colors)
            seenset = set()
            cur_orbit = set()
            cur_point = graph.vertices()[0]
            orbits = []
            while True:
                if cur_point in cur_orbit:
                    orbits += [cur_orbit.copy()]
                    cur_orbit = set()
                    # see if there are any more vertices
                    remaining = vset - seenset
                    if len(remaining) == 0:
                        break
                    cur_point = remaining.pop()
                seenset.add(cur_point)
                cur_orbit.add(cur_point)
                # apply a generator to this point
                for gen in gens:
                    # a generator is a list of cycles, so we go one level deeper
                    for cyc in gen:
                        try:
                            i = cyc.index(cur_point)
                            cur_point = cyc[(i + 1) % len(cyc)]
                            break
                        except:
                            continue

            # if we need the full automorphism group, we can compute it here.
            # A, orbits = graph.automorphism_group(partition=c + fixcolors, orbits=True,
            #                                      algorithm="bliss")
            # expand this node and add it to the queue
            for o in orbits:
                # print("Considering orbit %s" % o)
                e = o.pop() # pick an element in the orbit arbitrarily
                # check if this orbit is already true, or if it is any of the fixed colors
                if e in c[0] or e in [y for x in fixcolors for y in x]:
                    continue
                # we know this is an uncolored vertex
                newcolor = [c[0] + [e], [x for x in c[1] if x != e]]
                queue.append(newcolor)
    return acc

def genrep_bfs(G):
    colors = G.vertices()
    # folding function
    def add_vertex(acc, g, color):
        if len(color[0]) < len(g.vertices()) / 2.0:
            return acc + [color] + [color[::-1]]

        return acc + [color]
    return bfs_foldrep(G, colors, [], [], add_vertex)


def genrep_bfs_factor(G, varnodes, colornodes):
    colors = varnodes
    # folding function
    def add_vertex(acc, g, color):
        if len(color[0]) < len(varnodes) / 2.0:
            return acc + [color] + [color[::-1]]

        return acc + [color]
    return bfs_foldrep(G, colors, colornodes, [], add_vertex)

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
    pr = cProfile.Profile()
    pr.enable()

    # G = graphs.EmptyGraph()
    # for v in range(0, 300):
    #     G.add_vertex(v)
    # G = gen_complete_pairwise_factor(100)
    # G = graphs.CompleteGraph(5)
    # G = gen_friends_smokers(10)
    # G = graphs.CycleGraph(20)
    # G = gen_complete_extra(20)
    G = gen_friends_smokers_factor(7)
    # G = G.complement()
    # r = genrep_bfs(G)
    r = genrep_bfs_factor(G[0].complement(), G[1][0], [G[1][1]])
    for x in r:
        print(x)
    print("Found configurations: %d" % len(r))
    # print("Expected configurations: %d" % count_num_distinct(G))

    pr.disable()
    s = StringIO.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print s.getvalue()


if __name__ == "__main__":
    find_representatives()
