### orbit generation for Markov networks
### all cliques are assumed to be symmetric and identical
### this file is for experimenting; it is not what was used in the paper

from sage.all import *
from my_graphs import *
import cProfile, pstats, StringIO
from collections import deque
import my_bliss

### generate all orbits of `elements` under `gens`
def generate_orbits(elements, gens):
    vset = set(elements)
    seenset = set()
    orbits = []
    def step_point(point):
        res = []
        for gen in gens:
            # a generator is a list of cycles, so we go one level deeper
            for cyc in gen:
                try:
                    i = cyc.index(point)
                    res.append(cyc[(i + 1) % len(cyc)])
                except:
                    continue
        return res

    def dfs_find_orbit(gens, point, cur_orbit):
        frontier = step_point(point)
        cur_orbit.update(frontier)
        # iterate until fixed point
        for new_point in cur_orbit - set(frontier):
            cur_orbit.update(dfs_find_orbit(gens, new_point, cur_orbit))
        return cur_orbit

    while vset != seenset:
        diff = vset - seenset
        cur_point = diff.pop()
        new_orbit = dfs_find_orbit(gens, cur_point, set([cur_point]))
        seenset.update(new_orbit)
        orbits.append(new_orbit)

    return orbits

### bfs_foldrep
### visits all possible isomorphic graph colorings in `colors`. The colors in
### `fixcolors` are fixed throughout.
### if `orders` is True, provide `f` with the orders of the colored and
### uncolored graph automorphism group
def bfs_foldrep(graph, colors, fixcolors, acc, f, orders=False):
    # queue of colorings yet to be considered
    queue = deque()
    queue.append([[], colors])
    g_order = None
    if orders:
        g_order = graph.automorphism_group().order()
    # set of representative graph colorings
    reps = set()
    while len(queue) != 0:
        # print("queue: %s" % queue)
        c = queue.popleft()
        # TODO: turn this into a single GI call by having it return both the
        # automorphism group and the canonical fora
        part = c + fixcolors
        gcanon, cert = my_bliss.canonical_form(graph, partition=part, certificate=True,
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

        A, orbits = graph.automorphism_group(partition=part, orbits=True,
                                             algorithm="bliss")
        if orders:
            acc = f(acc, graph, c, g_order, A.order())
        else:
            acc = f(acc, graph, c)
        reps.add((gcanon, (tuple(c_canon[0]), tuple(c_canon[1]))))
        # print(len(reps))
        if len(c_canon[0]) + 1 <= len(colors) / 2.0:
            # if we can add more colors, try to

            # gens = my_bliss.raw_automorphism_generators(graph, partition=c + fixcolors)
            # orbits = generate_orbits(colors, gens)
            # if we need the full automorphism group, we can compute it here.
            # expand this node and add it to the queue
            for o in orbits:
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
    # folding function
    def add_vertex(acc, g, color):
        if len(color[0]) < len(varnodes) / 2.0:
            return acc + [color] + [color[::-1]]

        return acc + [color]
    return bfs_foldrep(G, varnodes, colornodes, [], add_vertex)


# given a graph g, count the number of distinct 2-colorings using Polya's
# theorem.
# TODO: This can in be done more efficiently without expanding the cycle index
# polynomial on 2 variables, but I can't figure out how.
def count_num_distinct(g):
    return g.automorphism_group().cycle_index().expand(2)((1,1))

def find_representatives():
    # n=3
    pr = cProfile.Profile()
    pr.enable()

    # G = graphs.EmptyGraph()
    # for v in range(0, 300):
    #     G.add_vertex(v)
    G = gen_complete_pairwise_factor(30)
    # G = graphs.CompleteGraph(5)
    # G = gen_pigeonhole(5,5)
    # G = graphs.CycleGraph(20)
    # G = gen_complete_extra(20)
    # G = gen_friends_smokers_factor(6)
    # G = G.complement()
    # r = genrep_bfs(G)
    r = genrep_bfs_factor(G[0], G[1][0], [G[1][1]])
    for x in r:
        print(x)
    print("Found configurations: %d" % len(r))
    # print("Expected configurations: %d" % count_num_distinct(G))

    pr.disable()
    s = StringIO.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    # print s.getvalue()


if __name__ == "__main__":
    set_gap_memory_pool_size(900000000)
    find_representatives()
