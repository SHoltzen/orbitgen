from sage.all import *
import cProfile, pstats, StringIO
import itertools
from sage.groups.perm_gps.partn_ref.refinement_graphs import isomorphic

def findsubsets(S,m):
    return set(itertools.combinations(S, m))

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
    A, orbits = graph.automorphism_group(partition=colors, orbits=True, algorithm="bliss")
    for o in orbits:
        e = o[0] # pick the first element in the orbit arbitrarily
        # check if this orbit is already true
        if e in colors[0]:
            continue
        # we know this is an uncolored vertex
        Dcolor = [colors[0], [c for c in colors[1] if c != e], [e]]
        D = graph.canonical_label(partition=Dcolor)
        Gpcolor = [colors[0] + [e], [c for c in colors[1] if c != e]]
        Gprime, c = graph.canonical_label(partition=Gpcolor,
                                          certificate=True, algorithm="bliss")

        # c is a map from old vertices to new vertices; we need to invert it
        inv_c = {v: k for k, v in c.iteritems()}

        # find lexically first vertex whose value is true in the canonical Gprime;
        # default to vertex 0
        first_t = next((v for v in Gprime.vertices() if v not in Gpcolor[1]), Gprime.vertices()[0])

        # need to convert these colorings into colorings of original vertices
        canoncolor = [[x for x in Gpcolor[0] if x != inv_c[first_t]],
                      Gpcolor[1],
                      [inv_c[first_t]]]
        Dcanon = graph.canonical_label(partition=canoncolor, algorithm="bliss")
        if Dcanon == D:
            acc = foldrep(graph, Gpcolor, acc, f)
    return acc

### genrep: generates a representative of each orbit class of a graph
def genrep(G):
    colors = [[], [i for i in G.vertices()]]
    # folding function
    def add_vertex(acc, g, color):
        # check if inversion is isomorphic to the current graph
        if len(color[0]) < len(g.vertices()) / 2.0:
            return acc + [color] + [color[::-1]]

        # we know it's a 50/50 color split, check if this is an isomorphic
        # inversion or not
        g2 = Graph(g)
        g.canonical_label(partition=color)
        g2.canonical_label(partition=color[::-1])
        if g == g2:
            # inversions are isomorphic, return one
            return acc + [color]
        else:
            # inversion are non-isomorphic, return both
           print("non isomorphic inversion: %s, %s" % (color, color[::-1]))
           return acc + [color] + [color[::-1]]
    return foldrep(G, colors, [], add_vertex)

### generates a friends and smokers graph with n people
def gen_friends_smokers(n):
    g = Graph(sparse=True)
    # make n smoker vertices
    smokers = [x for x in range(0,n)]
    # connect all the smokers
    smokeredges = findsubsets(smokers, 2)
    # make friends
    friends = []
    friendedges = []
    count = n
    for (s1,s2) in findsubsets(smokers, 2):
        friends += [count]
        friendedges += [(s1, count), (s2, count)]
        count += 1

    g.add_vertices(smokers)
    g.add_vertices(friends)
    g.add_edges(friendedges)
    g.add_edges(smokeredges)
    return g

# generates a graph which is fully-connected in m and connected across n
def gen_pigeonhole(n,m):
    g = Graph()
    v = []
    e = []
    # generate vertices
    for x in range(0,n):
        for y in range(0,m):
            v += [x*m + y]

    # generate edges
    # generate fully connected graph in n
    for x in range(0,n):
        for y in findsubsets(range(0, m), 2):
            e += [(x*m + y[0], x*m + y[1])]

    # connect between n and m
    for x in range(0,n):
        for y in range (0,m):
            e += [(x*m + y, ((x + 1) % n)*m + y)]

    g.add_vertices(v)
    g.add_edges(e)
    return g

# given a graph g, count the number of distinct 2-colorings using Polya's
# theorem.
# TODO: This can in be done more efficiently without expanding the cycle index
# polynomial on 2 variables, but I can't figure out how.
def count_num_distinct(g):
    return g.automorphism_group().cycle_index().expand(2)((1,1))

def main():
    # n=3
    # pr = cProfile.Profile()
    # pr.enable()

    # G = gen_friends_smokers(3)
    G = graphs.CycleGraph(8)
    num_vert = len(G.vertices())
    r = genrep(G)
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
    main()
