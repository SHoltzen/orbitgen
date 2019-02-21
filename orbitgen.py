from sage.all import *
import cProfile, pstats, StringIO

# graph: a sage graph
# colors: a list of lists which partitions the vertices in G into colors
#   colors[0]: True vertices
#   colors[1]: False vertices
#   colors[3] is reserved for the algorithm
def genrep(graph, colors):
    ret = [(graph, colors)]
    if len(colors[0]) >= (len(graph.vertices()) / 2):
        return ret
    # print("----\nCall with colors: %s" % colors)
    # print("graph: %s" % graph.edges())
    A, orbits = graph.automorphism_group(partition=colors, orbits=True, algorithm="bliss")
    # orbits = A.orbits()
    # print("A: %s, orbits: %s" % (A, orbits))
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
        # print("canon color: %s" % canoncolor)
        Dcanon = graph.canonical_label(partition=canoncolor, algorithm="bliss")
        # print("")
        # print("   Dcolor: %s" % Dcolor)
        # print("   Canonical deletion: %s" % canoncolor)
        # print("   Considering coloring %s" % (Gpcolor))
        if Dcanon == D:
            # print("   accepted")
            ret += genrep(graph, Gpcolor)
    return ret

# given a graph g, count the number of distinct 2-colorings
def count_num_distinct(g):
    return g.automorphism_group().cycle_index().expand(2)((1,1))

def main():
    # n=3
    # pr = cProfile.Profile()
    # pr.enable()

    G = graphs.CompleteGraph(10)
    colors = [[], [i for i in G.vertices()]]
    num_vert = len(G.vertices())
    r = genrep(G, colors)
    tot = 0
    for k in r:
        if len(k[1][0]) == num_vert/2.0:
            tot += 1
        else:
            tot+= 2
        print(k)
    print("Found configurations: %d" % tot)
    print("Expected configurations: %d" % count_num_distinct(G))
    # print("size: %d" % len(r))

    # pr.disable()
    # s = StringIO.StringIO()
    # sortby = 'cumulative'
    # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    # ps.print_stats()
    # print s.getvalue()


if __name__ == "__main__":
    main()
