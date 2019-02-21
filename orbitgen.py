from sage.all import *
import cProfile, pstats, StringIO

### foldrep
def foldrep(graph, colors, acc, f):
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
### graph: a sage graph
### colors: a list of lists which partitions the vertices in G into colors
###   colors[0]: True vertices
###   colors[1]: False vertices
###   colors[3] is reserved for the algorithm
### fold: apply a fold method to each orbit
def genrep(G):
    colors = [[], [i for i in G.vertices()]]
    return foldrep(G, colors, [], lambda acc, g, color: acc + [(g, color)])

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

    G = graphs.CycleGraph(4)
    num_vert = len(G.vertices())
    r = genrep(G)
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
