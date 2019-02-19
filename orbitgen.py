from sage.all import *

# graph: a sage graph
# colors: a list of lists which partitions the vertices in G into colors
#   colors[0]: True vertices
#   colors[1]: False vertices
#   colors[3] is reserved for the algorithm
def genrep(graph, colors):
    ret = [(graph, colors)]
    if graph.vertices() == colors[0].sort():
        return ret
    # print("colors: %s" % colors)
    A = graph.automorphism_group(partition=colors)
    orbits = A.orbits()
    # print("orbits: %s" % orbits)
    for o in orbits:
        # print("current orbit: %s" % o)
        e = o[0] # pick the first element in the orbit arbitrarily
        # check if this orbit is already true
        if e in colors[0]:
            continue
        # we know this is an uncolored vertex
        Dcolor = [colors[0], [c for c in colors[1] if c != e], [e]]
        D = graph.canonical_label(partition=Dcolor)
        Gpcolor = [colors[0] + [e], [c for c in colors[1] if c != e]]
        Gprime, c = graph.canonical_label(partition=Gpcolor, certificate=True)

        # find lexically first vertex whose value is true in the canonical Gprime;
        # default to vertex 0
        first_t = next((v for v in Gprime.vertices() if v not in Gpcolor[1]), Gprime.vertices()[0])
        canoncolor = [[c for c in Gpcolor[0] if c != first_t], Gpcolor[1], [first_t]]
        # print("canon color: %s" % canoncolor)
        Dcanon = graph.canonical_label(partition=canoncolor)

        if Dcanon == D:
            ret += genrep(Gprime, Gpcolor)
    return ret

def main():
    n=4
    G = graphs.CycleGraph(n)
    # print(G.vertices())
    colors = [[], [i for i in range(0,n)]]
    print(genrep(G, colors))

if __name__ == "__main__":
    main()
