from sage.all import *

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
    A = graph.automorphism_group(partition=colors)
    orbits = A.orbits()
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
        Gprime, c = graph.canonical_label(partition=Gpcolor, certificate=True)

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
        Dcanon = graph.canonical_label(partition=canoncolor)
        # print("")
        # print("   Dcolor: %s" % Dcolor)
        # print("   Canonical deletion: %s" % canoncolor)
        # print("   Considering coloring %s" % (Gpcolor))
        if Dcanon == D:
            # print("   accepted")
            ret += genrep(graph, Gpcolor)
    return ret

def main():
    n=4
    G = graphs.CycleGraph(n)
    colors = [[], [i for i in range(0,n)]]
    r = genrep(G, colors)
    for k in r:
        print(k)
    print(len([x for x in r if len(x[0]) >= 3]))
    print("size: %d" % len(r))

if __name__ == "__main__":
    main()
