from sage.all import *
import itertools

def findsubsets(S,m):
    return set(itertools.combinations(S, m))

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

# generates a complete graph with n vertices with some extra nodes on each
# vertex
def gen_complete_extra(n):
    g = Graph()
    v = []
    e = []
    # generate complete graph
    for x in range(0,n):
        v += [x]
    for x in findsubsets(range(0, n), 2):
        e += [(x[0], x[1])]

    # now generate an extra vertex for each point and add it
    for x in range(0,n):
        v += [n + x]
        e += [(x, n + x)]

    g.add_vertices(v)
    g.add_edges(e)
    return g
