### utility file for generating useful factor graphs, and graphs used
### in the experiments

from sage.all import *
import itertools

def findsubsets(S,m):
    return set(itertools.combinations(S, m))

### generates a complete pairwise factor graph
### returns a the graph and a tuple whose first element are the variables
### and second element are the factors
def gen_complete_pairwise_factor(n):
    g = Graph(sparse=True)
    v = [x for x in range(0,n)]
    factors = []
    edges = []
    count = n
    for (s1,s2) in findsubsets(v, 2):
        factors += [count]
        edges += [(s1, count), (s2, count)]
        count += 1

    g.add_vertices(v)
    g.add_vertices(factors)
    g.add_edges(edges)
    return (g, (v, factors))


### generates a factor graph with one big factor and n variables connected to it
### returns a the graph and a tuple whose first element are the variables
### and second element are the factors
def gen_single_big_factor(n):
    g = Graph(sparse=True)
    # make n smoker vertices
    v = [x for x in range(0,n)]
    # connect all the smoker to the factors
    # make friends
    factors = [n]
    edges = []
    # friends = []
    # friendedges = []
    count = n
    for va in v:
        edges += [(n, va)]

    g.add_vertices(v)
    g.add_vertices(factors)
    g.add_edges(edges)
    return (g, (v, factors))

### generates a complete pairwise factor graph
### returns a the graph and a tuple whose first element are the variables
### and second element are the factors
def gen_friends_smokers_factor(n):
    g = Graph(sparse=True)
    # make n smoker vertices
    smokers = [x for x in range(0,n)]
    # connect all the smoker to the factors
    # make friends
    factors = []
    edges = []
    count = n
    friends = []
    for (s1,s2) in findsubsets(smokers, 2):
        cur_factor = count
        factors += [cur_factor]
        count += 1
        edges += [(s1, count), (s2, count)]
        friends += [count]
        edges += [(cur_factor, count)]
        count += 1

    g.add_vertices(smokers)
    g.add_vertices(friends)
    g.add_vertices(factors)
    g.add_edges(edges)
    return (g, (smokers + friends, factors))


def gen_pigeonhole(n,m):
    g = Graph()
    v = []
    e = []
    # generate vertices
    for x in range(0,n):
        for y in range(0,m):
            v += [(x,y)]

    # generate edges
    # generate fully connected graph in n
    for x in range(0,n):
        for y in findsubsets(range(0, m), 2):
            t = tuple([(x, y[0]), (x, y[1])])
            e.append(t)

    # connect between n and m
    for x in range(0,n):
        for y in range (0,m):
            t = tuple([(x, y), (((x + 1) % n), y)])
            e.append(t)
    g.add_vertices(v)
    g.add_edges(e)
    return g


def gen_pigeonhole_fg(n,m):
    # n = holes, m = pigeons
    g = Graph()
    v = []
    e = []
    pigeon_factors = []
    hole_factors = []
    # generate vertices
    for x in range(0,n):
        for y in range(0,m):
            v += [(x,y)]

    count = 0
    # generate edges
    # generate fully connected graph in n
    for x in range(0,n):
        for y in findsubsets(range(0, m), 2):
            count += 1
            pigeon_factors.append(count)
            e.append(tuple([(x, y[0]), count]))
            e.append(tuple([count, (x, y[1])]))
    # connect between n and m
    for x in range(0,n-1):
        for y in range (0,m):
            count += 1
            hole_factors.append(count)
            e.append(tuple([(x, y), count]))
            e.append(tuple([count, (((x + 1) % n), y)]))
    g.add_vertices(v)
    g.add_vertices(pigeon_factors)
    g.add_vertices(hole_factors)
    g.add_edges(e)
    return (g, (v, [hole_factors, pigeon_factors]))




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

