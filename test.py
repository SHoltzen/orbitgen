### defines a test suite for testing basic functionality

from sage.all import *
import orbitgen
from my_graphs import *
import unittest
from factor import *

# generate a complete pairwise factor graph with half the factors different colors
def gen_complete_pairwise_factorgraph_half(n):
    (g, (v, factors)) = gen_complete_pairwise_factor(n)
    def potential(state):
        p = 0.0
        for v in state.itervalues():
            if v:
                p += 1
        return p
    # make half the factors different colors
    return FactorGraph(g, v, [factors[:len(factors)/2],factors[len(factors)/2:]], potential)

def gen_complete_pairwise_factorgraph(n):
    (g, (v, factors)) = gen_complete_pairwise_factor(n)
    def potential(state):
        p = 0.0
        for v in state.itervalues():
            if v:
                p += 1
        return p
    return FactorGraph(g, v, [factors], potential)


def count_num_generated(G):
    colors = [[], [i for i in G.vertices()]]
    num_vert = len(G.vertices())
    l = orbitgen.genrep_bfs(G)
    return len(l)

class TestGen(unittest.TestCase):
    def test_complete(self):
        G = graphs.CompleteGraph(10)
        self.assertEqual(orbitgen.count_num_distinct(G),
                         count_num_generated(G))

    def test_complete2(self):
        G = graphs.CompleteGraph(11)
        self.assertEqual(orbitgen.count_num_distinct(G),
                         count_num_generated(G))

    def test_cyclic(self):
        G = graphs.CycleGraph(11)
        self.assertEqual(orbitgen.count_num_distinct(G),
                         count_num_generated(G))

    def test_cyclic2(self):
        G = graphs.CycleGraph(10)
        self.assertEqual(orbitgen.count_num_distinct(G),
                         count_num_generated(G))

    def test_petersen(self):
        G = graphs.PetersenGraph()
        self.assertEqual(orbitgen.count_num_distinct(G),
                         count_num_generated(G))

    def test_circlelader(self):
        G = graphs.CircularLadderGraph(4)
        self.assertEqual(orbitgen.count_num_distinct(G),
                         count_num_generated(G))

    def test_star(self):
        G = graphs.StarGraph(7)
        self.assertEqual(orbitgen.count_num_distinct(G),
                         count_num_generated(G))

    def test_aug_complete(self):
        G = gen_complete_extra(8)
        self.assertEqual(orbitgen.count_num_distinct(G),
                         count_num_generated(G))

    def test_pairwise(self):
        fg = gen_complete_pairwise_factorgraph(10)
        self.assertEqual(g.partition(), fg.brute_force_partition)

    def test_half_pairwise(self):
        fg = gen_complete_pairwise_factorgraph_half(10)
        self.assertEqual(g.partition(), fg.brute_force_partition)




if __name__ == '__main__':
    unittest.main()
