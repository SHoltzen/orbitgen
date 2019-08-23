### defines a test suite for testing basic functionality
### invoke `sage test.py`

from sage.all import *
import orbitgen
from my_graphs import *
import unittest
from factor import *

# holes pigeons, m holes
def mk_pigeonhole_fg(n, m, order=True):
    w1 = 10000000
    w2 = 100000
    (g, (variables, factors)) = gen_pigeonhole_fg(n, m)
    def potential(state):
        total = 0.0
        # to see every pigeon in exactly one hole
        for p in range(0, n):
            # check the holes for the pigeons
            in_hole = False
            for h in range(0, m):
                if state[(p, h)]:
                    if in_hole:
                        return 0.000000001
                    else:
                        in_hole = True
            if in_hole:
                total += w1


        # check to see no no hole has 2 pigeons
        for h in range(0, m):
            for (p1, p2) in findsubsets(range(0, n), 2):
                if not state[(p1, h)] or not state[(p2, h)]:
                    total += w2

        return total
    return FactorGraph(g, variables, factors, potential)


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
        self.assertAlmostEqual(fg.partition(), fg.brute_force_partition())

    def test_half_pairwise(self):
        fg = gen_complete_pairwise_factorgraph_half(10)
        self.assertAlmostEqual(fg.partition(), fg.brute_force_partition())

    def test_pigeonhole(self):
        fg = mk_pigeonhole_fg(2,4)
        self.assertAlmostEqual(fg.partition(), fg.brute_force_partition())

    def test_pigeonhole_2(self):
        fg = mk_pigeonhole_fg(3,6)
        self.assertEqual(int(fg.partition()), int(fg.brute_force_partition()))


if __name__ == '__main__':
    unittest.main()
