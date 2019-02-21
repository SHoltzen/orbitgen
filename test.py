from sage.all import *
import orbitgen
import unittest

def count_num_generated(G):
    colors = [[], [i for i in G.vertices()]]
    num_vert = len(G.vertices())
    l = orbitgen.genrep(G)
    tot = 0
    for k in l:
        if len(k[1][0]) == num_vert/2.0:
            tot += 1
        else:
            tot+= 2
    return tot

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

if __name__ == '__main__':
    unittest.main()
