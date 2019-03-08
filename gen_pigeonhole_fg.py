import sys
import itertools

def findsubsets(S,m):
    return set(itertools.combinations(S, m))

# num pigeons
n = 20
print("MARKOV")
print(n * 2)
for i in range(0,n * 2):
    sys.stdout.write("2 " )
print("")
# print factor description
subsets = findsubsets(range(0,n), 2)
print(len(subsets) * 2 + n)
for (a,b) in subsets:
    print("2 %s %s " % (a,b))
    print("2 %s %s " % (a + n, b + n))
# print holes
for a in range(0, n):
    print("2 %s %s" % (a, a+n))
print("")
# print tables
for i in range(0, len(subsets) * 2 + n):
    print("4")
    print("  1 1000")
    print("  1000 1")
