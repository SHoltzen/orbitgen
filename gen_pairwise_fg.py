import sys
import itertools

def findsubsets(S,m):
    return set(itertools.combinations(S, m))


n = 24
print("MARKOV")
print(n)
for i in range(0,n):
    sys.stdout.write("2 " )
print("")
# print factor description
subsets = findsubsets(range(0,n), 2)
print(len(subsets))
for (a,b) in subsets:
    print("2 %s %s " % (a,b))
print("")
# print tables
for i in range(0, len(subsets)):
    print("4")
    print("  0.2 0.3")
    print("  0.3 0.2")
