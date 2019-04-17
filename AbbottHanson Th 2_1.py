P=[[1,4],[2,3]]
Q=[[1,4],[2,3]]

p=len(P)
q=len(Q)
m=max([max(l) for l in P])

n=max([max(l) for l in Q])

N=2*m+1
R=[[] for i in range(p+q)]


Sn=4

for b in range(n+1):
    for c in range(1,m+1):
        for i in range(p):
            if c in P[i] : R[i].append(N*b+c)

for b in range(1,n+1):
    for c in range(0,m+1):
        for i in range(q):
            if b in Q[i]: R[i+p].append(N * b - c)

print(R)
print(max(max(l) for l in R))

