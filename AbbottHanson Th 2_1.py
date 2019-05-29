### Pour 10

P=[[4,5,15,16,22,28,29,39,40,41,42,48,49,59],
   [2,3,8,14,19,20,24,25,36,46,47,51,62,73],
   [7,9,11,12,13,17,27,31,32,33,35,37,53,56,57,61,79],
   [1,6,10,18,21,23,26,30,34,38,43,45,50,54,65,74],
   [44,52,55,58,60,63,64,66,67,68,69,70,71,72,75,76,77,78,80]]
   
Q=[[4,5,15,16,22,28,29,39,40,41,42,48,49,59],
   [2,3,8,14,19,20,24,25,36,46,47,51,62,73],
   [7,9,11,12,13,17,27,31,32,33,35,37,53,56,57,61,79],
   [1,6,10,18,21,23,26,30,34,38,43,45,50,54,65,74],
   [44,52,55,58,60,63,64,66,67,68,69,70,71,72,75,76,77,78,80]]

for E in [P,Q]:
    for part in E:
        n = len(part)
        for i in range(n):
            part.append(161-part[i])
            
"""P=[[1,2,4,8],
   [3,5,6,7]]

Q=[[1,2,4,8],
   [3,5,6,7]]"""

p=len(P)
q=len(Q)
m=max([max(l) for l in P])

n=max([max(l) for l in Q])

N=2*m+1
R=[[] for i in range(p+q)]

for b in range(n+1):
    for c in range(1,m+1):
        for i in range(p):
            if c in P[i] : R[i].append(N*b+c)

for b in range(1,n+1):
    for c in range(0,m+1):
        for i in range(q):
            if b in Q[i]: R[i+p].append(N * b - c)

for P in R:
    P.sort()

def P10():
   global R
   return R

"""print(R)
print(len(R))
print(max(max(l) for l in R))"""