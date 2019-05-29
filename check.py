import numpy as np

### Pour 10

P1=[[4,5,15,16,22,28,29,39,40,41,42,48,49,59],
   [2,3,8,14,19,20,24,25,36,46,47,51,62,73],
   [7,9,11,12,13,17,27,31,32,33,35,37,53,56,57,61,79],
   [1,6,10,18,21,23,26,30,34,38,43,45,50,54,65,74],
   [44,52,55,58,60,63,64,66,67,68,69,70,71,72,75,76,77,78,80]]
   
Q=[[4,5,15,16,22,28,29,39,40,41,42,48,49,59],
   [2,3,8,14,19,20,24,25,36,46,47,51,62,73],
   [7,9,11,12,13,17,27,31,32,33,35,37,53,56,57,61,79],
   [1,6,10,18,21,23,26,30,34,38,43,45,50,54,65,74],
   [44,52,55,58,60,63,64,66,67,68,69,70,71,72,75,76,77,78,80]]

for E in [P1,Q]:
    for part in E:
        n = len(part)
        for i in range(n):
            part.append(161-part[i])
            
"""P1=[[1,4],
   [2,3]]

Q=[[1,4],
   [2,3]]"""

p=len(P1)
q=len(Q)
m=max([max(l) for l in P1])

n=max([max(l) for l in Q])

N=2*m+1
R=[[] for i in range(p+q)]

for b in range(n+1):
    for c in range(1,m+1):
        for i in range(p):
            if c in P1[i] : R[i].append(N*b+c)

for b in range(1,n+1):
    for c in range(0,m+1):
        for i in range(q):
            if b in Q[i]: R[i+p].append(N * b - c)

for P in R:
    P.sort()

def list_to_bin(liste):
    k = len(liste)
    n = max(max(liste[i]) for i in range(k))
    output = np.zeros((k,n+1),int)
    for i in range(k):
        for j in liste[i]:
            output[i][j] = 1
    return [L[1:] for L in output]

def bin_to_list(part):
    k = len(part)
    n = len(part[0])
    output = []
    for i in range(k):
        output.append([])
        for j in range(n):
            if part[i][j]:
                output[i].append(j)
    return output

def can_you_add_new_nb(part,weak=False):
    n = len(part)
    if n <= 1:
        return True
    if weak:
        return not (True in np.logical_and(part[:(n//2)],part[n-1:(n-1)//2:-1]))
    else:
        return not (True in np.logical_and(part[:((n+1)//2)],part[n-1:n//2 - 1:-1]))

def is_sum_free(part,weak=False):
    n = len(part[0])
    for P in part:
        for i in range(n):
            if P[i] and (not can_you_add_new_nb(P[:i],weak)):
                print(P[:i])
                return False
    return True

print(R)
print(max(max(P) for P in R))
print(is_sum_free(list_to_bin(R)))