import numpy as np

def partition_reader(n):
    with open("partitions/partition_" + str(n) + ".txt") as file:
        lines = file.readlines()
        output = []
        for line in lines:
            line_list = line.split()
            output.append(list(map(int,line_list)))
    return output

### Pour 10

k = 12

P1 = partition_reader(6)
   
Q = partition_reader(6)

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
    output = np.zeros((k,n),int)
    for i in range(k):
        for j in liste[i]:
            output[i][j-1] = 1
    return output

def bin_to_list(part):
    k = len(part)
    n = len(part[0])
    output = []
    for i in range(k):
        output.append([])
        for j in range(n):
            if part[i][j]:
                output[i].append(j+1)
    return output

def can_you_add_new_nb(part,weak=False):
    n = len(part)
    if n <= 1:
        return True
    if weak:
        if not (True in np.logical_and(part[:(n//2)],part[n-1:(n-1)//2:-1])):
            return True
        else:
            print(list(np.logical_and(part[:(n//2)],part[n-1:(n-1)//2:-1])).index(True))
            return False
    else:
        if not (True in np.logical_and(part[:((n+1)//2)],part[n-1:n//2 - 1:-1])):
            return True
        else:
            print(list(np.logical_and(part[:((n+1)//2)],part[n-1:n//2 - 1:-1])).index(True))
            return False

def is_sum_free(part,weak=False):
    n = len(part[0])
    for P in part:
        for i in range(n):
            if P[i] and (not can_you_add_new_nb(P[:i],weak)):
                output = []
                for j in range(i+1):
                    if P[j]:
                        output.append(j+1)
                print(output)
                return False
    return True

print(is_sum_free(list_to_bin(R)))

with open("partitions/partition_" + str(k) + ".txt","w") as file:
    for part in R:
        for i in part:
            file.write(str(i)+" ")
        file.write("\n")