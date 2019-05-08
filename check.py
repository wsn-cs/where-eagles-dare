import numpy as np

def list_to_bin(liste):
    k = len(liste)
    n = max(max(liste[i]) for i in range(k))
    output = np.zeros((k,n))
    for i in range(k):
        for j in liste[i]:
            output[i][j] = 1
    return output

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
                return False
    return True