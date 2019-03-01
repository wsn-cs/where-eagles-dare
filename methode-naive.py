
k = 2

def is_wsf(l) : 
    for i in range(len(l)):
        for j in range(i+1,len(l)):
            if l[i] != l[j] != (l[i]+l[j]) and  (l[i]+l[j]) in l : 
                return False
    return True

#tarité = k
def dfs(part,l,r) :

    if r == len(l)  : 
        wsf = True
        for i in range(k) :
            wsf = wsf and is_wsf(part[i])
        print(part, wsf)
        if wsf : 
            print(part)

        return wsf

    else :
        has_wsf_part = False

        for i in range(k) :
            part[i].append(l[r])
            has_wsf_part = has_wsf_part or dfs(part, l, r+1 )
            part[i].remove(l[r])

        return has_wsf_part

def naif():
    has_wsf_part = True
    n = 0
    while has_wsf_part :
        n +=1 
        l = [i+1 for i in range(n)]
        part = [[] for i in range(k)]
        has_wsf_part = dfs(part, l,0)
    print(n-1)

#WS(2) = 8 
#WS(3) = déjà trop long
naif()