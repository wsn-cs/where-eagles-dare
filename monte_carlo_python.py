import random 

def possible_inserer(element, liste):
    for x in liste: 
        for y in liste : 
            if x+y == element : 
                return False
    return True

def copie_part(partition) : 
    p = []
    for i in range(len(partition)) : 
        p.append(list(partition[i]))
    return p

def level0(debut, partition, nb_boites):
    p = copie_part(partition)
    k = random.randint(0,nb_boites-1)
    n = debut
    sum_free = possible_inserer(debut, p[k])
    while sum_free : 
        p[k].append(n)
        n+=1
        k = random.randint(0,nb_boites-1)
        sum_free = possible_inserer(n, p[k])
            
    return n, p

def repetition_level0(debut, partition, nb_boites):
    n_max = 0
    partition_max = None
    for _ in range(50) : 
        n, candidat = level0(debut, partition, nb_boites)
        if n > n_max : 
            n_max = n 
            partition_max = candidat
    return n_max, partition_max

def nested(level, debut, partition, nb_boites):
    if level == 0 : 
        return repetition_level0(debut, partition, nb_boites)
    else : 
        score_max = 0 
        boite_max = None
        for boite in range(0, nb_boites-1) : 

            if possible_inserer(debut, partition[boite]) : 
                p = copie_part(partition)
                p[boite].append(debut)
                score, _ = nested(level-1, debut+1, p, nb_boites)
                if score > score_max : 
                    score_max, boite_max,  = score, boite
        if boite_max == None : 
            return debut, partition
        
        partition[boite_max].append(debut)
        return nested(level, debut+1, partition, nb_boites)
        

