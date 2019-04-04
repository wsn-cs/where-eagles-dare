import random 

def possible_inserer(element, liste):
    for x in liste: 
        for y in liste : 
            if x+y == element : 
                return False
    return True

def nested(level, debut, partition, nb_boites): 

    if level == 0 : 
        sum_free = True
        a_ajouter = debut
        while sum_free : 
            k = random.randint(0, nb_boites-1)
            sum_free = possible_inserer(a_ajouter, partition[k])
            if sum_free : 
                partition[k].append(a_ajouter)
                a_ajouter+=1
        return [a_ajouter-1, partition]
    else : 
        config_max = [0,partition]
        config_fille = [0,partition]
        for k in range(nb_boites):
            sum_free = possible_inserer(debut, partition[k])
            if sum_free : 
                partition_k = partition.copy()
                partition_k[k].append(debut)
                config_fille = nested(level-1, debut+1, partition_k, nb_boites)
                if config_fille[0] > config_max[0] : 
                    config_max = config_fille
        
        print(config_max)
        nested(level, debut+1, config_max[1], nb_boites)

print(nested(1,1,[[] for k in range(3)], 3))