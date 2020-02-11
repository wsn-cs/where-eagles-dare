//
//  schurNumberIO.c
//  schurNumberRecursive
//
//  Created by rubis on 28/01/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberIO.h"

#define ADD_POINT(set, x) set[(x) / mp_bits_per_limb] |= ((unsigned long)1 << ((x) % mp_bits_per_limb))

#define GET_POINT(set, x) (set[(x) / mp_bits_per_limb] & ((unsigned long)1 << ((x) % mp_bits_per_limb)))

unsigned long schurNumberGetSetMaximum(char *str) {
    char *ptr0 = str;
    char *ptr1 = str;
    
    unsigned long nmax = 0;
    while (*ptr1 != '\0') {
        unsigned long n = strtoul(ptr0, &ptr1, 10);
        if (n > nmax) {
            nmax = n;
        }
        ptr0 = ptr1;
    }
    return nmax;
}

void schurNumberGetSet(char *str, mp_limb_t *set) {
    /*Remplit l'ensemble d'entiers set grâce à la chaîne de caractères str et renvoie le plus grand élément trouvé. La variable set doit pointée vers un tableau de taille suffisante pour contenir l'ensemble.*/
    char *ptr0 = str;
    char *ptr1 = str;
    
    while (*ptr1 != '\0') {
        unsigned long k = strtoul(ptr0, &ptr1, 10);
        ADD_POINT(set, k);
        ptr0 = ptr1;
    }
}

unsigned long schurNumberGetPartition(unsigned long p, char **str, mp_limb_t **partition) {
    /*Remplit la variable partition à partir du tableau de p chaînes de caractères *str. La variable partition doit pointée vers un tableau à p entrées.*/
    unsigned long n = 0;
    
    for (unsigned long i = 0; i < p; i++) {
        unsigned long m = schurNumberGetSetMaximum(str[i]);
        if (m > n) {
            n = m;
        }
    }
    
    mp_size_t limbsize = (n>>6) + 1;
    for (unsigned long i = 0; i < p; i++) {
        partition[i] = calloc(sizeof(mp_limb_t), limbsize);
        schurNumberGetSet(str[i], partition[i]);
    }
    
    return n;
}

void schurNumberPrintSet(unsigned long n, mp_limb_t *set) {
    /*Affiche le contenu de l'ensemble d'entiers set inclus dans l'intervalle [1, n].*/
    
    for (unsigned long i = 1; i <= n; i++) {
        if (GET_POINT(set, i)) {
            printf(" %lu", i);
        }
    }
}

void schurNumberPrintPartition(unsigned long p, unsigned long n, mp_limb_t **partition) {
    /*Affiche une partition à p ensembles de [1, n].*/

    printf("Partition:\n");
    for (unsigned long i=0; i<p; i++) {
        printf("\t");
        schurNumberPrintSet(n, partition[i]);
        printf("\n");
    }
}

void schurNumberActionAlloc(schur_number_action_t *action, unsigned long p, void (*func)(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action)) {
    action->p = p;
    action->count = 0;
    action->nmax = 0;
    size_t size = p << (p+1);
    action->size = size;
    action->limbsizes = calloc(sizeof(mp_size_t), size);
    action->partitions = calloc(sizeof(mp_limb_t *), size);
    action->func = func;
}

void schurNumberActionDealloc(schur_number_action_t *action) {
    for (size_t i = 0; i < action->count; i++) {
        mp_limb_t *partition = action->partitions[i];
        free(partition);
    }
    free(action->limbsizes);
}

void schurNumberSaveBestPartition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action) {
    /*Compare n avec nmax. Si n ≥ nmax, le nmax est mis à jour et la partition est ajoutée aux partitions.*/
    unsigned long  p = action->p;
    
    if (n > action->nmax) {
        /*Vider partitions*/
        for (unsigned long k = 0; k < action->count; k++) {
            mp_limb_t *partition_flat = action->partitions[k];
            free(partition_flat);
            action->limbsizes[k] = 0;
            action->partitions[k] = NULL;
        }
        action->nmax = n;
        action->count = 0;
    }
    
    if (n == action->nmax) {
        /*Ajouter la partition.*/
        unsigned long k = action->count;
        mp_size_t limbsize = (n>>6) + 1;
        
        if (k >= action->size) {
            action->size *= 2;
            mp_limb_t **partitions = realloc(action->partitions, action->size);
            if (partition != NULL) {
                action->partitions = partitions;
            }
            action->limbsizes = realloc(action->limbsizes, action->size);
        }

        mp_limb_t *set = calloc(sizeof(mp_limb_t), p * limbsize);
        action->partitions[k] = set;
        action->limbsizes[k] = limbsize;
        
        for (unsigned long j = 0; j < p; j++) {
            mpn_copyd(set, partition[j], limbsize);
            set += limbsize;
        }
        action->count ++;
    }
}

void schurNumberSaveAllPartition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action) {
    /*La partition est ajoutée aux autres.*/
    unsigned long  p = action->p;
    mp_size_t limbsize = ((action->nmax)>>6) + 1;
    
    if (n > action->nmax) {
        action->nmax = n;
    }
    
    unsigned long k = action->count;
    
    if (k >= action->size) {
        action->size *= 2;
        mp_limb_t **partitions = realloc(action->partitions, action->size);
        mp_size_t *limbsizes = realloc(action->limbsizes, action->size);
        if (limbsizes && partitions) {
            action->partitions = partitions;
            action->limbsizes = limbsizes;
        } else {
            printf("Echec.");
        }
    }
    
    mp_limb_t *set = calloc(sizeof(mp_limb_t), p * limbsize);
    action->partitions[k] = set;
    action->limbsizes[k] = limbsize;
    
    for (unsigned long j = 0; j < p; j++) {
        mpn_copyd(set, partition[j], limbsize);
        set += limbsize;
    }
    
    action->count ++;
}

size_t schurNumberPrintPartitions(struct schurNumberIOAction *action) {
    /*Affiche toutes les partitions.*/
    unsigned long p = action->p;
    mp_limb_t **partition = calloc(sizeof(mp_limb_t *), p);
    
    for (unsigned long k = 0; k < action->count; k++) {
        mp_limb_t *set = action->partitions[k];
        mp_size_t limbsize = action->limbsizes[k];
        for (unsigned long j = 0; j < p; j++) {
            partition[j] = set;
            set += limbsize;
        }
        schurNumberPrintPartition(p, limbsize * mp_bits_per_limb - 1, partition);
    }
    
    free(partition);
    
    return action->count;
}
