//
//  schurNumberIO.h
//  schurNumberRecursive
//
//  Created by rubis on 28/01/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#ifndef schurNumberIO_h
#define schurNumberIO_h

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

struct schurNumberIOAction {
    unsigned long p;    // Nombre d'ensembles par partition
    
    size_t count;       // Nombre de partitions
    unsigned long nmax; // Taille maximale des partitions
    size_t size;        // Nombre maximale de partitions alloué
    
    mp_size_t *limbsizes;    // Tableau du nombre de limbes par ensembles de la partition
    mp_limb_t **partitions;  // Tableau des partitions
    
    void (*func)(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);
};

typedef struct schurNumberIOAction schur_number_action_t;

void schurNumberGetSet(char *str, mp_limb_t *set);
unsigned long schurNumberGetPartition(unsigned long p, char **str, mp_limb_t **partition);

void schurNumberPrintSet(unsigned long n, mp_limb_t *set);
void schurNumberPrintPartition(unsigned long p, unsigned long n, mp_limb_t **partition);

void schurNumberActionAlloc(schur_number_action_t *action, unsigned long p, void (*func)(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action));
void schurNumberActionDealloc(schur_number_action_t *action);

void schurNumberSaveBestPartition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);
void schurNumberSaveAllPartition(mp_limb_t **partition, unsigned long n, struct schurNumberIOAction *action);
size_t schurNumberPrintPartitions(struct schurNumberIOAction *action);

#endif /* schurNumberIO_h */
