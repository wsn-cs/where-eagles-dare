//
//  schurNumberNestedMonteCarloHeader.h
//  SchurNumber
//
//  Created by Gabriel Merlin on 22/04/2019.
//
//  Ce fichier est un header déclarant toutes les fonctions de calcul.
//

#ifndef schurNumberNestedMonteCarloHeader
#define schurNumberNestedMonteCarloHeader

#include <stdlib.h>
#include <gmp.h>

struct schur_partition_struct {
    unsigned int pmax;  // Nombre maximal de huches
    unsigned int p;     // Nombre de huches non vides
    
    unsigned long n;    // Entier courant
    mp_size_t limballoc;// Nombre de limbes alloués à chaque huche
    mp_size_t limbsize; // Nombre de limbes utilisés par chaque huche
    
    mp_limb_t **partition;
    mp_limb_t **partitioninvert;
    
    mp_limb_t *work0;
    mp_limb_t *work1;
};

typedef struct schur_partition_struct partition_t;

/*Fonction ajoutant les entiers un par un*/
unsigned long schurNumberSimpleNestedMonteCarlo(unsigned int p, unsigned long *narray, unsigned int level,
                                                unsigned int simulnum, unsigned int iternum, mp_limb_t **sfpartitionbestglobal);

#endif
