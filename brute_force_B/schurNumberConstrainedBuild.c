//
//  schurNumberConstrainedBuild.c
//  schurNumberRecursive
//
//  Created by rubis on 08/01/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#include "schurNumberConstrainedBuild.h"

unsigned long schurNumberConstrainedBuild(schur_number_partition_t *partitionstruc, mp_limb_t **constraint_partition, mp_size_t constraint_size, struct schurNumberIOAction *action) {
    /*Cette fonction trouve tous les ensembles B maximaux tels que B - B ne rencontre pas constraint_set et les écrit dans le fichier donné par fd. Elle renvoie le nombre de B trouvés.*/
    
    /*Définition des variables*/
    long number_set = 0;
    
    /*Initialisation de la partition*/
    mp_limb_t **partition = partitionstruc->partition;
    mp_limb_t **partitioninvert = partitionstruc->partitioninvert;
    mp_size_t limballoc = partitionstruc->limballoc;    // Nombre de limbes alloué à chaque ensemble de sfpartition
    mp_size_t limbsize = partitionstruc->limbsize;      // Nombre de limbes utilisés par les ensembles de sfpartition
    unsigned long nsize = mp_bits_per_limb * limbsize;      // Plus grand entier pouvant être contenu dans limbsize limbes
    unsigned long nalloc = mp_bits_per_limb * limballoc;    // Plus grand entier pouvant être contenu dans limballoc limbes
    
    mp_limb_t *work0;
    mp_limb_t *work1 = calloc(sizeof(mp_limb_t), limballoc);
    mp_limb_t *work2 = calloc(sizeof(mp_limb_t), constraint_size);
    
    // Initialisation des variables
    unsigned long n0 = partitionstruc->n;   // Taille de la partition initiale
    unsigned long n = n0 + 1;               // Taille de l'intervalle+1
    unsigned long nbest = n0;               // Taille de la plus grande partition trouvée
    unsigned long i = 0;                    // Huche où placer l'entier suivant
    unsigned long p = partitionstruc->p;    // Nombre de huches non vides
    
    unsigned char notSumFree = 1;
    unsigned char is_new_branch = 1;
    
    while (n > n0) {
        while (notSumFree && i < p) {
            /*Regarder si n peut être adjoint à l'ensemble i*/
            
            /*Calcul de set1 = (n+1) - partition[i] = partitioninvert[i] - (nsize - n).*/
            unsigned long nrem = nsize - n;
            work0 = partitioninvert[i] + (limballoc - limbsize);
            
            while (nrem > 0) {
                unsigned int shift = nrem % mp_bits_per_limb;
                if (!shift) {
                    shift = mp_bits_per_limb - 1;
                }
                
                if (nrem == nsize - n) {
                    mpn_rshift(work1, work0, limbsize, shift);     // work1 = work0 - (N - n)
                } else {
                    mpn_rshift(work1, work1, limbsize, shift);     // work1 -= shift
                }
                
                nrem -= shift;
            }
            
            /*Intersection*/
            mpn_and_n(work2, work1, constraint_partition[i], constraint_size);
            notSumFree = !mpn_zero_p(work2, constraint_size);
            i += notSumFree;
        }
        
        if (notSumFree) {
            /*n ne peut être placé nul part*/
            if (is_new_branch) {
                /*Ecrire la partition*/
                if (n > nbest) {
                    schurNumberPrintPartition(p, n, partition);
                    nbest = n;
                }
            }
            /*Retirer n-1 de la partition*/
            n--;
            i = 0;
            while (i < p && !GET_POINT(partition[i], n)) {
                i++;
            }
            DELETE_POINT(partition[i], n);
            DELETE_POINT(partitioninvert[i], nalloc - n);
            if (nsize >= n + mp_bits_per_limb) {
                limbsize--;
                nsize -= mp_bits_per_limb;
            }
            i++;
            is_new_branch = 0;
        } else {
            /*n est placé dans l'ensemble i*/
            ADD_POINT(partition[i], n);
            ADD_POINT(partitioninvert[i], nalloc - n);
            n++;
            if (nsize < n) {
                limbsize++;
                nsize += mp_bits_per_limb;
            }
            i = 0;
            notSumFree = 1;
            is_new_branch = 1;
        }
    }
    
    // Nettoyage
    free(work1);
    free(work2);
    
    return number_set;
}
