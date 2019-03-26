//
//  schurNumberSimpleNestedMonteCarlo.c
//  SchurNumber
//
//  Created by Gabriel Merlin on 22/03/2019.
//

#include <stdlib.h>
#include <gmp.h>
#include <stdio.h>

#if GMP_NUMB_BITS == 64
#define GMP_2EXP_NUMB_BITS 6
#elif GMP_NUMB_BITS == 128
#define GMP_2EXP_NUMB_BITS 7
#elif GMP_NUMB_BITS == 32
#define GMP_2EXP_NUMB_BITS 5
#endif

#define REALLOC(ptr, reptr, size, resize) \
do { \
reptr = realloc(ptr, resize);\
while (reptr == NULL && resize > size) { \
resize -= resize >> 1;\
reptr = realloc(ptr, resize);\
} \
if (resize != size) { \
 ptr = reptr;\
} \
} while(0)

struct schur_partition_struct {
    unsigned int pmax;  // Nombre maximal de huches
    unsigned int p;     // Nombre de huches non vides
    
    unsigned long n;    // Entier courant
    unsigned long nbest;// Meilleur entier dans la récursion
    mp_size_t limballoc;// Nombre de limbes alloués à chaque huche
    mp_size_t limbsize; // Nombre de limbes utilisés par chaque huche
    
    mp_limb_t **partition;
    mp_limb_t **partitioninvert;
    mp_limb_t **partitionbest;
    
    mp_limb_t *work0;
    mp_limb_t *work1;
};

typedef struct schur_partition_struct partition_t;

void partition_realloc(partition_t *partitionstruc) {
    unsigned int i, p;
    mp_size_t limballoc, limbrealloc;
    mp_limb_t *allocptr, *reallocptr;
    mp_limb_t **partition, **partitioninvert, **partitionbest;
    
    limballoc = partitionstruc->limballoc;
    limbrealloc = 2*limballoc;
    p = partitionstruc->pmax;
    
    allocptr = partitionstruc->work0;
    REALLOC(allocptr, reallocptr, limballoc, limbrealloc);
    allocptr = partitionstruc->work1;
    REALLOC(allocptr, reallocptr, limballoc, limbrealloc);
    
    partition = partitionstruc->partition;
    partitioninvert = partitionstruc->partitioninvert;
    partitionbest = partitionstruc->partitionbest;
    for (i=0; i<p; i++) {
        allocptr = *partition;
        REALLOC(allocptr, reallocptr, limballoc, limbrealloc);
        allocptr = *partitioninvert;
        REALLOC(allocptr, reallocptr, limballoc, limbrealloc);
        allocptr = *partitionbest;
        REALLOC(allocptr, reallocptr, limballoc, limbrealloc);
        partition++;
        partitioninvert++;
        partitionbest++;
    }
    for (i=0; i<p; i++) {
        partitioninvert--;
        mpn_lshift(*partitioninvert, *partitioninvert, limbrealloc, limbrealloc - limballoc);
    }
}

unsigned long schurNumberSimpleMonteCarloLevelIteration(partition_t *sfpartitionstruc, unsigned int level) {
    /*Effectue une simulation de niveau l sur sfpartition et renvoie le score du meilleur résultat.*/
    unsigned long n, nmax, nsimulated, nbest;
    unsigned int i, ibest;
    unsigned int p, prand, pmax;
    char isSumFree;
    mp_limb_t *work0;
    mp_limb_t *work1;
    mp_limb_t *set;
    mp_limb_t **sfpartition;
    mp_limb_t **sfpartitioninvert;
    mp_limb_t **sfpartitionbest;
    mp_size_t limballoc;    // Nombre de limbes alloué à chaque ensemble de sfpartition et à work0 et work1
    mp_size_t limbsize;     // Nombre de limbes utilisés par les ensembles de sfpartition
    mp_size_t wlimbsize;    // Nombre de limbes utilisés par work0
    unsigned int nmodbpl;
    unsigned int shift;
    mp_limb_t mask1, mask2;
    
    /*Mise en place de la partition*/
    n = sfpartitionstruc->n;
    nmax = mp_bits_per_limb * sfpartitionstruc->limballoc;
    pmax = sfpartitionstruc->pmax;
    p = sfpartitionstruc->p;
    limballoc = sfpartitionstruc->limballoc;
    limbsize = sfpartitionstruc->limbsize;
    sfpartition = sfpartitionstruc->partition;     // Tableau contenant les ensembles de la partition
    sfpartitioninvert = sfpartitionstruc->partitioninvert; // Tableau contenant les inverses des ensembles de la partition
    sfpartitionbest = sfpartitionstruc->partitionbest;// Tableau contenant la plus grande partition
    work0 = sfpartitionstruc->work0;
    work1 = sfpartitionstruc->work1;
    
    /*Initialisation des variables auxiliaires*/
    nmodbpl = n%mp_bits_per_limb;
    shift = mp_bits_per_limb - nmodbpl;
    wlimbsize = ((n+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
    
    if (!level) {
        /*Simulation de niveau 0*/
        isSumFree = 1;
        if (p == pmax) {
            prand = p;  // Le tirage de la partition est fait parmi 0,…,prand-1
        } else {
            prand = p + 1;
        }
        while (isSumFree) {
            /* Sélectionner une huche aléatoirement */
            i = arc4random_uniform(prand);
            
            if (i < p) {
                /*La huche sélectionnée n'est pas déjà vide. Il faut tester si elle reste sans-somme*/
                set = sfpartitioninvert[i];
                mpn_copyd(work0, &set[limballoc - wlimbsize], wlimbsize);
                mpn_rshift(work1, work0, wlimbsize, shift);
                set = sfpartition[i];
                mpn_and_n(work0, set, work1, wlimbsize);
                isSumFree = mpn_zero_p(work0, wlimbsize);
            } else {
                /*La huche sélectionnée est vide*/
                isSumFree = 1;
                p++;
                if (p < pmax) {
                    prand++;
                }
            }
            
            if (isSumFree) {
                /* Ajouter n+1 à la huche i */
                mask1 = (mp_limb_t)1<<nmodbpl;
                sfpartition[i][limbsize -1] |= mask1;
                mask2 = (mp_limb_t)1<<(shift - 1);
                sfpartitioninvert[i][limballoc - limbsize] |= mask2;
                /* Incrémenter n */
                n++;
                nmodbpl = n%mp_bits_per_limb;
                shift = mp_bits_per_limb - nmodbpl;
                wlimbsize = ((n+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
                
                if (n >= nmax) {
                    partition_realloc(sfpartitionstruc);
                    nmax = mp_bits_per_limb * sfpartitionstruc->limballoc;
                }
                
                if (!nmodbpl) {
                    /* mp_bits_per_limb divise n, donc il faut commencer à remplir un nouveau limbe. */
                    for (i=0; i<pmax; i++) {
                        sfpartition[i][limbsize] = (mp_limb_t)0;
                        sfpartitioninvert[i][limballoc - limbsize -1] = (mp_limb_t)0;
                    }
                    limbsize++;
                }
            }
        }
        
        if (n > (sfpartitionstruc->nbest)) {
            for (i=0; i<p; i++) {
                mpn_copyi(sfpartitionbest[i], sfpartition[i], limbsize);
            }
            sfpartitionstruc->nbest = n;
        }
        nbest = n;
        
        /*Revenir à la partition initiale*/
        wlimbsize = sfpartitionstruc->limbsize;
        shift = mp_bits_per_limb - (sfpartitionstruc->n)%mp_bits_per_limb;
        mask1 = (mp_limb_t)GMP_NUMB_MAX>>shift;
        mask2 = (mp_limb_t)GMP_NUMB_MAX<<shift;
        for (i=0; i<p; i++) {
            mpn_zero(sfpartition[i] + wlimbsize, limbsize - wlimbsize);
            mpn_zero(sfpartitioninvert[i] + limballoc - limbsize -1, limbsize - wlimbsize);
            sfpartition[i][wlimbsize-1] |= mask1;
            sfpartitioninvert[i][limballoc - wlimbsize] |= mask2;
        }
        
    } else {
    _simulbegin:
        /*Simulation de niveau > 0*/
        ibest = 0;
        nbest = n;
        mask1 = (mp_limb_t)1<<nmodbpl;
        mask2 = (mp_limb_t)1<<(shift - 1);
        
        if (n >= nmax) {
            partition_realloc(sfpartitionstruc);
            nmax = mp_bits_per_limb * sfpartitionstruc->limballoc;
        }
        
        for (i=0; i<p; i++) {
            /*Tester si il est possible d'ajouter n+1 à la huche i*/
            set = sfpartitioninvert[i];
            mpn_copyd(work0, set + limballoc - wlimbsize, wlimbsize);
            mpn_rshift(work1, work0, wlimbsize, shift);
            set = sfpartition[i];
            mpn_and_n(work0, set, work1, wlimbsize);
            isSumFree = mpn_zero_p(work0, wlimbsize);
            if (isSumFree) {
                /*Ajouter n+1 à la huche i*/
                sfpartition[i][limbsize -1] |= mask1;
                sfpartitioninvert[i][limballoc - limbsize] |= mask2;
                sfpartitionstruc->n = n+1;
                /*Lancer la simulation au niveau l-1*/
                nsimulated = schurNumberSimpleMonteCarloLevelIteration(sfpartitionstruc, level-1);
                if (nsimulated > nbest) {
                    ibest = i;
                    nbest = nsimulated;
                }
                /*Retirer n+1 de la huche i*/
                sfpartition[i][limbsize -1] ^= mask1;
                sfpartitioninvert[i][limballoc - limbsize] ^= mask2;
                sfpartitionstruc->n = n;
            }
        }
        
        if (p < pmax) {
            /*La partition contient une huche vide. Ajouter n+1 dans cette huche.*/
            sfpartition[p][limbsize -1] |= mask1;
            sfpartitioninvert[p][limballoc - limbsize] |= mask2;
            sfpartitionstruc->n = n+1;
            sfpartitionstruc->p = p+1;
            /*Lancer la simulation au niveau l-1*/
            nsimulated = schurNumberSimpleMonteCarloLevelIteration(sfpartitionstruc, level-1);
            if (nsimulated > nbest) {
                ibest = p;
                nbest = nsimulated;
            }
            /*Retirer n+1 de la huche i*/
            sfpartition[p][limbsize -1] ^= mask1;
            sfpartitioninvert[p][limballoc - limbsize] ^= mask2;
            sfpartitionstruc->n = n;
            sfpartitionstruc->p = p;
        }
        
        if (nbest > n) {
            /* Ajouter définitivement n+1 à la partition */
            sfpartition[ibest][limbsize -1] |= mask1;
            sfpartitioninvert[ibest][limballoc - limbsize] |= mask2;
            if (ibest == p) {
                p++;
            }
            /* Incrémenter n */
            n++;
            sfpartitionstruc->n = n+1;
            nmodbpl = n%mp_bits_per_limb;
            shift = mp_bits_per_limb - nmodbpl;
            wlimbsize = ((n+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
            goto _simulbegin;
        }
    }
    
    return sfpartitionstruc->nbest;
}

unsigned long schurNumberSimpleNestedMonteCarlo(unsigned int p, unsigned int level) {
    unsigned long nbest;
    unsigned int i;
    mp_limb_t *work0;
    mp_limb_t *work1;
    partition_t sfpartitionstruc;
    mp_limb_t **sfpartition;
    mp_limb_t **sfpartitioninvert;
    mp_limb_t **sfpartitionbest;
    
    /*Initialisation*/
    work0 = calloc(p, sizeof(mp_limb_t));
    work1 = calloc(p, sizeof(mp_limb_t));
    sfpartition = calloc(p, sizeof(mp_limb_t *));  //Tableau contenant la partition
    sfpartitioninvert = calloc(p, sizeof(mp_limb_t *));    //Tableau contenant les ensembles "inverses" de la partition
    sfpartitionbest = calloc(p, sizeof(mp_limb_t *));
    for (i=0; i<p; i++) {
        sfpartition[i] = calloc(p, sizeof(mp_limb_t));
        sfpartitioninvert[i] = calloc(p, sizeof(mp_limb_t));
        sfpartitionbest[i] = calloc(p, sizeof(mp_limb_t));
    }
    sfpartitionstruc.limballoc = p;
    sfpartitionstruc.limbsize = 1;
    sfpartitionstruc.work0 = work0;
    sfpartitionstruc.work1 = work1;
    sfpartitionstruc.partition = sfpartition;
    sfpartitionstruc.partitioninvert = sfpartitioninvert;
    sfpartitionstruc.partitionbest = sfpartitionbest;
    sfpartitionstruc.pmax = p;
    sfpartitionstruc.nbest = 1;
    
    /*Création de la première partition {{1}}*/
    sfpartitionstruc.n = 1;  // Taille de l'intervalle
    sfpartitionstruc.p = 1;  // Nombre de huches non vides
    *sfpartition[0] = (mp_limb_t)1;
    sfpartitioninvert[0][p - 1] = (mp_limb_t)1<<(mp_bits_per_limb-1);
    
    /*Lancement de la simulation de niveau level*/
    nbest = schurNumberSimpleMonteCarloLevelIteration(&sfpartitionstruc, level);
    
    /*Nettoyage*/
    for (i=0; i<p; i++) {
        free(sfpartition[i]);
        free(sfpartitioninvert[i]);
        free(sfpartitionbest[i]);
    }
    free(sfpartition);
    free(sfpartitioninvert);
    free(sfpartitionbest);
    free(work0);
    free(work1);
    
    return nbest;
}
