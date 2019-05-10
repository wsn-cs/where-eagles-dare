//
//  schurNumberSimpleMonteCarlo.c
//  SchurNumber
//
//  Created by Gabriel Merlin on 22/03/2019.
//

#include "schurNumberNestedMonteCarloHeader.h"
#include "random_uniform.h"

unsigned long schurNumberSimpleMonteCarloLevelIteration(partition_t *sfpartitionstruc, unsigned int level, mp_limb_t **sfpartitionbest,
                                                        unsigned int *pbestptr, unsigned int simulnum0) {
    /*Effectue une simulation de niveau l sur sfpartition et renvoie le score du meilleur résultat.
     La partition correspondant est conservée dans sfpartitionbest.*/
    unsigned long n, nmax, nsimulated, nbest;
    unsigned int i, j, simul;
    unsigned int p, prand, pmax, pbestrec;
    char isSumFree;
    mp_limb_t *work0;
    mp_limb_t *work1;
    mp_limb_t **sfpartition;
    mp_limb_t **sfpartitioninvert;
    mp_limb_t **sfpartitionbestrec;//Tableau servant au niveau en-dessous
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
    work0 = sfpartitionstruc->work0;
    work1 = sfpartitionstruc->work1;
    
    /*Initialisation des variables auxiliaires*/
    nmodbpl = n%mp_bits_per_limb;
    shift = mp_bits_per_limb - nmodbpl;
    wlimbsize = ((n+1) / (2 * GMP_NUMB_BITS)) + 1;
    
    if (!level) {
        /*Simulation de niveau 0*/
        nbest = n;
        if (p == pmax) {
            prand = p;  // Le tirage de la partition est fait parmi 0,…,prand-1
        } else {
            prand = p + 1;
        }
        isSumFree = 1;
     
        while (isSumFree) {
            
            if (n >= nmax) {
                /*Augmentation de la taille de la partition.*/
                partition_realloc(sfpartitionstruc, sfpartitionbest);
                nmax = mp_bits_per_limb * sfpartitionstruc->limballoc;
                if (n == nmax) {
                    prand = 0;
                }
            }
         
            if (!nmodbpl) {
                /* mp_bits_per_limb divise n, donc il faut commencer à remplir un nouveau limbe. */
                for (i=0; i<pmax; i++) {
                    sfpartition[i][limbsize] = (mp_limb_t)0;
                    sfpartitioninvert[i][limballoc - limbsize -1] = (mp_limb_t)0;
                }
                limbsize++;
            }
            
            /* Sélectionner une huche aléatoirement parmi 0,…,prand-1*/
            i = RANDOM_UNIFORM(prand);
            
            /*Tester si il est possible de mettre n+1 dans la huche i.*/
            mpn_rshift(work1, sfpartitioninvert[i] + limballoc - limbsize, limbsize, shift);
            mpn_and_n(work0, sfpartition[i], work1, wlimbsize);
            isSumFree = mpn_zero_p(work0, wlimbsize);
            
            if (isSumFree) {
                /* Ajouter n+1 à la huche i */
                mask1 = (mp_limb_t)1<<nmodbpl;
                sfpartition[i][limbsize -1] |= mask1;
                mask2 = (mp_limb_t)1<<(shift - 1);
                sfpartitioninvert[i][limballoc - limbsize] |= mask2;
                if (i == p) {
                    /*Une nouvelle huche est commencée.*/
                    p++;
                    if (p == pmax) {
                        prand = p;  // Le tirage de la partition est fait parmi 0,…,prand-1
                    } else {
                        prand = p + 1;
                    }
                }
                /* Incrémenter n */
                n++;
                nmodbpl = n%mp_bits_per_limb;
                shift = mp_bits_per_limb - nmodbpl;
                wlimbsize = ((n+1) / (2 * GMP_NUMB_BITS)) + 1;
            }
        }
        
        /*Sauvegarder la partition obtenue dans sfpartitionbest.*/
        for (i=0; i<p; i++) {
            mpn_copyi(sfpartitionbest[i], sfpartition[i], limbsize);
        }
        nbest = n-1;
        *pbestptr = p;
        
        /*Revenir à la partition initiale*/
        wlimbsize = sfpartitionstruc->limbsize;
        shift = mp_bits_per_limb - (sfpartitionstruc->n)%mp_bits_per_limb;
        mask1 = (mp_limb_t)GMP_NUMB_MAX>>shift;
        mask2 = (mp_limb_t)GMP_NUMB_MAX<<shift;
        for (i=0; i<p; i++) {
            mpn_zero(sfpartition[i] + wlimbsize, limbsize - wlimbsize);
            mpn_zero(sfpartitioninvert[i] + limballoc - limbsize, limbsize - wlimbsize);
            sfpartition[i][wlimbsize-1] &= mask1;
            sfpartitioninvert[i][limballoc - wlimbsize] &= mask2;
        }
        
    } else {
        /*Simulation de niveau > 0*/
        /*Création de la meilleure partition utilisée au niveau level-1*/
        sfpartitionbestrec = calloc(pmax, sizeof(mp_limb_t *));
        for (i=0; i<pmax; i++) {
            sfpartitionbestrec[i] = calloc(limballoc, sizeof(mp_limb_t));
        }
        nbest = 0;
    _simulbegin:
        mask1 = (mp_limb_t)1<<nmodbpl;
        mask2 = (mp_limb_t)1<<(shift - 1);
        
        if (n >= nmax) {
            /*Augmentation de la taille de la partition.*/
            partition_realloc(sfpartitionstruc, sfpartitionbest);
            limballoc = sfpartitionstruc->limballoc;
            nmax = mp_bits_per_limb * limballoc;
            if (n == nmax) {
                return nbest;
            }
        }
        
        if (!nmodbpl) {
            /* mp_bits_per_limb divise n, donc il faut commencer à remplir un nouveau limbe. */
            for (i=0; i<pmax; i++) {
                sfpartition[i][limbsize] = (mp_limb_t)0;
                sfpartitioninvert[i][limballoc - limbsize -1] = (mp_limb_t)0;
            }
            limbsize++;
            sfpartitionstruc->limbsize = limbsize;
        }
        
        for (i=0; i<p; i++) {
            /*Tester si il est possible d'ajouter n+1 à la huche i*/
            mpn_rshift(work1, sfpartitioninvert[i] + limballoc - limbsize, limbsize, shift);
            mpn_and_n(work0, sfpartition[i], work1, wlimbsize);
            isSumFree = mpn_zero_p(work0, wlimbsize);
            if (isSumFree) {
                /*Ajouter n+1 à la huche i*/
                sfpartition[i][limbsize -1] |= mask1;
                sfpartitioninvert[i][limballoc - limbsize] |= mask2;
                sfpartitionstruc->n = n+1;
                /*Lancer la simulation au niveau l-1*/
                if (level == 1) {
                    for (simul=0; simul<simulnum0; simul++) {
                        nsimulated = schurNumberSimpleMonteCarloLevelIteration(sfpartitionstruc, 0, sfpartitionbestrec, &pbestrec, simulnum0);
                        if (nsimulated > nbest) {
                            /*La partition obtenue est meilleure que celle précédente.*/
                            nbest = nsimulated;
                            *pbestptr = pbestrec;
                            for (j=0; j<pbestrec; j++) {
                                mpn_copyi(sfpartitionbest[j], sfpartitionbestrec[j], 1 + (nbest / GMP_NUMB_BITS));
                            }
                        }
                    }
                } else {
                    nsimulated = schurNumberSimpleMonteCarloLevelIteration(sfpartitionstruc, level-1, sfpartitionbestrec, &pbestrec, simulnum0);
                    if (nsimulated > nbest) {
                        /*La partition obtenue est meilleure que celle précédente.*/
                        nbest = nsimulated;
                        *pbestptr = pbestrec;
                        for (j=0; j<pbestrec; j++) {
                            mpn_copyi(sfpartitionbest[j], sfpartitionbestrec[j], 1 + (nbest / GMP_NUMB_BITS));
                        }
                    }
                }
                /*Revenir à la partition initiale.*/
                for (j=0; j < pmax; j++) {
                    mpn_zero(sfpartition[j] + limbsize, limballoc - limbsize);
                    mpn_zero(sfpartitioninvert[j], limballoc - limbsize - 1);
                    sfpartition[j][limbsize-1] &= GMP_NUMB_MAX>>shift;
                    sfpartitioninvert[j][limballoc - limbsize] &= GMP_NUMB_MAX<<shift;
                }
                sfpartitionstruc->n = n;
                sfpartitionstruc->p = p;
                sfpartitionstruc->limbsize = limbsize;
            }
        }
        
        if (p < pmax) {
            /*La partition contient une huche vide. Ajouter n+1 dans cette huche.*/
            sfpartition[p][limbsize -1] |= mask1;
            sfpartitioninvert[p][limballoc - limbsize] |= mask2;
            sfpartitionstruc->n = n+1;
            sfpartitionstruc->p = p+1;
            /*Lancer la simulation au niveau l-1*/
            if (level == 1) {
                for (simul=0; simul<simulnum0; simul++) {
                    nsimulated = schurNumberSimpleMonteCarloLevelIteration(sfpartitionstruc, 0, sfpartitionbestrec, &pbestrec, simulnum0);
                    if (nsimulated > nbest) {
                        /*La partition obtenue est meilleure que celle précédente.*/
                        nbest = nsimulated;
                        *pbestptr = pbestrec;
                        for (j=0; j<pbestrec; j++) {
                            mpn_copyi(sfpartitionbest[j], sfpartitionbestrec[j], 1 + (nbest / GMP_NUMB_BITS));
                        }
                    }
                }
            } else {
                nsimulated = schurNumberSimpleMonteCarloLevelIteration(sfpartitionstruc, level-1, sfpartitionbestrec, &pbestrec, simulnum0);
                if (nsimulated > nbest) {
                    /*La partition obtenue est meilleure que celle précédente.*/
                    nbest = nsimulated;
                    *pbestptr = pbestrec;
                    for (j=0; j<pbestrec; j++) {
                        mpn_copyi(sfpartitionbest[j], sfpartitionbestrec[j], 1 + (nbest / GMP_NUMB_BITS));
                    }
                }
            }
            /*Revenir à la partition initiale.*/
            for (j=0; j < pmax; j++) {
                mpn_zero(sfpartition[j] + limbsize, limballoc - limbsize);
                mpn_zero(sfpartitioninvert[j], limballoc - limbsize -1);
                sfpartition[j][limbsize-1] &= GMP_NUMB_MAX>>shift;
                sfpartitioninvert[j][limballoc - limbsize] &= GMP_NUMB_MAX<<shift;
            }
            sfpartitionstruc->n = n;
            sfpartitionstruc->p = p;
            sfpartitionstruc->limbsize = limbsize;
        }
        
        if (nbest > n) {
            /*Déterminer i correspondant à la meilleure partition obtenue à ce niveau.*/
            i = 0;
            while (!(sfpartitionbest[i][limbsize - 1] & ((mp_limb_t)1<<nmodbpl))) {
                i++;
            }
            /* Ajouter définitivement n+1 à la partition */
            sfpartition[i][limbsize -1] |= mask1;
            sfpartitioninvert[i][limballoc - limbsize] |= mask2;
            if (i == p) {
                /*Une nouvelle huche est commencée.*/
                p++;
                sfpartitionstruc->p = p;
            }
            /* Incrémenter n */
            n++;
            sfpartitionstruc->n = n+1;
            nmodbpl = n%mp_bits_per_limb;
            shift = mp_bits_per_limb - nmodbpl;
            wlimbsize = ((n+1) / (2 * GMP_NUMB_BITS)) + 1;
            goto _simulbegin;
        } else {
            for (i=0; i<pmax; i++) {
                free(sfpartitionbestrec[i]);
            }
            free(sfpartitionbestrec);
        }
    }
    
    return nbest;
}

unsigned long schurNumberSimpleNestedMonteCarlo(unsigned int p, unsigned long *narray, unsigned int level, unsigned int simulnum, unsigned int iternum, mp_limb_t **sfpartitionbestglobal, partition_t *sfpartitionbeginstruc) {
    /*Fonction à lancer pour trouver une borne inférieure au nombre de Schur S(p).
     Elle alloue la mémoire nécessaire aux simulations et lance simulnum simulations de niveau level.
     Les simulations de niveau 0 seront effectués simulnum0 fois.*/
    unsigned long nbest, nsimulated;
    unsigned int i, j, pbest;
    mp_size_t limballoc;
    partition_t sfpartitionstruc;
    mp_limb_t **sfpartitionbest;
    
    /*Initialisation*/
    limballoc = (1<<(2*p))>>mp_bits_per_limb;
    partition_init(p, mp_bits_per_limb * limballoc, &sfpartitionstruc);
    sfpartitionbest = calloc(p, sizeof(mp_limb_t *));       //Tableau servant au niveau en-dessous
    for (i=0; i<p; i++) {
        sfpartitionbest[i] = calloc(limballoc, sizeof(mp_limb_t));
    }
    
    nbest = 1;
    for (i=0; i<simulnum; i++) {
        /*Création de la première partition à partir de sfpartitionbegin*/
        partition_copy(&sfpartitionstruc, sfpartitionbeginstruc);
        
        /*Lancement de la simulation de niveau level*/
        nsimulated = schurNumberSimpleMonteCarloLevelIteration(&sfpartitionstruc, level, sfpartitionbest, &pbest, iternum);
        narray[i] = nsimulated + 1;
        if (nsimulated > nbest) {
            nbest = nsimulated;
            for (j=0; j<pbest; j++) {
                mpn_copyd(sfpartitionbestglobal[j], sfpartitionbest[j], 1 + (nbest / GMP_NUMB_BITS));
            }
        }
    }
    
    /*Nettoyage*/
    partition_unalloc(&sfpartitionstruc);
    for (i=0; i<p; i++) {
        free(sfpartitionbest[i]);
    }
    free(sfpartitionbest);
    
    return nbest + 1;
}
