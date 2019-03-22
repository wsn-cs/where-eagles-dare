//
//  schurNumberSimpleMonteCarlo.c
//  schurNumberMonteCarlo
//
//  Created by Gabriel Merlin on 21/03/2019.
//

#include <stdlib.h>
#include <gmp.h>

#define GMP_2EXP_NUMB_BITS 6

unsigned long schurNumberSimpleNestedMonteCarlo(unsigned int pmax, unsigned int nestinglevel) {
    /*
      Cette fonction trouve une borne inférieure du nombre de Schur de pmax grâce à la méthode de Monte Carlo imbriqué.
      Les actions se réduisent au choix d'une huche où placer l'entier suivant.
      L'algorithme est écrit récursivement.
      Il y a quelques problèmes au niveau des réallocations.
     */
    unsigned long nglobal, nnest, ntest;
    unsigned long nbest, nbestlocal;
    unsigned long inest, itest, index, indexbestlocal;
    unsigned long pglobal, pnest, prand, ptest; // Nombre de huches non vides de sfpartition aux différents niveaux
    unsigned long l;
    char isSumFree;
    mp_limb_t *work0;
    mp_limb_t *work1;
    mp_limb_t *set;
    mp_limb_t **sfpartition;
    mp_limb_t **sfpartitioninvert;
    mp_limb_t **sfpartitionbest;
    mp_size_t limballoc;    // Nombre de limbes alloué à chaque ensemble de sfpartition et à work0 et work1
    mp_size_t limbsizeglobal, limbsizenest, limbsizetest;     // Nombre de limbes utilisés par les ensembles de sfpartition
    mp_size_t wlimbsize;    // Nombre de limbes utilisés par work0 et work1
    unsigned int nmodbpl;
    unsigned int shift;
    mp_limb_t mask;
    
    /*Initialisation*/
    limballoc = pmax;
    work0 = calloc(limballoc, sizeof(mp_limb_t));
    work1 = calloc(limballoc, sizeof(mp_limb_t));
    sfpartition = calloc(pmax, sizeof(mp_limb_t *));  //Tableau contenant la partition
    sfpartitioninvert = calloc(pmax, sizeof(mp_limb_t *));    //Tableau contenant les ensembles "inverses" de la partition
    sfpartitionbest = calloc(pmax, sizeof(mp_limb_t *));  //Tableau contenant la meilleure partition trouvée
    for (i=0; i<pmax; i++) {
        sfpartition[i] = calloc(limballoc, sizeof(mp_limb_t));
        sfpartitioninvert[i] = calloc(limballoc, sizeof(mp_limb_t));
        sfpartitionbest[i] = calloc(limballoc, sizeof(mp_limb_t));
    }
    
    /*Création de la première partition {{1}}*/
    nglobal = 1;  // Taille de l'intervalle
    nbest = 1;    // Meilleur n trouvé
    i = 0;        // Huche où placer l'entier suivant
    pglobal = 1;  // Nombre de huches non vides
    prand = 2;    // Nombre de huches à sélectionner
    *sfpartition[0] = (mp_limb_t)1;
    sfpartitioninvert[0][limballoc - 1] = (mp_limb_t)1<<(mp_bits_per_limb-1);
    limbsizeglobal = 1;
    wlimbsize = 1;    // Taille en limbe de la seconde moitié
    nmodbpl = 1;
    shift = mp_bits_per_limb - nmodbpl;
    
    while (1) {
        /*Trouver les suites de placement de longueur nestinglevel valides*/
        l = 1;    // Niveau d'imbrication actuelle
        nnest = nglobal;
        pnest = pglobal;
        nbestlocal = nglobal;
        limbsizenest = limbsizeglobal;
        inest = 0;
        index = 0;
        
        while (l > 0 && l <= nestinglevel) {
            if (inest < pnest) {
                // Tester si l'ensemble obtenu en ajoutant nnest = n+l à la huche i est sans-somme
                set = sfpartitioninvert[inest];
                mpn_copyd(work0, &set[limballoc - wlimbsize], wlimbsize);
                mpn_rshift(work1, work0, wlimbsize, shift);
                set = sfpartition[inest];
                mpn_and_n(work0, set, work1, wlimbsize);
                isSumFree = mpn_zero_p(work0, wlimbsize);
                if (isSumFree) {
                    // Ajouter nnest à la huche inest
                    mask = (mp_limb_t)1<<nmodbpl;
                    sfpartition[inest][limbsizenest -1] |= mask;
                    mask = (mp_limb_t)1<<(shift - 1);
                    sfpartitioninvert[inest][limballoc - limbsizenest] |= mask;
                    if (l == 1) {
                        index = inest;
                    }
                    l++;
                    nnest++;
                    nmodbpl = (nnest)%mp_bits_per_limb;
                    shift = mp_bits_per_limb - nmodbpl;
                    wlimbsize = ((nnest+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
                    
                    if (nnest > nbest) {
                        if (nnest >= mp_bits_per_limb*limballoc) {
                            limballoc *= 2;
                            work0 = realloc(work0, limballoc);
                            work1 = realloc(work1, limballoc);
                            for (inest=0; inest<pmax; inest++) {
                                sfpartition[inest] = realloc(sfpartition[inest], limballoc);
                                sfpartitioninvert[inest] = realloc(sfpartitioninvert[inest], limballoc);
                                mpn_lshift(sfpartitioninvert[inest], sfpartitioninvert[inest], limballoc, limballoc>>1);
                                sfpartitionbest[inest] = realloc(sfpartitionbest[inest], limballoc);
                            }
                        }
                    }
                    if (!nmodbpl) {
                        for (i=0; i<pmax; i++) {
                            sfpartition[inest][limbsizenest] = (mp_limb_t)0;
                            sfpartitioninvert[inest][limballoc - limbsizenest -1] = (mp_limb_t)0;
                        }
                        limbsizenest++;
                    }
                    inest = 0;
                } else {
                    inest++;
                }
            } else if (inest == pnest && pnest < pmax) {
                // Cas particulier où il faut remplir une nouvelle huche
                mask = (mp_limb_t)1<<nmodbpl;
                sfpartition[inest][limbsizenest -1] |= mask; // Ajoute nnest à la huche inest
                mask = (mp_limb_t)1<<(shift - 1);
                sfpartitioninvert[inest][limballoc - limbsizenest] |= mask;
                if (l == 1) {
                    index = inest;
                }
                l++;
                nnest++;
                nmodbpl = nnest%mp_bits_per_limb;
                shift = mp_bits_per_limb - nmodbpl;
                wlimbsize = ((nnest+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
                if (nnest > nbest) {
                    if (nnest >= 64*limballoc -1) {
                        limballoc *= 2;
                        work0 = realloc(work0, limballoc);
                        work1 = realloc(work1, limballoc);
                        for (inest=0; inest<pmax; i++) {
                            sfpartition[inest] = realloc(sfpartition[inest], limballoc);
                            sfpartitioninvert[inest] = realloc(sfpartitioninvert[inest], limballoc);
                            mpn_lshift(sfpartitioninvert[inest], sfpartitioninvert[inest], limballoc, limballoc>>1);
                            sfpartitionbest[inest] = realloc(sfpartitionbest[inest], limballoc);
                        }
                    }
                }
                if (!nmodbpl) {
                    for (inest=0; inest<pmax; i++) {
                        sfpartition[inest][limbsizenest] = (mp_limb_t)0;
                        sfpartitioninvert[inest][limballoc - limbsizenest -1] = (mp_limb_t)0;
                    }
                    limbsizenest++;
                }
                inest = 0;
                pnest++;
            } else {
                // Dépiler
                if (!nmodbpl) {
                    limbsizenest--;
                }
                l--;
                nnest--;
                nmodbpl = nnest%mp_bits_per_limb;
                shift = mp_bits_per_limb - nmodbpl;
                wlimbsize = ((nnest+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
                mask = (mp_limb_t)1<<nmodbpl;
                inest = 0;
                while (inest < pnest && !(sfpartition[inest][limbsizenest-1] & mask)) {
                    inest++;
                }
                sfpartition[inest][limbsizenest-1] ^= mask;
                mask = (mp_limb_t)1<<(shift - 1);
                sfpartitioninvert[inest][limballoc - limbsizenest] ^= mask;
                if (inest == pnest-1) {
                    if (mpn_zero_p(sfpartition[inest], limbsizenest)) {
                        pnest--;
                    }
                }
                inest++;
            }
            
            if (l == nestinglevel) {
                /*Réaliser la descente aléatoire de niveau 0*/
                indexbestlocal = 0;
                ntest = nnest;
                ptest = pnest;
                
                // Définition de prand
                if (ptest == pmax) {
                    prand = ptest;
                } else {
                    prand = ptest + 1;
                }
                limbsizetest = limbsizenest;
                isSumFree = 1;
                while (isSumFree) {
                    // Sélectionner une huche aléatoirement
                    itest = arc4random_uniform(prand);
                    if (itest < ptest) {
                        // Quand la huche n'est pas vide, tester si l'ensemble obtenu en ajoutant ntest+1 à la huche i est sans-somme
                        mpz_add_ui(iternum, iternum, 1);
                        set = sfpartitioninvert[itest];
                        mpn_copyd(work0, &set[limballoc - wlimbsize], wlimbsize);
                        mpn_rshift(work1, work0, wlimbsize, shift);
                        set = sfpartition[itest];
                        mpn_and_n(work0, set, work1, wlimbsize);
                        isSumFree = mpn_zero_p(work0, wlimbsize);
                    } else {
                        // Placer ntest+1 dans une huche vide
                        isSumFree = 1;
                        ptest++;
                        if (ptest < pmax) {
                            prand++;
                        }
                    }
                    if (isSumFree) {
                        // Ajouter ntest+1 à la huche i
                        mask = (mp_limb_t)1<<nmodbpl;
                        sfpartition[itest][limbsizetest -1] |= mask;
                        mask = (mp_limb_t)1<<(shift - 1);
                        sfpartitioninvert[itest][limballoc - limbsizetest] |= mask;
                        ntest++;
                        nmodbpl = ntest%mp_bits_per_limb;
                        shift = mp_bits_per_limb - nmodbpl;
                        wlimbsize = ((nnest+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
                        if (ntest > nbest) {
                            if (ntest >= mp_bits_per_limb*limballoc) {
                                limballoc *= 2;
                                work0 = realloc(work0, limballoc);
                                work1 = realloc(work1, limballoc);
                                for (itest=0; itest<pmax; itest++) {
                                    sfpartition[itest] = realloc(sfpartition[itest], limballoc);
                                    sfpartitioninvert[itest] = realloc(sfpartitioninvert[itest], limballoc);
                                    mpn_lshift(sfpartitioninvert[itest], sfpartitioninvert[itest], limballoc, limballoc>>1);
                                    sfpartitionbest[itest] = realloc(sfpartitionbest[itest], limballoc);
                                }
                            }
                        }
                        if (!nmodbpl) {
                            for (itest=0; itest<pmax; itest++) {
                                sfpartition[itest][limbsizetest] = (mp_limb_t)0;
                                sfpartitioninvert[itest][limballoc - limbsizetest -1] = (mp_limb_t)0;
                            }
                            limbsizetest++;
                        }
                    }
                }
                if (ntest > nbestlocal) {
                    indexbestlocal = index;
                    nbestlocal = ntest;
                }
                if (ntest > nbest) {
                    // Copier sfpartition dans sfpartitionbest
                    for (itest=0; itest<ptest; itest++) {
                        mpn_copyi(sfpartitionbest[itest], sfpartition[itest], limbsizetest);
                    }
                    ntest = nbest;
                }
                
                // Vider les huches nouvellement remplies
                for (itest=pnest; itest<ptest; itest++) {
                    mpn_zero_p(sfpartition[itest], limbsizetest);
                    mpn_zero_p(sfpartitioninvert[itest], limbsizetest);
                }
                // Nettoyer les autres et réinitialiser ntest
                ntest = nnest;
                nmodbpl = ntest%mp_bits_per_limb;
                shift = mp_bits_per_limb - nmodbpl;
                wlimbsize = ((nnest+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
                mask = (mp_limb_t)GMP_NUMB_MAX>>(shift);
                for (itest=0; itest<pnest; itest++) {
                    sfpartition[itest][limbsizenest -1] &= mask;
                }
                mask = (mp_limb_t)GMP_NUMB_MAX<<(nmodbpl);
                for (itest=0; itest<pnest; itest++) {
                    sfpartitioninvert[itest][limballoc - limbsizenest] &= mask;
                }
                
                // Dépiler
                if (!(nnest%mp_bits_per_limb)) {
                    limbsizenest--;
                }
                l--;
                nnest--;
                nmodbpl = nnest%mp_bits_per_limb;
                shift = mp_bits_per_limb - nmodbpl;
                wlimbsize = ((nnest+1)>>(GMP_2EXP_NUMB_BITS + 1)) + 1;
                mask = (mp_limb_t)1<<nmodbpl;
                i = 0;
                while (i < pnest && !(sfpartition[i][limbsizenest-1] & mask)) {
                    i++;
                }
                sfpartition[i][limbsizenest-1] ^= mask;
                mask = (mp_limb_t)1<<(shift - 1);
                sfpartitioninvert[i][limballoc - limbsizenest] ^= mask;
                if (i == pnest-1) {
                    if (mpn_zero_p(sfpartition[i], limbsizenest)) {
                        pnest--;
                    }
                }
                i++;
            }
            
        }
        
        if (nbestlocal == nglobal) {
            // Aucune partition plus grande n'a été trouvée
            break;
        }
        
        // Placer n+1 dans la huche indexbestlocal
        nmodbpl = nglobal%mp_bits_per_limb;
        shift = mp_bits_per_limb - nmodbpl;
        mask = (mp_limb_t)1<<nmodbpl;
        sfpartition[indexbestlocal][limbsizeglobal - 1] |= mask;
        mask = (mp_limb_t)1<<(shift - 1);
        sfpartitioninvert[indexbestlocal][limballoc - limbsizeglobal] |= mask;
        nglobal++;
        if (index == pglobal) {
            pglobal++;
        }
        if (!nmodbpl) {
            limbsizeglobal++;
        }
    }
    
    /*Nettoyage*/
    for (i=0; i<pmax; i++) {
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
