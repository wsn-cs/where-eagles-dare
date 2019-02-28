//
//  schurNumberIterative.c
//  SchurNumberIterative
//
//  Created by Gabriel Merlin on 25/02/2019.
//

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

void printPartition(unsigned long p, mpz_t *partition) {
    /*Affiche une partition.*/
    unsigned long i;
    unsigned long j;
    printf("Partition:\n");
    for (i=0; i<p; i++) {
        printf("\t");
        for (j=0; j < 64*(partition[i]->_mp_size); j++) {
            if (mpz_tstbit(partition[i], j)) {
                printf(" %lu", j+1);
            }
        }
        printf("\n");
    }
}

unsigned long schurNumberIterative(unsigned long pmax, unsigned long *nbests) {
    /*
     Cette fonction calcule successivement les nombres de Schur S(p) pour p<= pmax.
     Elle remplit le tableau nbests.
     
     La partition est représentée comme un tableau de grands entiers sfpartition.
     Chaque grand entier représente un ensemble grâce à ses bits: si le bit k vaut 1,
     alors k appartient à l'ensemble; sinon il n'y appartient pas.
     
     Pour tester si il est possible d'ajouter n+1 à l'ensemble P sans-somme,
     on effectue un & entre les [(n+1)/2] premiers bits avec les [(n+1)/2] derniers.
     Si le résultat est 0, il est possible d'ajouter n+1.
     */
    unsigned long n;
    unsigned long i;
    unsigned long p;
    unsigned long limbnum;
    unsigned long rem;
    unsigned long limbmax;
    char isSumFree;
    mp_limb_t *work0;
    mp_limb_t *work1;
    mp_limb_t *set;
    mpz_t *sfpartition;
    mpz_t *sfpartitioninvert;
    
    /*Initialisation*/
    limbmax = pmax; //Taille de l'entier work0 et work1
    work0 = calloc(limbmax, sizeof(mp_limb_t));
    work1 = calloc(limbmax, sizeof(mp_limb_t));
    sfpartition = calloc(pmax, sizeof(mpz_t));  //Tableau contenant la partition
    sfpartitioninvert = calloc(pmax, sizeof(mpz_t));    //Tableau contenant les ensembles "inverses" de la partition
    for (i=0; i<pmax; i++) {
        mpz_init2(sfpartition[i], 64*pmax);
        mpz_init2(sfpartitioninvert[i], 64*pmax);
        nbests[i] = i+1;
    }
    
    /*Création de la première partition {{1}}*/
    n = 1;  // Taille de l'intervalle
    i = 0;  // Position où placer l'entier suivant
    p = 1;  // Nombre de huches non vides
    mpz_set_ui(sfpartition[0], 1);
    mpz_setbit(sfpartitioninvert[0], 64*(sfpartitioninvert[0]->_mp_alloc) -1);
    
    /*Itération jusqu'à énumérer toutes les partions sum-free à au plus pmax huches*/
    while (n>0) {
        // Placer n+1 dans une des huches en conservant la sans-sommité
        limbnum = ((n+1)>>7)+1;
        rem = 64 - (n%64);
        printf("%lu %lu %lu %lu\n", n, rem, rem%0x10000000, ((unsigned long)rem>>16));
        if (i < p) {
            // Tester si l'ensemble obtenu en ajoutant n à la huche i est sans-somme
            set = sfpartitioninvert[i]->_mp_d;
            mpn_copyd(work1, &set[sfpartitioninvert[i]->_mp_alloc - limbnum], limbnum);
            mpn_rshift(work0, work1, limbnum, rem%0x10000000);
            mpn_rshift(work1, work0, limbnum, ((unsigned long)rem>>16));
            set = sfpartition[i]->_mp_d;
            printf("%lu %lu\n", *set, *work1);
            mpn_and_n(work0, set, work1, limbnum);
            isSumFree = mpn_zero_p(work0, limbnum);
            if (isSumFree) {
                mpz_setbit(sfpartition[i], n); // Ajoute n+1 à la huche i
                mpz_setbit(sfpartitioninvert[i], 64*(sfpartitioninvert[i]->_mp_alloc) - n - 1);
                n++;
                limbnum = ((n+1)>>7)+1;
                rem = 64 - (n%64);
                if (n > nbests[p-1]) {
                    nbests[p-1] = n;
                    if (n >= 64*limbmax -1) {
                        work0 = realloc(work0, 2*limbmax);
                        work1 = realloc(work1, 2*limbmax);
                    }
                }
                i = 0;
            } else {
                i++;
            }
        } else if (i == p && p < pmax) {
            // Cas particulier où il faut remplir une nouvelle huche
            set = sfpartitioninvert[i]->_mp_d;
            mpn_copyd(work1, &set[sfpartitioninvert[i]->_mp_alloc - limbnum], limbnum);
            mpn_rshift(work0, work1, limbnum, rem%0x10000000);
            mpn_rshift(work1, work0, limbnum, ((unsigned long)rem>>16));
            set = sfpartition[i]->_mp_d;
            printf("%lu %lu\n", *set, *work1);
            mpn_and_n(work0, set, work1, limbnum);
            isSumFree = mpn_zero_p(work0, limbnum);
            if (isSumFree) {
                mpz_setbit(sfpartition[i], n); // Ajoute n+1 à la huche i
                mpz_setbit(sfpartitioninvert[i], 64*(sfpartitioninvert[i]->_mp_alloc) - n - 1);
                n++;
                if (n > nbests[p]) {
                    nbests[p] = n;
                    if (n >= 64*limbmax -1) {
                        limbmax *= 2;
                        work0 = realloc(work0, limbmax);
                        work1 = realloc(work1, limbmax);
                    }
                }
                i = 0;
                p++;
            }
        } else {
            // Dépiler
            n --;
            i = 0;
            while (i < p && !mpz_tstbit(sfpartition[i], n)) {
                i++;
            }
            mpz_clrbit(sfpartition[i], n);
            mpz_clrbit(sfpartitioninvert[i], 64*(sfpartitioninvert[i]->_mp_alloc) - n -1);
            if (i == p-1) {
                if (!mpz_size(sfpartition[i])) {
                    p--;
                }
            }
            i++;
        }
    }
    
    /*Nettoyage*/
    for (i=0; i<pmax; i++) {
        mpz_clear(sfpartition[i]);
        mpz_clear(sfpartitioninvert[i]);
    }
    free(sfpartition);
    free(sfpartitioninvert);
    free(work0);
    free(work1);
    return nbests[pmax-1];
}

int main(int argc, const char * argv[]) {
    int i;
    unsigned long schurNumbers[3];
    schurNumberIterative(3, schurNumbers);
    for (i=0; i<3; i++) {
        printf("Schur Number of %i : %lu\n", i+1, schurNumbers[i]);
    }
    return 0;
}

