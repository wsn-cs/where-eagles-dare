//
//  schurNumberSymmetricBranch&Bound.c
//  schurNumberSymmetric
//
//  Created by wsn-cs on 30/05/2019.
//

#include "schurNumberSymmetricHeader.h"
#include "rscan1.c"

unsigned long schurNumberSymmetricBranch(unsigned long n, partition_t *partitionstruc, unsigned long *depths, char completesearch) {
    /*
     Cette fonction trouve les partitions sans-sommes symmétriques de {1,…,n} à p huches,
     de profondeur équivalente > mindepth.
     
     La partition est représentée comme un tableau de grands entiers sfpartition.
     Chaque grand entier représente un ensemble grâce à ses bits: si le bit k vaut 1,
     alors k appartient à l'ensemble; sinon il n'y appartient pas.
     A cause de la symétrie, seuls les [(n+1)/2] premiers bits sont représentés.
     
     Pour tester si il est possible d'ajouter m+1 et n-m à l'ensemble P sans-somme,
     on effectue un & entre les [(m+1)/2] premiers bits avec les [(m+1)/2] derniers.
     Si le résultat est 0, il est possible d'ajouter m+1.
     
     Lorsqu'un entier m n'a pu être inséré dans aucune huches, on revient au sup sur les huches
     du plus petit n'≥[(m+1)/2] tel que n' et m-n' appartiennent simultanément à cette huche.
     */
    unsigned long m;
    unsigned long mblocking;
    unsigned long mblockingmax;
    unsigned long mmed;
    unsigned int p, q;
    unsigned long i;
    unsigned long iblocking;
    char isSumFree;
    char isAppendable;  //Vrai si il est possible de former une partition sans-somme contenant n
    mp_limb_t *work0;
    mp_limb_t *work1;
    mp_limb_t **sfpartition;
    mp_limb_t **sfpartitioninvert;
    mp_size_t limballoc;    // Nombre de limbes alloué à chaque ensemble de sfpartition et à work0 et work1
    mp_size_t limbsize;     // Nombre de limbes utilisés par les ensembles de sfpartition
    mp_size_t wlimbsize;    // Nombre de limbes utilisés par work0
    unsigned int mmodbpl;
    unsigned int shift;
    unsigned long nmed;
    unsigned long maxdepth;
    char partitionfound;
    mp_limb_t mask;
    mp_limb_t mask2;
    
    /*Initialisation*/
    p = partitionstruc->pmax;
    nmed = (n+1) >> 1;
    limballoc = partitionstruc->limballoc;
    sfpartition = partitionstruc->partition;
    sfpartitioninvert = partitionstruc->partitioninvert;
    work0 = partitionstruc->work0;
    work1 = partitionstruc->work1;
    
    /*Création de la première partition {{1}}*/
    m = partitionstruc->n;  // Taille de l'intervalle
    q = partitionstruc->p;  // Nombre de huches non vides
    i = 0;                  // Huche où placer l'entier suivant
    limbsize = partitionstruc->limbsize;
    wlimbsize = 1;          // Taille en limbe de la seconde moitié
    mmodbpl = 1;
    shift = mp_bits_per_limb - mmodbpl;
    mmed = 0;
    mblockingmax = 0;
    isAppendable = 0;
    maxdepth = 0;
    partitionfound = 0;
    
    /*Itération jusqu'à énumérer toutes les partions sans-somme à au plus p huches*/
    while (m > 0 && m < nmed) {
        // Placer m+1 dans une des huches en conservant la sans-sommité
        while (i < q) {
            // Tester si l'ensemble obtenu en ajoutant m+1 à la huche i est sans-somme
            mpn_rshift(work1, sfpartitioninvert[i] + limballoc - limbsize, limbsize, shift);
            mpn_and_n(work0, sfpartition[i], work1, wlimbsize);
            isSumFree = mpn_zero_p(work0, wlimbsize);
            
            // Tester si l'ensemble obtenu en ajoutant n-m à la huche i est sans-somme
            if (isSumFree) {
                // Il suffit de vérifier si il existe deux entiers ≤ (n+1)/2 de somme n-m, en amenant n-m-1 sur 1
                mpn_rshift(work1, sfpartitioninvert[i] + limballoc - limbsize, limbsize, mp_bits_per_limb - ((n-m-1) % mp_bits_per_limb));
                mask = (mp_limb_t)1<<mmodbpl;
                sfpartition[i][limbsize -1] |= mask;
                mpn_and_n(work0, sfpartition[i], work1, limbsize);
                isSumFree = mpn_zero_p(work0, limbsize);
                
                if (isSumFree) {
                    // Ajouter n+1 à la huche i
                    mask2 = (mp_limb_t)1<<(shift - 1);
                    sfpartitioninvert[i][limballoc - limbsize] |= mask2;
                    m++;
                    mmodbpl = m % mp_bits_per_limb;
                    shift = mp_bits_per_limb - mmodbpl;
                    mmed = (m+1)>>1;
                    wlimbsize = (mmed / mp_bits_per_limb) + 1;
                    mmed--;
                    mblockingmax = 0;
                    if (!mmodbpl) {
                        for (i=0; i<p; i++) {
                            sfpartition[i][limbsize] = (mp_limb_t)0;
                            sfpartitioninvert[i][limballoc - limbsize -1] = (mp_limb_t)0;
                        }
                        limbsize++;
                    }
                    i = 0;
                    isAppendable = 0;
                    if (m >= nmed) {
                        partitionfound = 1;
                        if (!completesearch) {
                            return maxdepth;
                        }
                        // Retirer m
                        sfpartition[i][limbsize -1] ^= mask;
                        sfpartitioninvert[i][limballoc - limbsize] ^= mask2;
                        i++;
                    }
                } else {
                    // Retirer m
                    sfpartition[i][limbsize -1] ^= mask;
                    // Trouver le plus petit m' tel que n-(m+m') appartienne à la huche i
                    mblocking = mpn_scan1(work0, mmed);
                    if (mblocking > mblockingmax) {
                        mblockingmax = mblocking;
                        iblocking = i;
                    }
                    i++;
                }
            } else {
                // Trouver le plus petit m' tel que m-m' appartienne à la huche i
                mblocking = m - mpn_rscan1(work0, mmed);
                if (mblocking > mblockingmax) {
                    mblockingmax = mblocking;
                    iblocking = i;
                }
                i++;
            }
        }
        if (i == q && q < p) {
            // Cas particulier où il faut remplir une nouvelle huche
            mask = (mp_limb_t)1<<mmodbpl;
            sfpartition[i][limbsize -1] |= mask; // Ajoute m+1 à la huche i
            mask2 = (mp_limb_t)1<<(shift - 1);
            sfpartitioninvert[i][limballoc - limbsize] |= mask2;
            m++;
            mmodbpl = m % mp_bits_per_limb;
            shift = mp_bits_per_limb - mmodbpl;
            mmed = (m+1)>>1;
            wlimbsize = (mmed / mp_bits_per_limb) + 1;
            mmed--;
            mblockingmax = 0;
            if (!mmodbpl) {
                for (i=0; i<p; i++) {
                    sfpartition[i][limbsize] = (mp_limb_t)0;
                    sfpartitioninvert[i][limballoc - limbsize -1] = (mp_limb_t)0;
                }
                limbsize++;
            }
            
            if (m < depths[q]) {
                // Si la profondeur est ≤ mindepth, alors il faut revenir en arrière
                i = p;
                isAppendable = 1;
            } else {
                i = 0;
                isAppendable = 0;
                if (q == p-1 && m > maxdepth) {
                    maxdepth = m;
                }
            }
            q++;
        } else if(isAppendable) {
            // Dépiler d'un entier
            if (!mmodbpl) {
                // L'entier m commençait un limbe
                limbsize--;
            }
            m--;
            mmodbpl = m % mp_bits_per_limb;
            shift = mp_bits_per_limb - mmodbpl;
            mmed = (m+1)>>1;
            wlimbsize = (mmed / mp_bits_per_limb) + 1;
            mmed--;
            mask = (mp_limb_t)1 << mmodbpl;
            // Trouver la huche i contenant l'entier à dépiler
            i = 0;
            while (i < q && !(sfpartition[i][limbsize-1] & mask)) {
                i++;
            }
            sfpartition[i][limbsize-1] ^= mask;
            mask = (mp_limb_t)1<<(shift - 1);
            sfpartitioninvert[i][limballoc - limbsize] ^= mask;
            if (i == q - 1) {
                if (mpn_zero_p(sfpartition[i], limbsize)) {
                    // Cette huche est désormais vide
                    q--;
                }
            }
            i++;
            isAppendable = 1;
        } else {
            // Dépiler jusqu'à revenir à nblockingmax
            m = mblockingmax;
            limbsize = (m / mp_bits_per_limb) + 1;
            mmodbpl = m % mp_bits_per_limb;
            shift = mp_bits_per_limb - mmodbpl;
            mmed = (m+1)>>1;
            wlimbsize = (mmed / mp_bits_per_limb) + 1;
            mmed--;
            mask = (mp_limb_t)GMP_NUMB_MAX >> shift;
            mask2 = (mp_limb_t)GMP_NUMB_MAX << (shift);
            for (i=0; i<p; i++) {
                sfpartition[i][limbsize-1] &= mask;
                sfpartitioninvert[i][limballoc - limbsize] &= mask2;
            }
            while (mpn_zero_p(sfpartition[q-1], limbsize)) {
                q--;
            }
            i = iblocking + 1;
            isAppendable = 1;
        }
    }
    
    if (partitionfound) {
        return maxdepth;
    }
    
    return 0;
}
