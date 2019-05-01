//
//  schurNumberScanPartition.c
//  schurNumberMonteCarlo
//
//  Created by Gabriel Merlin on 01/05/2019.
//

#include <stdio.h>
#include "schurNumberNestedMonteCarloHeader.h"

void schurNumberScanPartitionFromFile(char *filename, partition_t *partitionstruc) {
    /*Crée une partition à partir d'un fichier texte.
     Des virgules séparent les entiers au sein d'un même ensemble
     et des points séparent les ensembles entre eux.*/
    FILE *fp;
    int c;
    char *intstr;
    unsigned long m, mmodbpl, mrem;
    size_t len, maxlen;// Taille de la chaîne de caractères codant un entier
    unsigned int i;
    unsigned int p;  // Nombre d'ensembles de la partition
    unsigned long n; // Nombre d'entiers
    mp_size_t limballoc;
    mp_limb_t **partition, **partitioninvert;
    
    
    /*Ouverture du flux*/
    fp = fopen(filename, "r");
    
    /*Lecture préliminaire comptant les virgules et les points*/
    p = 1;
    n = 0;
    c = getc(fp);
    len = 0;
    maxlen = 0;
    while (c != EOF) {
        switch (c) {
            case ',':
                n++;
                if (len > maxlen) {
                    maxlen = len;
                }
                len = 0;
                break;
            
            case '.':
                n++;
                p++;
                if (len > maxlen) {
                    maxlen = len;
                }
                len = 0;
                break;
                
            default:
                len++;
                break;
        }
        c = fgetc(fp);
    }
    rewind(fp);
    
    /*Allocation de la partition*/
    limballoc = 1 + n / mp_bits_per_limb;
    partition = calloc(p, sizeof(mp_limb_t *));
    partitioninvert = calloc(p, sizeof(mp_limb_t *));
    for (i=0; i<p; i++) {
        partition[i] = calloc(limballoc, sizeof(mp_limb_t));
        partitioninvert[i] = calloc(limballoc, sizeof(mp_limb_t));
    }
    
    /*Mise en place des éléments*/
    intstr = calloc(maxlen+1, sizeof(char));
    c = fgetc(fp);
    len = 0;
    p = 0;
    while (c != EOF) {
        switch (c) {
            case ',':
                intstr[len] = '\0';
                m = atol(intstr)-1;
                /*Placer m dans la bonne partition*/
                mmodbpl = m / mp_bits_per_limb;
                mrem = m%mp_bits_per_limb;
                partition[p][mmodbpl] |= (mp_limb_t)1<<mrem;
                partitioninvert[p][limballoc - mmodbpl - 1] |= (mp_limb_t)1 << (mp_bits_per_limb - mrem - 1);
                len = 0;
                break;
                
            case '.':
                intstr[len] = '\0';
                m = atol(intstr)-1;
                /*Placer m dans la bonne partition*/
                mmodbpl = m / mp_bits_per_limb;
                mrem = m%mp_bits_per_limb;
                partition[p][mmodbpl] |= (mp_limb_t)1<<mrem;
                partitioninvert[p][limballoc - mmodbpl - 1] |= (mp_limb_t)1 << (mp_bits_per_limb - mrem - 1);
                p++;
                len = 0;
                break;
                
            default:
                intstr[len] = c;
                len++;
                break;
        }
        c = fgetc(fp);
    }
    free(intstr);
    
    /*Sauvegarde dans partition_t*/
    partitionstruc->limballoc = limballoc;
    partitionstruc->limbsize = limballoc;
    partitionstruc->partition = partition;
    partitionstruc->partitioninvert = partitioninvert;
    partitionstruc->pmax = p;
    partitionstruc->p = p;
    partitionstruc->n = n;
    
    /*Fermeture du flux*/
    fclose(fp);
}
