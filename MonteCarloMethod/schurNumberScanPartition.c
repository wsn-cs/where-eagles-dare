//
//  schurNumberScanPartition.c
//  schurNumberMonteCarlo
//
//  Created by Gabriel Merlin on 01/05/2019.
//

#include <stdio.h>
#include "schurNumberNestedMonteCarloHeader.h"

unsigned int schurNumberScanPartitionFromFile(char *filename, partition_t *partitionstruc) {
    /*Crée une partition à partir d'un fichier texte.
     Des virgules séparent les entiers au sein d'un même ensemble
     et des points séparent les ensembles entre eux.
     La fonction renvoie le nombre d'ensembles p de la partition, ou 0 si un problème survient.*/
    FILE *fp;
    int c;
    char *intstr;
    unsigned long m, mmodbpl, mdivbpl;
    size_t len, maxlen;// Taille de la chaîne de caractères codant un entier
    unsigned int p;  // Nombre d'ensembles de la partition
    unsigned long n; // Nombre d'entiers
    mp_size_t limballoc;
    mp_limb_t **partition, **partitioninvert;
    
    
    /*Ouverture du flux*/
    fp = fopen(filename, "r");
    
    if (!fp) {
        /*Le fichier ne peut être ouvert.*/
        return 0;
    }
    
    /*Lecture préliminaire comptant les virgules et les points*/
    p = 0;
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
    
    /*Initialisation de la partition*/
    partition_init(p, n, partitionstruc);
    limballoc = partitionstruc->limballoc;
    
    /*Mise en place des éléments*/
    intstr = calloc(maxlen+1, sizeof(char));
    c = fgetc(fp);
    len = 0;
    partition = partitionstruc->partition;
    partitioninvert = partitionstruc->partitioninvert;
    while (c != EOF) {
        switch (c) {
            case ',':
                intstr[len] = '\0';
                m = atol(intstr)-1;
                /*Placer m dans la bonne partition*/
                mdivbpl = m / mp_bits_per_limb;
                mmodbpl = m % mp_bits_per_limb;
                (*partition)[mdivbpl] |= (mp_limb_t)1<<mmodbpl;
                (*partitioninvert)[limballoc - mdivbpl - 1] |= (mp_limb_t)1 << (mp_bits_per_limb - mmodbpl - 1);
                len = 0;
                break;
                
            case '.':
                intstr[len] = '\0';
                m = atol(intstr)-1;
                /*Placer m dans la bonne partition*/
                mdivbpl = m / mp_bits_per_limb;
                mmodbpl = m % mp_bits_per_limb;
                (*partition)[mdivbpl] |= (mp_limb_t)1<<mmodbpl;
                (*partitioninvert)[limballoc - mdivbpl - 1] |= (mp_limb_t)1 << (mp_bits_per_limb - mmodbpl - 1);
                partition++;
                partitioninvert++;
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
    partitionstruc->limbsize = limballoc;
    partitionstruc->p = p;
    partitionstruc->n = n;
    
    /*Fermeture du flux*/
    fclose(fp);
    
    return p;
}
