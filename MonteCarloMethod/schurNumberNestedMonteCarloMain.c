//
//  schurNumberNestedMonteCarloMain.c
//  schurNumberMonteCarlo
//
//  Created by Weak Schur Number CS on 22/04/2019.
//  Ce fichier contient la fonction main qui permet de construire un exécutable.
//

#include "asprintf.h"
#include <stdio.h>
#include <unistd.h>
#include <libgen.h>
#include <string.h>
#include <math.h>
#include "schurNumberNestedMonteCarloHeader.h"

void usage(char *cmdname) {
    fprintf(stderr,
            "usage: %s [-h] [-pv] [-l level] [-s simulnum] [-i iternum] partnumber\n"\
            "\t-b filename : Specify a filename for an initial partition\n"\
            "\t-o filename : Create a file named filename where to save the generated numbers.\n"\
            "\t-p : Print the found partition to stdout.\n"\
            "\t-v : Print mean and standard deviation to stdout.\n"\
            "\t-l level: Specify the nesting level. By default level = 1.\n"\
            "\t-s simulnum: Specify the number of simulation at the top level. By default simulnum = 1.\n"\
            "\t-i iternum: Specify the number of iteration at level 0. By default iternum = 64.\n"\
            "\t-h: Print usage message.\n",
            basename(cmdname));
}

void printPartition(unsigned int p, unsigned long n, mp_limb_t **partition) {
    /*Affiche une partition.*/
    unsigned long limbn;
    mp_size_t limbsize;
    unsigned int i;
    mp_bitcnt_t j;
    mp_limb_t *set;
    mp_limb_t limb;
    limbsize = (n>>6) + 1;
    printf("Partition:\n");
    for (i=0; i<p; i++) {
        printf("\t");
        set = partition[i];
        for (limbn=0; limbn<limbsize; limbn++) {
            limb = set[limbn];
            for (j=0; j<mp_bits_per_limb; j++) {
                if (limb & ((mp_limb_t)1<<j)) {
                    printf(" %lu", limbn*mp_bits_per_limb + j + 1);
                }
            }
        }
        printf("\n");
    }
}

unsigned int makeInitialPartition(char *filename, partition_t *partitionbeginstruc) {
    /*Si filename != NULL, lit le fichier filename et initie partitionbeginstruc à la partition
        contenue dans le fichier. Si filename == NULL, partitionbeginstruc = {{1}}.
     Renvoie la taille de la partition en l'absence de problème, 0 sinon.*/
    unsigned int p;
    
    if (filename) {
        /*Débuter à partir d'une partition contenue dans le fichier filename.*/
        p = schurNumberScanPartitionFromFile(filename, partitionbeginstruc);
    } else {
        /*Prendre pour partition initiale {{1}}*/
        p = 1;
        partition_init(1, 1, partitionbeginstruc);
        partitionbeginstruc->p = 1;
        partitionbeginstruc->limbsize = 1;
        partitionbeginstruc->n = 1;
        **(partitionbeginstruc->partition) = (mp_limb_t)1;
        **(partitionbeginstruc->partitioninvert) = (mp_limb_t)1<<(mp_bits_per_limb - 1);
    }
    
    return p;
}

char saveDistribution(char *filename, unsigned long *distribution, size_t size, unsigned long nmax) {
    /*Enregistre les n entiers de distribution dans un fichier au format tsv.
       Le format est entier_n\tnombre_d'occurences\n.
     Renvoie 1 si tout s'est bien passé, 0 sinon.*/
    FILE *fp;
    unsigned long n;
    unsigned long i;
    unsigned long count, total;
    
    fp = fopen(filename, "w");
    if (fp) {
        n = nmax;
        total = 0;
        while (total < size) {
            count = 0;  // Compte le nombre de fois qu'une partition de [1, n] a été obtenue
            for (i=0; i<size; i++) {
                if (distribution[i] == n) {
                    count++;
                    total++;
                }
            }
            if (count) {
                /*Au moins une partition sans-somme de [1, n] a été obtenue.*/
                fprintf(fp, "%lu\t%lu\r\n", n, count);
            }
            n--;
        }
        fclose(fp);
        return 1;
    }
    
    return 0;
}

int main(int argc, const char * argv[]) {
    char optc;
    unsigned int p, i;
    unsigned int level;
    unsigned int simulnum;
    unsigned int iternum;
    char printpartition;
    char statistics;
    char fileproblem;
    unsigned long *narray;
    unsigned long nmax;
    double nmean, nvar, delta;
    char *bfilename;
    char *ofilename;
    mp_limb_t **sfpartitionbestglobal;
    partition_t partitionbeginstruc;
    
    /*Set variables to default*/
    level = 1;
    simulnum = 1;
    iternum = 64;
    printpartition = 0;
    statistics = 0;
    fileproblem = 0;
    bfilename = NULL;
    ofilename = NULL;
    
    while ((optc = getopt(argc, argv, "hpvl:s:i:b:o:")) != -1) {
        /*Parse arguments*/
        switch (optc) {
                
            case 'b':
                asprintf(&bfilename, "%s", optarg);
                break;
                
            case 'o':
                asprintf(&ofilename, "%s", optarg);
                break;
                
            case 'h':
                usage(argv[0]);
                break;
                
            case 'p':
                printpartition = 1;
                break;
            
            case 'v':
                statistics = 1;
                break;
                
            case 'l':
                level = atoi(optarg);
                break;
                
            case 's':
                simulnum = atoi(optarg);
                break;
                
            case 'i':
                iternum = atoi(optarg);
                break;
                
            default:
                fprintf(stderr, "schurNumber: -%c: invalid option\n", optc);
                usage(argv[0]);
                return 0;
                break;
        }
    }
    
    if (argc - optind <= 0) {
        fprintf(stderr, "schurNumber: %s: invalid argument\n", argv[optind]);
        usage(argv[0]);
        return 0;
    }
    p = atoi(argv[optind]);
    if (p == 0) {
        fprintf(stderr, "schurNumber: %s: invalid argument\n", argv[optind]);
        usage(argv[0]);
        return 0;
    }
    
    i = makeInitialPartition(bfilename, &partitionbeginstruc);
    if (!i) {
        /*Problème dans l'ouverture du fichier bfilename*/
        fprintf(stderr, "%s: unable to open file %s.\n", basename(argv[0]), bfilename);
        fileproblem = 1;
    }
    if (i > p) {
        fprintf(stderr, "%s: the partition in file %s contains more sets (%u sets) than specified (%u sets). Unable to do any computations.\n", basename(argv[0]), bfilename, i, p);
        fileproblem = 1;
    }
    if (bfilename) {
        free(bfilename);
    }
    if (ofilename) {
        i = saveDistribution(ofilename, NULL, 0, 0);
        if (!i) {
            fprintf(stderr, "%s: unable to open %s for writing.\n", basename(argv[0]), ofilename);
            free(ofilename);
            fileproblem = 1;
        }
    }
    
    if (!fileproblem) {
        /*Allocation des variables.*/
        narray = calloc(simulnum, sizeof(unsigned long));
        sfpartitionbestglobal = calloc(p, sizeof(mp_limb_t *)); //Contiendra la plus grande partition sans-somme trouvée.
        for (i=0; i<p; i++) {
            sfpartitionbestglobal[i] = calloc(p, sizeof(mp_limb_t));
        }
        
        nmax = schurNumberSimpleNestedMonteCarlo(p, narray, level, simulnum, iternum, sfpartitionbestglobal, &partitionbeginstruc);
        
        printf("Schur Number S(%u) ≥ %lu\n", p, nmax);
        
        if (printpartition) {
            /*Afficher la meilleure partition sans-somme trouvée.*/
            printPartition(p, nmax, sfpartitionbestglobal);
        }
        
        if (statistics) {
            /*Calculer la valeur moyenne et l'écart-type grâce à la méthode de Welford.*/
            nmean = 0;
            nvar = 0;
            for (i=0; i < simulnum; i++) {
                delta = (double)narray[i] - nmean;
                nmean += delta/(i+1);
                nvar += delta * ((double)narray[i] - nmean);
            }
            nvar = nvar / (simulnum-1);
            printf("Moyenne : %f\nEcart-type : %f\n", nmean, sqrt(nvar));
        }
        
        if (ofilename) {
            i = saveDistribution(ofilename, narray, simulnum, nmax);
            if (!i) {
                fprintf(stderr, "%s: unable to write in file %s.\n", basename(argv[0]), ofilename);
            }
            free(ofilename);
        }
        
        /*Déallocation des variables.*/
        free(narray);
        for (i=0; i<p; i++) {
            free(sfpartitionbestglobal[i]);
        }
        free(sfpartitionbestglobal);
        partition_unalloc(&partitionbeginstruc);
    }
    
    return 0;
}
