//
//  schurNumberNestedMonteCarloMain.c
//  schurNumberMonteCarlo
//
//  Created by Weak Schur Number CS on 22/04/2019.
//  Ce fichier contient la fonction main qui permet de construire un exécutable.
//

#include <stdio.h>
#include <stdarg.h>
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

void printPartition(unsigned long p, unsigned long n, mp_limb_t **partition) {
    /*Affiche une partition.*/
    unsigned long limbn;
    mp_size_t limbsize;
    unsigned long i;
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
                if (limb & ((unsigned long)1<<j)) {
                    printf(" %lu", limbn*mp_bits_per_limb + j + 1);
                }
            }
        }
        printf("\n");
    }
}

/*Dans le cas où la fonction asprintf n'est pas implémentée, par exemple sous Windows.*/

#ifndef _vscprintf
/* For some reason, MSVC fails to honour this #ifndef. */
/* Hence function renamed to _vscprintf_so(). */
int _vscprintf_so(const char * format, va_list pargs) {
    int retval;
    va_list argcopy;
    va_copy(argcopy, pargs);
    retval = vsnprintf(NULL, 0, format, argcopy);
    va_end(argcopy);
    return retval;}
#endif // _vscprintf

#ifndef vasprintf
int vasprintf(char **strp, const char *fmt, va_list ap) {
    int len = _vscprintf_so(fmt, ap);
    if (len == -1) return -1;
    char *str = malloc((size_t) len + 1);
    if (!str) return -1;
    int r = vsnprintf(str, len + 1, fmt, ap); /* "secure" version of vsprintf */
    if (r == -1) return free(str), -1;
    *strp = str;
    return r;}
#endif // vasprintf

#ifndef asprintf
int asprintf(char *strp[], const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    int r = vasprintf(strp, fmt, ap);
    va_end(ap);
    return r;}
#endif // asprintf

int main(int argc, const char * argv[]) {
    char optc;
    unsigned int p, i;
    unsigned int level;
    unsigned int simulnum;
    unsigned int iternum;
    char printpartition;
    char statistics;
    unsigned long *narray;
    unsigned long n, nmax;
    unsigned long count, total;
    double nmean, nvar, delta;
    char *bfilename;
    char *ofilename;
    FILE *ofileptr;
    mp_limb_t **sfpartitionbestglobal;
    partition_t partitionbeginstruc;
    partition_t *partitionbeginptr;
    
    /*Set variables to default*/
    level = 1;
    simulnum = 1;
    iternum = 64;
    printpartition = 0;
    statistics = 0;
    bfilename = NULL;
    ofilename = NULL;
    partitionbeginptr = NULL;
    
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
    
    
    if (bfilename) {
        /*Débuter à partir d'une partition contenue dans le fichier filename.*/
        schurNumberScanPartitionFromFile(bfilename, &partitionbeginstruc);
        if (!i) {
            fprintf(stderr, "%s: unable to open file %s.\n", basename(argv[0]), bfilename);
            free(bfilename);
            return 0;
        }
        if (i < p) {
            fprintf(stderr, "%s: the partition in file %s contains more sets (%u sets) than specified (%u sets). Unable to do any computations.\n", basename(argv[0]), bfilename, i, p);
            free(bfilename);
            return 0;
        }
        partitionbeginptr = &partitionbeginstruc;
    }
    
    /*Allocation des variables.*/
    narray = calloc(simulnum, sizeof(unsigned long));
    sfpartitionbestglobal = calloc(p, sizeof(mp_limb_t *)); //Contiendra la plus grande partition sans-somme trouvée.
    for (i=0; i<p; i++) {
        sfpartitionbestglobal[i] = calloc(p, sizeof(mp_limb_t));
    }
    
    nmax = schurNumberSimpleNestedMonteCarlo(p, narray, level, simulnum, iternum, sfpartitionbestglobal, partitionbeginptr);
    
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
        /*Placer les entiers de narray dans un fichier au format tsv.
         Le format est entier_n\tnombre_d'occurences\n.*/
        ofileptr = fopen(ofilename, "w");
        if (ofileptr) {
            n = nmax;
            total = 0;
            while (total < simulnum) {
                count = 0;  // Compte le nombre de fois qu'une partition de [1, n] a été obtenue
                for (i=0; i<simulnum; i++) {
                    if (narray[i] == n) {
                        count++;
                        total++;
                    }
                }
                if (count) {
                    /*Au moins une partition sans-somme de [1, n] a été obtenue.*/
                    fprintf(ofileptr, "%lu\t%lu\r\n", n, count);
                }
                n--;
            }
            fclose(ofileptr);
        } else {
            fprintf(stderr, "%s: unable to write in %s.\n", basename(argv[0]), ofilename);
        }
        free(ofilename);
    }
    
    /*Déallocation des variables.*/
    free(narray);
    for (i=0; i<p; i++) {
        free(sfpartitionbestglobal[i]);
    }
    free(sfpartitionbestglobal);
    
    if (bfilename) {
        /*Déallocation de la partition initiale*/
        p = partitionbeginstruc.p;
        for (i=0; i<p; i++) {
            free(partitionbeginstruc.partition[i]);
            free(partitionbeginstruc.partitioninvert[i]);
        }
        free(partitionbeginstruc.partition);
        free(partitionbeginstruc.partitioninvert);
        free(bfilename);
    }
    
    return 0;
}

