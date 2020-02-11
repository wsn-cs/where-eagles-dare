//
//  main.c
//  schurNumberRecursive
//
//  Created by Gabriel Merlin on 01/12/2019.
//  Copyright © 2019 Gabriel Merlin. All rights reserved.
//

#include <stdio.h>
#include <libgen.h>
#include <unistd.h>
#include <fcntl.h>
#include "schurNumberConstrainedBuild.h"
#include "schurNumberIO.h"

void usage(char *cmdname) {
    fprintf(stderr,
            "usage: %s [-h] p set1 .. setp\n"\
            "\t-h : Print usage message\n",
            basename(cmdname));
}

int main(int argc, const char * argv[]) {
    char c;
    char should_execute = 1;
    
    while ((c = getopt(argc, argv, "h")) != -1) {
        switch (c) {
                
            case 'h':
                should_execute = 0;
                break;
                
            default:
                fprintf(stderr, "schurNumber: -%c: invalid option\n", c);
                should_execute = 0;
        }
    }

    /*Définition de la partition de contrainte*/
    unsigned long p = strtoul(argv[optind], NULL, 10);
    unsigned long pcons = argc - optind - 1;    // Nombre d'ensembles effectivement fournis

    if (p == 0) {
        fprintf(stderr, "schurNumber: no constraint partition\n");
        should_execute = 0;
    }
    
    if (p < pcons) {
        pcons = p;
    }
    
    if (should_execute) {
        mp_limb_t **constraint_partition = calloc(sizeof(mp_limb_t *), p);
        unsigned long n = schurNumberGetPartition(pcons, argv + optind + 1, constraint_partition);
        mp_size_t limbsize = ((unsigned long)n >> 6) + 1;
        if (pcons < p) {
            for (unsigned long j = pcons; j < p; j++) {
                constraint_partition[j] = calloc(sizeof(mp_limb_t), limbsize);
            }
        }
        
        /*Création de la partition*/
        schur_number_partition_t partition_s;
        schur_number_partition_alloc(&partition_s, 4 * limbsize, p);
        ADD_POINT(partition_s.partition[0], 0);
        
        schurNumberConstrainedBuild(&partition_s, constraint_partition, limbsize);
        
        schur_number_partition_dealloc(&partition_s);
        free(constraint_partition);
    } else {
        usage(argv[0]);
    }
    
    return 0;
}
