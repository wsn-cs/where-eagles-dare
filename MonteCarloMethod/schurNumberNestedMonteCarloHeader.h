//
//  schurNumberNestedMonteCarloHeader.h
//  SchurNumber
//
//  Created by Gabriel Merlin on 22/04/2019.
//
//  Ce fichier est un header déclarant toutes les fonctions de calcul.
//

#ifndef schurNumberNestedMonteCarloHeader
#define schurNumberNestedMonteCarloHeader

#include <stdlib.h>
#include <gmp.h>

/*Fonction ajoutant les entiers un par un*/
unsigned long schurNumberSimpleMonteCarloLevelIteration(partition_t *sfpartitionstruc, unsigned int level,
  mp_limb_t **sfpartitionbest, unsigned int *pbestptr, unsigned int iternum);

#endif
