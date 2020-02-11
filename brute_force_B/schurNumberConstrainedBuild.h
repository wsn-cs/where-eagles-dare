//
//  schurNumberConstrainedBuild.h
//  schurNumberRecursive
//
//  Created by rubis on 08/01/2020.
//  Copyright © 2020 rubis. All rights reserved.
//

#ifndef schurNumberConstrainedBuild_h
#define schurNumberConstrainedBuild_h

#include <stdio.h>
#include "schurNumberIO.h"
#include "schurNumberPartitionStruc.h"

unsigned long schurNumberConstrainedBuild(schur_number_partition_t *partitionstruc, mp_limb_t **constraint_partition, mp_size_t constraint_size, struct schurNumberIOAction *action);

#endif /* schurNumberConstrainedBuild_h */
