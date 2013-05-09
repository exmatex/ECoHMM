#ifndef __DEFORMATION_H_
#define __DEFORMATION_H_

#include "pmdTypes.h"

extern void forwardDeformation(SimFlat *s);
extern void reverseDeformation(SimFlat *s);
extern void matVec3(real_t *mat, real_t *vec);

#endif
