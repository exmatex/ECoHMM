#ifndef _REDISTRIBUTE_H_
#define _REDISTRIBUTE_H_

#include "pmdTypes.h"

// Drive local and remote redistribution of atoms
extern void redistribute(SimFlat* sim);

// Return list of boxes in halo
int* getHaloBoxIds(SimFlat* sim);

#endif
