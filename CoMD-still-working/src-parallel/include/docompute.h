
#ifndef __DOCOMPUTE_H
#define __DOCOMPUTE_H

#include "pmdTypes.h"
#include "mytype.h"

void printArray(real_t* array, int n, char* name);
struct SimFlat* initSimFromCmdLine(int argc, char** argv);
void* doComputeWork(SimThreadData *data);

#endif
