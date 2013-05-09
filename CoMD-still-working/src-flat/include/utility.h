
#ifndef __UTILITY_H
#define __UTILITY_H

#include "pmdTypes.h"

double timeNow();
void breakinme();
char *dupString(char *str);
void simulationAbort(int ineCode, char *inmsg);
struct SimFlat *blankSimulation(struct BasePotential *pot);
void destroySimulation(SimFlat **ps);

#endif
