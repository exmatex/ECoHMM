#ifndef __PMD_H_
#define __PMD_H_

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "mytype.h"
#include "mycommand.h"
#include "pmdTypes.h"
#include "eamTypes.h"
#include "ljTypes.h"
#include "domains.h"
#include "constants.h"

#define DEBUGLEVEL 0
#define PMDDEBUGPRINTF(xxx,...) {if(xxx>DEBUGLEVEL) printf(__VA_ARGS__);}
#define fPMDDEBUGPRINTF(xxx,...) {if(xxx>DEBUGLEVEL) fprintf(__VA_ARGS__);}

extern SimFlat *blankSimulation(struct BasePotential *pot);
extern void destroySimulation(SimFlat **s);

extern void writeClsman(SimFlat *s, char *fname);
extern void printIt(SimFlat *s,FILE *fp);

extern double nTimeSteps(int n, SimFlat *s, real_t dt);

/* utility routines */
extern void breakinme();
extern void simulationAbort(int ineCode, char *inmsg);
extern double timeNow();
extern int computeForce(SimFlat *s);
extern void *doComputeWork(SimFlat *sim);
extern SimFlat *initSimFromCmdLine(int argc, char **argv, int *eamFlag, int *gpuFlag);

#endif
