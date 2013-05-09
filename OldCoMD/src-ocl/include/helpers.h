#ifndef HELPER_H
#define HELPER_H

#include "cl_utils.h"
#include "pmdOCL.h"
#include "pmdTypes.h"
#include "domains.h"
#include "constants.h"
#include "quaternion.h"

#define DIAG_LEVEL 0

typedef struct HostVec {
    cl_real* x;
    cl_real* y;
    cl_real* z;
} HostVec;

typedef struct DevVec {
    cl_mem x;
    cl_mem y;
    cl_mem z;
} DevVec;

typedef struct HostGridSoa {
    HostVec rBox;
    cl_int* neighborList;
    cl_int* nNeighbors;
    cl_int* nAtoms;
    cl_real* bounds;
    int nBoxesIntSize;
    int nBoxesNeighborSize;
    int nBoxesRealSize;
} HostGridSoa;

typedef struct HostGridAos {
    cl_real4* rBox;
    cl_int* neighborList;
    cl_int* nNeighbors;
    cl_int* nAtoms;
    cl_real* bounds;
    int nBoxesIntSize;
    int nBoxesNeighborSize;
    int nBoxesRealSize;
} HostGridAos;

typedef struct DevGridSoa {
    DevVec rBox;
    cl_mem neighborList;
    cl_mem nNeighbors;
    cl_mem nAtoms;
    cl_mem bounds;
} DevGridSoa;

typedef struct DevGridAos {
    cl_mem rBox;
    cl_mem neighborList;
    cl_mem nNeighbors;
    cl_mem nAtoms;
    cl_mem bounds;
} DevGridAos;

typedef struct HostEamPot {
    cl_real* rho;
    cl_real* phi;
    cl_real* F;
    cl_int* nValues;
    cl_real cutoff;
    int rhoPotSize;
    int phiPotSize;
    int fPotSize;
} HostEamPot;

typedef struct DevEamPot {
    cl_mem rho;
    cl_mem phi;
    cl_mem F;
    cl_mem nValues;
    cl_real cutoff;
} DevEamPot;

typedef struct HostEamCh {
    cl_real* rho;
    cl_real* phi;
    cl_real* F;
    cl_int* nValues;
    cl_real cutoff;
    int rhoChebSize;
    int phiChebSize;
    int fChebSize;
} HostEamCh;

typedef struct DevEamCh {
    cl_mem rho;
    cl_mem phi;
    cl_mem F;
    cl_mem nValues;
    cl_real cutoff;
} DevEamCh;

typedef struct HostLjPot {
    cl_real cutoff;
    cl_real sigma;
    cl_real epsilon;
} HostLjPot;

typedef struct DevLjPot {
    cl_real cutoff;
    cl_real sigma;
    cl_real epsilon;
} DevLjPot;

typedef struct HostSimSoa {
    // grid array values
    HostVec r;
    HostVec p;
    HostVec f;
    cl_real* e;
    cl_real* m;
    cl_real* fi;
    cl_real* rho;
    // grid data 
    HostGridSoa grid;
    // real scalars
    cl_real dt;
    cl_real rMass;
    cl_real cfac;
    cl_real energy;
    // integer values
    int arraySize;
    int nBoxesIntSize;
    int nBoxesRealSize;
    int nBoxesNeighborSize;
    int nCells;
    int nTot;
    // eam flag
    int eamFlag;
    HostEamPot eamPot;
    HostEamCh eamCh;
    HostLjPot ljPot;
} HostSimSoa;

typedef struct HostSimAos {
    // grid array values
    cl_real4* r;
    cl_real4* p;
    cl_real4* f;
    cl_real* e;
    cl_real* m;
    cl_real* fi;
    cl_real* rho;
    // grid data 
    HostGridAos grid;
    // real scalars
    cl_real dt;
    cl_real rMass;
    cl_real cfac;
    cl_real energy;
    // integer values
    int arraySize;
    int nBoxesIntSize;
    int nBoxesRealSize;
    int nBoxesNeighborSize;
    int nCells;
    int nTot;
    // eam flag
    int eamFlag;
    HostEamPot eamPot;
    HostEamCh eamCh;
    HostLjPot ljPot;
} HostSimAos;

typedef struct DevSimSoa {
    // grid array values
    DevVec r;
    DevVec p;
    DevVec f;
    cl_mem e;
    cl_mem m;
    cl_mem fi;
    cl_mem rho;
    // grid data 
    DevGridSoa grid;
    // note all scalars can be passed directly to the kernels
    // real scalars
    cl_real dt;
    cl_real rMass;
    cl_real cfac;
    cl_real energy;
    // integer values
    int arraySize;
    int nBoxesIntSize;
    int nBoxesRealSize;
    int nBoxesNeighborSize;
    int nCells;
    DevEamPot eamPot;
    DevEamCh eamCh;
    DevLjPot ljPot;
} DevSimSoa;

typedef struct DevSimAos {
    // grid array values
    cl_mem r;
    cl_mem p;
    cl_mem f;
    cl_mem e;
    cl_mem m;
    cl_mem fi;
    cl_mem rho;
    // grid data 
    DevGridAos grid;
    // note all scalars can be passed directly to the kernels
    // real scalars
    cl_real dt;
    cl_real rMass;
    cl_real cfac;
    cl_real energy;
    // integer values
    int arraySize;
    int nBoxesIntSize;
    int nBoxesRealSize;
    int nBoxesNeighborSize;
    int nCells;
    DevEamPot eamPot;
    DevEamCh eamCh;
    DevLjPot ljPot;
} DevSimAos;
    

/* General helper utils */

void printArray(real_t* array, int n, char *name);

void printSim(SimFlat *s,FILE *fp);

void createDevVec(DevVec *a_D, int arraySize);

void getVector(cl_mem ax_D, cl_mem ay_D, cl_mem az_D,
	cl_real* ax_H, cl_real* ay_H, cl_real* az_H,
	int arraySize);

void putVec(HostVec a_H, DevVec a_D, int arraySize);

void putVector(cl_real* ax_H, cl_real* ay_H, cl_real* az_H, cl_mem ax_D, cl_mem ay_D, cl_mem az_D, int arraySize);

void putEamPot(HostEamPot eamPotH, DevEamPot eamPotD);

void oclRunKernel(cl_kernel kernel, cl_event *event, size_t* nGlobal, size_t* nLocal);

void getElapsedTime(cl_event event, cl_real* elapsed_time, cl_real* enqueuedTime);

void computeForceOCL(cl_kernel* forceKernels, cl_event* forceEvent, size_t* nGlobal, size_t* nLocal, int eamFlag, int nTot, cl_real* t_kern);

/* SoA (default) variants) */

void getPrintState(DevSimSoa simDevSoa, HostSimSoa simHostSoa);

void computePrintEnergy(DevSimSoa simDevSoa, HostSimSoa simHostSoa);

void printState(HostSimSoa simHostSoa, int nCells);

void createDevGrid(DevGridSoa *gridDev, int nBoxesRealSize, int nBoxesNeighborSize, int nBoxesIntSize);

void getVec(DevVec a_D, HostVec a_H, int arraySize);

void buildModulesSoa(cl_kernel *forceKernels, cl_kernel *advancePosition, cl_kernel *AdvanceVelocity, cl_kernel *Viz, 
        HostSimSoa simHostSoa, size_t *nLocal, size_t *nGlobal);

void initHostSim (HostSimSoa *simHostSoa, SimFlat *sim);

void initDevSim(DevSimSoa *simDevSoa, HostSimSoa *simHostSoa);

void putSim(HostSimSoa simHostSoa, DevSimSoa simDevSoa);

void putGrid(HostGridSoa gridHost, DevGridSoa gridDev);

void setLJArgs(cl_kernel ljForce, DevSimSoa simDevSoa);

void setEAMArgs(cl_kernel *forceKernels, DevSimSoa simDevSoa);

void setAVArgs(cl_kernel AdvanceVelocity, DevSimSoa simDevSoa, cl_real dt);

void setAPArgs(cl_kernel advancePosition, DevSimSoa simDevSoa, cl_real dt);

void FreeSims(HostSimSoa simHostSoa, DevSimSoa simDevSoa);

/* AoS variants */

void getPrintStateAoS(DevSimAos simDevSoa, HostSimAos simHostSoa);

void computePrintEnergyAoS(DevSimAos simDevSoa, HostSimAos simHostSoa);

void printStateAoS(HostSimAos simHostSoa, int nCells);

void createDevGridAoS(DevGridAos *gridDev, int nBoxesRealSize, int nBoxesNeighborSize, int nBoxesIntSize);

void getVecAoS(cl_mem a_D, cl_real4* a_H, int arraySize);

void buildModulesAoS(cl_kernel *forceKernels, cl_kernel *advancePosition, cl_kernel *AdvanceVelocity, cl_kernel *Viz, 
        HostSimAos simHostSoa, size_t *nLocal, size_t *nGlobal);

void initHostSimAoS (HostSimAos *simHostSoa, SimFlat *sim);

void initDevSimAoS(DevSimAos *simDevSoa, HostSimAos *simHostSoa);

void putSimAoS(HostSimAos simHostSoa, DevSimAos simDevSoa);

void putGridAoS(HostGridAos gridHost, DevGridAos gridDev);

void setLJArgsAoS(cl_kernel ljForce, DevSimAos simDevSoa);

void setEAMArgsAoS(cl_kernel *forceKernels, DevSimAos simDevSoa);

void setAVArgsAoS(cl_kernel AdvanceVelocity, DevSimAos simDevSoa, cl_real dt);

void setAPArgsAoS(cl_kernel advancePosition, DevSimAos simDevSoa, cl_real dt);

void FreeSimsAoS(HostSimAos simHostSoa, DevSimAos simDevSoa);

/* Graphics kernels */

void oclGraphics(cl_kernel vizKernel, DevSimSoa simDevSoa, size_t* nGlobal, size_t* nLocal);

void oclRender();

void oclInitInterop(int ncells);

#endif
