
#include "helpers.h"

#ifdef INTEROP_VIZ

#if defined (__APPLE__) || defined(MACOSX)
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#include <OpenGL/OpenGL.h>
#include "OpenCL/cl.h"
#include "OpenCL/cl_gl.h"
#include "OpenCL/cl_gl_ext.h"
#include "OpenCL/cl_ext.h"
#include "OpenGL/CGLDevice.h"
#define GL_SHARING_EXTENSION "cl_APPLE_gl_sharing"
#else
#include <GL/gl.h>
#include <GL/glx.h>
#include "CL/cl_gl.h"
#endif

GLuint vboBuffers[3];
cl_mem vboResources[3];
int g_ncells;
Quaternion q;
float cameraFOV;
float centerX, centerY, centerZ;

#endif

// Make sure PASS_1, PASS_2, PASS_3 flags are correctly set!
// Make sure PASS_2 matches that is eam_kernels.c

#define PASS_1 1
#define PASS_2 0
#define PASS_3 1

/** This chunk of code exists only because older OpenCL implementations
  on Apple didn't treat vector types correctly. If you see errors during build time
  about subscripted values, then you have a newer version of OpenCL, and these 'if'
  clauses are no longer needed. They appear a couple more times in this file.
  You can now set manually toggle these chunks of code using the APPLE_OCL_10 flag
  below. If you have OpenCL 1.1 or newer, set it to 0
 **/
#define APPLE_OCL_10 1

void dummyTest()
{
    cl_real4 dummy;

#if (defined (__APPLE__) || defined(MACOSX)) && (APPLE_OCL_10)
    dummy[0] = 1.0;
    dummy[1] = 1.0;
#else
    dummy.s[0] = 1.0;
    dummy.s[1] = 1.0;
#endif
}

void printArray(real_t* array, int n, char *name)
{
    printf("%s;\n", name);
    for (int i=0;i<n;i++)
    {
        printf("%d, %17.9e\n", i, array[i]);
    }
}

void computePrintEnergy( DevSimSoa simDevSoa, HostSimSoa simHostSoa)
{

    /**
      Copy the array of particle energies from device to host. 
      The total energy is summed and returned in the simHostSoa.energy variable
     **/

    double eLocal;

    oclCopyToHost(simDevSoa.e, simHostSoa.e, simHostSoa.arraySize);

    eLocal = 0.0;
    for (int iBox=0;iBox<simHostSoa.nCells;iBox++)
    {
        for (int iAtom=0;iAtom<simHostSoa.grid.nAtoms[iBox];iAtom++)
	{
            eLocal += simHostSoa.e[iBox*MAXATOMS + iAtom];
        }
    }
    simHostSoa.energy = (cl_real)eLocal;

    //printf("System energy = %30.20f\n", simHostSoa.energy);
    printf(" %30.20e", simHostSoa.energy);

}

void computePrintEnergyAoS( DevSimAos simDevSoa, HostSimAos simHostSoa)
{

    /**
      Copy the array of particle energies from device to host. 
      The total energy is summed and returned in the simHostSoa.energy variable
     **/

    double eLocal;

    oclCopyToHost(simDevSoa.e, simHostSoa.e, simHostSoa.arraySize);

    eLocal = 0.0;
    for (int iBox=0;iBox<simHostSoa.nCells;iBox++)
    {
        for (int iAtom=0;iAtom<simHostSoa.grid.nAtoms[iBox];iAtom++)
	{
            eLocal += simHostSoa.e[iBox*MAXATOMS + iAtom];
        }
    }
    simHostSoa.energy = (cl_real)eLocal;

    //printf("System energy = %30.20f\n", simHostSoa.energy);
    printf(" %30.20e", simHostSoa.energy);

}

void printState( HostSimSoa simHostSoa, int nToPrint)
{

    /**
      Print the box index, atom index, position, momentum and force for the 
      first nCells boxes of the simulation
     **/

    int i, iBox;
    int atomCount = 0;
    iBox = 0;
    printf("System state:\n");
    while (atomCount < nToPrint)
    {
        for (int iAtom=0;iAtom<simHostSoa.grid.nAtoms[iBox];iAtom++)
	{

            i = iBox*MAXATOMS + iAtom;

            printf("%02d, %02d, "
                    "X=(%+020.12e %+020.12e %+020.12e) 1 "
                    "P=(%+020.12e %+020.12e %+020.12e) "
                    "F=(%+020.12e %+020.12e %+020.12e)\n",
                    iBox, iAtom, 
                    simHostSoa.r.x[i],simHostSoa.r.y[i],simHostSoa.r.z[i],
                    simHostSoa.p.x[i],simHostSoa.p.y[i],simHostSoa.p.z[i],
                    simHostSoa.f.x[i],simHostSoa.f.y[i],simHostSoa.f.z[i]);
            atomCount ++;

        }
        iBox ++;
    }

}

void printStateAoS( HostSimAos simHostSoa, int nToPrint)
{

    /**
      Print the box index, atom index, position, momentum and force for the 
      first nCells boxes of the simulation
     **/

    int i, iBox;
    int atomCount = 0;
    iBox = 0;
    printf("System state:\n");
    while (atomCount < nToPrint)
    {
        for (int iAtom=0;iAtom<simHostSoa.grid.nAtoms[iBox];iAtom++)
	{

            i = iBox*MAXATOMS + iAtom;
#if (defined (__APPLE__) || defined(MACOSX)) && (APPLE_OCL_10)
            printf("%02d, %02d, "
                    "X=(%+020.12e %+020.12e %+020.12e) 1 "
                    "P=(%+020.12e %+020.12e %+020.12e) "
                    "F=(%+020.12e %+020.12e %+020.12e)\n",
                    iBox, iAtom, 
                    simHostSoa.r[i][0],simHostSoa.r[i][1],simHostSoa.r[i][2],
                    simHostSoa.p[i][0],simHostSoa.p[i][1],simHostSoa.p[i][2],
                    simHostSoa.f[i][0],simHostSoa.f[i][1],simHostSoa.f[i][2]);
#else
            printf("%02d, %02d, "
                    "X=(%+020.12e %+020.12e %+020.12e) 1 "
                    "P=(%+020.12e %+020.12e %+020.12e) "
                    "F=(%+020.12e %+020.12e %+020.12e)\n",
                    iBox, iAtom, 
                    simHostSoa.r[i].s[0],simHostSoa.r[i].s[1],simHostSoa.r[i].s[2],
                    simHostSoa.p[i].s[0],simHostSoa.p[i].s[1],simHostSoa.p[i].s[2],
                    simHostSoa.f[i].s[0],simHostSoa.f[i].s[1],simHostSoa.f[i].s[2]);
#endif
            atomCount ++;
        }
        iBox ++;
    }

}

void createDevGrid(DevGridSoa *gridDev, int nBoxesRealSize, int nBoxesNeighborSize, int nBoxesIntSize) 
{
    /** Create the device buffers to hold the grid data:
      box centers, neighbor lists info and bounds
     **/

    oclCreateReadWriteBuffer(&gridDev->rBox.x, nBoxesRealSize);
    oclCreateReadWriteBuffer(&gridDev->rBox.y, nBoxesRealSize);
    oclCreateReadWriteBuffer(&gridDev->rBox.z, nBoxesRealSize);

    oclCreateReadWriteBuffer(&gridDev->neighborList, nBoxesNeighborSize);
    oclCreateReadWriteBuffer(&gridDev->nNeighbors, nBoxesIntSize);
    oclCreateReadWriteBuffer(&gridDev->nAtoms, nBoxesIntSize);

    oclCreateReadWriteBuffer(&gridDev->bounds, sizeof(cl_real)*3);
}

void createDevGridAoS(DevGridAos *gridDev, int nBoxesRealSize, int nBoxesNeighborSize, int nBoxesIntSize) 
{
    /** Create the device buffers to hold the grid data:
      box centers, neighbor lists info and bounds
     **/

    oclCreateReadWriteBuffer(&gridDev->rBox, nBoxesRealSize*r3);

    oclCreateReadWriteBuffer(&gridDev->neighborList, nBoxesNeighborSize);
    oclCreateReadWriteBuffer(&gridDev->nNeighbors, nBoxesIntSize);
    oclCreateReadWriteBuffer(&gridDev->nAtoms, nBoxesIntSize);

    oclCreateReadWriteBuffer(&gridDev->bounds, sizeof(cl_real4));
}

void createDevVec(DevVec *aDev, int arraySize) 
{
    oclCreateReadWriteBuffer(&aDev->x, arraySize);
    oclCreateReadWriteBuffer(&aDev->y, arraySize);
    oclCreateReadWriteBuffer(&aDev->z, arraySize);
}

void getVecAoS(cl_mem aDev,
        cl_real4* aHost,
        int arraySize)
{
    oclCopyToHost(aDev, aHost, r3*arraySize);
}

void getVec(DevVec aDev,
        HostVec aHost,
        int arraySize)
{
    oclCopyToHost(aDev.x, aHost.x, arraySize);
    oclCopyToHost(aDev.y, aHost.y, arraySize);
    oclCopyToHost(aDev.z, aHost.z, arraySize);
}

void getVector(cl_mem axDev, cl_mem ayDev, cl_mem azDev,
        cl_real* axHost, cl_real* ayHost, cl_real* azHost,
        int arraySize)
{
    oclCopyToHost(axDev, axHost, arraySize);
    oclCopyToHost(ayDev, ayHost, arraySize);
    oclCopyToHost(azDev, azHost, arraySize);
}

void putVector(
        cl_real* axHost, cl_real* ayHost, cl_real* azHost,
        cl_mem axDev, cl_mem ayDev, cl_mem azDev,
        int arraySize)
{
    oclCopyToDevice(axHost, axDev, arraySize);
    oclCopyToDevice(ayHost, ayDev, arraySize);
    oclCopyToDevice(azHost, azDev, arraySize);
}

void putVec(
        HostVec aHost,
        DevVec aDev,
        int arraySize)
{
    oclCopyToDevice(aHost.x, aDev.x, arraySize);
    oclCopyToDevice(aHost.y, aDev.y, arraySize);
    oclCopyToDevice(aHost.z, aDev.z, arraySize);
}

void putGrid(HostGridSoa gridHost, DevGridSoa gridDev)
{
    putVec(gridHost.rBox, gridDev.rBox, gridHost.nBoxesRealSize);

    oclCopyToDevice(gridHost.nNeighbors, gridDev.nNeighbors, gridHost.nBoxesIntSize);
    oclCopyToDevice(gridHost.nAtoms, gridDev.nAtoms, gridHost.nBoxesIntSize);
    oclCopyToDevice(gridHost.neighborList, gridDev.neighborList, gridHost.nBoxesNeighborSize);

    oclCopyToDevice(gridHost.bounds, gridDev.bounds, sizeof(cl_real)*3);
}

void putGridAoS(HostGridAos gridHost, DevGridAos gridDev)
{
    oclCopyToDevice(gridHost.rBox, gridDev.rBox, gridHost.nBoxesRealSize*r3);

    oclCopyToDevice(gridHost.nNeighbors, gridDev.nNeighbors, gridHost.nBoxesIntSize);
    oclCopyToDevice(gridHost.nAtoms, gridDev.nAtoms, gridHost.nBoxesIntSize);
    oclCopyToDevice(gridHost.neighborList, gridDev.neighborList, gridHost.nBoxesNeighborSize);

    oclCopyToDevice(gridHost.bounds, gridDev.bounds, sizeof(cl_real4));
}

void putEamPot(HostEamPot eamPotH, DevEamPot eamPotD)
{
    oclCopyToDevice(eamPotH.rho, eamPotD.rho, eamPotH.rhoPotSize);
    oclCopyToDevice(eamPotH.phi, eamPotD.phi, eamPotH.phiPotSize);
    oclCopyToDevice(eamPotH.F, eamPotD.F, eamPotH.fPotSize);

    oclCopyToDevice(eamPotH.nValues, eamPotD.nValues, sizeof(cl_int)*3);
}

void putEamCh(HostEamCh eamChH, DevEamCh eamChD)
{
    oclCopyToDevice(eamChH.rho, eamChD.rho, eamChH.rhoChebSize);
    oclCopyToDevice(eamChH.phi, eamChD.phi, eamChH.phiChebSize);
    oclCopyToDevice(eamChH.F, eamChD.F, eamChH.fChebSize);

    oclCopyToDevice(eamChH.nValues, eamChD.nValues, sizeof(cl_int)*3);
}


void oclRunKernel(cl_kernel kernel, cl_event *event, size_t* nGlobal, size_t* nLocal)
{
    int err = clEnqueueNDRangeKernel(commandq, kernel, 2, NULL, nGlobal, nLocal, 0, NULL, event);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to enqueue kernel! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    }
    clWaitForEvents(1, event);
}

void printSim(SimFlat *s,FILE *fp) 
{
    /** Print the base simulation data
      Note this is in a slightly different order than the OpenCL code returns
     **/

    for (int i=0; i<s->nBoxes; i++)
    {
        int j;
        int *id;

        for (int iOff=i*MAXATOMS,j=0; j<s->nAtoms[i]; j++,iOff++)
	{
            if ( s->id[iOff] < 10)
	    {
                fprintf(fp,
                        "%02d %02d "
                        "X=(%+020.12e %+020.12e %+020.12e) 1 "
                        "P=(%+020.12e %+020.12e %+020.12e) "
                        "F=(%+020.12e %+020.12e %+020.12e)\n",
                        i,
                        s->id[iOff]+1,
                        s->r[iOff][0],s->r[iOff][1],s->r[iOff][2],
                        s->p[iOff][0],s->p[iOff][1],s->p[iOff][2],
                        s->f[iOff][0],s->f[iOff][1],s->f[iOff][2]
                       );
            }
        }
    }
    return;
}

void setEAMArgs( cl_kernel *forceKernels, DevSimSoa simDevSoa)
{ 
    /** Set the kernel arguments for the three EAM force computation kernels
     **/

#if (USE_CHEBY) 
    DevEamCh localPot;
    localPot = simDevSoa.eamCh;
#else 
    DevEamPot localPot;
    localPot = simDevSoa.eamPot;
#endif

    printf("Setting EAM kernel arguments\n");
    printf("Kernel 1...");
    fflush(stdout);
    // set kernel arguments for EAM_Force_1
    int err = 0;
    int nArg = 0;
    err  = clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.r.x);
    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.r.y);
    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.r.z);

    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.f.x);
    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.f.y);
    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.f.z);

    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.e);
    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.rho);
    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.fi);

    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.grid.rBox.x);
    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.grid.rBox.y);
    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.grid.rBox.z);

    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.grid.bounds);
    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.grid.neighborList);
    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.grid.nNeighbors);
    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.grid.nAtoms);

    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &localPot.nValues);

    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &localPot.phi);
    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &localPot.rho);
    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &localPot.F);

    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_real), &localPot.cutoff);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set EAM_Force_1 arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        fflush(stdout);
        exit(1);
    }
    else
    {
        printf("done\n");
        fflush(stdout);
    }

    // set kernel arguments for EAM_Force_2
    printf("Kernel 2...");
    fflush(stdout);
    err = 0;
    nArg = 0;
    err  = clSetKernelArg(forceKernels[1], nArg++, sizeof(cl_mem), &simDevSoa.fi);
    err |= clSetKernelArg(forceKernels[1], nArg++, sizeof(cl_mem), &simDevSoa.e);
    err |= clSetKernelArg(forceKernels[1], nArg++, sizeof(cl_mem), &simDevSoa.rho);

    err |= clSetKernelArg(forceKernels[1], nArg++, sizeof(cl_mem), &simDevSoa.grid.nAtoms);

    err |= clSetKernelArg(forceKernels[1], nArg++, sizeof(cl_mem), &localPot.F);
    err |= clSetKernelArg(forceKernels[1], nArg++, sizeof(cl_mem), &localPot.nValues);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set EAM_Force_2 arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        fflush(stdout);
        exit(1);
    }
    else
    {
        printf("done\n");
        fflush(stdout);
    }

    // set kernel arguments for EAM_Force_3
    printf("Kernel 3...");
    err = 0;
    nArg = 0;
    err  = clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.r.x);
    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.r.y);
    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.r.z);

    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.f.x);
    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.f.y);
    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.f.z);

    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.fi);

    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.grid.rBox.x);
    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.grid.rBox.y);
    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.grid.rBox.z);

    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.grid.bounds);
    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.grid.neighborList);
    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.grid.nNeighbors);
    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.grid.nAtoms);

    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &localPot.nValues);
    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &localPot.rho);
    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_real), &localPot.cutoff);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set EAM_Force_3 arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    }
    else
    {
        printf("done\n");
    }

}
void setEAMArgsAoS(
        cl_kernel *forceKernels,
        DevSimAos simDevSoa)
{ 
    /** Set the kernel arguments for the three EAM force computation kernels
     **/

    int err, nArg;

    printf("Setting EAM kernel arguments\n");
    printf("Kernel 1...");
    fflush(stdout);
    // set kernel arguments for EAM_Force_1
    err = 0;
    nArg = 0;
    err  = clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.r);

    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.f);

    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.e);
    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.rho);
    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.fi);

    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.grid.rBox);

    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.grid.bounds);
    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.grid.neighborList);
    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.grid.nNeighbors);
    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.grid.nAtoms);

    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.eamPot.nValues);

    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.eamPot.phi);
    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.eamPot.rho);
    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_mem), &simDevSoa.eamPot.F);

    err |= clSetKernelArg(forceKernels[0], nArg++, sizeof(cl_real), &simDevSoa.eamPot.cutoff);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set EAM_Force_1 arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
        fflush(stdout);
    }
    else
    {
        printf("done\n");
        fflush(stdout);
    }

    // set kernel arguments for EAM_Force_2
    printf("Kernel 2...");
    err = 0;
    nArg = 0;
    err  = clSetKernelArg(forceKernels[1], nArg++, sizeof(cl_mem), &simDevSoa.fi);
    err |= clSetKernelArg(forceKernels[1], nArg++, sizeof(cl_mem), &simDevSoa.e);
    err |= clSetKernelArg(forceKernels[1], nArg++, sizeof(cl_mem), &simDevSoa.rho);

    err |= clSetKernelArg(forceKernels[1], nArg++, sizeof(cl_mem), &simDevSoa.grid.nAtoms);

    err |= clSetKernelArg(forceKernels[1], nArg++, sizeof(cl_mem), &simDevSoa.eamPot.F);
    err |= clSetKernelArg(forceKernels[1], nArg++, sizeof(cl_mem), &simDevSoa.eamPot.nValues);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set EAM_Force_2 arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
        fflush(stdout);
    }
    else
    {
        printf("done\n");
        fflush(stdout);
    }

    // set kernel arguments for EAM_Force_3
    printf("Kernel 3...");
    err = 0;
    nArg = 0;
    // field arrays
    err  = clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.r);
    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.f);
    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.fi);
    // grid data
    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.grid.rBox);
    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.grid.bounds);
    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.grid.neighborList);
    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.grid.nNeighbors);
    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.grid.nAtoms);
    //potential data
    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.eamPot.nValues);
    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_mem), &simDevSoa.eamPot.rho);
    err |= clSetKernelArg(forceKernels[2], nArg++, sizeof(cl_real), &simDevSoa.eamPot.cutoff);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set EAM_Force_3 arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        fflush(stdout);
        exit(1);
    }
    else
    {
        printf("done\n");
        fflush(stdout);
    }

}

void setLJArgs(cl_kernel ljForce,
        DevSimSoa simDevSoa)
{
    /** Set the kernel arguments for the LJ force kernel **/

    printf("Setting LJ kernel arguments\n");
    // set kernel arguments for ljForce
    int err = 0;
    int nArg = 0;
    // field arrays
    err  = clSetKernelArg(ljForce, nArg++, sizeof(cl_mem), &simDevSoa.r.x);
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_mem), &simDevSoa.r.y);
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_mem), &simDevSoa.r.z);
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_mem), &simDevSoa.f.x);
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_mem), &simDevSoa.f.y);
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_mem), &simDevSoa.f.z);
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_mem), &simDevSoa.e);
    // grid data
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_mem), &simDevSoa.grid.rBox.x);
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_mem), &simDevSoa.grid.rBox.y);
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_mem), &simDevSoa.grid.rBox.z);
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_mem), &simDevSoa.grid.bounds);
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_mem), &simDevSoa.grid.neighborList);
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_mem), &simDevSoa.grid.nNeighbors);
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_mem), &simDevSoa.grid.nAtoms);
    // potential data
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_real), &simDevSoa.ljPot.sigma);
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_real), &simDevSoa.ljPot.epsilon);
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_real), &simDevSoa.ljPot.cutoff);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set ljForce arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    }
    else
    {
        printf("ljForce arguments set\n");
    }
}

void setLJArgsAoS(cl_kernel ljForce,
        DevSimAos simDevSoa)
{
    /** Set the kernel arguments for the LJ force kernel **/

    printf("Setting LJ kernel arguments\n");
    // set kernel arguments for ljForce
    int err = 0;
    int nArg = 0;
    // field arrays
    err  = clSetKernelArg(ljForce, nArg++, sizeof(cl_mem), &simDevSoa.r);
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_mem), &simDevSoa.f);
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_mem), &simDevSoa.e);
    // grid data
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_mem), &simDevSoa.grid.rBox);
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_mem), &simDevSoa.grid.bounds);
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_mem), &simDevSoa.grid.neighborList);
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_mem), &simDevSoa.grid.nNeighbors);
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_mem), &simDevSoa.grid.nAtoms);
    // potential data
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_real), &simDevSoa.ljPot.sigma);
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_real), &simDevSoa.ljPot.epsilon);
    err |= clSetKernelArg(ljForce, nArg++, sizeof(cl_real), &simDevSoa.ljPot.cutoff);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set ljForceAos arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    }
    else
    {
        printf("ljForceAos arguments set\n");
    }
}

void setAVArgs(cl_kernel advanceVelocity,
        DevSimSoa simDevSoa,
        cl_real dt)
{
    /** Set the arguments for the advanceVelocity kernel.
      Because of the Verlet timestepping scheme we keep the timestep as a separate argument
     **/

    int err = 0;
    int nArg = 0;
    // field arrays
    err  = clSetKernelArg(advanceVelocity, nArg++, sizeof(cl_mem), &simDevSoa.p.x);
    err |= clSetKernelArg(advanceVelocity, nArg++, sizeof(cl_mem), &simDevSoa.p.y);
    err |= clSetKernelArg(advanceVelocity, nArg++, sizeof(cl_mem), &simDevSoa.p.z);
    err |= clSetKernelArg(advanceVelocity, nArg++, sizeof(cl_mem), &simDevSoa.f.x);
    err |= clSetKernelArg(advanceVelocity, nArg++, sizeof(cl_mem), &simDevSoa.f.y);
    err |= clSetKernelArg(advanceVelocity, nArg++, sizeof(cl_mem), &simDevSoa.f.z);
    // grid data
    err |= clSetKernelArg(advanceVelocity, nArg++, sizeof(cl_mem), &simDevSoa.grid.nAtoms);
    // timestep
    err |= clSetKernelArg(advanceVelocity, nArg++, sizeof(cl_real), &dt);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set advanceVelocity arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    }
    else
    {
        printf("advanceVelocity arguments set\n");
        printf("dt = %e\n", dt);
    }
}

void setAPArgs(cl_kernel advancePosition,
        DevSimSoa simDevSoa,
        cl_real dt)
{
    /** Set the arguments for the advancePosition kernel.
      Because of the Verlet timestepping scheme we keep the timestep as a separate argument
     **/

    int err = 0;
    int nArg = 0;
    // field arrays
    err  = clSetKernelArg(advancePosition, nArg++, sizeof(cl_mem), &simDevSoa.p.x);
    err |= clSetKernelArg(advancePosition, nArg++, sizeof(cl_mem), &simDevSoa.p.y);
    err |= clSetKernelArg(advancePosition, nArg++, sizeof(cl_mem), &simDevSoa.p.z);
    err |= clSetKernelArg(advancePosition, nArg++, sizeof(cl_mem), &simDevSoa.r.x);
    err |= clSetKernelArg(advancePosition, nArg++, sizeof(cl_mem), &simDevSoa.r.y);
    err |= clSetKernelArg(advancePosition, nArg++, sizeof(cl_mem), &simDevSoa.r.z);
    //err |= clSetKernelArg(advancePosition, nArg++, sizeof(cl_mem), &simDevSoa.m);
    // grid data
    err |= clSetKernelArg(advancePosition, nArg++, sizeof(cl_mem), &simDevSoa.grid.nAtoms);
    // timestep
    err |= clSetKernelArg(advancePosition, nArg++, sizeof(cl_real), &simDevSoa.rMass);
    err |= clSetKernelArg(advancePosition, nArg++, sizeof(cl_real), &dt);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set advancePosition arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    }
    else
    {
        printf("advancePosition arguments set\n");
        printf("dt = %e\n", dt);
    }
}

void setAVArgsAoS(cl_kernel advanceVelocity,
        DevSimAos simDevSoa,
        cl_real dt)
{
    /** Set the arguments for the advanceVelocity kernel.
      Because of the Verlet timestepping scheme we keep the timestep as a separate argument
     **/

    int err = 0;
    int nArg = 0;
    // field arrays
    err  = clSetKernelArg(advanceVelocity, nArg++, sizeof(cl_mem), &simDevSoa.p);
    err |= clSetKernelArg(advanceVelocity, nArg++, sizeof(cl_mem), &simDevSoa.f);
    err |= clSetKernelArg(advanceVelocity, nArg++, sizeof(cl_mem), &simDevSoa.grid.nAtoms);
    err |= clSetKernelArg(advanceVelocity, nArg++, sizeof(cl_real), &dt);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set advanceVelocityAos arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    }
    else
    {
        printf("advanceVelocityAos arguments set\n");
        printf("dt = %e\n", dt);
    }
}

void setAPArgsAoS(cl_kernel advancePosition,
        DevSimAos simDevSoa,
        cl_real dt)
{
    /** Set the arguments for the advancePosition kernel.
      Because of the Verlet timestepping scheme we keep the timestep as a separate argument
     **/

    int err = 0;
    int nArg = 0;
    // field arrays
    err  = clSetKernelArg(advancePosition, nArg++, sizeof(cl_mem), &simDevSoa.p);
    err |= clSetKernelArg(advancePosition, nArg++, sizeof(cl_mem), &simDevSoa.r);
    //err |= clSetKernelArg(advancePosition, nArg++, sizeof(cl_mem), &simDevSoa.m);
    err |= clSetKernelArg(advancePosition, nArg++, sizeof(cl_mem), &simDevSoa.grid.nAtoms);
    err |= clSetKernelArg(advancePosition, nArg++, sizeof(cl_real), &simDevSoa.rMass);
    err |= clSetKernelArg(advancePosition, nArg++, sizeof(cl_real), &dt);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set advancePositionAos arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    }
    else
    {
        printf("advancePositionAos arguments set\n");
        printf("dt = %e\n", dt);
    }
}

void getElapsedTime(cl_event event, cl_real* elapsed_time, cl_real* enqueuedTime)
{
    /** Helper routine to return the start-to-finish time (elapsed_time) 
      and the time from enqueueing to finish (enqueuedTime)
     **/

    cl_ulong t_start, t_end, t_enqueue;
    int err;
    size_t param_size;


    err = clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &t_end, &param_size);
    if (err != CL_SUCCESS)
    {
        printf("Error: %s\n", print_cl_errstring(err));
        printf("t_end = %llu\n", t_end);
        //exit(1);
    }

    err = clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &t_start, &param_size);
    if (err != CL_SUCCESS)
    {
        printf("Error: %s\n", print_cl_errstring(err));
        printf("t_start = %llu\n", t_start);
        //exit(1);
    }

    err = clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &t_enqueue, &param_size);

    *elapsed_time = (t_end - t_start)*1.0e-9 ;
    *enqueuedTime = (t_end - t_enqueue)*1.0e-9 ;
}

void MakeHostArrays(
        HostVec *r_H,
        HostVec *v_H,
        HostVec *f_H,
        HostGridSoa *gridHost,
        int arraySize,
        int nBoxesRealSize,
        int nBoxesNeighborSize,
        int nBoxesIntSize
        )
{
    /** Allocate all the host arrays needed for the base simulation data **/

    // host memory
    // location
    r_H->x = malloc(arraySize);
    r_H->y = malloc(arraySize);
    r_H->z = malloc(arraySize);

    // momenta
    v_H->x = malloc(arraySize);
    v_H->y = malloc(arraySize);
    v_H->z = malloc(arraySize);

    // forces
    f_H->x = malloc(arraySize);
    f_H->y = malloc(arraySize);
    f_H->z = malloc(arraySize);

    // box locations
    gridHost->rBox.x = malloc(nBoxesRealSize);
    gridHost->rBox.y = malloc(nBoxesRealSize);
    gridHost->rBox.z = malloc(nBoxesRealSize);

    gridHost->neighborList = malloc(nBoxesNeighborSize);
    gridHost->nNeighbors = malloc(nBoxesIntSize);
    gridHost->nAtoms = malloc(nBoxesIntSize);
    gridHost->bounds = malloc(3*sizeof(cl_real));

    gridHost->nBoxesIntSize = nBoxesIntSize;
    gridHost->nBoxesNeighborSize = nBoxesNeighborSize;
    gridHost->nBoxesRealSize = nBoxesRealSize;

}

void getPrintState(DevSimSoa simDevSoa, HostSimSoa simHostSoa) 
{
    getVec(simDevSoa.r, simHostSoa.r, simHostSoa.arraySize);
    getVec(simDevSoa.p, simHostSoa.p, simHostSoa.arraySize);
    getVec(simDevSoa.f, simHostSoa.f, simHostSoa.arraySize);
    printState(simHostSoa, 2);
}

void getPrintStateAoS(DevSimAos simDevSoa, HostSimAos simHostSoa) 
{
    getVecAoS(simDevSoa.r, simHostSoa.r, simHostSoa.arraySize);
    getVecAoS(simDevSoa.p, simHostSoa.p, simHostSoa.arraySize);
    getVecAoS(simDevSoa.f, simHostSoa.f, simHostSoa.arraySize);
    printStateAoS(simHostSoa, 2);
}

void initHostEAM(HostEamPot *eamPotH, SimFlat *sim) 
{
    /** Allocate and initialize all the EAM potential data needed **/

    int i;
    int n_v_rho;
    int n_v_phi;
    int n_v_F;

    // assign eam potential values
    printf("Using eam potential\n");
    EamPotential *new_pot;
    new_pot = (EamPotential*) sim->pot;
    eamPotH->cutoff  = new_pot->cutoff;
    printf("Cutoff = %e\n", eamPotH->cutoff);

    n_v_rho = new_pot->rho->n;
    eamPotH->rhoPotSize = (6 + new_pot->rho->n)*sizeof(cl_real);
    printf("rho potential size = %d\n", eamPotH->rhoPotSize);

    n_v_phi = new_pot->phi->n;
    eamPotH->phiPotSize = (6 + new_pot->phi->n)*sizeof(cl_real);
    printf("phi potential size = %d\n", eamPotH->phiPotSize);

    n_v_F = new_pot->f->n;
    eamPotH->fPotSize = (6 + new_pot->f->n)*sizeof(cl_real);
    printf("F potential size = %d\n", eamPotH->fPotSize);

    eamPotH->rho = malloc(eamPotH->rhoPotSize);
    eamPotH->phi = malloc(eamPotH->phiPotSize);
    eamPotH->F = malloc(eamPotH->fPotSize);
    eamPotH->nValues = malloc(3*sizeof(cl_int));

    // Note the EAM array has 3 extra values to account for over/under flow
    // We also add another 3 values to store x0, xn, invDx
    eamPotH->rho[n_v_rho+3] = new_pot->rho->x0;
    eamPotH->rho[n_v_rho+4] = new_pot->rho->xn;
    eamPotH->rho[n_v_rho+5] = new_pot->rho->invDx;

    for (i=0;i<n_v_rho+3;i++)
    {
        eamPotH->rho[i] = new_pot->rho->values[i-1];
    }

    eamPotH->phi[n_v_phi+3] = new_pot->phi->x0;
    eamPotH->phi[n_v_phi+4] = new_pot->phi->xn;
    eamPotH->phi[n_v_phi+5] = new_pot->phi->invDx;

    for (i=0;i<n_v_phi+3;i++)
    {
        eamPotH->phi[i] = new_pot->phi->values[i-1];
    }

    eamPotH->F[n_v_F+3] = new_pot->f->x0;
    eamPotH->F[n_v_F+4] = new_pot->f->xn;
    eamPotH->F[n_v_F+5] = new_pot->f->invDx;

    for (i=0;i<n_v_F+3;i++)
    {
        eamPotH->F[i] = new_pot->f->values[i-1];
    }

    eamPotH->nValues[0] = n_v_phi;
    eamPotH->nValues[1] = n_v_rho;
    eamPotH->nValues[2] = n_v_F;
}

void initHostEAMCh(HostEamCh *eamChH, SimFlat *sim) 
{
    /** Allocate and initialize all the EAM potential data needed **/

    int i;
    int n_v_rho;
    int n_v_phi;
    int n_v_F;

    // assign eam potential values
    printf("Using eam potential\n");
    EamPotential *new_pot;
    new_pot = (EamPotential*) sim->pot;
    EamCheby *chPot = sim->chPot;
    eamChH->cutoff  = new_pot->cutoff;
    printf("Cutoff = %e\n", eamChH->cutoff);

    /** The Chebyshev arrays are stored as:
      The interpolants for the potential
      The interpolants for the derivative
      The 3 values r0, rN and a dummy float to match the table 
     * **/

    n_v_rho = chPot->rho->n;
    eamChH->rhoChebSize = (3 + 2*chPot->rho->n)*sizeof(cl_real);
    printf("rho potential size = %d\n", eamChH->rhoChebSize);

    n_v_phi = chPot->phi->n;
    eamChH->phiChebSize = (3 + 2*chPot->phi->n)*sizeof(cl_real);
    printf("phi potential size = %d\n", eamChH->phiChebSize);

    n_v_F = chPot->f->n;
    eamChH->fChebSize = (3 + 2*chPot->f->n)*sizeof(cl_real);
    printf("F potential size = %d\n", eamChH->fChebSize);

    eamChH->rho = malloc(eamChH->rhoChebSize);
    eamChH->phi = malloc(eamChH->phiChebSize);
    eamChH->F = malloc(eamChH->fChebSize);
    eamChH->nValues = malloc(3*sizeof(cl_int));

    eamChH->rho[2*n_v_rho+0] = chPot->rho->x0;
    eamChH->rho[2*n_v_rho+1] = chPot->rho->xn;
    eamChH->rho[2*n_v_rho+2] = chPot->rho->invDx;

    for (i=0;i<n_v_rho;i++)
    {
        eamChH->rho[i] = chPot->rho->values[i];
        eamChH->rho[n_v_rho+i] = chPot->drho->values[i];
    }

    eamChH->phi[2*n_v_phi+0] = chPot->phi->x0;
    eamChH->phi[2*n_v_phi+1] = chPot->phi->xn;
    eamChH->phi[2*n_v_phi+2] = chPot->phi->invDx;

    for (i=0;i<n_v_phi;i++)
    {
        eamChH->phi[i] = chPot->phi->values[i];
        eamChH->phi[n_v_phi+i] = chPot->dphi->values[i];
    }

    eamChH->F[2*n_v_F+0] = chPot->f->x0;
    eamChH->F[2*n_v_F+1] = chPot->f->xn;
    eamChH->F[2*n_v_F+2] = chPot->f->invDx;

    for (i=0;i<n_v_F;i++)
    {
        eamChH->F[i] = chPot->f->values[i];
        eamChH->F[n_v_F+i] = chPot->df->values[i];
    }

    eamChH->nValues[0] = n_v_phi;
    eamChH->nValues[1] = n_v_rho;
    eamChH->nValues[2] = n_v_F;
}

void initHostLJ(HostLjPot *ljPot_H, SimFlat *sim)
{
    /** Allocate and initialize all the LJ potential data needed **/

    LjPotential *new_pot;
    new_pot = (LjPotential*) sim->pot;
    ljPot_H->sigma  = new_pot->sigma;
    ljPot_H->epsilon  = new_pot->epsilon;
    ljPot_H->cutoff  = new_pot->cutoff;

    printf("Using lj potential\n");
    printf("Sigma = %e\n", ljPot_H->sigma);
    printf("Epsilon = %e\n", ljPot_H->epsilon);
    printf("Cutoff = %e\n", ljPot_H->cutoff);

}

void initHostSim (HostSimSoa *simHostSoa, SimFlat *sim)
{
    /** Allocate and initialize all the host-side simulation data needed, 
      including the appropriate potential data 
     **/

    int iBox, iAtom, iOff;

    simHostSoa->nCells = sim->nBoxes;
    simHostSoa->cfac = bohr_per_atu_to_A_per_s;
    simHostSoa->rMass = (cl_real)((double)1822.83);
    simHostSoa->dt = 1.0e-15*simHostSoa->cfac;
    printf("dt = %e\n", simHostSoa->dt);

    int nMaxInBox = 0;
    int boxWithMax = 0;
    for (iBox=0;iBox<simHostSoa->nCells;iBox++)
    {
        if (sim->nAtoms[iBox] > nMaxInBox)
	{
            nMaxInBox = sim->nAtoms[iBox];
            boxWithMax = iBox;
        }
    }
    printf("Max atom count %d in box %d\n", nMaxInBox, boxWithMax);

    simHostSoa->arraySize = MAXATOMS*sim->nBoxes*sizeof(cl_real);
    simHostSoa->nBoxesNeighborSize = sim->nBoxes*NUMNEIGHBORS*sizeof(cl_int);
    simHostSoa->nBoxesIntSize = sim->nBoxes*sizeof(cl_int);
    simHostSoa->nBoxesRealSize = sim->nBoxes*sizeof(cl_real);
    simHostSoa->nTot = sim->nTot;

    // location
    simHostSoa->r.x = malloc(simHostSoa->arraySize);
    simHostSoa->r.y = malloc(simHostSoa->arraySize);
    simHostSoa->r.z = malloc(simHostSoa->arraySize);

    // momenta
    simHostSoa->p.x = malloc(simHostSoa->arraySize);
    simHostSoa->p.y = malloc(simHostSoa->arraySize);
    simHostSoa->p.z = malloc(simHostSoa->arraySize);

    // forces
    simHostSoa->f.x = malloc(simHostSoa->arraySize);
    simHostSoa->f.y = malloc(simHostSoa->arraySize);
    simHostSoa->f.z = malloc(simHostSoa->arraySize);

    // box locations
    simHostSoa->grid.rBox.x = malloc(simHostSoa->nBoxesRealSize);
    simHostSoa->grid.rBox.y = malloc(simHostSoa->nBoxesRealSize);
    simHostSoa->grid.rBox.z = malloc(simHostSoa->nBoxesRealSize);

    simHostSoa->grid.neighborList = malloc(simHostSoa->nBoxesNeighborSize);
    simHostSoa->grid.nNeighbors = malloc(simHostSoa->nBoxesIntSize);
    simHostSoa->grid.nAtoms = malloc(simHostSoa->nBoxesIntSize);
    simHostSoa->grid.bounds = malloc(3*sizeof(cl_real));

    simHostSoa->grid.nBoxesIntSize = simHostSoa->nBoxesIntSize;
    simHostSoa->grid.nBoxesNeighborSize = simHostSoa->nBoxesNeighborSize;
    simHostSoa->grid.nBoxesRealSize = simHostSoa->nBoxesRealSize;

    // mass, energy
    simHostSoa->m = malloc(simHostSoa->arraySize);
    simHostSoa->e = malloc(simHostSoa->arraySize);

    if(simHostSoa->eamFlag)
    {
        simHostSoa->fi = malloc(simHostSoa->arraySize);
        simHostSoa->rho = malloc(simHostSoa->arraySize);

        initHostEAM(&simHostSoa->eamPot, sim);
        initHostEAMCh(&simHostSoa->eamCh, sim);
    }
    else
    {
        initHostLJ(&simHostSoa->ljPot, sim);
    }

    for (iBox=0;iBox<sim->nBoxes;iBox++)
    {

        int* nbrBoxes;
        nbrBoxes = getNeighborBoxes(sim,iBox);

        simHostSoa->grid.nAtoms[iBox] = sim->nAtoms[iBox];

        simHostSoa->grid.rBox.x[iBox] = sim->dCenter[iBox][0];
        simHostSoa->grid.rBox.y[iBox] = sim->dCenter[iBox][1];
        simHostSoa->grid.rBox.z[iBox] = sim->dCenter[iBox][2];

        int j;
        simHostSoa->grid.nNeighbors[iBox] = nbrBoxes[-1];
        for (j=0;j<simHostSoa->grid.nNeighbors[iBox];j++)
	{
            simHostSoa->grid.neighborList[NUMNEIGHBORS*iBox + j] = nbrBoxes[j];
        }

        for (iAtom=0;iAtom<sim->nAtoms[iBox];iAtom++)
	{

            iOff = iBox*MAXATOMS + iAtom;

            simHostSoa->r.x[iOff] = sim->r[iOff][0];
            simHostSoa->r.y[iOff] = sim->r[iOff][1];
            simHostSoa->r.z[iOff] = sim->r[iOff][2];

            simHostSoa->p.x[iOff] = sim->p[iOff][0];
            simHostSoa->p.y[iOff] = sim->p[iOff][1];
            simHostSoa->p.z[iOff] = sim->p[iOff][2];

            simHostSoa->f.x[iOff] = sim->f[iOff][0];
            simHostSoa->f.y[iOff] = sim->f[iOff][1];
            simHostSoa->f.z[iOff] = sim->f[iOff][2];

        }

        for (j=0;j<3;j++)
	{
            simHostSoa->grid.bounds[j] = sim->bounds[j];
        }
    }

#ifdef INTEROP_VIZ
    centerX = sim->bounds[0]/2.0f;
    centerY = sim->bounds[1]/2.0f;
    centerZ = sim->bounds[2]/2.0f;
    printf("Center: %f %f %f\n", centerX, centerY, centerZ);
#endif

}

void initHostSimAoS (HostSimAos *simHostSoa, SimFlat *sim)
{
    /** Allocate and initialize all the host-side simulation data needed, 
      including the appropriate potential data 
     **/

    int iBox, iAtom, iOff;

    simHostSoa->nCells = sim->nBoxes;
    simHostSoa->cfac = bohr_per_atu_to_A_per_s;
    simHostSoa->rMass = (cl_real)((double)1822.83);
    simHostSoa->dt = 1.0e-15*simHostSoa->cfac;
    printf("dt = %e\n", simHostSoa->dt);

    int nMaxInBox = 0;
    int boxWithMax = 0;
    for (iBox=0;iBox<simHostSoa->nCells;iBox++)
    {
        if (sim->nAtoms[iBox] > nMaxInBox)
	{
            nMaxInBox = sim->nAtoms[iBox];
            boxWithMax = iBox;
        }
    }
    printf("Max atom count %d in box %d\n", nMaxInBox, boxWithMax);

    simHostSoa->arraySize = MAXATOMS*sim->nBoxes*sizeof(cl_real);
    simHostSoa->nBoxesNeighborSize = sim->nBoxes*NUMNEIGHBORS*sizeof(cl_int);
    simHostSoa->nBoxesIntSize = sim->nBoxes*sizeof(cl_int);
    simHostSoa->nBoxesRealSize = sim->nBoxes*sizeof(cl_real);

    // location
    simHostSoa->r = malloc(simHostSoa->arraySize*r3);

    // momenta
    simHostSoa->p = malloc(simHostSoa->arraySize*r3);

    // forces
    simHostSoa->f = malloc(simHostSoa->arraySize*r3);

    // box locations
    simHostSoa->grid.rBox       = malloc(simHostSoa->nBoxesRealSize*r3);

    simHostSoa->grid.neighborList = malloc(simHostSoa->nBoxesNeighborSize);
    simHostSoa->grid.nNeighbors = malloc(simHostSoa->nBoxesIntSize);
    simHostSoa->grid.nAtoms     = malloc(simHostSoa->nBoxesIntSize);
    simHostSoa->grid.bounds      = malloc(sizeof(cl_real4));

    simHostSoa->grid.nBoxesIntSize    = simHostSoa->nBoxesIntSize;
    simHostSoa->grid.nBoxesNeighborSize   = simHostSoa->nBoxesNeighborSize;
    simHostSoa->grid.nBoxesRealSize    = simHostSoa->nBoxesRealSize;
    simHostSoa->nTot             = sim->nTot;

    // mass, energy
    simHostSoa->m = malloc(simHostSoa->arraySize);
    simHostSoa->e = malloc(simHostSoa->arraySize);

    if(simHostSoa->eamFlag)
    {
        simHostSoa->fi = malloc(simHostSoa->arraySize);
        simHostSoa->rho = malloc(simHostSoa->arraySize);

        initHostEAM(&simHostSoa->eamPot, sim);
        initHostEAMCh(&simHostSoa->eamCh, sim);
    }
    else
    {
        initHostLJ(&simHostSoa->ljPot, sim);
    }

    for (iBox=0;iBox<sim->nBoxes;iBox++)
    {

        int* nbrBoxes;
        nbrBoxes = getNeighborBoxes(sim,iBox);

        simHostSoa->grid.nAtoms[iBox] = sim->nAtoms[iBox];

#if (defined (__APPLE__) || defined(MACOSX)) && (APPLE_OCL_10)
        simHostSoa->grid.rBox[iBox][0] = sim->dCenter[iBox][0];
        simHostSoa->grid.rBox[iBox][1] = sim->dCenter[iBox][1];
        simHostSoa->grid.rBox[iBox][2] = sim->dCenter[iBox][2];
#else
        simHostSoa->grid.rBox[iBox].s[0] = sim->dCenter[iBox][0];
        simHostSoa->grid.rBox[iBox].s[1] = sim->dCenter[iBox][1];
        simHostSoa->grid.rBox[iBox].s[2] = sim->dCenter[iBox][2];
#endif

        int j;
        simHostSoa->grid.nNeighbors[iBox] = nbrBoxes[-1];
        for (j=0;j<simHostSoa->grid.nNeighbors[iBox];j++)
	{
            simHostSoa->grid.neighborList[NUMNEIGHBORS*iBox + j] = nbrBoxes[j];
        }

        for (iAtom=0;iAtom<sim->nAtoms[iBox];iAtom++)
	{

            iOff = iBox*MAXATOMS + iAtom;
#if (defined (__APPLE__) || defined(MACOSX)) && (APPLE_OCL_10)
            simHostSoa->r[iOff][0] = sim->r[iOff][0];
            simHostSoa->r[iOff][1] = sim->r[iOff][1];
            simHostSoa->r[iOff][2] = sim->r[iOff][2];

            simHostSoa->p[iOff][0] = sim->p[iOff][0];
            simHostSoa->p[iOff][1] = sim->p[iOff][1];
            simHostSoa->p[iOff][2] = sim->p[iOff][2];

            simHostSoa->f[iOff][0] = sim->f[iOff][0];
            simHostSoa->f[iOff][1] = sim->f[iOff][1];
            simHostSoa->f[iOff][2] = sim->f[iOff][2];
#else
            simHostSoa->r[iOff].s[0] = sim->r[iOff][0];
            simHostSoa->r[iOff].s[1] = sim->r[iOff][1];
            simHostSoa->r[iOff].s[2] = sim->r[iOff][2];

            simHostSoa->p[iOff].s[0] = sim->p[iOff][0];
            simHostSoa->p[iOff].s[1] = sim->p[iOff][1];
            simHostSoa->p[iOff].s[2] = sim->p[iOff][2];

            simHostSoa->f[iOff].s[0] = sim->f[iOff][0];
            simHostSoa->f[iOff].s[1] = sim->f[iOff][1];
            simHostSoa->f[iOff].s[2] = sim->f[iOff][2];
#endif
        }

        for (j=0;j<3;j++)
	{
            simHostSoa->grid.bounds[j] = sim->bounds[j];
        }
    }

}

void initDevSim(DevSimSoa *simDevSoa, HostSimSoa *simHostSoa)
{
    /** Allocate all the device-side arrays needed for the simulation **/

    // allocate memory buffer on device
    printf("Allocating device memory...");

    // positions
    createDevVec(&simDevSoa->r, simHostSoa->arraySize);
    createDevVec(&simDevSoa->p, simHostSoa->arraySize);
    createDevVec(&simDevSoa->f, simHostSoa->arraySize);

    // particle mass
    oclCreateReadWriteBuffer(&simDevSoa->m, simHostSoa->arraySize);

    // particle energy
    oclCreateReadWriteBuffer(&simDevSoa->e, simHostSoa->arraySize);

    createDevGrid(&simDevSoa->grid, simHostSoa->nBoxesRealSize, simHostSoa->nBoxesNeighborSize, simHostSoa->nBoxesIntSize);

    if (simHostSoa->eamFlag)
    {
        oclCreateReadWriteBuffer(&simDevSoa->fi, simHostSoa->arraySize);
        oclCreateReadWriteBuffer(&simDevSoa->rho, simHostSoa->arraySize);

        //************************************************************************
        // EAM table data
        oclCreateReadWriteBuffer(&simDevSoa->eamPot.rho, simHostSoa->eamPot.rhoPotSize);
        oclCreateReadWriteBuffer(&simDevSoa->eamPot.phi, simHostSoa->eamPot.phiPotSize);
        oclCreateReadWriteBuffer(&simDevSoa->eamPot.F, simHostSoa->eamPot.fPotSize);

        oclCreateReadWriteBuffer(&simDevSoa->eamPot.nValues, sizeof(cl_int)*3);

        // add this here to make passing arguments to kernels easier
        simDevSoa->eamPot.cutoff = simHostSoa->eamPot.cutoff;

        //************************************************************************
        // EAM Chebychev coefficient data
        oclCreateReadWriteBuffer(&simDevSoa->eamCh.rho, simHostSoa->eamCh.rhoChebSize);
        oclCreateReadWriteBuffer(&simDevSoa->eamCh.phi, simHostSoa->eamCh.phiChebSize);
        oclCreateReadWriteBuffer(&simDevSoa->eamCh.F, simHostSoa->eamCh.fChebSize);

        oclCreateReadWriteBuffer(&simDevSoa->eamCh.nValues, sizeof(cl_int)*3);

        // add this here to make passing arguments to kernels easier
        simDevSoa->eamCh.cutoff = simHostSoa->eamCh.cutoff;
    }
    else
    {
        simDevSoa->ljPot.cutoff = simHostSoa->ljPot.cutoff;
        simDevSoa->ljPot.sigma = simHostSoa->ljPot.sigma;
        simDevSoa->ljPot.epsilon = simHostSoa->ljPot.epsilon;
    }
    printf("device memory allocated\n");

}

void initDevSimAoS(DevSimAos *simDevSoa, HostSimAos *simHostSoa)
{
    /** Allocate all the device-side arrays needed for the simulation **/

    // allocate memory buffer on device
    printf("Allocating device memory (AoS)...");

    // positions, momenta, force
    oclCreateReadWriteBuffer(&simDevSoa->r, simHostSoa->arraySize*r3);
    oclCreateReadWriteBuffer(&simDevSoa->p, simHostSoa->arraySize*r3);
    oclCreateReadWriteBuffer(&simDevSoa->f, simHostSoa->arraySize*r3);

    // particle mass
    oclCreateReadWriteBuffer(&simDevSoa->m, simHostSoa->arraySize);

    // particle energy
    oclCreateReadWriteBuffer(&simDevSoa->e, simHostSoa->arraySize);

    createDevGridAoS(&simDevSoa->grid, simHostSoa->nBoxesRealSize, simHostSoa->nBoxesNeighborSize, simHostSoa->nBoxesIntSize);

    if (simHostSoa->eamFlag)
    {
        oclCreateReadWriteBuffer(&simDevSoa->fi, simHostSoa->arraySize);
        oclCreateReadWriteBuffer(&simDevSoa->rho, simHostSoa->arraySize);

        //************************************************************************
        // EAM table data
        oclCreateReadWriteBuffer(&simDevSoa->eamPot.rho, simHostSoa->eamPot.rhoPotSize);
        oclCreateReadWriteBuffer(&simDevSoa->eamPot.phi, simHostSoa->eamPot.phiPotSize);
        oclCreateReadWriteBuffer(&simDevSoa->eamPot.F, simHostSoa->eamPot.fPotSize);

        oclCreateReadWriteBuffer(&simDevSoa->eamPot.nValues, sizeof(cl_int)*3);

        // add this here to make passing arguments to kernels easier
        simDevSoa->eamPot.cutoff = simHostSoa->eamPot.cutoff;

        //************************************************************************
        // EAM Chebychev coefficient data
        oclCreateReadWriteBuffer(&simDevSoa->eamCh.rho, simHostSoa->eamCh.rhoChebSize);
        oclCreateReadWriteBuffer(&simDevSoa->eamCh.phi, simHostSoa->eamCh.phiChebSize);
        oclCreateReadWriteBuffer(&simDevSoa->eamCh.F, simHostSoa->eamCh.fChebSize);

        oclCreateReadWriteBuffer(&simDevSoa->eamCh.nValues, sizeof(cl_int)*3);

        // add this here to make passing arguments to kernels easier
        simDevSoa->eamCh.cutoff = simHostSoa->eamCh.cutoff;
    }
    else
    {
        simDevSoa->ljPot.cutoff = simHostSoa->ljPot.cutoff;
        simDevSoa->ljPot.sigma = simHostSoa->ljPot.sigma;
        simDevSoa->ljPot.epsilon = simHostSoa->ljPot.epsilon;
    }
    printf("device memory allocated\n");

}

void putSim(HostSimSoa simHostSoa, DevSimSoa simDevSoa)
{
    /** Copy all the host-side simulation data to the corresponding device arrays **/

    printf("Copying data to device...");

    // copy the input arrays to the device
    // positions

    putVec(simHostSoa.r, simDevSoa.r, simHostSoa.arraySize);
    putVec(simHostSoa.p, simDevSoa.p, simHostSoa.arraySize);
    putVec(simHostSoa.f, simDevSoa.f, simHostSoa.arraySize);

    // mass
    oclCopyToDevice(simHostSoa.m, simDevSoa.m, simHostSoa.arraySize);

    // simulation data
    putGrid(simHostSoa.grid, simDevSoa.grid);

    if(simHostSoa.eamFlag)
    {
        putEamPot(simHostSoa.eamPot, simDevSoa.eamPot);
        putEamCh(simHostSoa.eamCh, simDevSoa.eamCh);
    }
    printf("data copied\n");

}

void putSimAoS(HostSimAos simHostSoa, DevSimAos simDevSoa)
{
    /** Copy all the host-side simulation data to the corresponding device arrays **/

    printf("Copying data to device (AoS)...");
    fflush(stdout);

    // copy the input arrays to the device
    // positions

    oclCopyToDevice(simHostSoa.r, simDevSoa.r, simHostSoa.arraySize*r3);
    oclCopyToDevice(simHostSoa.p, simDevSoa.p, simHostSoa.arraySize*r3);
    oclCopyToDevice(simHostSoa.f, simDevSoa.f, simHostSoa.arraySize*r3);

    printf("positions...");
    fflush(stdout);

    // mass
    oclCopyToDevice(simHostSoa.m, simDevSoa.m, simHostSoa.arraySize);

    printf("mass...");
    fflush(stdout);

    // simulation data
    putGridAoS(simHostSoa.grid, simDevSoa.grid);

    printf("grid...");
    fflush(stdout);

    if(simHostSoa.eamFlag)
    {
        putEamPot(simHostSoa.eamPot, simDevSoa.eamPot);
    }
    printf("data copied\n");
    fflush(stdout);

}

void FreeSims(HostSimSoa simHostSoa, DevSimSoa simDevSoa)
{
    /** clean up all the host and device memory objects **/

    free(simHostSoa.r.x);
    free(simHostSoa.r.y);
    free(simHostSoa.r.z);

    free(simHostSoa.p.x);
    free(simHostSoa.p.y);
    free(simHostSoa.p.z);

    free(simHostSoa.f.x);
    free(simHostSoa.f.y);
    free(simHostSoa.f.z);

    free(simHostSoa.m);
    free(simHostSoa.e);

    clReleaseMemObject(simDevSoa.r.x);
    clReleaseMemObject(simDevSoa.r.y);
    clReleaseMemObject(simDevSoa.r.z);

    clReleaseMemObject(simDevSoa.p.x);
    clReleaseMemObject(simDevSoa.p.y);
    clReleaseMemObject(simDevSoa.p.z);

    clReleaseMemObject(simDevSoa.f.x);
    clReleaseMemObject(simDevSoa.f.y);
    clReleaseMemObject(simDevSoa.f.z);

    clReleaseMemObject(simDevSoa.m);
    clReleaseMemObject(simDevSoa.e);

    if(simHostSoa.eamFlag)
    {
        free(simHostSoa.rho);
        free(simHostSoa.fi);

        clReleaseMemObject(simDevSoa.rho);
        clReleaseMemObject(simDevSoa.fi);

    }
}

void FreeSimsAoS(HostSimAos simHostSoa, DevSimAos simDevSoa)
{
    /** clean up all the host and device memory objects **/

    free(simHostSoa.r);
    free(simHostSoa.p);
    free(simHostSoa.f);
    free(simHostSoa.m);
    free(simHostSoa.e);

    clReleaseMemObject(simDevSoa.r);
    clReleaseMemObject(simDevSoa.p);
    clReleaseMemObject(simDevSoa.f);
    clReleaseMemObject(simDevSoa.m);
    clReleaseMemObject(simDevSoa.e);

    if(simHostSoa.eamFlag)
    {
        free(simHostSoa.rho);
        free(simHostSoa.fi);

        clReleaseMemObject(simDevSoa.rho);
        clReleaseMemObject(simDevSoa.fi);

    }
}

void buildModulesSoa(cl_kernel *forceKernels, 
        cl_kernel *advancePosition, 
        cl_kernel *advanceVelocity,
        cl_kernel *Viz, 
        HostSimSoa simHostSoa, 
        size_t *nLocal, 
        size_t *nGlobal)
{
    /** Build the kernels to compute force, and advance position and velocity.
      Return the appropriate global and local sizes
     **/

    cl_program timestepModule;
    cl_program ljModule;
    cl_program eamModule;
    cl_program vizModule;

    int err;

    // build the program from the kernel source file

    buildProgramFromFile(&timestepModule, "./src-ocl/timestep_kernels.c", context, deviceId);
    // only build the modules needed for the potential chosen
    if(simHostSoa.eamFlag)
    {
        buildProgramFromFile(&eamModule, "./src-ocl/eam_kernels.c", context, deviceId);
    }
    else
    {
        buildProgramFromFile(&ljModule, "./src-ocl/lj_kernels.c", context, deviceId);
    }
#ifdef INTEROP_VIZ
    buildProgramFromFile(&vizModule, "./src-ocl/viz_kernels.c", context, deviceId);
#endif

    printf("Program built\n");

    if(simHostSoa.eamFlag)
    {
        // create the EAM_Force_x kernels from the program
        forceKernels[0] = clCreateKernel(eamModule, "EAM_Force_1", &err);
        if (!forceKernels[0] || err != CL_SUCCESS)
        {
            printf("Error: Failed to create compute kernel EAM_Force_1!\n");
            exit(1);
        }
	else
	{
            printf("Kernel EAM_Force_1 built\n");
        }
        forceKernels[1] = clCreateKernel(eamModule, "EAM_Force_2", &err);
        if (!forceKernels[1] || err != CL_SUCCESS)
        {
            printf("Error: Failed to create compute kernel EAM_Force_2!\n");
            exit(1);
        }
	else
	{
            printf("Kernel EAM_Force_2 built\n");
        }
        forceKernels[2] = clCreateKernel(eamModule, "EAM_Force_3", &err);
        if (!forceKernels[2] || err != CL_SUCCESS)
        {
            printf("Error: Failed to create compute kernel EAM_Force_3!\n");
            exit(1);
        }
	else
	{
            printf("Kernel EAM_Force_3 built\n");
        }

    }
    else
    {
        // create the ljForce kernel from the program
        forceKernels[0] = clCreateKernel(ljModule, "ljForce", &err);
        if (!forceKernels[0] || err != CL_SUCCESS)
        {
            printf("Error: Failed to create compute kernel ljForce!\n");
            exit(1);
        }
	else
	{
            printf("Kernel ljForce built\n");
        }
    }
    // create the advanceVelocity kernel from the program
    *advanceVelocity = clCreateKernel(timestepModule, "advanceVelocitySoa", &err);
    if (!*advanceVelocity || err != CL_SUCCESS)
    {
        printf("Error: Failed to create compute kernel advanceVelocity!\n");
        exit(1);
    }
    else
    {
        printf("Kernel advanceVelocity built\n");
    }
    // create the advancePosition kernel from the program
    *advancePosition = clCreateKernel(timestepModule, "advancePositionSoa", &err);
    if (!*advancePosition || err != CL_SUCCESS)
    {
        printf("Error: Failed to create compute kernel advancePosition!\n");
        exit(1);
    }
    else
    {
        printf("Kernel advancePosition built\n");
    }

#ifdef INTEROP_VIZ
    // create the Viz kernel from the program
    *Viz = clCreateKernel(vizModule, "Viz", &err);
    if (!*Viz || err != CL_SUCCESS)
    {
        printf("Error: Failed to create compute kernel Viz!\n");
        exit(1);
    }
    else
    {
        printf("Kernel Viz built\n");
    }
#endif

    // determine allowable local work sizes for the device we chose
    if(simHostSoa.eamFlag)
    {
    }
    else
    {
        err = clGetKernelWorkGroupInfo(forceKernels[0], deviceId, CL_KERNEL_WORK_GROUP_SIZE, sizeof(nLocal), nLocal, NULL);
        if (err != CL_SUCCESS)
        {
            printf("Error: Failed to retrieve ljForce work group info! %d\n", err);
            printf("Error: %s\n", print_cl_errstring(err));
            exit(1);
        }
        printf("Maximum local size for ljForce is (%lu, %lu)\n", nLocal[0], nLocal[1]);
    }

    err = clGetKernelWorkGroupInfo(*advanceVelocity, deviceId, CL_KERNEL_WORK_GROUP_SIZE, sizeof(nLocal), nLocal, NULL);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to retrieve advanceVelocity work group info! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    }
    printf("Maximum local size for advanceVelocity is (%lu, %lu)\n", nLocal[0], nLocal[1]);

    err = clGetKernelWorkGroupInfo(*advancePosition, deviceId, CL_KERNEL_WORK_GROUP_SIZE, sizeof(nLocal), nLocal, NULL);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to retrieve advancePosition work group info! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    }
    printf("Maximum local size for advancePosition is (%lu, %lu)\n", nLocal[0], nLocal[1]);

    nLocal[1] += 1;

    // set the global size equal to the numer of atoms
    nGlobal[0] = MAXATOMS;
    nGlobal[1] = simHostSoa.nCells;

    // fudge the loop unrollling here
    nGlobal[0] = nGlobal[0]*UNROLL;
    nGlobal[1] = nGlobal[1]/UNROLL;

    // if the local size is greater than the total size, set local size to global size
    int i;
    for (i=0;i<2;i++)
    {
        if (nGlobal[i] < nLocal[i])
	{
            nLocal[i] = nGlobal[i];
        }
    }
    printf("Global and local sizes are (%lu, %lu), (%lu, %lu)\n", nGlobal[0], nGlobal[1], nLocal[0], nLocal[1]);

}
void buildModulesAoS(cl_kernel *forceKernels, 
        cl_kernel *advancePosition, 
        cl_kernel *advanceVelocity,
        cl_kernel *Viz, 
        HostSimAos simHostSoa, 
        size_t *nLocal, 
        size_t *nGlobal)
{
    /** Build the kernels to compute force, and advance position and velocity.
      Return the appropriate global and local sizes
     **/

    cl_program timestepModule;
    cl_program ljModule;
    cl_program eamModule;
    cl_program vizModule;

    int err;

    // build the program from the kernel source file

    buildProgramFromFile(&timestepModule, "./src-ocl/timestep_kernels.c", context, deviceId);
    // only build the modules needed for the potential chosen
    if(simHostSoa.eamFlag)
    {
        buildProgramFromFile(&eamModule, "./src-ocl/eam_kernels.c", context, deviceId);
    }
    else
    {
        buildProgramFromFile(&ljModule, "./src-ocl/lj_kernels_aos.c", context, deviceId);
    }
#ifdef INTEROP_VIZ
    buildProgramFromFile(&vizModule, "./src-ocl/viz_kernels.c", context, deviceId);
#endif

    printf("Program built\n");

    if(simHostSoa.eamFlag)
    {
        // create the EAM_Force_x kernels from the program
        forceKernels[0] = clCreateKernel(eamModule, "EAM_Force_1_AoS", &err);
        if (!forceKernels[0] || err != CL_SUCCESS)
        {
            printf("Error: Failed to create compute kernel EAM_Force_1!\n");
            exit(1);
        }
	else
	{
            printf("Kernel EAM_Force_1 built\n");
        }
        forceKernels[1] = clCreateKernel(eamModule, "EAM_Force_2_AoS", &err);
        if (!forceKernels[1] || err != CL_SUCCESS)
        {
            printf("Error: Failed to create compute kernel EAM_Force_2!\n");
            exit(1);
        }
	else
	{
            printf("Kernel EAM_Force_2 built\n");
        }
        forceKernels[2] = clCreateKernel(eamModule, "EAM_Force_3_AoS", &err);
        if (!forceKernels[2] || err != CL_SUCCESS)
        {
            printf("Error: Failed to create compute kernel EAM_Force_3!\n");
            exit(1);
        }
	else
	{
            printf("Kernel EAM_Force_3 built\n");
        }

    }
    else
    {
        // create the ljForce kernel from the program
        forceKernels[0] = clCreateKernel(ljModule, "ljForceAos", &err);
        if (!forceKernels[0] || err != CL_SUCCESS)
        {
            printf("Error: Failed to create compute kernel ljForceAos!\n");
            exit(1);
        }
	else
	{
            printf("Kernel ljForceAos built\n");
        }
    }
    // create the advanceVelocity kernel from the program
    *advanceVelocity = clCreateKernel(timestepModule, "advanceVelocityAos", &err);
    if (!*advanceVelocity || err != CL_SUCCESS)
    {
        printf("Error: Failed to create compute kernel advanceVelocityAos!\n");
        exit(1);
    }
    else
    {
        printf("Kernel advanceVelocity built\n");
    }
    // create the advancePosition kernel from the program
    *advancePosition = clCreateKernel(timestepModule, "advancePositionAos", &err);
    if (!*advancePosition || err != CL_SUCCESS)
    {
        printf("Error: Failed to create compute kernel advancePositionAos!\n");
        exit(1);
    }
    else
    {
        printf("Kernel advancePosition built\n");
    }

#ifdef INTEROP_VIZ
    // create the Viz kernel from the program
    *Viz = clCreateKernel(vizModule, "Viz", &err);
    if (!*Viz || err != CL_SUCCESS)
    {
        printf("Error: Failed to create compute kernel Viz!\n");
        exit(1);
    }
    else
    {
        printf("Kernel Viz built\n");
    }
#endif

    // determine allowable local work sizes for the device we chose
    if(simHostSoa.eamFlag)
    {
    }
    else
    {
        err = clGetKernelWorkGroupInfo(forceKernels[0], deviceId, CL_KERNEL_WORK_GROUP_SIZE, sizeof(nLocal), nLocal, NULL);
        if (err != CL_SUCCESS)
        {
            printf("Error: Failed to retrieve ljForce work group info! %d\n", err);
            printf("Error: %s\n", print_cl_errstring(err));
            exit(1);
        }
        printf("Maximum local size for ljForce is (%lu, %lu)\n", nLocal[0], nLocal[1]);
    }

    err = clGetKernelWorkGroupInfo(*advanceVelocity, deviceId, CL_KERNEL_WORK_GROUP_SIZE, sizeof(nLocal), nLocal, NULL);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to retrieve advanceVelocity work group info! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    }
    printf("Maximum local size for advanceVelocity is (%lu, %lu)\n", nLocal[0], nLocal[1]);

    err = clGetKernelWorkGroupInfo(*advancePosition, deviceId, CL_KERNEL_WORK_GROUP_SIZE, sizeof(nLocal), nLocal, NULL);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to retrieve advancePosition work group info! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    }
    printf("Maximum local size for advancePosition is (%lu, %lu)\n", nLocal[0], nLocal[1]);

    //nLocal[1] += 1;

    // set the global size equal to the numer of atoms
    nGlobal[0] = MAXATOMS;
    nGlobal[1] = simHostSoa.nCells;

    // if the local size is greater than the total size, set local size to global size
    int i;
    for (i=0;i<2;i++)
    {
        if (nGlobal[i] < nLocal[i])
	{
            nLocal[i] = nGlobal[i];
        }
    }
    printf("Global and local sizes are (%lu, %lu), (%lu, %lu)\n", nGlobal[0], nGlobal[1], nLocal[0], nLocal[1]);

}

void computeForceOCL(cl_kernel* forceKernels, 
        cl_event *forceEvent, 
        size_t* nGlobal, 
        size_t* nLocal, 
        int eamFlag, 
        int nTot,
        cl_real* tKern)
{
    /** Execute the appropriate force kernels **/

    int err;
    cl_real tSimple, tOverall;
    cl_real tTotal = 0.0;
    if (eamFlag)
    {
#if (PASS_1)
        if (DIAG_LEVEL > 1)
	{
            printf("Running EAM kernel 1..");
            fflush(stdout);
        }
        oclRunKernel(forceKernels[0], forceEvent, nGlobal, nLocal);
        err = clWaitForEvents(1, forceEvent);
        if (DIAG_LEVEL > 1)
	{
            printf("done\n");
            fflush(stdout);
        }
        getElapsedTime(*forceEvent, &tSimple, &tOverall);
        tTotal += tSimple;
#endif
#if (PASS_2)
        if (DIAG_LEVEL > 1)
	{
            printf("Running EAM kernel 2..");
            fflush(stdout);
        }
        oclRunKernel(forceKernels[1], forceEvent, nGlobal, nLocal);
        err = clWaitForEvents(1, forceEvent);
        if (DIAG_LEVEL > 1)
	{
            printf("done\n");
            fflush(stdout);
        }
        getElapsedTime(*forceEvent, &tSimple, &tOverall);
        tTotal += tSimple;
#endif
#if (PASS_3)
        if (DIAG_LEVEL > 1)
	{
            printf("Running EAM kernel 3..");
            fflush(stdout);
        }
        oclRunKernel(forceKernels[2], forceEvent, nGlobal, nLocal);
        err = clWaitForEvents(1, forceEvent);
        if (DIAG_LEVEL > 1)
	{
            printf("done\n");
            fflush(stdout);
        }
        getElapsedTime(*forceEvent, &tSimple, &tOverall);
        tTotal += tSimple;
#endif
        *tKern = tTotal;
        if (DIAG_LEVEL > 0)
            printf("Kernel EAM_Force executed in %.3e secs. (%e us/atom for %d atoms)\n", 
                    tTotal, 1.0e6*tTotal/nTot, nTot);

    }
    else
    {
        if (DIAG_LEVEL > 1)
	{
            printf("Running LJ kernel..");
        }
        oclRunKernel(forceKernels[0], forceEvent, nGlobal, nLocal);
        err = clWaitForEvents(1, forceEvent);
        if (DIAG_LEVEL > 1)
	{
            printf("done\n");
        }
        getElapsedTime(*forceEvent, &tSimple, &tOverall);
        if (DIAG_LEVEL > 0)
            printf("Kernel ljForce executed in %.3e secs. (%e us/atom for %d atoms)\n", 
                    tSimple, 1.0e6*tSimple/nTot, nTot);
        *tKern = tSimple;
    }

}

#ifdef INTEROP_VIZ
// write atom positions to VBO
void oclGraphics(cl_kernel vizKernel, DevSimSoa simDevSoa, size_t* nGlobal, size_t* nLocal)
{
    // Set the arguments of our kernel
    cl_int err;  unsigned int i;
    const n = g_ncells*MAXATOMS*3;    

    cl_mem cl_vertices = vboResources[0];   

    int nArg = 0;
    err  = clSetKernelArg(vizKernel, nArg++, sizeof(cl_mem), &simDevSoa.r.x);
    err |= clSetKernelArg(vizKernel, nArg++, sizeof(cl_mem), &simDevSoa.r.y);
    err |= clSetKernelArg(vizKernel, nArg++, sizeof(cl_mem), &simDevSoa.r.z);
    err |= clSetKernelArg(vizKernel, nArg++, sizeof(cl_mem), &simDevSoa.grid.nAtoms);
    err |= clSetKernelArg(vizKernel, nArg++, sizeof(cl_mem), &simDevSoa.grid.rBox.x);
    err |= clSetKernelArg(vizKernel, nArg++, sizeof(cl_mem), &simDevSoa.grid.rBox.y);
    err |= clSetKernelArg(vizKernel, nArg++, sizeof(cl_mem), &simDevSoa.grid.rBox.z);
    err |= clSetKernelArg(vizKernel, nArg++, sizeof(cl_mem), (void *) &cl_vertices);  
    clFinish(commandq);

    // Execute the kernel 
    err = clEnqueueNDRangeKernel(commandq, vizKernel, 2, NULL, nGlobal, NULL, 0, NULL, 0);  
    clFinish(commandq);

    // Read back and output results
    /*float* data = (float*)malloc(sizeof(float)*n); 
      err = clEnqueueReadBuffer(commandq, cl_vertices, CL_TRUE, 0, n*sizeof(float), data, 0, NULL, 0);
      clFinish(commandq);
      printf("Ran graphics kernel\n");
      for (i=0; i<n; i++) printf("%f ", data[i]);
      printf("\n");*/  
}

// render VBOs
void oclRender()
{
    cl_int err = clEnqueueReleaseGLObjects(commandq, 1, &vboResources[0], 0, 0, 0);
    clFinish(commandq);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(cameraFOV, 1.0f, 1.0f, 2000.0f);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0.0f, 0.0f, 150.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);

    glPushMatrix();
    //float halfGrid = 75.0f/2.0f;
    glTranslatef(-centerX, -centerY, -centerZ); //-halfGrid, -halfGrid, -halfGrid/2.0f);

    float rotationMatrix[16];
    QuaternionGetRotMat(rotationMatrix, q);
    glMultMatrixf(rotationMatrix);

    //float centerX = halfGrid;  float centerY = halfGrid;  float centerZ = halfGrid/2.0f;
    GLfloat matrix[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
    float offsetX = matrix[0]*centerX + matrix[1]*centerY + matrix[2]*centerZ; 
    float offsetY = matrix[4]*centerX + matrix[5]*centerY + matrix[6]*centerZ;
    float offsetZ = matrix[8]*centerX + matrix[9]*centerY + matrix[10]*centerZ;
    offsetX = centerX - offsetX; offsetY = centerY - offsetY; offsetZ = centerZ - offsetZ;
    glTranslatef(-offsetX, -offsetY, -offsetZ);

    //glPointSize(5.0f);
    glColor4f(1.0f, 0.0f, 0.0f, 1.0f);

    glEnableClientState(GL_VERTEX_ARRAY);
    //glEnableClientState(GL_COLOR_ARRAY);
    //glEnableClientState(GL_NORMAL_ARRAY);

    glBindBuffer(GL_ARRAY_BUFFER, vboBuffers[0]);
    glVertexPointer(3, GL_FLOAT, 0, 0);

    //glBindBuffer(GL_ARRAY_BUFFER, vboBuffers[1]);
    //glColorPointer(4, GL_FLOAT, 0, 0);

    //glBindBuffer(GL_ARRAY_BUFFER, vboBuffers[2]);
    //glNormalPointer(GL_FLOAT, 4*sizeof(float), 0);

    glDrawArrays(GL_POINTS, 0, g_ncells*MAXATOMS); 

    glPopMatrix();

    clEnqueueAcquireGLObjects(commandq, 1, &(vboResources[0]), 0, 0, NULL); 
    clFinish(commandq);
}


// initialize vertex buffer objects
void oclInitInterop(int ncells)
{
    g_ncells = ncells;
    int vboSize = ncells*MAXATOMS*3*sizeof(float);
    int numBuffers = 1; 
    int i;
    glGenBuffers(numBuffers, vboBuffers);
    for (i=0; i<numBuffers; i++)
    {
        unsigned int buffer_size = vboSize*sizeof(float);
        glBindBuffer(GL_ARRAY_BUFFER, vboBuffers[i]);
        glBufferData(GL_ARRAY_BUFFER, buffer_size, 0, GL_DYNAMIC_DRAW);
    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    for (i=0; i<numBuffers; i++) 
    {
        cl_int status = CL_SUCCESS;
        vboResources[i] = clCreateFromGLBuffer(context, CL_MEM_WRITE_ONLY, vboBuffers[i], &status);
    }

    clEnqueueAcquireGLObjects(commandq, 1, &(vboResources[0]), 0, 0, NULL); 
    clFinish(commandq);
}

#endif

