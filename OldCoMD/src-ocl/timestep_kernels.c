// timestep subroutines

/**
  Since OpenCL doesn't pick up #include properly, we need to manually switch real_t from 
  float to double in each kernel file individually.
  **/

#define N_MAX_NEIGHBORS 27
#define PERIODIC 1

#define KERN_DIAG 0

#if defined(cl_khr_fp64)  // Khronos extension available?
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)  // AMD extension available?
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

/* CL_REAL_T is set to single or double depending on compile time flags */
typedef CL_REAL_T real_t; 
typedef CL_REAL4_T cl_real4; 


__kernel void advanceVelocitySoa (
        __global real_t* px,
        __global real_t* py,
        __global real_t* pz,
        __global const real_t* fx,
        __global const real_t* fy,
        __global const real_t* fz,
        __global const int* nAtoms,
        const real_t dt
        )

{
    int iAtom = get_global_id(0);
    int iBox = get_global_id(1);

    int offset = get_global_size(0);
    int iOff = iAtom + offset*iBox;

    real_t dtLocal = dt;

#if(KERN_DIAG > 0) 
    if (iAtom == 0 && iBox == 0) printf(" AV dt = %e\n", dt);

    if (nAtoms[iBox] > 0) printf("%d, %d, %d\n", iBox, iAtom, nAtoms[iBox]);
#endif

    if (iAtom < nAtoms[iBox])
    {
        px[iOff] -= dtLocal*fx[iOff];
        py[iOff] -= dtLocal*fy[iOff];
        pz[iOff] -= dtLocal*fz[iOff];
    }

}

__kernel void advancePositionSoa (
        __global const real_t* px,
        __global const real_t* py,
        __global const real_t* pz,
        __global real_t* xPos,
        __global real_t* yPos,
        __global real_t* zPos,
        //__global const real_t* mass,
        __global const int* nAtoms,
        const real_t rMass,
        const real_t dt
        )

{
    int iAtom = get_global_id(0);
    int iBox = get_global_id(1);

    int offset = get_global_size(0);
    int iOff = iAtom + offset*iBox;

#if(KERN_DIAG > 0) 

    if (nAtoms[iBox] > 0) printf("%d, %d, %d\n", iBox, iAtom, nAtoms[iBox]);
#endif

    real_t dtLocal = dt/rMass; 

#if(KERN_DIAG > 0) 
    if (iAtom == 0 && iBox == 0) printf(" AP dt = %e\n", dtLocal);
    if (iAtom == 0 && iBox == 0) printf("%f, %f, %f\n", dtLocal, dt, rMass);
#endif

    if (iAtom < nAtoms[iBox])
    {
        xPos[iOff] += dtLocal*px[iOff];
        yPos[iOff] += dtLocal*py[iOff];
        zPos[iOff] += dtLocal*pz[iOff];

#if(KERN_DIAG > 0) 
    if (iAtom == 0 && iBox == 0) printf("%d, %d, %f, %f, %f\n", iBox, iAtom, xPos[iOff], yPos[iOff], zPos[iOff]);
#endif
    }

}

__kernel void advanceVelocityAos (
        __global cl_real4* p,
        __global const cl_real4* f,
        __global const int* nAtoms,
        const real_t dt
        )

{
    int iAtom = get_global_id(0);
    int iBox = get_global_id(1);

    int offset = get_global_size(0);
    int iOff = iAtom + offset*iBox;

    cl_real4 dtLocal = {dt, dt, dt, 0.0};
    //cl_real4 dtLocal = {0.0, 0.0, 0.0, 0.0};

#if(KERN_DIAG > 0) 
    if (iAtom == 0 && iBox == 0) printf("AoS AV dt = %e\n", dt);

    if (nAtoms[iBox] > 0) printf("%d, %d, %d\n", iBox, iAtom, nAtoms[iBox]);
#endif

    if (iAtom < nAtoms[iBox])
    {
        p[iOff] -= dtLocal*f[iOff];
    }

}

__kernel void advancePositionAos (
        __global const cl_real4* p,
        __global cl_real4* pos,
        //__global const real_t* mass,
        __global const int* nAtoms,
        const real_t rMass,
        const real_t dt
        )

{
    int iAtom = get_global_id(0);
    int iBox = get_global_id(1);

    int offset = get_global_size(0);
    int iOff = iAtom + offset*iBox;

#if(KERN_DIAG > 0) 
    if (iAtom == 0 && iBox == 0) printf("AoS AP dt = %e\n", dt);

    if (nAtoms[iBox] > 0) printf("%d, %d, %d\n", iBox, iAtom, nAtoms[iBox]);
#endif

    real_t rdt = dt/rMass;

    cl_real4 dtLocal = {rdt, rdt, rdt, 0.0};

#if(KERN_DIAG > 0) 
    if (iAtom == 0 && iBox == 0) printf("%e, %e, %e\n", dtLocal.x, dt, rMass);
#endif

    if (iAtom < nAtoms[iBox])
    {
        pos[iOff] += dtLocal*p[iOff];

#if(KERN_DIAG > 0) 
    if (iAtom == 0 && iBox == 0) printf("%d, %d, %f, %f, %f\n", iBox, iAtom, pos[iOff].x, pos[iOff].y, pos[iOff].z);
#endif
    }

}


