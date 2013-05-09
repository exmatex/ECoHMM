
#define N_MAX_ATOMS 64
#define N_MAX_NEIGHBORS 27
#define PERIODIC 1

#define ALLOW_PRINTF 0

#if defined(cl_khr_fp64)  // Khronos extension available?
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)  // AMD extension available?
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

typedef CL_REAL_T real_t; 

__kernel void Viz(
        __global real_t* xPos,
        __global real_t* yPos,
        __global real_t* zPos,
        const __global int* nAtoms,
        __global real_t* x_cen,
        __global real_t* y_cen,
        __global real_t* z_cen, 
        __global float* vertices) 
{ 
  int iAtom = get_global_id(0);
  int iBox = get_global_id(1);

  int offset = N_MAX_ATOMS; 

  if (iAtom < nAtoms[iBox]) 
  {
    vertices[iAtom*3 + offset*iBox*3 + 0] = x_cen[iBox] + xPos[iAtom + offset*iBox];
    vertices[iAtom*3 + offset*iBox*3 + 1] = y_cen[iBox] + yPos[iAtom + offset*iBox];
    vertices[iAtom*3 + offset*iBox*3 + 2] = z_cen[iBox] + zPos[iAtom + offset*iBox];
  }
  else
  {
    vertices[iAtom*3 + offset*iBox*3 + 0] = 0.0f;
    vertices[iAtom*3 + offset*iBox*3 + 1] = 0.0f;
    vertices[iAtom*3 + offset*iBox*3 + 2] = 0.0f;
  }
}
