//Initial implementation of the MD code

/** 
  Since OpenCL doesn't pick up #include properly, we need to manually switch real_t from 
  float to double in each kernel file individually.
 **/

#define UNROLL 4
#define N_MAX_NEIGHBORS 27
#define PERIODIC 1

// diagnostic flag to allow multiple levels of debug output (on CPU only)
#define KERN_DIAG 0

#if defined(cl_khr_fp64)  // Khronos extension available?
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)  // AMD extension available?
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

/* CL_REAL_T is set to single or double depending on compile time flags */
typedef CL_REAL_T real_t;

// Simple version without local blocking to check for correctness
__kernel void ljForce(
      __global const real_t* xPos,
      __global const real_t* yPos,
      __global const real_t* zPos,
      __global real_t* fx,
      __global real_t* fy,
      __global real_t* fz,
      __global real_t* energy,
      __global const real_t* dcx,
      __global const real_t* dcy,
      __global const real_t* dcz,
      __global const real_t* bounds,
      __global const int* neighborList,
      __global const int* nNeighbors,
      __global const int* nAtoms,
      const real_t sigma,
      const real_t epsilon,
      const real_t cutoff) 
{

   // no loop unrolling
   /*
   int iAtom = get_global_id(0);
   int iBox = get_global_id(1);
   int maxAtoms = get_global_size(0);
   */

   // loop unrolling
   int iLocal = get_global_id(0);
   int iLine = get_global_id(1);
   int lineLength = get_global_size(0);

   int maxAtoms = lineLength/UNROLL;
   int iBox = iLine*UNROLL + iLocal/maxAtoms;
   int iAtom = iLocal  % maxAtoms;

   real_t dx, dy, dz;
   real_t r2, r6;
   real_t fr, e;

   real_t dxbox, dybox, dzbox;

   // accumulate local force value
   real_t fxItem, fyItem, fzItem;

   real_t rCut = cutoff;
   real_t rCut2 = rCut*rCut;
   real_t s6 = sigma*sigma*sigma*sigma*sigma*sigma;

   int iOffset, jOffset;
   int iParticle, jParticle;

   // zero out forces on particles
   fxItem = 0.0;
   fyItem = 0.0;
   fzItem = 0.0;

   e = 0.0;

   iOffset = iBox*maxAtoms; //N_MAX_ATOMS;
   iParticle = iOffset + iAtom;

   if (iAtom < nAtoms[iBox])
   {// each thread executes on a single atom in the box

#if(KERN_DIAG > 1) 
      printf("i = %d, %f, %f, %f\n", iParticle, xPos[iParticle], yPos[iParticle], zPos[iParticle]);
#endif

#if(KERN_DIAG > 0) 
      printf("iBox = %d, nNeighbors = %d\n", iBox, nNeighbors[iBox]);
#endif

      for (int j = 0; j<nNeighbors[iBox]; j++)
      {// loop over neighbor cells
         int jBox = neighborList[iBox*N_MAX_NEIGHBORS + j];
         jOffset = jBox*maxAtoms; //N_MAX_ATOMS;

         // compute box center offsets
         dxbox = dcx[iBox] - dcx[jBox];
         dybox = dcy[iBox] - dcy[jBox];
         dzbox = dcz[iBox] - dcz[jBox];

         // correct for periodic 
         if(PERIODIC)
         {
            if (dxbox<-0.5*bounds[0]) dxbox += bounds[0];
            else if (dxbox > 0.5*bounds[0] ) dxbox -= bounds[0];
            if (dybox<-0.5*bounds[1]) dybox += bounds[1];
            else if (dybox > 0.5*bounds[1] ) dybox -= bounds[1];
            if (dzbox<-0.5*bounds[2]) dzbox += bounds[2];
            else if (dzbox > 0.5*bounds[2] ) dzbox -= bounds[2];
         }

         // printf("dxbox, dybox, dzbox = %f, %f, %f\n", dxbox, dybox, dzbox);

         for (int jAtom = 0; jAtom<nAtoms[jBox]; jAtom++)
         {// loop over all groups in neighbor cell 

            jParticle = jOffset + jAtom; // global offset of particle

            dx = xPos[iParticle] - xPos[jParticle] + dxbox;;
            dy = yPos[iParticle] - yPos[jParticle] + dybox;;
            dz = zPos[iParticle] - zPos[jParticle] + dzbox;;

#if(KERN_DIAG > 1) 
            //printf("dx, dy, dz = %f, %f, %f\n", dx, dy, dz);
            //printf("i = %d, j = %d, %f, %f, %f\n", iParticle, jParticle, xPos[jParticle], yPos[jParticle], zPos[jParticle]);
#endif

            r2 = dx*dx + dy*dy + dz*dz;

            if ( r2 <= rCut2 && r2 > 0.0)
            {// no divide by zero

#if(KERN_DIAG > 1) 
               printf("%d, %d, %f\n", iParticle, jParticle, r2);
               //printf("r2, rCut = %f, %f\n", r2, rCut);
#endif

               // reciprocal of r2 now
               r2 = (real_t)1.0/r2;

               r6 = r2*r2*r2;

               e += r6*(s6*r6 - 1.0);

#if(KERN_DIAG > 1) 
               //printf("%d, %d, %f\n", iParticle, jParticle, r2);
               //printf("iParticle = %d, jParticle = %d, i_b = %d, r6 = %f\n", iParticle, jParticle, i_b, r6);
#endif

               fr = -4.0*epsilon*s6*r2*r6*(12.0*r6*s6 - 6.0);

               fxItem += dx*fr;
               fyItem += dy*fr;
               fzItem += dz*fr;

            } 
            else 
            {
            }


         } // loop over all atoms
      } // loop over neighbor cells

      fx[iParticle] = fxItem;
      fy[iParticle] = fyItem;
      fz[iParticle] = fzItem;

      // since we loop over all particles, each particle contributes 1/2 the pair energy to the total
      energy[iParticle] = e*2.0*epsilon*s6;

   }
}

__kernel void ljForceLocal(
      __global real_t* xPos,
      __global real_t* yPos,
      __global real_t* zPos,
      __global real_t* fx,
      __global real_t* fy,
      __global real_t* fz,
      __global real_t* energy,
      __global real_t* dcx,
      __global real_t* dcy,
      __global real_t* dcz,
      __global real_t* bounds,
      __global int* neighborList,
      __global int* nNeighbors,
      __global int* nAtoms,
      __local real_t* x_ii,
      __local real_t* y_ii,
      __local real_t* z_ii,
      __local real_t* x_ij,
      __local real_t* y_ij,
      __local real_t* z_ij,
      const real_t sigma,
      const real_t epsilon,
      const real_t cutoff) 
{




   real_t dx, dy, dz;
   real_t r2, r6;
   real_t fr;
   real_t dxbox, dybox, dzbox;

   real_t rCut = 5.0*sigma;
   real_t rCut2 = rCut*rCut;
   real_t s6 = sigma*sigma*sigma*sigma*sigma*sigma;

   int i_p;

   int iOffset, jOffset;
   int jBox;
   int iParticle, jParticle;

   // zero out forces, energy on particles
   real_t fx_ii = 0.0;
   real_t fy_ii = 0.0;
   real_t fz_ii = 0.0;

   real_t e = 0.0;

   int iBox = get_global_id(1);
   int iLocal = get_local_id(0);
   int maxAtoms = get_global_size(0);
   int nGroups = get_num_groups(0);
   int nItems = get_local_size(0);
   int groupId = get_group_id(0);

   int group_offset = nItems*groupId;
   int iAtom = group_offset + iLocal;
   iOffset = iBox*maxAtoms; //N_MAX_ATOMS;
   iParticle = iOffset + iAtom;

#if(KERN_DIAG) 
   printf("Number of work groups: %d\n", nGroups);
   printf("Number of work items: %d\n", nItems);
#endif

   if (iAtom < nAtoms[iBox])
   {// each thread executes on a single atom in the box

#if(KERN_DIAG) 
      //printf("iParticle = %d\n", iParticle);
      //printf("i_global = %d, iLocal = %d, i_b = %d, nItems = %d\n", i_global, iLocal, i_b, nItems);
#endif

      // load particle data into local arrays
      x_ii[iLocal] = xPos[iParticle];
      y_ii[iLocal] = yPos[iParticle];
      z_ii[iLocal] = zPos[iParticle];

      barrier(CLK_LOCAL_MEM_FENCE);

#if(KERN_DIAG) 
      //printf("x_ii, y_ii, z_ii = %f, %f, %f\n", x_ii[iLocal], y_ii[iLocal], z_ii[iLocal]);
      printf("%d, %f, %f, %f\n", iParticle, x_ii[iLocal], y_ii[iLocal], z_ii[iLocal]);
#endif

      for (int j = 0; j<nNeighbors[iBox]; j++)
      {// loop over neighbor cells
         jBox = neighborList[iBox*N_MAX_NEIGHBORS + j];
         jOffset = jBox*maxAtoms; //N_MAX_ATOMS;

         // compute box center offsets
         dxbox = dcx[iBox] - dcx[jBox];
         dybox = dcy[iBox] - dcy[jBox];
         dzbox = dcz[iBox] - dcz[jBox];

         // correct for periodic 
         if(PERIODIC)
         {
            if (dxbox<-0.5*bounds[0]) dxbox += bounds[0];
            else if (dxbox > 0.5*bounds[0] ) dxbox -= bounds[0];
            if (dybox<-0.5*bounds[1]) dybox += bounds[1];
            else if (dybox > 0.5*bounds[1] ) dybox -= bounds[1];
            if (dzbox<-0.5*bounds[2]) dzbox += bounds[2];
            else if (dzbox > 0.5*bounds[2] ) dzbox -= bounds[2];
         }

         // printf("dxbox, dybox, dzbox = %f, %f, %f\n", dxbox, dybox, dzbox);

         for (int j_b = 0; j_b<nGroups; j_b++)
         {// loop over all groups in neighbor cell 

            // use iLocal to load data in blocks of size nItems
            x_ij[iLocal] = xPos[iLocal + j_b*nItems + jOffset];
            y_ij[iLocal] = yPos[iLocal + j_b*nItems + jOffset];
            z_ij[iLocal] = zPos[iLocal + j_b*nItems + jOffset];

            barrier(CLK_LOCAL_MEM_FENCE);

            for (int jLocal=0;jLocal < nItems; jLocal ++)
            {// loop over all atoms in group

               jParticle = jLocal+ j_b*nItems; // global offset of particle

               dx = x_ii[iLocal] - x_ij[jLocal] + dxbox;
               dy = y_ii[iLocal] - y_ij[jLocal] + dybox;
               dz = z_ii[iLocal] - z_ij[jLocal] + dzbox;

#if(KERN_DIAG) 
               //printf("dx, dy, dz = %f, %f, %f\n", dx, dy, dz);
               printf("%d, %f, %f, %f\n", jParticle, x_ij[jLocal], y_ij[jLocal], z_ij[jLocal]);
               printf("%d, %d, %f, %f, %f\n", iParticle, jParticle, dx, dy, dz);
#endif

               r2 = dx*dx + dy*dy + dz*dz;

#if(KERN_DIAG) 
               printf("%d, %d, %f\n", iParticle, jParticle, r2);
               //printf("r2, rCut = %f, %f\n", r2, rCut);
#endif

               if ( r2 <= rCut2 && r2 > 0.0)
               {// no divide by zero

                  // reciprocal of r2 now
                  r2 = (real_t)1.0/r2;

                  r6 = r2*r2*r2;

                  e += r6*(s6*r6 - 1.0);

#if(KERN_DIAG) 
                  //printf("%d, %d, %f\n", iParticle, jParticle, r2);
                  //printf("iParticle = %d, jParticle = %d, i_b = %d, r6 = %f\n", iParticle, jParticle, i_b, r6);
#endif

                  fr = 4.0*epsilon*s6*r2*r6*(12.0*r6*s6 - 6.0);

                  fx_ii += dx*fr;
                  fy_ii += dy*fr;
                  fz_ii += dz*fr;

               } 
               else 
               {
               }

            } // loop over all atoms in group

         } // loop over all groups in neighbor cell
      } // loop over neighbor cells

      fx[iParticle] = fx_ii;
      fy[iParticle] = fy_ii;
      fz[iParticle] = fz_ii;


      // since we loop over all particles, each particle contributes 1/2 the pair energy to the total
      energy[iParticle] = e*2.0*epsilon*s6;

      barrier(CLK_LOCAL_MEM_FENCE);
   }

}

