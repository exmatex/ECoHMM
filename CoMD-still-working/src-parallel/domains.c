/*

Copyright (c) 2011, Los Alamos National Security, LLC All rights
reserved.  Copyright 2011. Los Alamos National Security, LLC. This
software was produced under U.S. Government contract DE-AC52-06NA25396
for Los Alamos National Laboratory (LANL), which is operated by Los
Alamos National Security, LLC for the U.S. Department of Energy. The
U.S. Government has rights to use, reproduce, and distribute this
software.

NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY
WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF
THIS SOFTWARE.

If software is modified to produce derivative works, such modified
software should be clearly marked, so as not to confuse it with the
version available from LANL.

Additionally, redistribution and use in source and binary forms, with
or without modification, are permitted provided that the following
conditions are met:

· Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

· Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

· Neither the name of Los Alamos National Security, LLC, Los Alamos
  National Laboratory, LANL, the U.S. Government, nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS
ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

/**
 * this file will split space up into domains so that
 * we will have fast neighbor listing
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>

#include "pmd.h"
#include "parallel.h"
#include "geometry.h"
#include "memUtils.h"

#define freeMe(s,element) {if(s->element) suAlignedFree(s->element);  s->element = NULL;}

void getBoxIxyz(SimFlat* s, int iBox, int* iret)
{
   // given a box id return the ix,iy,iz coordinates
   iret[0] = iBox%s->lNbx[0];
   iret[1] = (iBox/s->lNbx[0])%s->lNbx[1];
   iret[2] = iBox/s->lNbx[0]/s->lNbx[1];

   return;
}

/** 
 * given x,y,z in world co-ordinates return the
 * box in which those coordinates fall
 */
int getBoxIDWorldCoords(SimFlat* s, const real2_t xyz[3], int* ibx)
{
   int iBox;

   ibx[0] = (int)(floor((xyz[0]-s->minLBound[0])/s->boxSize[0]));
   ibx[1] = (int)(floor((xyz[1]-s->minLBound[1])/s->boxSize[1]));
   ibx[2] = (int)(floor((xyz[2]-s->minLBound[2])/s->boxSize[2]));

   iBox = ibx[0]+s->lNbx[0]*ibx[1]+s->lNbx[0]*s->lNbx[1]*ibx[2];

   return iBox;
}

void destroyDomains(SimFlat* s)
{
   // release memory for particles
   freeMe(s,dCenter);
   freeMe(s,nAtoms);
   freeMe(s,id);
   freeMe(s,iType);
   freeMe(s,mass);
   freeMe(s,r);
   freeMe(s,p);
   freeMe(s,f);
   freeMe(s,fi);
   freeMe(s,rho);
   s->stateflag = SIM_NOSTATE;

   return;
}

int allocDomainArrays(SimFlat* s)
{
   int ierr;
   int i;

   ierr = 0;
   s->dCenter = (real3*)suAlignedCalloc(s->nTotalBoxes*sizeof(real3));
   s->nAtoms = (int*)suAlignedCalloc(s->nTotalBoxes*sizeof(int));
   s->id = (int*)suAlignedCalloc(MAXATOMS*s->nTotalBoxes*sizeof(int));
   s->iType = (int*)suAlignedCalloc(MAXATOMS*s->nTotalBoxes*sizeof(int));
   s->mass = (real_t*)suAlignedCalloc(MAXATOMS*s->nTotalBoxes*sizeof(real_t));
   s->r = (real3*)suAlignedCalloc(MAXATOMS*s->nTotalBoxes*sizeof(real3));
   s->p = (real3*)suAlignedCalloc(MAXATOMS*s->nTotalBoxes*sizeof(real3));
   s->f = (real4*)suAlignedCalloc(MAXATOMS*s->nTotalBoxes*sizeof(real4));
   s->fi = (real_t*)suAlignedCalloc(MAXATOMS*s->nTotalBoxes*sizeof(real_t));
   s->rho = (real_t*)suAlignedCalloc(MAXATOMS*s->nTotalBoxes*sizeof(real_t));

#ifdef USE_IN_SITU_VIZ
   s->centro = (real_t*)suAlignedCalloc(MAXATOMS*s->nTotalBoxes*sizeof(real_t));
#endif 

   ierr += ((s->dCenter)?0:1);
   ierr += ((s->nAtoms)?0:1);
   ierr += ((s->id)?0:1);
   ierr += ((s->iType)?0:1);
   ierr += ((s->mass)?0:1);
   ierr += ((s->r)?0:1);
   ierr += ((s->p)?0:1);
   ierr += ((s->f)?0:1);
   ierr += ((s->fi)?0:1);
   ierr += ((s->rho)?0:1);

   if(ierr) destroyDomains(s);
   s->stateflag |= SIM_ALLOCED;

   return ierr;

}

int copyDomainArrays(SimFlat* destination, SimFlat* source)
{
   memcpy(destination->dCenter, source->dCenter, source->nTotalBoxes*sizeof(real3)); 
   memcpy(destination->nAtoms, source->nAtoms, source->nTotalBoxes*sizeof(int)); 
   memcpy(destination->id, source->id, MAXATOMS*source->nTotalBoxes*sizeof(int)); 
   memcpy(destination->iType, source->iType, MAXATOMS*source->nTotalBoxes*sizeof(int)); 
   memcpy(destination->mass, source->mass, MAXATOMS*source->nTotalBoxes*sizeof(real_t)); 
   memcpy(destination->r, source->r, MAXATOMS*source->nTotalBoxes*sizeof(real3)); 
   memcpy(destination->p, source->p, MAXATOMS*source->nTotalBoxes*sizeof(real3)); 
   memcpy(destination->f, source->f, MAXATOMS*source->nTotalBoxes*sizeof(real4)); 
   memcpy(destination->fi, source->fi, MAXATOMS*source->nTotalBoxes*sizeof(real_t)); 
   memcpy(destination->rho, source->rho, MAXATOMS*source->nTotalBoxes*sizeof(real_t)); 

   return 0;
}

void allocDomains(SimFlat* s)
{
   if ( ! s ) suAbort(1,(char *) "s unallocated in allocDomains()");
   //if ( ! s->pot)  suAbort(2,(char *) "allocDomains() called without s->pot()");

   // decide how many boxes are needed
   s->nBoxes = 1;

   real_t strain[3];
   // dummy strain field
   strain[0] = s->defGrad;
   // temporary kluge to make Viz work
   if (strain[0] == 0.0) strain[0]=1.0;
   strain[1] = 1.0;
   strain[2] = 1.0;

   // box factor
   real_t boxFactor = 1.0;

   // calculate number of boxes and box size
   for(int j=0; j<3; j++)
   {
      s->gNbx[j] = (int)floor(s->bounds[j]/(s->pot->cutoff*boxFactor*strain[j]));
      s->lNbx[j] = s->gNbx[j] / s->nproc[j];
      subdivideGeometry(s);

      if (s->lNbx[j] < 1) 
         suAbort(3,(char *) "Need at least one cutoff wide domain in allocDomains()");

      s->nBoxes = s->nBoxes*s->lNbx[j];
      s->boxSize[j] = s->lSize[j]/(real_t)s->lNbx[j];
   }

/* from EAM
2 Boxes wide
NBufCell = 4 * ((Xcells+4) * (Ycells + Zcells + 4) + Ycells*Zcells);

1 box wide
NBufCell = 2 * ((Xcells+2) * (Ycells + Zcells + 2) + Ycells*Zcells);
*/
#ifdef DOUBLE_WIDE
   s->nHaloBoxes = 4 * ((s->lNbx[0] + 4) *
                         (s->lNbx[1] + s->lNbx[2] + 4) +
                         (s->lNbx[1] * s->lNbx[2]));
#else
   s->nHaloBoxes = 2 * ((s->lNbx[0] + 2) *
                         (s->lNbx[1] + s->lNbx[2] + 2) +
                         (s->lNbx[1] * s->lNbx[2]));
#endif

   s->nTotalBoxes = s->nBoxes + s->nHaloBoxes;

   if(allocDomainArrays(s)) suAbort(4,(char *) "Unable to allocate domains in allocDomains()");

   // calculate box centers
   for(int iBox=0;iBox<s->nBoxes; iBox++)
   {
      int ib3[3];

      getBoxIxyz(s,iBox,ib3);

      for(int j=0; j<3; j++)
      {
         s->dCenter[iBox][j] = s->minLBound[j] + s->boxSize[j]*(real_t)ib3[j];
      }
   }

   return;
}

void copyDomains(SimFlat* destination, SimFlat* source)
{
    destination->nTot = source->nTot;

    copyDomainArrays(destination, source);
}

void copyAtom(SimFlat* s, int iatom, int iBox, int jAtom, int jBox)
{
   // copy atom iatom in domain iBox to atom jAtom in box jBox
   const int iOff = MAXATOMS*iBox+iatom;
   const int jOff = MAXATOMS*jBox+jAtom;
   s->id[jOff] = s->id[iOff];
   s->iType[jOff] = s->iType[iOff];
   s->mass[jOff] = s->mass[iOff];
   memcpy(s->r[jOff],s->r[iOff],sizeof(real3));
   memcpy(s->f[jOff],s->f[iOff],sizeof(real4));
   memcpy(s->p[jOff],s->p[iOff],sizeof(real3));
   s->fi[jOff] = s->fi[iOff];
   s->rho[jOff] = s->rho[iOff];

   return;
}

void moveAtom(SimFlat* s, int iId, int iBox, int jBox)
{
   int nj,ni;

   nj = s->nAtoms[jBox];
   copyAtom(s, iId, iBox, nj, jBox);
   s->nAtoms[jBox]++;
   if(s->nAtoms[jBox]>= MAXATOMS)
   {
     simulationAbort(-11,"Increase maxatoms");
   }

   s->nAtoms[iBox]--;
   ni = s->nAtoms[iBox];
   if(ni) copyAtom(s,ni,iBox,iId,iBox);
  
   return;
}

void removeAtom(SimFlat* s, int iId, int iBox)
{
    int ni;

    s->nAtoms[iBox]--;
    ni = s->nAtoms[iBox];
    if(ni) copyAtom(s,ni,iBox,iId,iBox);

    return;
}

/**
 * finds an appropriate box for an atom based on the
 * spatial coordinates and puts it in there.
 *
 * reallocates memory if needed.
 */
void putAtomInBox(SimFlat* s,
        const int id, const char flagMove, const int iType,
        const real_t mass,
        const real2_t x,const real2_t y,const real2_t z,
        const real2_t px,const real2_t py,const real2_t pz,
        const real2_t fx,const real2_t fy,const real2_t fz)
{
   int iBox;
   int iOff;
   int ibx[3];
   real2_t xyz[3] = {x,y,z};

   if ( ! s->stateflag ) suAbort(-10,(char *) "\n\n  ERROR: s not allocated in putAtomInBox()\n\n\n");

   // push atom into primary period
   for(int m=0; m<3; m++)
   {
      if(xyz[m] < 0.0 ) xyz[m] += s->maxGBound[m];
      else if (xyz[m] >= s->maxGBound[m] ) xyz[m] -= s->maxGBound[m];
   }

   // Find correct box.
   iBox = getBoxIDWorldCoords(s,xyz,ibx);
   iOff = iBox*MAXATOMS;
   iOff += s->nAtoms[iBox];

   // assign values to array elements
   s->nTot++;
   s->nAtoms[iBox]++;
   s->id[iOff] = id;
   s->iType[iOff] = iType;
   s->mass[iOff] = mass;

   for(int m=0; m<3; m++)
   {
      s->r[iOff][m] = (real_t)(xyz[m]-s->dCenter[iBox][m]);
   }

   s->p[iOff][0] = (real_t)px;
   s->p[iOff][1] = (real_t)py;
   s->p[iOff][2] = (real_t)pz;

   s->f[iOff][0] = (real_t)fx;
   s->f[iOff][1] = (real_t)fy;
   s->f[iOff][2] = (real_t)fz;

   return;
}

/**
 * box functions
 */

int getIBoxFromIxyz3(SimFlat* s, int* ixyz)
{
   int ibox = 0;
   int ix, iy, iz;

   ix = ixyz[0]; iy = ixyz[1]; iz = ixyz[2];


  
   // Halo in Z+
   if (iz == s->lNbx[2])
   {
      ibox = s->nBoxes + 2*s->lNbx[2]*s->lNbx[1] + 2*s->lNbx[2]*(s->lNbx[0]+2) +
         (s->lNbx[0]+2)*(s->lNbx[1]+2) + (s->lNbx[0]+2)*(iy+1) + (ix+1);
   }

   // Halo in Z-
   else if (iz == -1)
   {
      ibox = s->nBoxes + 2*s->lNbx[2]*s->lNbx[1] + 2*s->lNbx[2]*(s->lNbx[0]+2) +
         (s->lNbx[0]+2)*(iy+1) + (ix+1);
   }

   // Halo in Y+
   else if (iy == s->lNbx[1])
   {
      ibox = s->nBoxes + 2*s->lNbx[2]*s->lNbx[1] + s->lNbx[2]*(s->lNbx[0]+2) +
         (s->lNbx[0]+2)*iz + (ix+1);
   }

   // Halo in Y-
   else if (iy == -1)
   {
      ibox = s->nBoxes + 2*s->lNbx[2]*s->lNbx[1] + iz*(s->lNbx[0]+2) + (ix+1);
   }

   // Halo in X+
   else if (ix == s->lNbx[0])
   {
      ibox = s->nBoxes + s->lNbx[1]*s->lNbx[2] + iz*s->lNbx[1] + iy;
   }

   // Halo in X-
   else if (ix == -1)
   {
      ibox = s->nBoxes + iz*s->lNbx[1] + iy;
   }

   // No halos
   else
   {
      ibox = ix + s->lNbx[0]*iy + s->lNbx[0]*s->lNbx[1]*iz;
   }

   return ibox;
}

int getIBoxFromIxyz3NP(SimFlat* s, int* ixyz)
{
   // for non-periodic boundary conditions, decide whether we
   // are at a boundary or not **/

   if (! PERIODIC)
   {
      for(int j=0; j<3; j++)
      {
         if((ixyz[j]<0)||(ixyz[j]>=s->lNbx[j])) return -1; 
      }
   }

   return getIBoxFromIxyz3(s,ixyz);
}

int getIBoxFromIxyz(SimFlat* s, int ix, int iy, int iz)
{
   int iBox;
   int ixyz[3] = {ix,iy,iz};

   return getIBoxFromIxyz3NP(s,ixyz);
}

void getIxyz3(SimFlat* s, int ibox, int* i3)
{

   int ix, iy, iz;

   // If a local box
   if(ibox < s->nBoxes)
   {
      ix = (ibox % s->lNbx[0]);
      iy = (ibox / s->lNbx[0]) % (s->lNbx[1]);
      iz = (ibox / s->lNbx[0] / s->lNbx[1]) % (s->lNbx[2]);
   }
   // It's a halo box
   else 
   {
      int ink;
      ink = ibox - s->nBoxes;
      if (ink < 2*s->lNbx[1]*s->lNbx[2])
      {
         if (ink < s->lNbx[1]*s->lNbx[2]) 
         {
            ix = 0;
         }
         else 
         {
            ink -= s->lNbx[1]*s->lNbx[2];
            ix = s->lNbx[0] + 1;
         }
         iy = 1 + ink % s->lNbx[1];
         iz = 1 + ink / s->lNbx[1];
      }
      else if (ink < (2 * s->lNbx[2] * (s->lNbx[1] + s->lNbx[0] + 2))) 
      {
         ink -= 2 * s->lNbx[2] * s->lNbx[1];
         if (ink < ((s->lNbx[0] + 2) *s->lNbx[2])) 
         {
            iy = 0;
         }
         else 
         {
            ink -= (s->lNbx[0] + 2) * s->lNbx[2];
            iy = s->lNbx[1] + 1;
         }
            ix = ink % (s->lNbx[0] + 2);
            iz = 1 + ink / (s->lNbx[0] + 2);
      }
      else 
      {
         ink -= 2 * s->lNbx[2] * (s->lNbx[1] + s->lNbx[0] + 2);
         if (ink < ((s->lNbx[0] + 2) * (s->lNbx[1] + 2))) 
         {
            iz = 0;
         }
         else 
         {
            ink -= (s->lNbx[0] + 2) * (s->lNbx[1] + 2);
            iz = s->lNbx[2] + 1;
         }
         ix = ink % (s->lNbx[0] + 2);
         iy = ink / (s->lNbx[0] + 2);
      }

      // Calculated as off by 1
      ix--;
      iy--;
      iz--;
   }

   i3[0] = ix;
   i3[1] = iy;
   i3[2] = iz;
}

void getIxyz(SimFlat* s, int iBox, int* ix, int* iy, int* iz)
{
   int ixyz[3];

   getIxyz3(s,iBox,ixyz);

   *ix = ixyz[0];
   *iy = ixyz[1];
   *iz = ixyz[2];

   return;
}

void reBoxAll(SimFlat* s)
{
   for(int iBox=0; iBox<s->nBoxes; iBox++)
   {
      int i;
      int ixold[3];
      int iOff;
      getIxyz3(s, iBox, ixold);
      i=0;
      iOff = iBox*MAXATOMS;

      while(i<s->nAtoms[iBox])
      {
         int ixnew[3];
         int jBox;
         real3 rNew;
         memcpy(rNew,s->r[iOff],sizeof(rNew));

         for(int k=0; k<3; k++)
         {
	         if(s->r[iOff][k] < 0.0)
            {
	            ixnew[k] = ixold[k]-1;
	            rNew[k] += s->boxSize[k];
	         }
	         else if(s->r[iOff][k] >= s->boxSize[k])
            {
	            ixnew[k] = ixold[k]+1;
	            rNew[k] -= s->boxSize[k];
	         }
	         else
            {
	            ixnew[k] = ixold[k];
	         }
         }

         jBox = getIBoxFromIxyz3NP(s,ixnew);
         if((jBox==iBox)||jBox<0||jBox>=s->nBoxes)
         {
	         // do nothing if same box or non-periodic boundary
	         i++;
	         iOff++;
         }
         else
         {
	         // move atom to new box
	         memcpy(s->r[iOff],rNew,sizeof(rNew));
	         moveAtom(s,i,iBox,jBox);
         }
      }
   }

   return;
}

void reBoxAll2(SimFlat* sim)
{
    int iOff, iBox;

    for (int iz = 0; iz < sim->lNbx[2]; iz++)
    {
        for (int iy = 0; iy < sim->lNbx[1]; iy++)
        {
            for (int ix = 0; ix < sim->lNbx[0]; ix++)
            {
                if (ix == 0 || ix == sim->lNbx[0]-1 ||
                    iy == 0 || iy == sim->lNbx[1]-1 ||
                    iz == 0 || iz == sim->lNbx[2]-1)
                {
                    iBox = getIBoxFromIxyz(sim, ix, iy, iz);

                    iOff = MAXATOMS * iBox;

                    for (int ii = 0; ii < sim->nAtoms[iBox]; ii++, iOff++)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            if (sim->r[iOff][k] < 0.0 ||
                               sim->r[iOff][k] >= sim->boxSize[k])
                            {
                                printf("%d remove atom %d from box %d ixyz %d %d %d pr %f %f %f dCenter %f %f %f\n", 
				getMyParallel(), ii, iBox, ix, iy, iz, 
				sim->r[iOff][0], sim->r[iOff][1], sim->r[iOff][2], 
				sim->dCenter[iBox][0], sim->dCenter[iBox][1], sim->dCenter[iBox][2]);

                                removeAtom(sim,ii,iBox);
                            }
                        }
                    }
                }
            }
        }
    }
}

/**
 * returns an array whose first element is number of neighbors
 * followed by the ids of those neighbors.
 * the neighbor list is terminated by -1.
 */
int *getNeighborBoxes(SimFlat* s, int iBox)
{
   static int actualNbrs[1+NUMNEIGHBORS];
   int *nbrs;
   int ix, iy, iz;
   int in;

   memset(actualNbrs,-1,1+NUMNEIGHBORS*sizeof(int));
   nbrs = actualNbrs+1;

   if(s->nBoxes == 1)
   {
      nbrs[-1] = 1;
      nbrs[0] = 0;

      return nbrs;
   }

   getIxyz(s, iBox, &ix, &iy, &iz);

   in = 0;

   // we now get ids from (ix-1,iy-1,iz-1) to (ix+1,iy+1,iz+1)
   // which includes the current box as neighbor 13

   for(int i=ix-1; i<=ix+1; i++)
   {
      real3 tcenter;
      int tx, ty, tz;

      for(int j=iy-1; j<=iy+1; j++)
      {
         for(int k=iz-1; k<=iz+1; k++)
         {
               tcenter[0] = s->minLBound[0] + i * s->boxSize[0];
               tcenter[1] = s->minLBound[1] + j * s->boxSize[1];
               tcenter[2] = s->minLBound[2] + k * s->boxSize[2];

	         nbrs[in] = getIBoxFromIxyz(s,i,j,k);

            getIxyz(s, nbrs[in], &tx, &ty, &tz);
/*
            printf("%d iBox %d xyz %d %d %d dc %f %f %f neighbor %d xyz %d %d %d dc %f %f %f \n", 
              getMyParallel(), iBox, ix, iy, iz,
              s->dCenter[iBox][0], s->dCenter[iBox][1], s->dCenter[iBox][2],
              nbrs[in], i, j, k,
              s->dCenter[nbrs[in]][0], s->dCenter[nbrs[in]][1], s->dCenter[nbrs[in]][2]) ;
*/

	         if(nbrs[in] >= 0 ) in++;
         }
      }
   }
   
   nbrs[-1] = in;

   return nbrs;
}
