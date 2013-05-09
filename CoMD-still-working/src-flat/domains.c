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
 * we will have fast neighbor listing **/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include "pmd.h"
#include "memUtils.h"


#define freeMe(s,element) {if(s->element) suAlignedFree(s->element);  s->element = NULL;}


void getBoxIxyz(SimFlat *s, int iBox, int *iret)
{
   /* given a box id return the ix,iy,iz coordinates */
   iret[0] = iBox%s->nbx[0];
   iret[1] = (iBox/s->nbx[0])%s->nbx[1];
   iret[2] = iBox/s->nbx[0]/s->nbx[1];
   return;
}
static int getBoxIDWorldCoords(SimFlat *s, const real2_t xyz[3], int *ibx)
{
   int iBox;
   /* given x,y,z in world co-ordinates return the
      box in which those coordinates fall */

   ibx[0] = (int)(floor(xyz[0]/s->boxSize[0]));
   ibx[1] = (int)(floor(xyz[1]/s->boxSize[1]));
   ibx[2] = (int)(floor(xyz[2]/s->boxSize[2]));

   iBox = ibx[0]+s->nbx[0]*ibx[1]+s->nbx[0]*s->nbx[1]*ibx[2];
   return iBox;
}

void destroyDomains(SimFlat *s)
{
   /* release memory for particles */
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

static int allocDomainArrays(SimFlat *s)
{
   int ierr;
   int i;

   ierr = 0;
   s->dCenter = (real3*)suAlignedCalloc(s->nBoxes*sizeof(real3));
   s->nAtoms = (int*)suAlignedCalloc(s->nBoxes*sizeof(int));
   s->id = (int*)suAlignedCalloc(MAXATOMS*s->nBoxes*sizeof(int));
   s->iType = (int*)suAlignedCalloc(MAXATOMS*s->nBoxes*sizeof(int));
   s->mass = (real_t*)suAlignedCalloc(MAXATOMS*s->nBoxes*sizeof(real_t));
   s->r = (real3*)suAlignedCalloc(MAXATOMS*s->nBoxes*sizeof(real3));
   s->p = (real3*)suAlignedCalloc(MAXATOMS*s->nBoxes*sizeof(real3));
   s->f = (real4*)suAlignedCalloc(MAXATOMS*s->nBoxes*sizeof(real4));
   s->fi = (real_t*)suAlignedCalloc(MAXATOMS*s->nBoxes*sizeof(real_t));
   s->rho = (real_t*)suAlignedCalloc(MAXATOMS*s->nBoxes*sizeof(real_t));

#ifdef USE_IN_SITU_VIZ
   s->centro = (real_t*)suAlignedCalloc(MAXATOMS*s->nBoxes*sizeof(real_t));
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

static int copyDomainArrays(SimFlat *destination, SimFlat* source)
{

   memcpy(destination->dCenter, source->dCenter, source->nBoxes*sizeof(real3)); 
   memcpy(destination->nAtoms, source->nAtoms, source->nBoxes*sizeof(int)); 
   memcpy(destination->id, source->id, MAXATOMS*source->nBoxes*sizeof(int)); 
   memcpy(destination->iType, source->iType, MAXATOMS*source->nBoxes*sizeof(int)); 
   memcpy(destination->mass, source->mass, MAXATOMS*source->nBoxes*sizeof(real_t)); 
   memcpy(destination->r, source->r, MAXATOMS*source->nBoxes*sizeof(real3)); 
   memcpy(destination->p, source->p, MAXATOMS*source->nBoxes*sizeof(real3)); 
   memcpy(destination->f, source->f, MAXATOMS*source->nBoxes*sizeof(real4)); 
   memcpy(destination->fi, source->fi, MAXATOMS*source->nBoxes*sizeof(real_t)); 
   memcpy(destination->rho, source->rho, MAXATOMS*source->nBoxes*sizeof(real_t)); 

   return 0;
}

void allocDomains(SimFlat *s)
{
   if ( ! s ) suAbort(1,(char *) "s unallocated in allocDomains()");
   if ( ! s->pot)  suAbort(2,(char *) "allocDomains() called without s->pot()");

   /* decide how many boxes are needed */
   s->nBoxes = 1;
   real_t strain[3];
   // dummy strain field
   strain[0] = s->defGrad;
   // temporary kluge to make Viz work
   if (strain[0] == 0.0) strain[0]=1.0;
   strain[1] = 1.0;
   strain[2] = 1.0;
   // box factor

   real_t boxFactor = s->bf;
   if (boxFactor == 0.0) boxFactor = 1.0;

   for (int j=0; j<3; j++)
   {
      //printf("bounds = %e, cutoff = %e, box factor = %e, strain = %e\n",
	    //s->bounds[j], s->pot->cutoff, boxFactor, strain[j]);
      s->nbx[j] = (int)floor(s->bounds[j]/(s->pot->cutoff*boxFactor*strain[j]));
      //printf("nbx(%d) = %d\n", j, s->nbx[j]);
      if (s->nbx[j] < 1) suAbort(3,(char *) "Need at least one cutoff wide domain in allocDomains()");
      s->nBoxes = s->nBoxes*s->nbx[j];
      s->boxSize[j] = s->bounds[j]/(real_t)s->nbx[j];
   }

   if(allocDomainArrays(s)) suAbort(4,(char *) "Unable to allocate domains in allocDomains()");

   //printf("____________size is [%d,%d,%d] [%.2f %.2f %.2f]\n",
	 //s->nbx[0],s->nbx[1],s->nbx[2],
	 //s->bounds[0],s->bounds[1],s->bounds[2]
	 //);
   for (int iBox=0;iBox<s->nBoxes; iBox++)
   {
      int ib3[3],j;
      getBoxIxyz(s,iBox,ib3);
      for (int j=0; j<3; j++)
      {
	 s->dCenter[iBox][j] = s->boxSize[j]*(real_t)ib3[j];
      }
   }

   return;
}

void copyDomains(SimFlat *destination, SimFlat *source)
{
   destination->nTot = source->nTot;

   copyDomainArrays(destination, source);
}

void copyAtom(SimFlat *s, int iAtom, int iBox, int jAtom, int jBox)
{
   /* copy atom iAtom in domain iBox to atom jAtom in box jBox */
   const int iOff = MAXATOMS*iBox+iAtom;
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

void moveAtom(SimFlat *s, int iId, int iBox, int jBox)
{
   int nj,ni;
   nj = s->nAtoms[jBox];
   copyAtom(s,iId, iBox, nj, jBox);
   s->nAtoms[jBox]++;
   if(s->nAtoms[jBox]>= MAXATOMS) simulationAbort(-11,"Increase maxatoms");

   s->nAtoms[iBox]--;
   ni = s->nAtoms[iBox];
   if(ni) copyAtom(s,ni,iBox,iId,iBox);

   return;
}

void putAtomInBox(SimFlat *s,
      const int id, const char flagMove, const int iType,
      const real_t mass,
      const real2_t x,const real2_t y,const real2_t z,
      const real2_t px,const real2_t py,const real2_t pz,
      const real2_t fx,const real2_t fy,const real2_t fz)
{

   /**
    * finds an appropriate box for an atom based on the
    * spatial cooridnates and puts it in there.
    *
    * reallocates memory if needed.
    **/
   int iBox;
   int i;
   int iOff;
   int ibx[3];
   real2_t xyz[3] = {x,y,z};

   if ( ! s->stateflag ) suAbort(-10,(char *) "\n\n  ERROR: s not allocated in putAtomInBox()\n\n\n");

   /**
    * push atom into primary period **/
   for (int m=0; m<3; m++)
   {
      if(xyz[m] < 0.0 ) xyz[m] += s->bounds[m];
      else if (xyz[m] >= s->bounds[m] ) xyz[m] -= s->bounds[m];
   }

   /**
    * Find correct box.
    * for now, only one box **/
   iBox = getBoxIDWorldCoords(s,xyz,ibx);
   iOff = iBox*MAXATOMS;
   iOff += s->nAtoms[iBox];

   /**
    * assign values to array elements **/
   s->nTot++;
   s->nAtoms[iBox]++;
   s->id[iOff] = id;
   s->iType[iOff] = iType;
   s->mass[iOff] = mass;
   for (int m=0; m<3; m++)
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
 * static box functions **/
static int getIBoxFromIxyz3(SimFlat *s, int *ixyz)
{
   int iBox=0;

   for (int j=0; j<3; j++)
   {
      if(ixyz[j]<0) ixyz[j] += s->nbx[j];
      else if ( ixyz[j] >= s->nbx[j] ) ixyz[j] -= s->nbx[j];
   }

   iBox = ixyz[0] + ixyz[1]*s->nbx[0] + ixyz[2]*s->nbx[0]*s->nbx[1];

   return iBox;

}
static int getIBoxFromIxyz3NP(SimFlat *s, int *ixyz)
{

   /**
    * for non-periodic boundary conditions, decide whether we
    * are at a boundary or not **/

   if( ! PERIODIC)
   {
      for (int j=0; j<3; j++)
      {
	 if((ixyz[j]<0)||(ixyz[j]>=s->nbx[j])) return -1; 
      }
   }
   return getIBoxFromIxyz3(s,ixyz);

}
static int getIBoxFromIxyz(SimFlat *s, int ix, int iy, int iz)
{
   int iBox;
   int ixyz[3] = {ix,iy,iz};
   return getIBoxFromIxyz3NP(s,ixyz);
}

static void getIxyz3(SimFlat *s, int iBox, int *i3)
{
   i3[0] = (iBox%s->nbx[0]);
   i3[1] = (iBox/s->nbx[0])%(s->nbx[1]);
   i3[2] = (iBox/s->nbx[0]/s->nbx[1])%(s->nbx[2]);
   return;
}
static void getIxyz(SimFlat *s, int iBox, int *ix, int *iy, int *iz)
{
   int ixyz[3];
   getIxyz3(s,iBox,ixyz);
   *ix = ixyz[0];
   *iy = ixyz[1];
   *iz = ixyz[2];
   return;
}

void reBoxAll(SimFlat *s)
{
   for (int iBox=0; iBox<s->nBoxes; iBox++)
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
	 for (int k=0; k<3; k++)
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
	 if((jBox<0)||(jBox==iBox))
	 {
	    /* do nothing if same box or non-periodic boundary */
	    i++;
	    iOff++;
	 }
	 else 
	 {
	    /* move atom to new box */
	    memcpy(s->r[iOff],rNew,sizeof(rNew));
	    moveAtom(s,i,iBox,jBox);
	 }
      }
   }
   return;
}

int *getNeighborBoxes(SimFlat *s, int iBox)
{
   /**
    * returns an array whose first element is number of neighbors
    * followed by the ids of those neighbors.
    * the neighbor list is terminated by -1.
    **/
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

   /* we now get ids from (ix-1,iy-1,iz-1) to (ix+1,iy+1,iz+1)
    * which includes the current box as neighbor 13 */

   for (int i=ix-1; i<=ix+1; i++)
   {
      for (int j=iy-1; j<=iy+1; j++)
      {
	 for (int k=iz-1; k<=iz+1; k++)
	 {
	    nbrs[in] = getIBoxFromIxyz(s,i,j,k);
	    if(nbrs[in] >= 0 ) in++;
	 }
      }
   }

   nbrs[-1] = in;
   return nbrs;
}

