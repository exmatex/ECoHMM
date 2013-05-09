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

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <float.h>
#include <time.h>
#include <math.h>

#include "pmd.h"
#include "utility.h"
#include "ic_fcc.h"
#include "geometry.h"
#include "parallel.h"

#define DOICTIMERS 4

// Helper utilities for random numbers: taken from numerical recipes for now

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

real_t ran1(long* idum)
{
   long k;
   int j;
   static long iy=0;
   static long iv[NTAB];
   real_t temp;

   if (*idum <= 0 || !iy) 
   {
      if (-(*idum) < 1) *idum=1;
      else *idum = -(*idum);
      for (j=NTAB+7;j>=0;j--) 
      {
         k=(*idum)/IQ;
         *idum=IA*(*idum-k*IQ)-IR*k;
         if (*idum < 0) *idum += IM;
         if (j < NTAB) iv[j] = *idum;
      }
      iy=iv[0];
   }
   k=(*idum)/IQ;
   *idum=IA*(*idum-k*IQ)-IR*k;
   if (*idum < 0) *idum += IM;
   j=iy/NDIV;
   iy=iv[j];
   iv[j] = *idum;
   if ((temp=AM*iy) > RNMX) return RNMX;
   else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

real_t gasdev(long* idum)
{
   //real_t ran1(long* idum);
   static int iset=0;
   static real_t gset;
   real_t fac,rsq,v1,v2;

   if  (iset == 0) 
   {
      do 
      {
         v1=2.0*ran1(idum)-1.0;
         v2=2.0*ran1(idum)-1.0;
         rsq=v1*v1+v2*v2;
      } while (rsq >= 1.0 || rsq == 0.0);
      fac=sqrt(-2.0*log(rsq)/rsq);
      gset=v1*fac;
      iset=1;
      return v2*fac;
   } 
   else 
   {
      iset=0;
      return gset;
   }
}

// targetKineticEnergy = macroScaleEnergy - potentialEnergy
void applyKineticEnergy(double targetKineticEnergy, SimFlat* s) 
{
   int N = s->nTot;
   real_t v[N][3];
   real_t m = s->pot->mass;

   if (printParallel())
   {
      printf("In velocity allocation step\n");
      printf("Mass is %e\n", m);
      printf("N is %d\n", N);
   }

   long ranseed[3] = {-1.0, -2.0, -3.0};

   // randomize velocities
   for (int i = 0; i < N; i++)
   {
       v[i][0] = gasdev(&ranseed[0]);
       v[i][1] = gasdev(&ranseed[1]);
       v[i][2] = gasdev(&ranseed[2]);
   }
   //printf("Initial velocity: %e, %e, %e\n", v[N][0], v[N][1], v[N][2]);

   // calculate net momentum
   double netPx = 0.0;
   double netPy = 0.0;
   double netPz = 0.0;
   for (int i = 0; i < N; i++)
   {
       netPx += v[i][0] * m;
       netPy += v[i][1] * m;
       netPz += v[i][2] * m;
   }
   //printf("Net momenta: %e, %e, %e\n", netPx, netPy, netPz);

   // shift velocities so that net momentum is zero
   for (int i = 0; i < N; i++)
   {
      v[i][0] -= netPx * (1.0 / (m * N));
      v[i][1] -= netPy * (1.0 / (m * N));
      v[i][2] -= netPz * (1.0 / (m * N));
   }

   // calculate existing kinetic energy
   double ke = 0.0;
   for (int i = 0; i < N; i++)
   {
      ke += 0.5 * m * (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
   }
   if (printParallel())
      printf("Initial kinetic energy is %e\n", ke);

   // rescale velocities to match target kinetic energy
   double scale = sqrt(targetKineticEnergy / ke);
   if (printParallel())
      printf("Scale factor = %e\n", scale);
   for (int i = 0; i < N; i++)
   {
      v[i][0] *= scale;
      v[i][1] *= scale;
      v[i][2] *= scale;
   }
   // calculate existing kinetic energy
   ke = 0.0;
   for (int i = 0; i < N; i++)
   {
      ke += 0.5 * m * (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
   }
   if (printParallel())
      printf("Final kinetic energy is %e\n", ke);

   // assign the values to the simulation struct at the end
   int i = 0;
   int n_max_in_box = 0;
   int box_with_max = 0;
   for (int iBox=0;iBox<s->nBoxes;iBox++)
   {
      for(int iOff=MAXATOMS*iBox,ii=0; ii<s->nAtoms[iBox]; ii++,iOff++)
      {
         for(int k=0; k<3; k++) s->p[iOff][k] = v[i][k];
         i++;
         if (s->nAtoms[iBox] >= n_max_in_box)
         {
            n_max_in_box = s->nAtoms[iBox];
            box_with_max = iBox;
         }
      }
   }
   if (printParallel()) 
   {
      printf("End velocity allocation step\n");

      printf("On processor %d max atom count %d in box %d\n", getMyParallel(), n_max_in_box, box_with_max);
   }

}

void assignTKE(Command cmd, SimFlat* s) 
{

   // compute potential energy of system
   computeForce(s);
   real_t pe = s->e;

   int nx = cmd.nx;
   int ny = cmd.ny;
   int nz = cmd.nz;
   real_t lat = cmd.lat;

   //real_t te = (tke-pe)*nx*ny*nz*lat*lat*lat;
   real_t te = (cmd.tke)*nx*ny*nz*lat*lat*lat;

   real_t ke = te - pe;

   if (printParallel())
      printf("Target kinetic energy = %e\n", ke);
   // test for zero KE

   ke = 0.0;

   // Assign the velocities here
   applyKineticEnergy(ke, s);

}

/**
 * Creates an fcc lattice with nx * ny * nz unit cells and lattice constant lat
 *
 */
SimFlat* createFccLattice(Command cmd, struct BasePotential* pot)
{
   int nx = cmd.nx;
   int ny = cmd.ny;
   int nz = cmd.nz;
   real_t lat = cmd.lat;
   real_t defGrad = cmd.defGrad;
   real_t boxFactor = cmd.bf;
   real_t tke = cmd.tke;

   SimFlat* s = NULL;
   int i, j, k, n, itype, natoms;
   real_t  x, y, z, halflat;
   real_t x0, y0, z0;

#ifdef DOICTIMERS  
   clock_t start,old,count=0;
   real2_t fx,fy,fz, px,py,pz;
   fx=fy=fz=px=py=pz=0.0;
#endif

#ifdef DOICTIMERS  
   start = clock();
#endif

   if (printParallel()) 
   {
      printf("In createFccLattice:");
      printf("Potential cutoff %e\n", pot->cutoff);
   }

   s = blankSimulation(pot);

   s->pot = pot;

   if ( ! s ) simulationAbort(-50,(char *) "Unable to create Simulation data structure");

   if (printParallel())
   {
      printf("In createFccLattice:");
      printf("Potential cutoff %e\n", s->pot->cutoff);
   }

   // Optional simulation comment (blank for now)
   s->comment = (char*)suAlignedCalloc(1024*sizeof(char));

   natoms = 4*nx*ny*nz;
   halflat = lat / 2.0;

   // zero variables
   memset(s->bounds,0,sizeof(real3));
   memset(s->minLBound,0,sizeof(real3));
   memset(s->maxLBound,0,sizeof(real3));
   memset(s->minGBound,0,sizeof(real3));
   memset(s->maxGBound,0,sizeof(real3));
   memset(s->boxSize,0,sizeof(real3));
   memset(s->lSize,0,sizeof(real3));

   // stretch in x direction
   s->defGrad = defGrad;

   // initialialize global periodic bounds
   s->minGBound[0] = 0;
   s->minGBound[1] = 0;
   s->minGBound[2] = 0;

   s->maxGBound[0] = nx * lat * defGrad;
   s->maxGBound[1] = ny * lat;
   s->maxGBound[2] = nz * lat;

   // global bounds
   s->bounds[0] = s->maxGBound[0] - s->minGBound[0];
   s->bounds[1] = s->maxGBound[1] - s->minGBound[1];
   s->bounds[2] = s->maxGBound[2] - s->minGBound[2];

   // set proc counts in all directions
   memset(s->nproc,0,sizeof(real3));
   s->nproc[0] = cmd.xproc;
   s->nproc[1] = cmd.yproc;
   s->nproc[2] = cmd.zproc;

   if (printParallel())
   {
      printf("Potential cutoff %e\n", s->pot->cutoff);
   }

   // initialize geometry
   initializeGeometry(s, pot->cutoff);

   if (printParallel())
   {
      printf("Potential cutoff %e\n", s->pot->cutoff);
   }

#if DOICTIMERS  >2
   old = clock();
#endif

   s->pot = pot;

   // create boxes
   allocDomains(s);

#if DOICTIMERS  >2
   count = clock()-old+count;
#endif

   // create and place atoms

   // Skip ahead until we get to the box assigned to this processor
   x = y = z = 0.0;
   x0 = y0 = z0 = 0.0;
   i = j = k = n = 0;

   z = lat / 4.0;
   y = lat / 4.0;
   x = lat * defGrad / 4.0;
   while (x < s->minLBound[0])
   {
      x += halflat * defGrad;
      i++;
   }
   while (y < s->minLBound[1])
   {
      y += halflat;
      j++;
   }
   while (z < s->minLBound[2])
   {
      z += halflat;
      k++;
   }
   z0 = z; y0 = y; x0 = x;

   //printf( "x = %f y = %f z = %f i = %d j = %d k = %d n = %d\n", x, y, z, i, j, k, n);

   // create atoms on this processor
  
   z = z0;
   while (z < s->maxLBound[2])
   {
      y = y0;
      while (y < s->maxLBound[1])
      {
         x = x0;
         while (x < s->maxLBound[0])
         {
            if ((i+j+k) % 2)
               putAtomInBox(s,n++,1,1,s->pot->mass,x,y,z,px,py,pz,fx,fy,fz);
	         x += halflat * defGrad;
	         i++;
         }
         y += halflat;
         j++;
      }
      z += halflat;
      k++;
   }

   // set total atoms in simulation
   s->nTotalAtoms = addIntParallel(s->nTot);

#ifdef DOICTIMERS  
   old = clock();
#endif

   if(DOICTIMERS)
      if (printParallel())
         printf("\n    ---- FCC initial condition took %.2gs for %d atoms\n\n",
	         (float)(old-start)/(float)(CLOCKS_PER_SEC), s->nTotalAtoms);

   return s;
}
