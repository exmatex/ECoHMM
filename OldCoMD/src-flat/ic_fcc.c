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

#define DOICTIMERS 0

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



void matInv3x3 (real_t *in, real_t *out)
{
   out[0] = in[8]*in[4] - in[7]*in[5];
   out[1] = in[7]*in[2] - in[8]*in[1];
   out[2] = in[5]*in[1] - in[4]*in[2];
   out[3] = in[6]*in[5] - in[8]*in[3];
   out[4] = in[8]*in[0] - in[6]*in[2];
   out[5] = in[3]*in[2] - in[5]*in[0];
   out[6] = in[7]*in[3] - in[6]*in[4];
   out[7] = in[6]*in[1] - in[7]*in[0];
   out[8] = in[4]*in[0] - in[3]*in[1];

   real_t det = in[0]*(in[8]*in[4] - in[7]*in[5])
      - in[3]*(in[8]*in[1] - in[7]*in[2])
      + in[6]*(in[5]*in[1] - in[4]*in[2]);

   // add the determinant to the tensor 
   in[9] = det;

   real_t invDet = 1.0/det;

   //printf("Determinant, inverse; %g, %g\n", det, invDet);

   for (int m = 0; m<9; m++) 
   {
      out[m] = out[m]*invDet;
   }
   // add the determinant to the inverse
   out[9] = invDet;

   real_t test[9];

   for (int i = 0;i<3;i++) 
   {
      for (int j = 0;j<3;j++) 
      {
         int m = 3*i + j;
         test[m] = 0.0;
         for (int k = 0; k<3; k++) 
         {
            test[m] +=in[3*i + k]*out[3*k + j];
         }
      }

   }
   //printf("Input matrix:\n");
   //printf("[ %g %g %g] \n [%g %g %g] \n [%g %g %g] \n",
   //      in[0], in[1], in[2],
   //      in[3], in[4], in[5],
   //      in[6], in[7], in[8]);

   //printf("Output matrix:\n");
   //printf("[ %g %g %g] \n [%g %g %g] \n [%g %g %g] \n",
   //      out[0], out[1], out[2],
   //      out[3], out[4], out[5],
   //      out[6], out[7], out[8]);

   //printf("Test matrix:\n");
   //printf("[ %g %g %g] \n [%g %g %g] \n [%g %g %g] \n",
   //      test[0], test[1], test[2],
   //      test[3], test[4], test[5],
   //      test[6], test[7], test[8]);
}

real_t ran1(long *idum)
{
   int j;
   long k;
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

real_t gasdev(long *idum)
{
   //real_t ran1(long *idum);
   static int iset=0;
   static real_t gset;
   real_t fac,rsq,v1,v2;

   if  (iset == 0)
   {
      do {
         v1=2.0*ran1(idum)-1.0;
         v2=2.0*ran1(idum)-1.0;
         rsq=v1*v1+v2*v2;
      } while (rsq >= 1.0 || rsq == 0.0);
      fac=sqrt(-2.0*log(rsq)/rsq);
      gset=v1*fac;
      iset=1;
      return v2*fac;
   } else {
      iset=0;
      return gset;
   }
}

// targetKineticEnergy = macroScaleEnergy - potentialEnergy
void applyKineticEnergy(double targetKineticEnergy, SimFlat *s) 
{

   int N = s->nTot;
   real_t v[N][3];
   real_t m = s->pot->mass;

   //printf("In velocity allocation step\n");
   //printf("Mass is %e\n", m);
   //printf("N is %d\n", N);

   long ranseed[3] = {-1.0, -2.0, -3.0};

   // randomize velocities
   for (int i = 0; i < N; i++)
   {
      v[i][0] = gasdev(&ranseed[0]);
      v[i][1] = gasdev(&ranseed[1]);
      v[i][2] = gasdev(&ranseed[2]);
   }
   ////printf("Initial velocity: %e, %e, %e\n", v[N][0], v[N][1], v[N][2]);

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
   //printf("Initial kinetic energy is %e\n", ke);

   // rescale velocities to match target kinetic energy
   double scale = sqrt(targetKineticEnergy / ke);
   //printf("Scale factor = %e\n", scale);
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
   //printf("Final kinetic energy is %e\n", ke);

   // assign the values to the simulation struct at the end
   int i = 0;
   int n_max_in_box = 0;
   int box_with_max = 0;
   for (int iBox=0;iBox<s->nBoxes;iBox++)
   {
      int ii;
      for (int iOff=MAXATOMS*iBox,ii=0; ii<s->nAtoms[iBox]; ii++,iOff++)
      {
         for (int k=0; k<3; k++) s->p[iOff][k] = v[i][k];
         i++;
         if (s->nAtoms[iBox] >= n_max_in_box)
         {
            n_max_in_box = s->nAtoms[iBox];
            box_with_max = iBox;
         }
      }
   }
   //printf("End velocity allocation step\n");

   //printf("Max atom count %d in box %d\n", n_max_in_box, box_with_max);

}

void assignTKE(Command cmd, SimFlat *s) 
{

   // compute potential energy of system
   computeForce(s);
   real_t pe = s->e;
   //printf("Potential energy of the system: %e\n", pe);

   int nx = cmd.nx;
   int ny = cmd.ny;
   int nz = cmd.nz;
   real_t lat = cmd.lat;

   //real_t te = (tke-pe)*nx*ny*nz*lat*lat*lat;
   real_t te = (cmd.tke)*nx*ny*nz*lat*lat*lat;

   real_t ke = te - pe;

   //printf("Target kinetic energy = %e\n", ke);
   // test for zero KE

   ke = 0.0;

   // Assign the velocities here
   applyKineticEnergy(ke, s);

}

SimFlat *createFccLattice(Command cmd, struct BasePotential *pot) 
{
   /**
    * Creates an fcc lattice with nx * ny * nz unit cells and lattice constant lat
    *
    **/
   int nx = cmd.nx;
   int ny = cmd.ny;
   int nz = cmd.nz;
   real_t lat = cmd.lat;
   real_t defGrad = cmd.defGrad;
   real_t boxFactor = cmd.bf;
   real_t tke = cmd.tke;

   SimFlat *s = NULL;
   int     i, j, k, n, itype, nAtoms;
   real_t  x, y, z, halflat;
#ifdef DOICTIMERS  
   clock_t start,old,count=0;
   real2_t fx,fy,fz, px,py,pz;
   fx=fy=fz=px=py=pz=0.0;
#endif

#ifdef DOICTIMERS  
   start = clock();
#endif
   s = blankSimulation(pot);
   if ( ! s ) simulationAbort(-50,(char *) "Unable to create Simulation data structure");

   /* Optional simulation comment (blank for now) */
   s->comment = (char*)suAlignedCalloc(1024*sizeof(char));

   nAtoms = 4*nx*ny*nz;
   halflat = lat / 2.0;

   /* periodic  boundaries  */
   memset(s->bounds,0,sizeof(real3));
   memset(s->boxSize,0,sizeof(real3));
   // stretch in x direction
   s->defGrad = defGrad;
   s->bounds[0] = nx * lat * defGrad;
   s->bounds[1] = ny * lat;
   s->bounds[2] = nz * lat;

   s->bf = boxFactor;

#if DOICTIMERS  >2
   old = clock();
#endif
   allocDomains(s);

#if DOICTIMERS  >2
   count = clock()-old+count;
#endif
   i = j = k = n = 0;
   z = lat / 4.0;
   while (z < s->bounds[2])
   {
      y = lat / 4.0;
      while (y < s->bounds[1])
      {
         x = lat * defGrad / 4.0;
         while (x < s->bounds[0])
         {
            if ((i+j+k) % 2)
            {
               putAtomInBox(s,n++,1,1,s->pot->mass,x,y,z,px,py,pz,fx,fy,fz);
            }
            x += halflat * defGrad;
            i++;
         }
         y += halflat;
         j++;
      }
      z += halflat;
      k++;
   }
   //printf("%d atoms placed in the lattice\n", n);


#ifdef DOICTIMERS  
   old = clock();
#endif
   if(DOICTIMERS) printf("\n    ---- FCC initial condition took %.2gs for %d atoms\n\n",
         (float)(old-start)/(float)(CLOCKS_PER_SEC), n);

   return s;
}



