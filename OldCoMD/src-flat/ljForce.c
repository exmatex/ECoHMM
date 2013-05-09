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
 * The file for force subroutines.
 * Initially we will only have Lennard Jones potential here.
 *
 **/
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include "pmd.h"
#include "domains.h"

#define POT_SHIFT 0

extern int LJ(void *s);
extern int PERIODIC;

void ljDestroy(void **inppot)
{
   LjPotential *pot;
   if ( ! inppot ) return;
   pot = *(LjPotential **)inppot;
   if ( ! pot ) return;
   suAlignedFree(pot);
   *inppot = NULL;
   return;
}


LjPotential *getLjPot()
{
   /* return a new potential with LJ for copper */
   LjPotential *pot;
   pot = (LjPotential*)suAlignedMalloc(sizeof(LjPotential));
   pot->force = LJ;
   pot->destroy = ljDestroy;
   pot->sigma = 1.53;
   pot->epsilon = 0.0085;
   pot->cutoff = 3.0*pot->sigma;
   pot->mass = 1.0;
   return pot;
}

int LJ(void *inS)
{
   /**
    * calculates forces for the 12-6 Lennard Jones potential **/
   /**
    * Notes on LJ:
    *
    * http://en.wikipedia.org/wiki/Lennard_Jones_potential
    *
    * LJ is a simple potential of the form:
    *
    * e_lj(r) = 4*epsilon*((sigma/r)^12 - (sigma/r)^6)
    * F(r) = 4*epsilon*(12*sigma^12/r^13 - 6*sigma^6/r^7)
    *
    * epsilon and sigma are the adjustable parameters in the potential.
    *    epsilon = well depth
    *    sigma   = hard sphere diameter
    *
    * You can also adjust the 12 & 6, but few people do.
    *
    * Some typical values for epsilon and sigma:
    *
    *   material            epsilon            sigma   
    *     N                 36.20 K             3.299 A
    *     O                 44.06 K             2.956 A
    **/
   SimFlat *s = (SimFlat *) inS;
   int *nbrBoxes;
   LjPotential *pot;
   real_t sigma = -1;
   real_t epsilon = -1;
   real_t rCut=-1;
   real_t rCut2 = -1;
   double etot;

   pot = (LjPotential *) s->pot;
   sigma = pot->sigma;
   epsilon = pot->epsilon;
   rCut = pot->cutoff;
   rCut2 = rCut*rCut;

   /**
    * point to energy and force and zero them out **/
   etot = 0.0;
   s->e = (real_t) 0.0;

   /* zero forces / energy */
   memset(s->f,0,s->nBoxes*MAXATOMS*sizeof(real4));

   real_t s6 = sigma*sigma*sigma*sigma*sigma*sigma;

   //virial stress added here
   for (int m=0; m<9; m++) 
   {
      s->stress[m] = 0.0;
   }

#if (LATTICE_DIAG)
   int* numInCutoff = malloc(s->nTot*sizeof(int));
   for (int i=0;i<s->nTot;i++)
   {
      numInCutoff[i] = 0;
   }
#endif

#if POT_SHIFT
   real_t r6cut = 1.0/(rCut2*rCut2*rCut2);
   real_t eShift = r6cut*(s6*r6cut - 1.0);
#else 
   real_t eShift = 0.0;
#endif

   matVec3(s->strain, s->bounds);

   for (int iBox=0; iBox<s->nBoxes; iBox++)
   {/* loop over all boxes in system */
      int nIBox = s->nAtoms[iBox];
      if ( ! nIBox) continue;
      nbrBoxes = getNeighborBoxes(s,iBox);
      for (int jTmp=0; jTmp<nbrBoxes[-1]; jTmp++)
      {/* loop over neighbor boxes */
	 real3 drBox;
	 int jBox = nbrBoxes[jTmp];
	 if(jBox<0) break; /* done with neighbor boxes */
	 for (int j=0; j<3; j++)
	 {
	    drBox[j] = s->dCenter[iBox][j]-s->dCenter[jBox][j];
	    if(PERIODIC)
	    {
	       if(drBox[j]<-0.5*s->bounds[j]) drBox[j] += s->bounds[j];
	       else if (drBox[j] > 0.5*s->bounds[j] ) drBox[j] -= s->bounds[j];
	    }
	 }

	 int nJBox = s->nAtoms[jBox];

	 for (int iOff=iBox*MAXATOMS,ii=0; ii<nIBox; ii++,iOff++)
	 {/* loop over atoms in iBox */
	    // kinetic energy contribution to the virial stress
	    for (int ii=0; ii<3; ii++) 
	    {
	       for (int jj=0; jj<3; jj++) 
	       {
		  int m = 3*ii + jj;
		  s->stress[m] -= s->p[iOff][ii]*s->p[iOff][jj]/s->mass[iOff];
	       }
	    }
	    int i = s->id[iOff];  /* the ij-th atom in iBox */
	    for (int jOff=MAXATOMS*jBox,ij=0; ij<nJBox; ij++,jOff++)
	    {/* loop over atoms in iBox */
	       real_t dr[3];
	       int j = s->id[jOff];  /* the ij-th atom in iBox */
	       if ( j <= i ) continue;
	       real_t r2 = 0.0;

	       for (int m=0; m<3; m++)
	       {
		  dr[m] = drBox[m]+s->r[iOff][m]-s->r[jOff][m];
		  r2+=dr[m]*dr[m];
	       }

	       if ( r2 > rCut2) continue;

	       /**
		* Important note:
		*
		* from this point on r actually refers to 1.0/r
		*
		**/
	       r2 =(real_t)1.0/r2;
	       //r = sqrt(r2);
	       real_t r6 = (r2*r2*r2);
	       //r7 = r6*r;
	       s->f[iOff][3] += 0.5*r6*(s6*r6 - 1.0); 
	       s->f[jOff][3] += 0.5*r6*(s6*r6 - 1.0); 
	       etot += r6*(s6*r6 - 1.0)-eShift; 
	       //f = 4.0*epsilon*s6*r7*(12.0*s6*r6 - 6.0);
	       // different formulation to avoid sqrt computation
	       real_t fr = - 4.0*epsilon*s6*r6*r2*(12.0*s6*r6 - 6.0);

	       for (int m=0; m<3; m++)
	       {
		  s->f[iOff][m] += dr[m]*fr;
		  s->f[jOff][m] -= dr[m]*fr;
	       }
	       for (int ii=0; ii<3; ii++) 
	       {
		  for (int jj=0; jj<3; jj++) 
		  {
		     int m = 3*ii + jj;
		     s->stress[m] -= 1.0*fr*dr[ii]*dr[jj];
		  }
	       }

#if (LATTICE_DIAG >0)
	       numInCutoff[i]+=1;
	       numInCutoff[j]+=1;
#endif

	    } /* loop over atoms in jBox */
	 } /* loop over atoms in iBox */
      } /* loop over neighbor boxes */
   } /* loop over all boxes in system */

#if (LATTICE_DIAG)
   for (int i=0;i<s->nTot;i++)
   {
      //printf("Number of atoms within cutoff: %d, %d\n", i, numInCutoff[i]);
   }

   free(numInCutoff);
#endif

   etot = etot*4.0*epsilon*s6;
   s->e = (real_t) etot;

   matVec3(s->invStrain, s->bounds);

   // renormalize stress
   for (int m=0; m<9; m++) 
   {
      // divide by the undeformed volume
      s->stress[m] = s->stress[m]/(s->bounds[0]*s->bounds[1]*s->bounds[2]);
      // normalize the volume accounting for the strain
      s->stress[m] = s->stress[m]/s->strain[9];
   }

   return 0;
}

