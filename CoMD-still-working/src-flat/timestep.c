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
 * The file for timestep subroutines.
 * Initially we will only have verlet
 *
 **/
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include "pmd.h"

static void advanceVelocity(SimFlat *s, const int nStart, const int nTodo, real_t dt)
{
   for (int iBox=nStart; iBox<nStart+nTodo; iBox++)
   {
      /* loop over my boxes in system */
      for (int iOff=MAXATOMS*iBox,ii=0; ii<s->nAtoms[iBox]; ii++,iOff++)
      {
         for (int k=0; k<3; k++) s->p[iOff][k] -= dt*s->f[iOff][k];
      }
   }
}

static void advancePosition(SimFlat *s, const int nStart, const int nTodo, real_t dt)
{
   real_t rMass;  

   /**
    * All computations are done in atomic units.
    * (http://en.wikipedia.org/wiki/Atomic_units)
    * In these units the mass of an electrion is
    * assumed to be 1.0.  Since our mass is in
    * amu, we need to multiply it by 1822.83 to
    * get the mass in amu.
    **/

   /* convert mass to atomic units */
   rMass = (real_t)(amu_to_m_e*(double)s->pot->mass);

   for (int iBox=nStart; iBox<nStart+nTodo; iBox++)
   {
      /* loop over my boxes in system */
      if ( ! s->nAtoms[iBox]) continue;
      for (int iOff=MAXATOMS*iBox,ii=0; ii<s->nAtoms[iBox]; ii++,iOff++)
      {
         /* loop over atoms in iBox */
         for (int k=0; k<3; k++)
         {
            s->r[iOff][k] += dt*s->p[iOff][k]/rMass;
         }
      }
   }
   /* move atoms to primary period */
   reverseDeformation(s);
   reBoxAll(s);
   forwardDeformation(s);

}

int computeForce(SimFlat *s)
{
   return s->pot->force(s);
}

double nTimeSteps(int n, SimFlat *s, real_t dt)
{
   extern void printIt(SimFlat *sim,FILE *fp);
   /**
    * Standard verlet algorithm:
    *   1: advance positions half a timestep using current velocities
    *   2: compute forces
    *   3: advance velocities (momenta) a full timestep
    *   4: advance positions half a timestep to bring in sync with velocities.
    **/

   /* convert dt to atomic units */
   dt = dt * bohr_per_atu_to_A_per_s;

   for (int i=0; i<n; i++)
   {
      advancePosition(s,0,s->nBoxes,(dt/2.0));
      computeForce(s);
      advanceVelocity(s,0,s->nBoxes,dt); 
      advancePosition(s,0,s->nBoxes,(dt/2.0));
   }

   /* compute force to make consistent */
   computeForce(s);

   return s->e;

}


