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
 * This file contains functions for the geometry calls.
 */

#include <stdarg.h>
#include <stdio.h>
#include <strings.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include "pmdTypes.h"
#include "geometry.h"
#include "parallel.h"

#define   NOR       0x01000000

/**
 * compute the real processor number given the (i, j, k) coordinates
 */
int processorNum(SimFlat* s, int i, int j, int k)
{
   int adrx, adry, adrz, proc;

   adrz = (int) getMyParallel() / (s->nproc[0] * s->nproc[1]);
   adry = (int) (getMyParallel() - (adrz * s->nproc[0] * s->nproc[1])) / s->nproc[0];
   adrx = getMyParallel() - (s->nproc[0] * s->nproc[1]) * adrz - s->nproc[0] * adry;

   proc = s->nproc[0] * s->nproc[1] * ((adrz + k + s->nproc[2]) % s->nproc[2]) +
                        s->nproc[0] * ((adry + j + s->nproc[1]) % s->nproc[1]) +
                                      ((adrx + i + s->nproc[0]) % s->nproc[0]);

   proc = NOR | proc;

   return proc;
}

/**
 * initialize processor geometry
 */
int initProcessorGeometry(int xproc, int yproc, int zproc)
{
   int numproc;
   numproc = xproc * yproc * zproc;
   if (numproc != getTotalParallel())
   {
      printf("initProcessorGeometry(). Number of processors must match total size (%d)",
              getTotalParallel());
      exit(1);
   }

   return 0;
}

/**
 * initialize global geometry
 */
void initializeGeometry(SimFlat* s, real_t cutoff)
{
   // # of boxes in each dimension on processor
   for (int i = 0; i < 3; i++)
      s->lNbx[i] = ((s->maxGBound[i] - s->minGBound[i])/cutoff)/s->nproc[i];

   // # of boxes in each dimension system
   for (int i = 0; i < 3; i++)
      s->gNbx[i] = ((s->maxGBound[i] - s->minGBound[i])/cutoff);

   if ((s->lNbx[0] < 2) || (s->lNbx[1] < 2) || (s->lNbx[2] < 2))
   {
      double dx, mincutoff = 1000000.0;
      dx = (s->maxGBound[0] - s->minGBound[0])/(2*s->nproc[0]);
      if (dx < mincutoff) mincutoff = dx;
      dx = (s->maxGBound[1] - s->minGBound[1])/(2*s->nproc[1]);
      if (dx < mincutoff) mincutoff = dx;
      dx = (s->maxGBound[2] - s->minGBound[2])/(2*s->nproc[2]);
      if (dx < mincutoff) mincutoff = dx;
      printf("Must have at least 2x2x2 cells (cutoff < %g)",mincutoff);
      exit(1);
   }

   for (int i = 0; i < 3; i++)
      s->boxSize[i] = (s->maxGBound[i] - s->minGBound[i]) / ((double) s->lNbx[i] * s->nproc[i]),

   subdivideGeometry(s);
   printGeometry(s);

}

/**
 * subdivide geometry
 */
void subdivideGeometry(SimFlat* s)
{
   // calculate x,y,z size on a processor
   for (int i = 0; i < 3; i++)
      s->lSize[i] = (s->maxGBound[i] - s->minGBound[i]) / s->nproc[i];

   // calculate x,y,z size for a box
   for (int i = 0; i < 3; i++) 
      s->boxSize[i] = s->lSize[i] / s->lNbx[i];

   // calculate processor i,j,k (region of space) for this processor
   s->proc[2] = (int) (getMyParallel() / (s->nproc[0] * s->nproc[1]));
   s->proc[1] = (int) ((getMyParallel() - s->proc[2] *s->nproc[0] * s->nproc[1]) / s->nproc[0]);
   s->proc[0] = getMyParallel() - s->proc[2] * s->nproc[0] * s->nproc[1] - s->proc[1] * s->nproc[0];

   // calculate local bounds on this processor
   for (int i = 0; i < 3; i++)
      s->minLBound[i] = s->minGBound[i] + s->proc[i] * s->lSize[i];

   for (int i = 0; i < 3; i++)
      s->maxLBound[i] = s->minGBound[i] + (s->proc[i]+1)* s->lSize[i];

/*
   printf("%d Min Local Bounds : [ %0.17f, %0.17f, %0.17f ]\n",getMyParallel(), s->minLBound[0], s->minLBound[1], s->minLBound[2]);
   printf("%d Max Local Bounds : [ %0.17f, %0.17f, %0.17f ]\n",getMyParallel(), s->maxLBound[0], s->maxLBound[1], s->maxLBound[2]);
*/

   // Determine if this is a boundary processor
   for (int i = 0; i < 3; i++)
   {
      if ((s->proc[i] == 0) || (s->proc[i] == (s->nproc[i] - 1))) 
         s->boundProc[i] = 1;
      else
         s->boundProc[i] = 0;
   }

}

/**
 * update geometry
 */
void updateGeometry(SimFlat* s, real_t xmin, real_t ymin, real_t zmin, 
                                real_t xmax, real_t ymax, real_t zmax)
{
   if ((xmin >= xmax) && (ymin >= ymax) && (zmin >= zmax))
   {
      printf("updateGeometry().  Invalid range : min = (%g,%g,%g), max = (%g,%g,%g)",
              xmin,ymin,zmin,xmax,ymax,zmax);
      exit(1);
   }

   // save new global bounds
   s->minGBound[0] = xmin;
   s->minGBound[1] = ymin;
   s->minGBound[2] = zmin;

   s->maxGBound[0] = xmax;
   s->maxGBound[1] = ymax;
   s->maxGBound[2] = zmax;

   subdivideGeometry(s);
   printGeometry(s);
}

/**
 * print geometry
 */
void printGeometry(SimFlat* s)
{
   if (printParallel())
   {
      printf("\nGlobal Geometry Parameters\n");
      printf("==============================\n");
      printf("  Processors         : (%d,%d,%d)\n", s->nproc[0], s->nproc[1], s->nproc[2]);
      printf("  Min Global Bounds  : [ %0.17f, %0.17f, %0.17f ]\n",s->minGBound[0], s->minGBound[1], s->minGBound[2]);
      printf("  Max Global Bounds  : [ %0.17f, %0.17f, %0.17f ]\n",s->maxGBound[0], s->maxGBound[1], s->maxGBound[2]);
      printf("  Global Boxes       : (%d,%d,%d) = %d\n", s->gNbx[0], s->gNbx[1], s->gNbx[2], s->gNbx[0]*s->gNbx[1]*s->gNbx[2]);
      printf("  Local Boxes        : (%d,%d,%d) = %d\n", s->lNbx[0], s->lNbx[1], s->lNbx[2], s->lNbx[0]*s->lNbx[1]*s->lNbx[2]);
      printf("  Boxsize            : [ %0.17f, %0.17f, %0.17f ]\n", s->boxSize[0], s->boxSize[1], s->boxSize[2]);
   }
}
