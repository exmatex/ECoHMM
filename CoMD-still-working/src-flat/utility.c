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
 * The data for our simulation will be held in a
 * structure of arrays.  To make it go fast on the
 * CBE we have decided to separate out the three
 * positions and three momenta into separate
 * arrays.  This will allow for more overlap
 * between memory access and computation.
 **/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include <sys/time.h>

#include "pmdTypes.h"
#include "pmd.h"
#include "memUtils.h"
#include "utility.h"

#define offsetPtr(ptr,offset) (((char *)ptr) + offset)

/* some utility routines */
double timeNow()
{
   struct timeval tnow;
   double t;
   gettimeofday(&tnow,NULL);
   t = (double) tnow.tv_sec + ((double)tnow.tv_usec)*1.0e-6;
   return t;
}

void breakinme()
{
   /** a blank break**/
   return;
}

int suAbort(int status,char *message)
{
   fprintf(stderr,"\n\n\n   ABORT!: \n   %s\n\n\n",message);
   breakinme();				\
      exit(status);
}

char *dupString(char *str)
{
   char *s;
   s = (char*)suAlignedCalloc((strlen(str)+1)*sizeof(char));
   strcpy(s,str);
   return s;
}

void simulationAbort(int ineCode, char *inmsg)
{
   int eCode;					
   char *msg;					
   eCode = ineCode;				
   msg = inmsg;				
   if(! msg)
   {
      msg = (char *) "Unknown error";		
   }					
   fprintf(stderr,"\n\n"			
         "    Error code: %d\n"	
         "    %s\n\n",			
         eCode, msg);			
   exit(eCode);				
}

struct SimFlat *blankSimulation(struct BasePotential *pot)
{
   SimFlat *s;

   s = (SimFlat*)suAlignedMalloc((1*sizeof(SimFlat)));

   if(pot)
   {
      if ( pot->cutoff < 0.0) simulationAbort(-1,(char *) "Blanksimulation() got a negative cutofff");
      s->pot = pot;
   }
   else 
   {
      extern LjPotential *getLjPot();
      s->pot = (BasePotential *) getLjPot();
   }

   //printf("In blankSimulation: ");
   //printf("Potential cutoff %e\n", s->pot->cutoff);

   s->stateflag = 0;

   return s;
}


void destroySimulation(SimFlat **ps)
{
   /**
    * frees all data associated with *s and frees *s **/
   extern void destroyDomains(SimFlat *s);
   SimFlat *s;
   BasePotential *pot;

   if ( ! ps ) return;

   s = *ps;

   if ( ! s ) return;
   pot = s->pot;
   if ( pot) pot->destroy(&pot);
   destroyDomains(s);
   suAlignedFree(s);
   *ps = NULL;
   return;
}

