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

#include <stdio.h>
#include "pmd.h"

typedef struct myrl_t {
   int box;
   int id;
} myrl_t;

static myrl_t *getRL(SimFlat *s)
{
   int *iDone = NULL;
   myrl_t *rl = NULL;
   rl = (myrl_t*)calloc(s->nTot,sizeof(myrl_t));
   iDone = (int*)calloc(s->nTot,sizeof(int));
   for (int ib=0; ib<s->nBoxes; ib++)
   {
      for (int iOff=MAXATOMS*ib,i=0; i<s->nAtoms[ib]; i++,iOff++)
      {
         int id = s->id[iOff];
         rl[id].id = i;
         rl[id].box = ib;
         if(iDone[id])
         {
            fprintf(stderr,"_____regot id=%d box %d,atom %d\n",id,ib,i);
         }
         else 
         {
            iDone[id]++;
         }
      }
   }
   for (int i=0; i<s->nTot; i++)
   {
      if ( ! iDone[i])
      {
         fprintf(stderr,"_____didnt get id=%d \n",i);
      }
   }
   free(iDone);
   fflush(stderr);

   return rl;
}

void printIt(SimFlat *s,FILE *fp)
{
   for (int i=0; i<s->nBoxes; i++)
   {
      int *id;

      for (int iOff=i*MAXATOMS,j=0; j<s->nAtoms[i]; j++,iOff++)
      {
         if ( s->id[iOff] < 10)
         {
            fprintf(fp,
                  "%02d, %02d, X=(%+020.12e %+020.12e %+020.12e) 1 P=(%+020.12e %+020.12e %+020.12e) F=(%+020.12e %+020.12e %+020.12e)\n",
                  i,
                  s->id[iOff]+1,
                  s->r[iOff][0],s->r[iOff][1],s->r[iOff][2],
                  s->p[iOff][0],s->p[iOff][1],s->p[iOff][2],
                  s->f[iOff][0],s->f[iOff][1],s->f[iOff][2]
                  );
         }
      }
   }
   return;
}

void writeClsman(SimFlat *s, char *fname)
{
   myrl_t *rl;
   FILE *fp;
   int ti[16];
   double td[16];


   fp = fopen(fname,"wb");

   ti[0] = 3*sizeof(int);/*header*/
   ti[1]=s->nTot;/*nAtoms*/
   ti[2]=0;/*nmove*/
   ti[3]=0;/**/
   ti[4] = 3*sizeof(int);/*header*/
   ti[5]=0;/*header*/
   ti[6]=0;/*header*/
   fwrite(ti,7,sizeof(int),fp);

   rl = getRL(s);

   ti[0] = 3*sizeof(double);/*header*/
   fwrite(ti,1,sizeof(int),fp);

   /* bounds */
   td[0] = (double)s->bounds[0];
   td[1] = (double)s->bounds[1];
   td[2] = (double)s->bounds[2];
   fwrite(td,3,sizeof(double),fp); 

   ti[0] = 3*sizeof(double);/*header*/
   fwrite(ti,1,sizeof(int),fp);

   ti[0] = s->nTot*(6*sizeof(double) + 1*sizeof(int));
   fwrite(ti,1,sizeof(int),fp);

   for (int i=0; i<s->nTot; i++)
   {
      const int jtype=1;
      int id = MAXATOMS*rl[i].box+rl[i].id;
      double r[3],p[3];
      for (int j=0; j<3; j++)
      {
         td[j] = (double)s->dCenter[rl[i].box][j] + (double)s->r[id][j];
         td[j+3] = (double)s->p[id][j];
      }
      fwrite(td,3,sizeof(double),fp);
      fwrite(&jtype,1,sizeof(int),fp);
      fwrite(td+3,3,sizeof(double),fp);
   }
   ti[0] = s->nTot*(6*8 + 1*4);
   fwrite(ti,1,sizeof(int),fp);
   fflush(fp);
   fclose(fp);
   free(rl);
   return;
}

/* gzASCIIWrite Writer written by John Patchett 6/7/2012
 * to match the expectations of fromFileGzip()
 * This writes ascii text to a gzip file using zlib methods
 *                                                           */

#include "zlib.h"
void gzASCIIWrite(SimFlat *s, char *fname)
{
   myrl_t *rl;                  //From the science ... ???
   gzFile fp;                   //File pointer
   double td[16];               //temporary storage for x,y,z, px,py,pz per atom
   char stringToWrite[16348];   //String that will get written to file
   int bytesWritten;           //used to ensure each file write nominally worked

   fp = gzopen(fname,"wb");

   /* 1. number of atoms followed by nmove */
   sprintf(stringToWrite, "%d 0\n", s->nTot);
   bytesWritten = gzputs(fp, stringToWrite);
   if (bytesWritten < 1) exit(-99);

   /* 2. Comment Line                      */
   sprintf(stringToWrite, "%s", s->comment );
   bytesWritten = gzputs(fp, stringToWrite);
   if (bytesWritten < 1) exit(-98);


   /* 3. Write the bounds -- 3 reals on 1 line */
   sprintf(stringToWrite, "%.8lf %.8lf %.8lf\n",
         (double)s->bounds[0], (double)s->bounds[1], (double)s->bounds[2]);
   bytesWritten = gzputs(fp, stringToWrite);
   if (bytesWritten < 1) exit(-97);

   /* 4. for each atom x,y,z,itype,px,py,pz*/
   rl = getRL(s);
   for (int i=0; i<s->nTot; i++)
   {
      const int jtype=1;
      int id = MAXATOMS*rl[i].box+rl[i].id;
      int j;
      double r[3],p[3];
      for (int j=0; j<3; j++)
      {
         td[j] = (double)s->dCenter[rl[i].box][j] + (double)s->r[id][j];
         td[j+3] = (double)s->p[id][j];
      }
      sprintf(stringToWrite, "%.8lf %.8lf %.8lf %d %.8lf %.8lf %.8lf\n",
            td[0], td[1], td[2], jtype, td[3], td[4], td[5]);
      bytesWritten = gzputs(fp, stringToWrite);
      if (bytesWritten < 1) exit(-97);
   }
   gzclose(fp);
}


