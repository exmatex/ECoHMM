#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include "pmd.h"
#include "memUtils.h"


// form the matrix-vector product in place
void matVec3 (real_t *mat, real_t *vec)
{
   real_t prod[3];

   prod[0] = mat[0]*vec[0] + mat[1]*vec[1] + mat[2]*vec[2];
   prod[1] = mat[3]*vec[0] + mat[4]*vec[1] + mat[5]*vec[2];
   prod[2] = mat[6]*vec[0] + mat[7]*vec[1] + mat[8]*vec[2];

   vec[0] = prod[0];
   vec[1] = prod[1];
   vec[2] = prod[2];
}

void forwardDeformation(SimFlat *s)
{
   for (int iBox=0; iBox < s->nBoxes; iBox++)
   {
      // apply deformation to all box centers
      real_t boxCenter[3];
      for (int m = 0; m<3; m++)
      {
	 boxCenter[m] = s->dCenter[iBox][m];
      }
      matVec3(s->strain, boxCenter);
      for (int m = 0; m<3; m++)
      {
	 s->dCenter[iBox][m] = boxCenter[m];
      }
      // apply deformation to all particles
      for (int iOff=MAXATOMS*iBox,ii=0; ii<s->nAtoms[iBox]; ii++,iOff++)
      {
	 /* loop over atoms in iBox */
	 real_t atomPos[3];
	 for (int m=0; m<3; m++)
	 {
	    atomPos[m] = s->r[iOff][m];
	 }
	 matVec3(s->strain, atomPos);
	 for (int m=0; m<3; m++)
	 {
	    s->r[iOff][m] = atomPos[m];
	 }
      }
   }

}

void reverseDeformation(SimFlat *s)
{
   for (int iBox=0; iBox < s->nBoxes; iBox++)
   {
      // apply deformation to all box centers
      real_t boxCenter[3];
      for (int m = 0; m<3; m++)
      {
	 boxCenter[m] = s->dCenter[iBox][m];
      }
      matVec3(s->invStrain, boxCenter);
      for (int m = 0; m<3; m++)
      {
	 s->dCenter[iBox][m] = boxCenter[m];
      }
      // apply deformation to all particles
      for (int iOff=MAXATOMS*iBox,ii=0; ii<s->nAtoms[iBox]; ii++,iOff++)
      {
	 /* loop over atoms in iBox */
	 real_t atomPos[3];
	 for (int m=0; m<3; m++)
	 {
	    atomPos[m] = s->r[iOff][m];
	 }
	 matVec3(s->invStrain, atomPos);
	 for (int m=0; m<3; m++)
	 {
	    s->r[iOff][m] = atomPos[m];
	 }
      }
   }

}
