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
 * The include file for domains.c and
 * contains function declarations
 */

#ifndef __BOXIT_H_
#define __BOXIT_H_

#include "pmdTypes.h"

#define NUMNEIGHBORS 27

extern void zeroForcesInBox(SimFlat* s, int iBox);

extern void allocDomains(SimFlat* s);
extern void copyDomains(SimFlat* destination, SimFlat* source);
extern int *getNeighborBoxes(SimFlat* s, int iBox);
extern void putAtomInBox(SimFlat* s,
			 const int id, const char flagMove, const int iType,
			 const real_t mass,
			 const real2_t x,const real2_t y,const real2_t z,
			 const real2_t px,const real2_t py,const real2_t pz,
			 const real2_t fx,const real2_t fy,const real2_t fz);
extern void reBoxAll(SimFlat* s);
extern void reBoxAll2(SimFlat* s);
extern int getIBoxFromIxyz3(SimFlat* s, int* xyz);
extern int getIBoxFromIxyz(SimFlat* s, int x, int y, int z);
extern int getIBoxFromIxyz3NP(SimFlat* s, int* xyz);
extern void getIxyz3(SimFlat* s, int iBox, int* i3 );
extern void moveAtom(SimFlat* s, int j, int iBox, int jBox);
extern void removeAtom(SimFlat* s, int j, int iBox);

#endif
