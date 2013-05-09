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
 * File: eamTypes.h
 * Author Sriram Swaminarayan
 * \brief contains the definition of the EAM potential.
 **/

#ifndef __EAM_TYPES_H_
#define __EAM_TYPES_H_

//#include "pmd.h"
#include "mytype.h"

/**
 * EAM uses interpolation.  We store interpolation tables in a struct
 * PotentialArray.  The interpolation is supported on the range
 * \f$[x_0,x_n]\f$, has n values, and each interval has the width
 * \f$1/invDx\f$.
 **/
typedef struct PotentialArray {
  int n;          /**< the number of entries in the array **/
  real_t x0;      /**< the starting ordinate range **/
  real_t xn;      /**< the ending ordinate range **/
  real_t invDx;   /**< the inverse of the interval \f$=n/(x_n-x_0)\f$**/
  real_t *values; /**< the abscissa values **/
} PotentialArray;

/**
 * - All potentials are expected to conform to the
 * following units:
 *   -# atomic distances are in Angstroms
 *   -# atomic masses are in AMUs (atomic mass units)
 *   -# forces are returned in hartrees/angstrom
 *   -# energies are computed in hartrees
 *
 * - EAM potential adds to these:
 *   -# phi, the phi array
 *   -# rho, the rho array
 *   -# f,   the F array.
 *
 * Note that phi, rho, and f are of type PotentialArray.
 *
 **/
typedef struct EamPotential {
  real_t cutoff;       /**< the potential cutoff distance **/
  real_t mass;           /**< mass of atoms **/
  int (*force)(void *s); /**< the force function **/
  void (*destroy)(void **pot); /**< the deallocate function **/
  struct PotentialArray *phi; /**< the phi array **/
  struct PotentialArray *rho; /**< the rho array **/
  struct PotentialArray *f;   /**< the F array   **/
  real_t lat;           /**< lattice constant **/
} EamPotential;

// the Chebychev 'potential' is based on the EAM potential type
typedef struct EamCheby {
  real_t cutoff;       /**< the potential cutoff distance **/
  real_t mass;           /**< mass of atoms **/
  struct PotentialArray *phi; /**< the phi array **/
  struct PotentialArray *rho; /**< the rho array **/
  struct PotentialArray *f;   /**< the F array   **/
  struct PotentialArray *dphi; /**< the phi array **/
  struct PotentialArray *drho; /**< the rho array **/
  struct PotentialArray *df;   /**< the F array   **/
  real_t lat;           /**< lattice constant **/
} EamCheby;


extern struct BasePotential *setEamPot(char *dir, char *file);

#endif
