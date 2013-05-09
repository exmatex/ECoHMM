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

#ifndef __PMDTYPES_H_
#define __PMDTYPES_H_

#include "mytype.h"
#include "memUtils.h"
#include "eamTypes.h"

#define MAXATOMS 64 
#define SIM_NOSTATE 0
#define SIM_ALLOCED 1


/**
 * \brief the base form off of which all potentials will be set
 *
 * - All potentials are expected to conform to the
 * following units:
 *   -# atomic distances are in Angstroms
 *   -# atomic masses are in AMUs (atomic mass units)
 *   -# forces are returned in hartrees/angstrom
 *   -# energies are computed in hartrees
 * Note that phi, rho, and f are of type PotentialArray.
 *
 */
typedef struct BasePotential {
   real_t cutoff;         /**< potential cutoff distance in Angstroms **/
   real_t mass;           /**< mass of atoms in atomic mass units **/
   int (*force)(void *s); /**< the actual parameter is "struct simulation_t *s" **/
   void (*destroy)(void *pot); /**< destruction of the potential **/
} BasePotential;



#ifdef USE_IN_SITU_VIZ
typedef struct {   /**< Neighbor info used for centrosymmetry parameter **/
   real3    r;
   real_t   rsq;
   int      index;
} Neighbor;
typedef Neighbor Neighbor12[12];
#endif

typedef struct fileAtom {
   float  x, y, z;       
   float  bond_order;
   float  centrosymmetry;
} FileAtom;
  

/** def array offset by n boxes
 */
#define arrayAtBox(array,boxid) ((array)+MAXATOMS*boxid)

/**
 * The basic flat simulation data structure with MAXATOMS in every box
 */
typedef struct SimFlat {
   int stateflag; //!< unused for now

   // processor-specific data
   int nproc[3];     //!< number of processors in each dimension
   int proc[3];      //!< i,j,k for this processor
   int boundProc[3]; //!< boundary processor flag in each dimension

   // global bounds data
   real3 minGBound; //!< minimum global bounds
   real3 maxGBound; //!< maximum global bounds
   real3 bounds;    //!< periodic bounds - global periodic size
                    //!< change to gSize - maxGBound - minGBound

   // local bounds data
   real3 minLBound;  //!< minimum local bounds on processor
   real3 maxLBound;  //!< maximum local bounds on processor
   real3 lSize;      //!< size on one processor in each dimension
                     //!< maxLBound - minLBound

   // box data
   int gNbx[3];      //!< number of boxes in each dimension for simulation
   int lNbx[3];      //!< number of boxes in each dimension on processor
   int nBoxes;       //!< total number of local boxes on processor
   int nHaloBoxes;   //!< total number of remote halo/ghost boxes on processor
   int nTotalBoxes;  //!< total number of boxes on processor
                     //!< nboxes + nghostboxes
   real3 boxSize;    //!< size of box in each dimension
   real3* dCenter;   //!< box center n each dimension

   // atom-specific data
   int nTot;      //!< total number of atoms on this processor/
   int* nAtoms;   //!< total number of atoms in each box
   int* id;       //!< The original ID of the atom
   int* iType;    //!< the type of atoms
   real_t* mass;  //!< mass of the atoms
   real3* r;      //!< positions
   real3* p;      //!< momenta of atoms
   real4* f;      //!< fx, fy, fz, energy
   real_t* rho;   //!< rhosum for EAM potential
   real_t* fi;    //!< rhobar for EAM potential

   // halo-specific data
   int haloSizes[3];   //!< number of boxes in each dimension for halo
   int haloOffsets[3]; //!< offsets in each dimension for halo
   
#ifdef USE_IN_SITU_VIZ
   real_t* centro;
#endif

   // the total potential energy of the simulation
   real_t e;        //!< the total energy of the system
   int nTotalAtoms;  //!< The total number of atoms in the simulation

   BasePotential *pot; /**< the potential**/

   char* comment; //!< free form string that describes the simulation

   // new variables to store meso-micro data. Will be tensors eventually
   real_t stress;      //!< virial stress
   real_t defGrad;     //!< deformation gradient
   real_t bf;          //!< box factor
   EamCheby* chPot; //!< Chebychev coefficients

} SimFlat;

/**
  add-ons for visualization
*/

#ifdef USE_IN_SITU_VIZ
#include "AtomVisualize.h"

#else
typedef void AtomVisualize;
#endif

typedef struct simThreadData
{
   AtomVisualize* viz;
   SimFlat* sim;
   int niter;
   int nsteps;
   int eamFlag;
   int gpuFlag;
} SimThreadData;

#endif
