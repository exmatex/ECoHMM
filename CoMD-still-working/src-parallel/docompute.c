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

#include "pmd.h"
#include "cheby.h"
#include "ic_fcc.h"
#include "docompute.h"
#include "geometry.h"
#include "read.h"
#include "yamlOutput.h"
#include "mycommand.h"
#include "parallel.h"
#include "redistribute.h"

#define NUM_CHEBY 16

void printArray(real_t* array, int n, char *name)
{
   printf("%s;\n", name);
   for (int i = 0; i<n; i++)
   {
      printf("%d, %17.9e\n", i, array[i]);
   }
}

void printSimulationDataYaml(FILE* file,
      SimFlat* sim,
      Command cmd,
      BasePotential* pot)
{
   if (printParallel())
   {
      fprintf(file, "Simulation data:\n");
      if (cmd.doeam)
      {  
	      EamPotential *eam_pot = (EamPotential*) pot;
	      fprintf(file, "  EAM potential values:\n");
	      fprintf(file, "    cutoff: "EMT1"\n", eam_pot->cutoff);
	      fprintf(file, "    mass: "EMT1"\n", eam_pot->mass);
	      fprintf(file, "  phi potential:\n");
	      fprintf(file, "    n: %d\n", eam_pot->phi->n);
	      fprintf(file, "    x0: "EMT1"\n", eam_pot->phi->x0);
	      fprintf(file, "    xn: "EMT1"\n", eam_pot->phi->xn);
#if (USE_CHEBY)
	      fprintf(file, "  Used %d Chebyshev coefficients\n", NUM_CHEBY);
#endif
      }
      else
      {
	      LjPotential* lj_pot = (LjPotential*) pot;
	      fprintf(file, "  LJ potential values:\n");
	      fprintf(file, "    cutoff: "EMT1"\n", lj_pot->cutoff);
	      fprintf(file, "    mass: "EMT1"\n", lj_pot->mass);
      }
      fprintf(file, "  nAtoms: %d\n", sim->nTot);
   }
}

SimFlat* initSimFromCmdLine(int argc, char** argv) {
   Command cmd;
   BasePotential* pot;
   SimFlat* sim;
   SimThreadData* threadData;

   int hbox;

   if (printParallel())
   {
      printf("    Double Precision=%s\n", (sizeof(real_t)==4?"false":"true") );
      fprintf(yamlFile,"Double Precision: %s\n", (sizeof(real_t)==4?"false":"true") );
   }

   // get command line params
   parseCommandLine(&cmd, argc, argv);

   printCmd(&cmd);

   // initialize processor geometry
   initProcessorGeometry(cmd.xproc, cmd.yproc, cmd.zproc);

   // decide whether to get LJ or EAM potentials
   if(cmd.doeam) 
   {
      if (printParallel()) printf("Setting EAM potential\n");
      pot = setEamPot(cmd.potDir, cmd.potName);
   } 
   else 
   {
      if (printParallel()) printf("Setting LJ potential\n");
      pot = (BasePotential*) getLjPot();
   }

   if ( ! pot ) simulationAbort(-2,(char *) "Unable to initialize potential");

   if (printParallel())
      printf("Simulation data\n");
   if (cmd.doeam)
   {
      EamPotential* newPot;
      newPot = (EamPotential*) pot;
      if (printParallel())
      {
	      printf("EAM potential values:\n");
	      printf("cutoff = %e\n", newPot->cutoff);
	      printf("mass = %e\n", newPot->mass);
	      printf("phi potential:\n");
	      printf("n = %d\n", newPot->phi->n);
      }
      cmd.lat = newPot->lat;
      if (printParallel())
      {
	      printf("cmd.lat = "EMT1"\n", cmd.lat);
#if (DIAG_LEVEL > 1)
	      printf("x0 = "EMT1"\n", newPot->phi->x0);
	      printf("xn = "EMT1"\n", newPot->phi->xn);
	      printArray(newPot->phi->values, newPot->phi->n, "phi ref");
	      printArray(newPot->rho->values, newPot->rho->n, "rho ref");
	      printArray(newPot->f->values, newPot->f->n, "f ref");
#endif
      }
   } 
   else
   {
      LjPotential* newPot;
	   newPot = (LjPotential*) pot;
      if (printParallel())
      {
	      printf("LJ potential values:\n");
	      printf("cutoff = "EMT1"\n", newPot->cutoff);
	      printf("mass = "EMT1"\n", newPot->mass);
      }

      // Base lattice places nearest neighbor at zero-force distance
      //cmd.lat = 1.587401052*newPot->sigma ;      // This is for Lennard-Jones = 2^(2/3)*sigma
      // Empirically determined prefactor to minimize virial stress
      cmd.lat = 0.97415*1.587401052*newPot->sigma ;      // This is for Lennard-Jones = 2^(2/3)
   }

   if (printParallel()) 
   {
      printf("Before sim initialization: ");
      printf("Potential cutoff %e\n", pot->cutoff);
   }

   if(strcmp(cmd.filename, ""))
   {    
      // Read in file
      if (printParallel())
	      printf("\n\n    File is __%s__\n\n", cmd.filename);
      sim = fromFileASCII(cmd, pot);
      if ( ! sim ) simulationAbort(-3,(char *) "Input file does not exist");
   } 
   else
   {
      if (printParallel())
	      printf("cmd.lat = %e\n", cmd.lat);
      sim = createFccLattice(cmd, pot);
   }

   if (cmd.doeam)
   {
      EamPotential *newPot = (EamPotential*) pot;
      printf("%e, %e\n", newPot->phi->x0, newPot->phi->xn);
      sim->chPot = setChebPot(newPot, NUM_CHEBY);
#if (DIAG_LEVEL > 0)
      if (printParallel()) 
      {
	      printf("Chebychev coefficients:\n");
	      fflush(stdout);
	      printf("%d, %d, %d\n",
	         sim->chPot->phi->n,
	         sim->chPot->rho->n,
	         sim->chPot->f->n);
	      fflush(stdout);

	      printArray(sim->chPot->phi->values, sim->chPot->phi->n, "phi");
	      printArray(sim->chPot->rho->values, sim->chPot->rho->n, "rho");
	      printArray(sim->chPot->f->values,   sim->chPot->f->n,   "f");
	      printArray(sim->chPot->dphi->values, sim->chPot->dphi->n, "dphi");
	      printArray(sim->chPot->drho->values, sim->chPot->drho->n, "drho");
	      printArray(sim->chPot->df->values,   sim->chPot->df->n,   "df");
      }
#endif
   }

   if (printParallel())
   {
      printf("    total atoms is: %d\n", sim->nTotalAtoms);
      printf("box factor is ("EMT1", "EMT1", "EMT1")\n", 
	      sim->boxSize[0]/pot->cutoff,
	      sim->boxSize[1]/pot->cutoff,
	      sim->boxSize[2]/pot->cutoff);
   }

   // initial output for consistency check
   redistribute(sim);
   //reBoxAll(sim);

   // apply the desired kinetic energy to the system
   assignTKE(cmd, sim);

   if (printParallel())
   {
      printCmdYaml(&cmd, yamlFile);
      printSimulationDataYaml(yamlFile, sim, cmd, pot);
      printf("Initialization finished\n");
   }

   return sim;
}

void* doComputeWork(SimThreadData* data)
{
   double eone;
   double ts, te;
   int iter;
   int niter = 20;
   //int niter = 5;
   int nsteps = 10;
   //int nsteps = 1;
   double dt;
   real_t etotal, estress;

#if (DIAG_LEVEL > 1)
   if (printParallel())
   {
      printf("Chebychev coefficients:\n");
      fflush(stdout);
      printf("%d, %d, %d\n",
	      data->sim->chPot->phi->n,
	      data->sim->chPot->rho->n,
	      data->sim->chPot->f->n);
      fflush(stdout);
   }
#endif

   if (printParallel()) 
   {
      printf("Starting simulation\n");
      fflush(stdout);
   }

   ts = timeNow();
   computeForce(data->sim); 
   te = timeNow();

#if (DIAG_LEVEL > 1)
   if (printParallel())
      printIt(data->sim, stdout);
#endif

   // convert the timestep
   dt = 1.0e-15;

   if (printParallel())
   {
      printf("Starting timesteps\n");
      fflush(stdout);
   }

   for(iter=0; iter<niter; iter++)
   {
      double ns;
      ns = (iter==0?1.0:(double)(1+nsteps));
      /** note the 'nsteps+1' below.
       * this is because we have an extra
       * force call at the end of timesteps()
       */

      // sum atoms and energy across all processers
      data->sim->nTot = 0;
      for (int i = 0; i < data->sim->nBoxes; i++)
      {
	      data->sim->nTot += data->sim->nAtoms[i];
      }
      data->sim->nTotalAtoms = addIntParallel(data->sim->nTot);
      etotal = addRealParallel(data->sim->e);

      // sum stress across all processors
      estress = addRealParallel(data->sim->stress);

      if (printParallel())
      {
	      printf(" %3d %30.20g computed in %.3fs (%8.4f us/atom for %d atoms)\n",
	         iter*nsteps, etotal, (te-ts), 1.0e6*(te-ts)/(double)ns/(double)data->sim->nTotalAtoms, data->sim->nTotalAtoms);
	      printf("Virial stress = %g\n", estress);
      }

      ts = timeNow();
      eone = nTimeSteps(nsteps,data->sim,dt);
      te = timeNow();

#if (DIAG_LEVEL > 1)
      if (printParallel())
	      printIt(data->sim, stdout);
#endif

#ifdef USE_IN_SITU_VIZ
      data->viz->updateData(data->sim);
#endif
   }

   // sum atoms and energy across all processers
   data->sim->nTot = 0;
   for (int i = 0; i < data->sim->nBoxes; i++)
   {
      data->sim->nTot += data->sim->nAtoms[i];
   }
   data->sim->nTotalAtoms = addIntParallel(data->sim->nTot);
   etotal = addRealParallel(data->sim->e);

   if (printParallel()) 
   {
      printf(" %3d %30.20f computed in %.3fs (%8.4f us/atom for %d atoms)\n",
	      iter*nsteps, etotal, (te-ts), 1.0e6*(te-ts)/(double)(nsteps+1)/(double)data->sim->nTotalAtoms, data->sim->nTotalAtoms);
      fprintf(yamlFile, "MD Loop:\n"
         "  Force Time: %f microseconds/atom\n",
         1.0e6*(te-ts)/(double)(nsteps+1)/(double)data->sim->nTotalAtoms);

#if (DIAG_LEVEL > 1)
      printIt(data->sim, stdout);
#endif

   }

   return NULL;
}
