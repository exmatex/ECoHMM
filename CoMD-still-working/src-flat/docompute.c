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
#include "read.h"
#include "yamlOutput.h"
#include "domains.h"

void printArray(real_t* array, int n, char *name)
{
   //printf("%s;\n", name);
   for (int i=0;i<n;i++)
   {
      //printf("%d, %17.9e\n", i, array[i]);
   }
}

void setStrain(SimFlat *sim)
{
   real_t testMat[10];
   real_t testInv[10];
   for (int m = 0; m<10;m++) 
   {
      testMat[m] = 0.0;
   }
   testMat[0] = 1.0;
   testMat[4] = 1.0;
   testMat[8] = 1.0;
   testMat[1] = 0.0;
   testMat[3] = 0.0;

   matInv3x3(testMat, testInv);

   // transfer to sim structs, including determinant values
   for (int m = 0; m<10; m ++) 
   {
      sim->strain[m] = testMat[m];
      sim->invStrain[m] = testInv[m];
   }
}

void printSimulationDataYaml(FILE* file,
      SimFlat* sim,
      Command cmd,
      struct BasePotential* pot)
{
   fprintf(file, "Simulation data:\n");
   if (cmd.doeam)
   {
      EamPotential *eamPot = (EamPotential*) pot;
      fprintf(file, "  EAM potential values:\n");
      fprintf(file, "    cutoff: "EMT1"\n", eamPot->cutoff);
      fprintf(file, "    mass: "EMT1"\n", eamPot->mass);
      fprintf(file, "  phi potential:\n");
      fprintf(file, "    n: %d\n", eamPot->phi->n);
      fprintf(file, "    x0: "EMT1"\n", eamPot->phi->x0);
      fprintf(file, "    xn: "EMT1"\n", eamPot->phi->xn);
#if (USE_CHEBY)
      fprintf(file, "  Used %d Chebyshev coefficients\n", NUM_CHEBY);
#endif
   }
   else
   {
      LjPotential* ljPot = (LjPotential*) pot;
      fprintf(file, "  LJ potential values:\n");
      fprintf(file, "    cutoff: "EMT1"\n", ljPot->cutoff);
      fprintf(file, "    mass: "EMT1"\n", ljPot->mass);
   }
   fprintf(file, "  nTot: %d\n", sim->nTot);
}

SimFlat *initSimFromCmdLine(int argc, char **argv) 
{
   Command cmd;
   struct BasePotential *pot;
   SimFlat *sim;
   SimThreadData* threadData;

   //printf("    Double Precision=%s\n", (sizeof(real_t)==4?"false":"true") );
   //fprintf(yamlFile,"Double Precision: %s\n", (sizeof(real_t)==4?"false":"true") );

   /* get command line params */
   parseCommandLine(&cmd, argc, argv);

   printCmd(&cmd);

   /* decide whether to get LJ or EAM potentials */
   if(cmd.doeam) 
   {
      //printf("Setting eam potential\n");
      pot = setEamPot(cmd.potDir, cmd.potName);
   } 
   else 
   {
      pot = (BasePotential *) getLjPot();
   }

   if ( ! pot ) simulationAbort(-2,(char *) "Unable to initialize potential");

   //printf("Simulation data\n");
   if (cmd.doeam)
   {
      EamPotential *newPot;
      newPot = (EamPotential*) pot;
      //printf("EAM potential values:\n");
      //printf("cutoff = "EMT1"\n", newPot->cutoff);
      //printf("mass = "EMT1"\n", newPot->mass);
      //printf("phi potential:\n");
      //printf("n = %d\n", newPot->phi->n);
      cmd.lat = newPot->lat;
      //printf("cmd.lat = "EMT1"\n", cmd.lat);
      //printf("%e, %e\n", newPot->phi->x0, newPot->phi->xn);
#if (DIAG_LEVEL > 1)
      //printf("x0 = "EMT1"\n", newPot->phi->x0);
      //printf("xn = "EMT1"\n", newPot->phi->xn);
      printArray(newPot->phi->values, newPot->phi->n, "phi ref");
      printArray(newPot->rho->values, newPot->rho->n, "rho ref");
      printArray(newPot->f->values, newPot->f->n, "f ref");
#endif
   } 
   else 
   {
      LjPotential *newPot;
      newPot = (LjPotential*) pot;
      //printf("LJ potential values:\n");
      //printf("cutoff = "EMT1"\n", newPot->cutoff);
      //printf("mass = "EMT1"\n", newPot->mass);
      // Base lattice places nearest neighbor at zero-force distance
      cmd.lat = 1.587401052*newPot->sigma ;      // This is for Lennard-Jones = 2^(2/3)*sigma
      // Empirically determined prefactor to minimize virial stress
      //cmd.lat = 0.97415*1.587401052*newPot->sigma ;      // This is for Lennard-Jones = 2^(2/3)
   }

   if(strcmp(cmd.filename, ""))
   {
      /* Read in file */
      //printf("\n\n    File is __%s__\n\n", cmd.filename);
      sim = fromFileASCII(cmd, pot);
      if ( ! sim ) simulationAbort(-3,(char *) "Input file does not exist");
   } 
   else 
   {
      //printf("cmd.lat = "EMT1"\n", cmd.lat);
      sim = createFccLattice(cmd, pot);
   }

   setStrain(sim);

   if (cmd.doeam)
   {
      EamPotential *newPot = (EamPotential*) pot;
      //printf("%e, %e\n", newPot->phi->x0, newPot->phi->xn);
      sim->chPot = setChebPot(newPot, NUM_CHEBY);
#if (DIAG_LEVEL > 0)
      //printf("Chebychev coefficients:\n");
      fflush(stdout);
      //printf("%d, %d, %d\n", 
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
#endif
   }

   //printf("    total atoms is: %d\n", sim->nTot);
   //printf("box factor is ("EMT1", "EMT1", "EMT1")\n", 
	 //sim->boxSize[0]/pot->cutoff,
	 //sim->boxSize[1]/pot->cutoff,
	 //sim->boxSize[2]/pot->cutoff);

   /* initial output for consistency check */
   reBoxAll(sim);

   /* apply the desired kinetic energy to the system */
   assignTKE(cmd, sim);

   // must do deformation after reboxing
   forwardDeformation(sim);

   //printCmdYaml(&cmd, yamlFile);
   //printSimulationDataYaml(yamlFile, sim, cmd, pot);

   //printf("Initialization finished\n");
   return sim;
}


void *doComputeWork(SimThreadData *data)
{
   double eone;
   double ts, te;
   int niter = 1;
   int nsteps = 1;
   double dt;

#if (DIAG_LEVEL > 1)
   //printf("Chebychev coefficients:\n");
   fflush(stdout);
   //printf("%d, %d, %d\n", 
	 data->sim->chPot->phi->n,
	 data->sim->chPot->rho->n,
	 data->sim->chPot->f->n);
   fflush(stdout);
#endif

   //printf("Starting simulation\n");
   fflush(stdout);
   ts = timeNow();
   computeForce(data->sim); 
   te = timeNow();

#if (DIAG_LEVEL > 1)
   printIt(data->sim, stdout);
#endif

   /* convert the timestep */
   dt = 1.0e-15;

   //printf("Starting timesteps\n");
   fflush(stdout);
   for (int iter=0; iter<niter; iter++)
   {
      double ns;
      ns = (iter==0?1.0:(double)(1+nsteps));
      /* note the 'nsteps+1' below.
       * this is because we have an extra
       * force call at the end of timesteps()
       */
      //printf(" %3d %30.20e computed in %.3fs (%8.4f us/atom for %d atoms)\n",
	    //iter*nsteps, data->sim->e/data->sim->nTot,(te-ts),1.0e6*(te-ts)/(double)ns/(double)data->sim->nTot,data->sim->nTot);
      ////printf("Virial stress = %g\n", data->sim->stress[0]);

      //printf("Stress tensor :\n [%10e, %10e, %10e]\n [%10e, %10e, %10e]\n [%10e, %10e, %10e]\n",
	   // data->sim->stress[0], data->sim->stress[1], data->sim->stress[2],
	   // data->sim->stress[3], data->sim->stress[4], data->sim->stress[5],
	   // data->sim->stress[6], data->sim->stress[7], data->sim->stress[8]);

      //printf("Strain tensor :\n [%10e, %10e, %10e]\n [%10e, %10e, %10e]\n [%10e, %10e, %10e]\n",
	   // data->sim->strain[0], data->sim->strain[1], data->sim->strain[2],
	   // data->sim->strain[3], data->sim->strain[4], data->sim->strain[5],
	   // data->sim->strain[6], data->sim->strain[7], data->sim->strain[8]);

      ts = timeNow();
      eone = nTimeSteps(nsteps,data->sim,dt);
      te = timeNow();

#if (DIAG_LEVEL > 1)
      printIt(data->sim, stdout);
#endif

#ifdef USE_IN_SITU_VIZ
      data->viz->updateData(data->sim);
#endif
   }

   //printf(" %3d %30.20e computed in %.3fs (%8.4f us/atom for %d atoms)\n",
	 //niter*nsteps, data->sim->e/data->sim->nTot,(te-ts),1.0e6*(te-ts)/(double)(nsteps+1)/(double)data->sim->nTot,data->sim->nTot);
   //fprintf(yamlFile, "MD Loop:\n"
	 //"  Force Time: %f microseconds/atom\n",
	// 1.0e6*(te-ts)/(double)(nsteps+1)/(double)data->sim->nTot);

  fprintf(stdout, "%10e\n", data->sim->stress[0]);

/*
   FILE *sFile;
   sFile = fopen("stressTensor.txt","w");

   fprintf(sFile, "%10e, %10e, %10e %10e, %10e, %10e %10e, %10e, %10e",
	 data->sim->stress[0], data->sim->stress[1], data->sim->stress[2],
	 data->sim->stress[3], data->sim->stress[4], data->sim->stress[5],
	 data->sim->stress[6], data->sim->stress[7], data->sim->stress[8]);

   fclose(sFile);
*/
   return NULL;

#if (DIAG_LEVEL > 1)
   printIt(data->sim, stdout);
#endif

}
