#include "mytype.h"
#include "pmdOCL.h"
#include "cheby.h"
#include "ic_fcc.h"
#include "read.h"
#include "helpers.h"
#include "yamlOutput.h"

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

SimFlat *initSimFromCmdLine(int argc, char **argv, int *eamFlag, int *gpuFlag)
{
   Command cmd;
   struct BasePotential *pot;
   SimFlat *sim;
   SimThreadData* threadData;

   printf("    Double Precision=%s\n", (sizeof(real_t)==4?"false":"true") );
   fprintf(yamlFile,"Double Precision: %s\n", (sizeof(real_t)==4?"false":"true") );

   /* get command line params */
   parseCommandLine(&cmd, argc, argv);

   printCmd(&cmd);

   *eamFlag = cmd.doeam;
   *gpuFlag = cmd.usegpu;

   /* decide whether to get LJ or EAM potentials */
   if(cmd.doeam) 
   {
      pot = setEamPot(cmd.potDir, cmd.potName);
   }
   else
   {
      pot = (BasePotential *) getLjPot();
   }

   if ( ! pot ) simulationAbort(-2,(char *) "Unable to initialize potential");

   printf("Simulation data\n");
   if (cmd.doeam)
   {
      EamPotential *new_pot;
      new_pot = (EamPotential*) pot;
      printf("EAM potential values:\n");
      printf("cutoff = %e\n", new_pot->cutoff);
      printf("mass = %e\n", new_pot->mass);
      printf("phi potential:\n");
      printf("n = %d\n", new_pot->phi->n);
      cmd.lat = new_pot->lat;
      printf("cmd.lat = %e\n", cmd.lat);
   }
   else
   {
      LjPotential *new_pot;
      new_pot = (LjPotential*) pot;
      printf("LJ potential values:\n");
      printf("cutoff = %e\n", new_pot->cutoff);
      printf("mass = %e\n", new_pot->mass);
      // Base lattice places nearest neighbor at zero-force distance
      cmd.lat = 1.587401052*new_pot->sigma ;      // This is for Lennard-Jones = 2^(2/3)*sigma
      // Empirically determined prefactor to minimize virial stress
      //cmd.lat = 0.97415*1.587401052*new_pot->sigma ;      // This is for Lennard-Jones = 2^(2/3)
      // Low 0ccupancy cells for testing
      //cmd.lat = 2.5*new_pot->sigma;

      printf("cmd.lat = %e\n", cmd.lat);
   }

   if(strcmp(cmd.filename, ""))
   {
      /* Read in file */
      printf("\n\n    File is __%s__\n\n", cmd.filename);
      sim = fromFileASCII(cmd, pot);
      if ( ! sim ) simulationAbort(-3,(char *) "Input file does not exist");
   }
   else
   {
      sim = createFccLattice(cmd, pot);
   }

   setStrain(sim);

   printf("Initial condition set\n");

   if (cmd.doeam)
   {
      EamPotential *new_pot = (EamPotential*) pot;
      printf("%e, %e\n", new_pot->phi->x0, new_pot->phi->xn);
      sim->chPot = setChebPot(new_pot, 32);
#if (DIAG_LEVEL > 0)
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
#endif
   }

   printf("    total atoms is: %d\n", sim->nTot);
   printf("box factor is (%e, %e, %e)\n", 
         sim->boxSize[0]/pot->cutoff,
         sim->boxSize[1]/pot->cutoff,
         sim->boxSize[2]/pot->cutoff);

   /* initial output for consistency check */
   reBoxAll(sim);

   /* apply the desired kinetic energy to the system */
   assignTKE(cmd, sim);

   printCmdYaml(&cmd, yamlFile);
   printSimulationDataYaml(yamlFile, sim, cmd, pot);

   return sim;
}


void *doComputeWork(SimFlat* sim)
{
   int iter;
   double eone;
   double ts, te;
   int niter = 10;
   int nsteps = 10;
   double dt;

   ts = timeNow();
   computeForce(sim); 
   te = timeNow();

   printIt(sim, stdout);

   /* convert the timestep */
   dt = 1.0e-15;

   for (iter=0; iter<niter; iter++)
   {
      double ns;
      ns = (iter==0?1.0:(double)(1+nsteps));
      /* note the 'nsteps+1' below.
       * this is because we have an extra
       * force call at the end of timesteps()
       */
      printf(" %3d %30.20f computed in %.3es (%8.4f us/atom for %d atoms)\n",
            iter*nsteps, sim->e,(te-ts),1.0e6*(te-ts)/(double)ns/(double)sim->nTot,sim->nTot);
      printf("Virial stress = %g\n", sim->stress[0]);

      printf("Stress tensor :\n [%10e, %10e, %10e]\n [%10e, %10e, %10e]\n [%10e, %10e, %10e]\n",
            sim->stress[0], sim->stress[1], sim->stress[2],
            sim->stress[3], sim->stress[4], sim->stress[5],
            sim->stress[6], sim->stress[7], sim->stress[8]);


      ts = timeNow();
      eone = nTimeSteps(nsteps,sim,dt);
      te = timeNow();
   }

   printf(" %3d %30.20f computed in %.3es (%8.4f us/atom for %d atoms)\n",
         iter*nsteps, sim->e,(te-ts),1.0e6*(te-ts)/(double)(nsteps+1)/(double)sim->nTot,sim->nTot);
   fprintf(yamlFile, "MD Loop:\n"
         "  Reference Time: %f microseconds/atom\n",
         1.0e6*(te-ts)/(double)(nsteps+1)/(double)sim->nTot);

   return NULL;
}
