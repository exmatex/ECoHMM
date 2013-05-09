/**
 * a simple md simulator
 **/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <pthread.h>

#include "mytype.h"
#include "pmdOCL.h"
#include "helpers.h"
#include "yamlOutput.h"

#ifdef INTEROP_VIZ

#if defined (__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#else
#include <glut.h>
#endif

int mouse_old_x, mouse_old_y, mouse_buttons;
int g_threadStatus, g_clFinished;
extern Quaternion q;
extern float cameraFOV;
#endif

cl_kernel advanceVelocitySoa;
cl_kernel advancePositionSoa;
cl_kernel vizSoa;
cl_kernel *forceKernelsSoa;

cl_kernel advanceVelocityAos;
cl_kernel advancePositionAos;
cl_kernel vizAos;
cl_kernel *forceKernelsAos;

int err;
size_t nLocal[2];
size_t nGlobal[2];

cl_event forceEvent;
cl_event apEvent;
cl_event avEvent;

real_t tExec;
real_t tOverall;
real_t tLocal;
cl_real dt;

HostSimSoa simHostSoa;
DevSimSoa simDevSoa;

HostSimAos simHostAos;
DevSimAos simDevAos;

double ts, te;

int eamFlag;
int gpuFlag;
SimFlat *sim;
int iter, ns;
int nIter, nSteps;


void computeIterationSoA()
{
   real_t tKern, tEnq;
   real_t tAcc = 0.0;
   for (int ns = 0;ns < nSteps; ns++) 
   {
      // advance particle positions dt/2
      oclRunKernel(advancePositionSoa, &apEvent, nGlobal, nLocal);
      clWaitForEvents(1, &apEvent);
      getElapsedTime(apEvent, &tKern, &tEnq);
      tAcc += tKern;

      if (DIAG_LEVEL > 1)
	 getPrintState(simDevSoa, simHostSoa);

      // compute force
      computeForceOCL(forceKernelsSoa, &forceEvent, nGlobal, nLocal, eamFlag, simHostSoa.nTot, &tKern);
      tAcc += tKern;

      if (DIAG_LEVEL > 1)
	 getPrintState(simDevSoa, simHostSoa);

      // advance velocity a full timestep
      oclRunKernel(advanceVelocitySoa, &avEvent, nGlobal, nLocal);
      clWaitForEvents(1, &avEvent);
      getElapsedTime(avEvent, &tKern, &tEnq);
      tAcc += tKern;

      if (DIAG_LEVEL > 1)
	 getPrintState(simDevSoa, simHostSoa);

      // advance particle positions dt/2
      oclRunKernel(advancePositionSoa, &apEvent, nGlobal, nLocal);
      clWaitForEvents(1, &apEvent);
      getElapsedTime(apEvent, &tKern, &tEnq);
      tAcc += tKern;

      if (DIAG_LEVEL > 1)
	 getPrintState(simDevSoa, simHostSoa);
   }

#ifdef INTEROP_VIZ 
   oclGraphics(vizSoa, simDevSoa, nGlobal, nLocal);
#endif

   // compute force
   computeForceOCL(forceKernelsSoa, &forceEvent, nGlobal, nLocal, eamFlag, simHostSoa.nTot, &tKern);

   if (DIAG_LEVEL > 0)
      getPrintState(simDevSoa, simHostSoa);

   printf(" %4d ", nSteps*(iter + 1));
   computePrintEnergy(simDevSoa, simHostSoa);
   printf("    computed in  = %e (%e us/atom for %d atoms)\n", 
	 tAcc,1.0e6*tAcc/(double)(sim->nTot*nSteps), sim->nTot);
}

void computeIterationAoS()
{
   real_t tKern, tEnq;
   real_t tAcc = 0.0;
   for (int ns = 0;ns < nSteps; ns++) 
   {
      // advance particle positions dt/2
      oclRunKernel(advancePositionAos, &apEvent, nGlobal, nLocal);
      clWaitForEvents(1, &apEvent);
      getElapsedTime(apEvent, &tKern, &tEnq);
      tAcc += tKern;

      if (DIAG_LEVEL > 1)
	 getPrintStateAoS(simDevAos, simHostAos);

      // compute force
      computeForceOCL(forceKernelsAos, &forceEvent, nGlobal, nLocal, eamFlag, simHostSoa.nTot, &tKern);
      tAcc += tKern;

      if (DIAG_LEVEL > 1)
	 getPrintStateAoS(simDevAos, simHostAos);

      // advance velocity a full timestep
      oclRunKernel(advanceVelocityAos, &avEvent, nGlobal, nLocal);
      clWaitForEvents(1, &avEvent);
      getElapsedTime(avEvent, &tKern, &tEnq);
      tAcc += tKern;

      if (DIAG_LEVEL > 1)
	 getPrintStateAoS(simDevAos, simHostAos);

      // advance particle positions dt/2
      oclRunKernel(advancePositionAos, &apEvent, nGlobal, nLocal);
      clWaitForEvents(1, &apEvent);
      getElapsedTime(apEvent, &tKern, &tEnq);
      tAcc += tKern;

      if (DIAG_LEVEL > 1)
	 getPrintStateAoS(simDevAos, simHostAos);
   }

#ifdef INTEROP_VIZ 
   //oclGraphics(vizSoa, simDevSoa, nGlobal, nLocal);
#endif

   // compute force
   computeForceOCL(forceKernelsAos, &forceEvent, nGlobal, nLocal, eamFlag, simHostSoa.nTot, &tKern);

   if (DIAG_LEVEL > 0)
      getPrintStateAoS(simDevAos, simHostAos);

   //printf(" after %4d steps: ", nSteps*(iter + 1));
   printf(" %4d ", nSteps*(iter + 1));
   computePrintEnergyAoS(simDevAos, simHostAos);
   printf("    computed in  = %e (%e us/atom for %d atoms)\n", 
	 tAcc,1.0e6*tAcc/(double)(sim->nTot*nSteps), sim->nTot);
}

void finishOclSoA()
{
   getElapsedTime(apEvent, &tExec, &tOverall);
   printf("Kernel advancePosition executed in %e secs.\n", tExec);
   getElapsedTime(avEvent, &tExec, &tOverall);
   printf("Kernel advanceVelocitySoa executed in %e secs.\n", tExec);

   real_t tRef = (real_t)(te - ts);

   printf("Reference wallclock elapsed time = %e (%e us/atom for %d atoms)\n", 
	 tRef,1.0e6*tRef/(double)(sim->nTot*nIter*nSteps), sim->nTot);
   fprintf(yamlFile, "MD Loop:\n"
	 "  Overall SoA Time: %f microseconds/atom\n",
	 1.0e6*(tRef)/(double)((nSteps+1)*nIter*sim->nTot));

   // copy result back from device
   // note this utility forces synchronous copy so 
   // it does not return until the copy is complete
   getVec(simDevSoa.f, simHostSoa.f, simHostSoa.arraySize);
   oclCopyToHost(simDevSoa.e, simHostSoa.e, simHostSoa.arraySize);

   computePrintEnergy(simDevSoa, simHostSoa);
   printf("\n");

   // clean up the OpenCL memory 
   FreeSims(simHostSoa, simDevSoa);

}


void finishOclAoS()
{
   getElapsedTime(apEvent, &tExec, &tOverall);
   printf("Kernel advancePositionAos executed in %e secs.\n", tExec);
   getElapsedTime(avEvent, &tExec, &tOverall);
   printf("Kernel advanceVelocityAos executed in %e secs.\n", tExec);

   real_t tRef = (real_t)(te - ts);

   printf("Reference wallclock elapsed time = %e (%e us/atom for %d atoms)\n", 
	 tRef,1.0e6*tRef/(double)(sim->nTot*nIter*nSteps), sim->nTot);
   fprintf(yamlFile, "MD Loop:\n"
	 "  Overall AoS Time: %f microseconds/atom\n",
	 1.0e6*(tRef)/(double)((nSteps+1)*nIter*sim->nTot));

   // copy result back from device
   // note this utility forces synchronous copy so 
   // it does not return until the copy is complete
   getVecAoS(simDevAos.f, simHostAos.f, simHostAos.arraySize);
   oclCopyToHost(simDevAos.e, simHostAos.e, simHostAos.arraySize);

   computePrintEnergyAoS(simDevAos, simHostAos);
   printf("\n");

   // clean up the OpenCL memory 
   FreeSimsAoS(simHostAos, simDevAos);

   clReleaseEvent(apEvent);
   clReleaseEvent(avEvent);

   oclCleanup();
}

void runRef()
{
   printf("************************************************************************\n");
   printf("Running reference simulation\n");
   // write initial state 
   //writeClsman(sim,(char *) "init.bin");

   // do the computation 
   (void)doComputeWork(sim); 

   // write final configuration 
   //writeClsman(sim,(char *) "final.bin");

   // free memory 
   destroySimulation(&sim);
}

void computeInitSoA()
{
   printf("************************************************************************\n");
   printf("Initializinging SoA OpenCL simulation\n");

   simDevSoa.rMass = (cl_real)(amu_to_m_e)*(double)(sim->pot->mass);
   simHostSoa.eamFlag = eamFlag;

   if(eamFlag) 
   {
      forceKernelsSoa = malloc(sizeof(cl_kernel)*3);
   } 
   else 
   {
      forceKernelsSoa = malloc(sizeof(cl_kernel)*1);
   }

   // set up arrays for OpenCL

   initHostSim(&simHostSoa, sim);

   initDevSim(&simDevSoa, &simHostSoa);

   putSim(simHostSoa, simDevSoa);

   // build the program from the kernel source file

   buildModulesSoa(forceKernelsSoa, &advancePositionSoa, &advanceVelocitySoa, &vizSoa, simHostSoa, nLocal, nGlobal);

   cl_real dthalf = 0.5*simHostSoa.dt;
   cl_real dtminushalf = -0.5*simHostSoa.dt;

   if (eamFlag) 
   {
      // set the arguments for all 3 EAM_Force kernels
      setEAMArgs(forceKernelsSoa, simDevSoa);

   } 
   else 
   {
      // set kernel arguments for ljForce
      setLJArgs(forceKernelsSoa[0], simDevSoa);
   }

   // set kernel arguments for advanceVelocitySoa
   setAVArgs(advanceVelocitySoa, simDevSoa, simHostSoa.dt);

   // set kernel arguments for advancePosition
   setAPArgs(advancePositionSoa, simDevSoa, dthalf);

   // Start the simulation here;
   printf("************************************************************************\n");
   printf("Starting SoA OpenCL simulation\n");

   //ts = timeNow();

   cl_real tKern;
   computeForceOCL(forceKernelsSoa, &forceEvent, nGlobal, nLocal, eamFlag, simHostSoa.nTot, &tKern);

   getVec(simDevSoa.f, simHostSoa.f, simHostSoa.arraySize);

   printf(" Initial energy:");
   computePrintEnergy(simDevSoa, simHostSoa);
   printf("\n");

   getPrintState(simDevSoa, simHostSoa);
}

void computeInitAoS()
{
   printf("************************************************************************\n");
   printf("Initializinging AoS OpenCL simulation\n");

   simDevAos.rMass = (cl_real)(amu_to_m_e)*(double)(sim->pot->mass);
   simHostAos.eamFlag = eamFlag;

   if(eamFlag) 
   {
      forceKernelsAos = malloc(sizeof(cl_kernel)*3);
   } 
   else 
   {
      forceKernelsAos = malloc(sizeof(cl_kernel)*1);
   }

   // set up arrays for OpenCL

   initHostSimAoS(&simHostAos, sim);

   initDevSimAoS(&simDevAos, &simHostAos);

   putSimAoS(simHostAos, simDevAos);

   // build the program from the kernel source file

   buildModulesAoS(forceKernelsAos, &advancePositionAos, &advanceVelocityAos, &vizAos, simHostAos, nLocal, nGlobal);

   cl_real dthalf = 0.5*simHostAos.dt;
   cl_real dtminushalf = -0.5*simHostAos.dt;

   // set kernel arguments for advanceVelocitySoa
   setAVArgsAoS(advanceVelocityAos, simDevAos, simHostAos.dt);

   // set kernel arguments for advancePosition
   setAPArgsAoS(advancePositionAos, simDevAos, dthalf);

   if (eamFlag) 
   {
      // set the arguments for all 3 EAM_Force kernels
      setEAMArgsAoS(forceKernelsAos, simDevAos);

   } 
   else 
   {
      // set kernel arguments for ljForce
      setLJArgsAoS(forceKernelsAos[0], simDevAos);
   }

   // Start the simulation here;
   printf("************************************************************************\n");
   printf("Starting AoS OpenCL simulation\n");

   //ts = timeNow();

   cl_real tKern;
   computeForceOCL(forceKernelsAos, &forceEvent, nGlobal, nLocal, eamFlag, simHostAos.nTot, &tKern);

   getVecAoS(simDevAos.f, simHostAos.f, simHostAos.arraySize);

   printf(" Initial energy:");
   computePrintEnergyAoS(simDevAos, simHostAos);
   printf("\n");

   getPrintStateAoS(simDevAos, simHostAos);
}


#ifdef INTEROP_VIZ
void keyboard(unsigned char key, int x, int y) 
{ 
   glutPostRedisplay(); 
}

void idle() 
{ 
   glutPostRedisplay(); 
}

void mouse(int button, int state, int x, int y) 
{
   if (state == GLUT_DOWN) mouse_buttons |= 1<<button;
   else if (state == GLUT_UP) mouse_buttons = 0;

   mouse_old_x = x;
   mouse_old_y = y;
   glutPostRedisplay();
}

void motion(int x, int y) 
{
   float dx = x - mouse_old_x;
   float dy = y - mouse_old_y;

   if (mouse_buttons == 1)
   {
      Quaternion newRotX;
      QuaternionSetEulerAngles(&newRotX, -0.2*dx*3.14159/180.0, 0.0, 0.0);
      QuaternionMul(&q, q, newRotX);

      Quaternion newRotY;
      QuaternionSetEulerAngles(&newRotY, 0.0, 0.0, -0.2*dy*3.14159/180.0);
      QuaternionMul(&q, q, newRotY);
   }
   else if (mouse_buttons == 4)
   {
      cameraFOV += dy/25.0f;
   }

   mouse_old_x = x;
   mouse_old_y = y;
   glutPostRedisplay();
}

void renderData() 
{
   if (iter == 0) ts = timeNow();
   if (iter++ < nIter) computeIterationSoA();
   if (iter == nIter) 
   { 
      te = timeNow();
      finishOclSoA();
   }

   glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   oclRender();
   glutSwapBuffers();
}
#endif


int main(int argc, char **argv) 
{

#ifdef INTEROP_VIZ
   q.x = q.y = q.z = 0.0f;  q.w = 1.0f;  
   mouse_old_x = 0; mouse_old_y = 0; mouse_buttons = 0; cameraFOV = 60.0f;
   g_threadStatus = 0;  g_clFinished = 0;
   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
   glutInitWindowSize(1024, 1024);
   glutCreateWindow("CoMD"); 
   glutDisplayFunc(renderData);
   glutKeyboardFunc(keyboard);
   glutMouseFunc(mouse);
   glutMotionFunc(motion);
   glutIdleFunc(idle);

   glewInit();
#endif

   yamlBegin();
   yamlAppInfo(stdout);

   /* get sim_flat from cmdline */
   sim = initSimFromCmdLine(argc, argv, &eamFlag, &gpuFlag);

   oclInit(gpuFlag);

#ifdef INTEROP_VIZ
   oclInitInterop(sim->nBoxes);
   computeInitSoA();
   iter = 0;  nIter = 10;  nSteps = 10;
   glutMainLoop();
#else
   nIter = 10;  nSteps = 10;
   // SoA test loop
   computeInitSoA();
   ts = timeNow();
   for (iter=0; iter<nIter; iter++) computeIterationSoA();
   te = timeNow();
   finishOclSoA();
   // AoS test loop
   computeInitAoS();
   ts = timeNow();
   for (iter=0; iter<nIter; iter++) computeIterationAoS();
   te = timeNow();
   finishOclAoS();
   runRef();
#endif
   yamlEnd();

}

