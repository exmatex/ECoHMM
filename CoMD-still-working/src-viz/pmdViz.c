/**
 * a simple md simulator
 **/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

#include "pmd.h"
#include "docompute.h"
#include "yamlOutput.h"
#include "utility.h"

#include <pthread.h>


void *computation_thread(void *ptr)
/* do the computation work */
{
  SimThreadData *data = (SimThreadData*)ptr;
  return doComputeWork(data);
}


int main(int argc, char **argv) {
  SimFlat *sim;
  SimThreadData* threadData;
printf("A\n");
  yamlBegin();
  yamlAppInfo(stdout);

  /* get sim_flat from cmdline */
  sim = initSimFromCmdLine(argc, argv);

  writeClsman(sim,(char *) "init.bin");
printf("B\n");
  /* start computation thread */
  if(1) {
    pthread_t computationThread;
    pthread_t visualizationThread;
    AtomVisualize* viz = new AtomVisualize();
printf("C\n");
    viz->initialize(sim);
    viz->updateData(sim);
    viz->renderData();
printf("D\n");      
    threadData = new SimThreadData;
    threadData->sim = sim;
    threadData->viz = viz;
printf("E\n");		  
    pthread_create(&computationThread, NULL, computation_thread, (void*) threadData);
    viz->interact();
  }

  writeClsman(sim,(char *) "final.bin");

  /* release sim_flat */
  destroySimulation(&sim); 	

  return 0;
}
