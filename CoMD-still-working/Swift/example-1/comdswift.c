
#include <assert.h>

// CoMD includes:
#include <docompute.h>
#include <pmd.h>
#include <utility.h>
#include <yamlOutput.h>

#include "comdswift.h"

static void COMDSWIFT_runSim(int argc, char** argv, double* virial_stress);

double
COMDSWIFT_runSimString(char* arg_string)
{
  int count = 1;
  int length = strlen(arg_string);
  char* end = arg_string+length;
  for (char* p = arg_string; p < end; p++)
    if (*p == ' ')
      count++;
  char** array = malloc(count*sizeof(char*));
  assert(array != NULL);
  int index = 0;
  for (char* p = arg_string; p < end; )
  {
    char* q;
    for (q = p+1; *q != ' ' && q < end; q++);
    int n = q-p;
    char* token = malloc((n+1)*sizeof(char));
    strncpy(token, p, n);
    token[n] = '\0';
    array[index++] = token;
    p += (n+1);
  }
  assert(index == count);
  double virial_stress;
  COMDSWIFT_runSim(count, array, &virial_stress);
  for (int i = 0; i < count; i++)
    free(array[i]);
  free(array);
  return virial_stress;
}

static void
COMDSWIFT_runSim(int argc, char** argv, double* virial_stress)
{
  yamlBegin();
  yamlAppInfo(stdout);
  SimFlat* sim = initSimFromCmdLine(argc, argv);
  SimThreadData* threadData =
    (SimThreadData*) suAlignedCalloc(sizeof(SimThreadData));
  threadData->sim = sim;
  threadData->viz = NULL;
  doComputeWork(threadData);
  *virial_stress = threadData->sim->stress;
  destroySimulation(&sim);
  free(threadData);
  yamlEnd();
}
