#include "yamlOutput.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "CoMD_info.h"

FILE* yamlFile = NULL;


static const char* CoMDVersion = "1.0";

static char* getTimeString()
{
   time_t rawtime;
   struct tm * timeinfo;
   time(&rawtime);
   timeinfo = localtime(&rawtime);

   char* timestring = malloc(32);
   sprintf(timestring,
           "%4d-%02i-%02d, %02d-%02d-%02d",
           timeinfo->tm_year+1900,
           timeinfo->tm_mon+1,
           timeinfo->tm_mday,
           timeinfo->tm_hour,
           timeinfo->tm_min,
           timeinfo->tm_sec);
   return timestring;
}


void yamlBegin()
{
   char filename[64];
   time_t rawtime;
   time (&rawtime);
   struct tm* ptm = localtime(&rawtime);
   char sdate[25];
  //use tm_mon+1 because tm_mon is 0 .. 11 instead of 1 .. 12
   sprintf (sdate,"%04d:%02d:%02d-%02d:%02d:%02d",
            ptm->tm_year + 1900, ptm->tm_mon+1,
            ptm->tm_mday, ptm->tm_hour, ptm->tm_min,ptm->tm_sec);
   sprintf(filename, "CoMDOCL.%s.yaml", sdate);
   yamlFile = fopen(filename, "w");
   yamlAppInfo(yamlFile);
}

void yamlEnd()
{
   fclose(yamlFile);
}



void yamlAppInfo(FILE* file)
{
   fprintf(file,"Mini-Application Name: CoMD\n");
   fprintf(file,"Mini-Application Version: %s\n", CoMDVersion);
   fprintf(file,"Platform:\n");
   fprintf(file,"  hostname: %s\n",         CoMD_HOSTNAME);
   fprintf(file,"  kernel name: %s\n",      CoMD_KERNEL_NAME);
   fprintf(file,"  kernel release: %s\n",   CoMD_KERNEL_RELEASE);
   fprintf(file,"  processor: %s\n",        CoMD_PROCESSOR);
   fprintf(file,"Build:\n");
   fprintf(file,"  CC: %s\n",               CoMD_COMPILER);
   fprintf(file,"  compiler version: %s\n", CoMD_COMPILER_VERSION);
   fprintf(file,"  CFLAGS: %s\n",           CoMD_CFLAGS);
   fprintf(file,"  using MPI: no\n");
   fprintf(file,"  Threading: none\n");
   char* timestring = getTimeString();
   fprintf(file,"Run Date/Time: %s\n", timestring);
   free(timestring);
   fprintf(file, "\n");
}
