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
 * mycommand.c
 *
 * A command line parser where arguments can be set.
 *
 * Uses cmdLineParser.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "pmd.h"
#include <time.h>
#include <cmdLineParser.h>

#include "parallel.h"
// the periodic boundary conditions flag
int PERIODIC = 1;

extern void freeArgs();

/**
 * given a struct Command, this function will print
 * the parameters of that struct.
 */
void printCmd(Command* cmd)
{
   if (printParallel())
      printf("\n"
	      "   File = %s\n"
	      "   doeam = %d\n"
	      "   potDir = %s\n"
	      "   potName = %s\n"
	      "   Periodic = %d\n"
	      "   Use gpu = %d\n"
	      "   Number of unit cells = %d x %d x %d\n"
	      "   Number of processors = %d x %d x %d\n"
	      "   Lattice constant = %g\n"
         "   Box factor = %g\n"
         "   Deformation gradient = %g\n",
	      cmd->filename,
	      cmd->doeam,
	      cmd->potDir,
	      cmd->potName,
	      cmd->periodic,
	      cmd->usegpu,
	      cmd->nx, cmd->ny, cmd->nz,
         cmd->xproc, cmd->yproc, cmd->zproc,
         cmd->lat,
         cmd->bf,
         cmd->tke,
         cmd->defGrad
         );

   return;
}

/**
 * given a struct Command, this function will print
 * the parameters of that struct.
 **/
void printCmdYaml(Command* cmd, FILE* file) {
   fprintf(file,
            "Global Run Parameters:\n"
            "  File: %s\n"
            "  doeam: %d\n"
            "  potDir: %s\n"
            "  potName: %s\n"
            "  Periodic: %d\n"
            "  Use gpu: %d\n"
            "  nx: %d\n"
            "  ny: %d\n"
            "  nz: %d\n"
            "  xproc: %d\n"
            "  yproc: %d\n"
            "  zproc: %d\n"
            "  Lattice constant: %g\n"
            "  Box factor: %g\n"
            "  Simulation target kinetic energy = %g\n"
            "  Deformation gradient: %g\n",
            cmd->filename,
            cmd->doeam,
            cmd->potDir,
            cmd->potName,
            cmd->periodic,
            cmd->usegpu,
            cmd->nx, cmd->ny, cmd->nz,
				cmd->xproc, cmd->yproc, cmd->zproc,
            cmd->lat,
            cmd->bf,
            cmd->tke,
            cmd->defGrad);
   return;
}

/**
 * Parse the command line parameters and fill in a
 * struct Command.
 *
 * @param cmd a pointer to struct Command that will be filled in
 * @param argc the number of arguments
 * @param argv the arguments array
 */
void parseCommandLine(Command* cmd, int argc, char** argv)
{
   char spus[1024];
   char debug[1024];
   int digit_optind = 0;
   int c;
   int noperiodic=0;
   int help=0;

   // fill up cmd with defaults
#ifdef GZIPSUPPORT
   strcpy(cmd->filename,"data/8k.inp.gz");
#else
   strcpy(cmd->filename,"data/8k.inp");
#endif
   strcpy(cmd->filename,"");
   strcpy(cmd->potDir,"pots");
   strcpy(cmd->potName,"ag");
   cmd->doeam = 0;
   cmd->usegpu = 0; // set default to use host
   cmd->periodic = 1;
   cmd->nx = 20;
   cmd->ny = 20;
   cmd->nz = 20;
   cmd->bf = 1.0;
   cmd->tke = 0.0;
   cmd->defGrad = 1.0;
   cmd->xproc = 1;
   cmd->yproc = 1;
   cmd->zproc = 1;

   // add arguments for processing
   addArg((char *) "help",       'h',  0,  'i',  &(help), 0,                            (char *) "print this message");
   addArg((char *) "infile",     'f',  1,  's',  cmd->filename,  sizeof(cmd->filename), (char *) "input file name");
   addArg((char *) "potDir",     'd',  1,  's',  cmd->potDir,    sizeof(cmd->potDir),   (char *) "potential directory");
   addArg((char *) "potName",    'p',  1,  's',  cmd->potName,   sizeof(cmd->potName),  (char *) "potential name");
   addArg((char *) "doeam",      'e',  0,  'i',  &(cmd->doeam),  0,                     (char *) "compute eam potentials");
   addArg((char *) "noperiodic", 'o',  0,  'i',  &noperiodic, 0,                        (char *) "do not use periodic bc");
   addArg((char *) "usegpu",     'g',  0,  'i',  &(cmd->usegpu), 0,                     (char *) "use a gpu for OpenCL");
   addArg((char *) "nx",         'x',  1,  'i',  &(cmd->nx), 0,                         (char *) "number of unit cells in x");
   addArg((char *) "ny",         'y',  1,  'i',  &(cmd->ny), 0,                         (char *) "number of unit cells in y");
   addArg((char *) "nz",         'z',  1,  'i',  &(cmd->nz), 0,                         (char *) "number of unit cells in z");
   addArg((char *) "bf",         'b',  1,  'd',  &(cmd->bf), 0,                         (char *) "box factor");
    addArg((char *) "defGrad",    's',  1,  'd',  &(cmd->defGrad), 0,                    (char *) "deformation gradient");
    addArg((char *) "tke",       't',  1,  'd',  &(cmd->tke), 0,                       (char *) "target kinetic energy");

   addArg((char *) "xproc", 'i', 1, 'i', &(cmd->xproc), 0, (char *) "processors in x direction");
   addArg((char *) "yproc", 'j', 1, 'i', &(cmd->yproc), 0, (char *) "processors in y direction");
   addArg((char *) "zproc", 'k', 1, 'i', &(cmd->zproc), 0, (char *) "processors in z direction");

   processArgs(argc,argv);

   if(help)
   {
      printArgs();
      freeArgs();
      exit(2);
   }
   freeArgs();
   cmd->periodic = !(noperiodic);
   PERIODIC=cmd->periodic;

   return;
}
