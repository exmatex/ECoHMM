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
 * CoMD: A Molecular Dynamics proxy applications suite
 *
 * CoMD is designed to be an extensible molecular dynamics
 * proxy applications suite.  The idea is that we will
 * capture the essense of different features of molecular
 * dynamics in this software and extend it to multiple
 * platforms.  The code will serve as a vehicle for co-design
 * by allowing others to extend it as needed to get better
 * performance on different architectures.
 *
 * -# A good introduction to molecular dynamics can be found
 *   in "Computer simulation of liquids" by M.P. Allen and
 *   D.J. Tildesley.
 *   ISBN-10: 0198556454 | ISBN-13: 978-0198556459
 * -# A description of how to run MD at large scales can be
 *   found in "369 Tflop/s molecular dynamics simulations
 *   on the Roadrunner general-purpose heterogeneous
 *   supercomputer" by S. Swaminarayan, K. Kadau, T.C. Germann,
 *   and G. Fossum published in Proceedings of SC'08.
 *
 * Our code allows the user to use either the Lennard-Jones
 * potential or the Embedded Atom Method Potential.  Both these
 * are discussed in #2 above.
 *
 *  
 **/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>

#include "utility.h"
#include "pmd.h"
#include "yamlOutput.h"
#include "docompute.h"

/**
 * All we do in the main subroutine is:
 * -# read in a data file
 * -# write out a file that can be read by
 *    clsman (http://www.t12.lanl.gov/home/afv/)
 * -# call doComputeWork() to compute on this data
 * -# write out a final configuration that
 *    is clsman compatible
 * -# free data
 *
 * command line parameters are:
 * -# -h            : print command line parameters
 * -# -f <filename> : the name of the inputfile
 * -# -e            : use EAM potentials
 * -# -p <potName>  : name of the EAM potential
 * -# -d <potDir>   : directory where EAM potential files reside
 * -# -z            : disable periodic boundary conditions
 *
 **/
int main(int argc, char **argv) 
{
    SimFlat *sim;
    SimThreadData* threadData;

    //yamlBegin();
    //yamlAppInfo(stdout);
    
    /* get sim_flat from cmdline */
    sim = initSimFromCmdLine(argc, argv);
    threadData = (SimThreadData*)suAlignedCalloc(sizeof(SimThreadData));
    threadData->sim = sim;
    threadData->viz = NULL;


    /* write initial state */
    //  writeClsman(sim,(char *) "init.bin");

    /* do the computation */
    //nanosleep((struct timespec[]){{0, 500000000}}, NULL);
    //sleep(20);
    (void)doComputeWork(threadData);

    /* write final configuration */
    // writeClsman(sim,(char *) "final.bin");

    /* free memory */
    destroySimulation(&sim);
    free(threadData);
    //yamlEnd();

    return(-22);
}


