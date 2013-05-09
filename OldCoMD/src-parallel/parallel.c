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
 * This file contains functions for the parallel calls.
 */

#include <stdarg.h>
#include <stdio.h>
#include <strings.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include "parallel.h"

// global parallel variables
int myParallel = -1;
int totalParallel = 0;

static  MPI_Status  status;
static  int rbuffer[65536]; //!< Temporary buffer for receiving

#ifdef DO_MPI
MPI_Comm myCommHosts = MPI_COMM_WORLD;
MPI_Comm myCommZplane = MPI_COMM_NULL;
#endif

int bytesReceived; //!< number of bytes received

/**
 * return total number of processors
 */
int getTotalParallel()
{
   return totalParallel;
}


/**
 * return local rank
 */
int getMyParallel()   
{
   return myParallel;
}

/**
 * return whether printing occurs from this rank
 */
int printParallel()
{
   if (myParallel == 0) return 1;
   return 0;
}

/**
 * initialize parallel system
 */
void initParallel(int* argc, char*** argv)
{
#ifdef DO_MPI
   MPI_Init(argc, argv);

   MPI_Comm_rank(myCommHosts, &myParallel);
   MPI_Comm_size(myCommHosts, &totalParallel);

   //printf("My rank = %d  Total ranks = %d\n", myParallel, totalParallel);
#endif
}

/**
 * finalize parallel system
 */
void destroyParallel()
{
#ifdef DO_MPI
   MPI_Finalize();
#endif
}

/**
 * barrier parallel
 */
void barrierParallel()
{
#ifdef DO_MPI
   MPI_Barrier(myCommHosts);
#endif
}

/**
 * send and receive
 */
void sendReceiveParallel(int sour, int stag, void* sbuffer, int slen,
                         int dest, int dtag, void* dbuffer, int dlen)
{
#ifdef DO_MPI
   int buf, bytes;

   MPI_Sendrecv(dbuffer, dlen, MPI_BYTE, dest, 0,
                rbuffer, slen, MPI_BYTE, sour, 0,
                myCommHosts, &status);
   MPI_Get_count(&status, MPI_BYTE, &bytesReceived);

//   printf("%d : Received %d bytes\n", myParallel, bytesReceived);
 
   bcopy(rbuffer, sbuffer, bytesReceived);

#endif
}

/**
 * blocking send
 */
void sendBlockingParallel(int dest, int tag, void* buffer, int len)
{
//   printf("send tag is %d\n", tag);

#ifdef DO_MPI
   MPI_Send(buffer, len, MPI_BYTE, dest, tag, myCommHosts);
#endif
}

/**
 * blocking receive
 */
void receiveBlockingParallel(int sour, int tag, void* buffer, int len)
{
#ifdef DO_MPI
   if ( tag < 0 )
   {
     printf("break at receive tag is %d\n", tag);
   }

//   printf("receive tag is %d\n", tag);

   MPI_Recv(buffer, len, MPI_BYTE, sour, tag, myCommHosts, &status);
   MPI_Get_count(&status, MPI_BYTE, &bytesReceived);
#endif
}

/**
 * non-blocking send
 */
void sendParallel(int dest, int tag, void* buffer, int len)
{
#ifdef DO_MPI
   MPI_Request req;
   int flag;

//   printf("isend tag is %d\n", tag);

   MPI_Isend(buffer, len, MPI_BYTE, dest, tag, myCommHosts, &req);
   MPI_Test(&req, &flag, MPI_STATUS_IGNORE);

   if (req != MPI_REQUEST_NULL)
       MPI_Request_free(&req);

  // Note that we ignore req and just overwrite

//   MPI_Isend(buffer, len, MPI_BYTE, dest, tag, myCommHosts, MPI_REQUEST_NULL);
#endif
}

/**
 * non-blocking receive
 */
void receiveParallel(int dest, int tag, void* buffer, int len)
{
#ifdef DO_MPI
   MPI_Request req;
   int flag;

   MPI_Irecv(buffer, len, MPI_BYTE, dest, tag, myCommHosts, &req);
   MPI_Test(&req, &flag, MPI_STATUS_IGNORE);

   if (req != MPI_REQUEST_NULL)
     MPI_Request_free(&req);
#endif
}

/**
 * parallel add int reduction
 */
int addIntParallel(int n)
{
#ifdef DO_MPI
   int result;

   MPI_Allreduce(&n, &result, 1, MPI_INT, MPI_SUM, myCommHosts);

   return result;
#endif
}

/**
 * parallel add real (real_t) reduction
 */
real_t addRealParallel(real_t a)
{
#ifdef DO_MPI
   real_t result;

#ifdef SINGLE
   MPI_Allreduce(&a, &result, 1, MPI_FLOAT, MPI_SUM, myCommHosts);
#else
   MPI_Allreduce(&a, &result, 1, MPI_DOUBLE, MPI_SUM, myCommHosts);
#endif

   return result;
#endif
}
