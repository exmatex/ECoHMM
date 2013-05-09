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

#ifndef _MEMUTILS_H_

#define _MEMUTILS_H_

#include <unistd.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>

#ifndef DMA_PAD 
#define DMA_PAD(iSize) (((iSize)+15) & ~15)
#endif

extern void breakinme();
extern int suAbort(int status, char *msg);
  
    
/**
 * use posix_memalign to alloc **/

/** \def ALIGNED_MALLOC(ptr,size)
 * \brief a macro that allocates size bites into ptr using posix_memalign
 **/
#define ALIGNED_MALLOC(ptr,size) {					\
    if(posix_memalign((void **) &ptr,128,DMA_PAD(size))){		\
      char t[1024];sprintf(t,"unable to dingo malloc %ld\n",(unsigned long) size); \
      breakinme();							\
      suAbort(195,t);							\
    }									\
    memset(ptr,0,size);							\
  }
    
/** \def ALIGNED_FREE(ptr)
 * \brief a macro that frees memory allocated by ALIGNED_MALLOC 
 **/
#define ALIGNED_FREE(ptr) {if(ptr) free(ptr);}
  
  
/** \def ALIGNED_CALLOC(ptr,size)
 * \brief a macro that allocates size bites into ptr using posix_memalign and zeroes them.
 **/
#define ALIGNED_CALLOC(ptr,size) {ALIGNED_MALLOC(ptr,size);memset(ptr,0,size);}

/**utility routine for allocaating aligned bytes and zeroing them**/
static void *suAlignedCalloc(size_t isize/**<number of bytes to alloc**/) {
  void *b;
  ALIGNED_CALLOC(b,isize);
  return b;
}

/**utility routine for allocaating aligned bytes**/
static void *suAlignedMalloc(size_t isize) {
  void *b;
  ALIGNED_MALLOC(b,isize);
  return b;
}

/** utility routine for freeing bytes allocated using either suAlignedMalloc() or suAlignedFree() **/
static void suAlignedFree(void *d) {
  ALIGNED_FREE(d);
}

#endif
