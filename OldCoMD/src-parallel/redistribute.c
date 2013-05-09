#include <unistd.h>
#include <stdlib.h>

#include "pmd.h"
#include "redistribute.h"
#include "geometry.h"
#include "parallel.h"
#include "domains.h"
#include "redistribute.h"

#define MAXVARS 13
#define IMAXVARS 2

#define FACE_X_MINUS 0
#define FACE_X_PLUS  1

#define FACE_Y_MINUS 2
#define FACE_Y_PLUS  3

#define FACE_Z_MINUS 4
#define FACE_Z_PLUS  5

#define   PROC_MASK  0x0000ffff

static int step = 0;

const int fLookup[3][2] = {{FACE_X_MINUS,FACE_X_PLUS},
                           {FACE_Y_MINUS,FACE_Y_PLUS},
                           {FACE_Z_MINUS,FACE_Z_PLUS}};

int shiftFlags[6] = {0};
int procAddress[6] = {0};
double shiftR[6] = {0.0};
double reflectR[6] = {0.0};
real_t *sendBuffer[6][MAXVARS] = {NULL};
real_t *receiveBuffer[6][4] = {NULL};
int *sendIbuffer[6][IMAXVARS] = {NULL};

int *sppd[3] = {NULL};
int *rppd[3] = {NULL};

real_t *scenter[3] = {NULL};
real_t *rcenter[3] = {NULL};

extern void adjustHaloCenters(SimFlat* sim);

/**
 * Determine halo box ids **/
int* getHaloBoxIds(SimFlat* sim) {

   int nx, ny, nz, nbx;
   int ioff;
   int iXtra;
   int* haloBoxIds=NULL;
   int* idMinus;
   int* idPlus;

   // Sim is empty
   if ( !sim ) 
   {
      if ( haloBoxIds ) suAlignedFree(haloBoxIds);
      haloBoxIds = NULL;
      return NULL;
   }

   // Collect halo box Ids
   if (!haloBoxIds) {

      nx = sim->lNbx[0];
      ny = sim->lNbx[1];
      nz = sim->lNbx[2];
      nbx = sim->nBoxes;
      iXtra = sim->nHaloBoxes;

      haloBoxIds = suAlignedCalloc(iXtra*sizeof(int));

      /**
      * XPlus and XMinus faces **/
      idMinus = haloBoxIds;
      idPlus = idMinus + sim->haloSizes[0];

      for (int iz = 0; iz < nz; iz++) 
      {
         for (int iy = 0; iy < ny; iy++, idMinus++, idPlus++) 
         {
            *idMinus = getIBoxFromIxyz(sim, 0, iy, iz);
            *idPlus = getIBoxFromIxyz(sim, nx-1, iy, iz);
         }
      }

      /**
       * YPlus and YMinus faces **/
      ioff = nx*ny-nx;
      idMinus = haloBoxIds + 2*sim->haloSizes[0];
      idPlus = idMinus + sim->haloSizes[1];
      for (int iz=0; iz<nz; iz++) 
      {

         /**
          * Contribution from XMinus face **/
         *idMinus = getIBoxFromIxyz(sim, -1, 0, iz);
         idMinus++;

         *idPlus = getIBoxFromIxyz(sim, -1, ny-1, iz);
         idPlus++;

         /**
            * Contribution from local faces **/
         for (int ix=0; ix<nx; ix++,idMinus++,idPlus++)
         {
            *idMinus = getIBoxFromIxyz(sim, ix, 0, iz);
	         *idPlus = getIBoxFromIxyz(sim, ix, ny-1, iz);
         }

         /**
            * Contribution from XPlus face **/
         *idMinus =  getIBoxFromIxyz(sim, nx, 0, iz); //iz*ny; 
         idMinus++;

         *idPlus =  getIBoxFromIxyz(sim, nx, ny-1, iz); //iz*ny + ny - 1; 
         idPlus++;

      }

      /**
       * ZPlus and ZMinus faces **/
      idMinus = haloBoxIds + 2*sim->haloSizes[0] + 2*sim->haloSizes[1];
      idPlus = idMinus + sim->haloSizes[2];

      /**
       * Contribution from YMinus face **/
      ioff = (nx+2)*(nz-1);
      for (int ix=-1; ix<nx+1; ix++, idMinus++, idPlus++)
      {
         *idMinus = getIBoxFromIxyz(sim, ix, -1, 0);
         *idPlus = getIBoxFromIxyz(sim, ix, -1, nz-1);
      }
      
      ioff = nbx - nx*ny;
      for (int iy=0; iy<ny; iy++) 
      {

         /**
          * Contribution from XMinus face **/
         *idMinus = getIBoxFromIxyz(sim, -1, iy, 0);
         idMinus++;

         *idPlus = getIBoxFromIxyz(sim, -1, iy, nz-1);
         idPlus++;

         /**
          * Contribution from local faces **/
         for (int ix=0; ix<nx; ix++, idMinus++, idPlus++)
         {
            *idMinus = getIBoxFromIxyz(sim, ix, iy, 0);
            *idPlus = getIBoxFromIxyz(sim, ix, iy, nz-1);
         }

         /**
          * Contribution from XPlus face **/
         *idMinus =  getIBoxFromIxyz(sim, nx, iy, 0);
         idMinus++;

         *idPlus =  getIBoxFromIxyz(sim, nx, iy, nz-1); 
         idPlus++;
      }

      /**
       * Contribution from YPlus face **/
      ioff = sim->haloOffsets[1] + sim->haloSizes[1];
      for (int ix=-1; ix<nx+1; ix++, idMinus++, idPlus++)
      {
         *idMinus =  getIBoxFromIxyz(sim, ix, ny, 0 );
         *idPlus =  getIBoxFromIxyz(sim, ix, ny, nz-1);
      }

      if(0) 
      {
         FILE* fp;
         int ioff, ioff2;
         int iface;
         int ix = 0;
         int iBox, hbox;
         char fname[80];

         sprintf(fname, "halo_boxes%d.txt", getMyParallel());
         fp = fopen(fname,"w");
         fprintf(fp, "%10d %10d %d %d (%d %d %d) \n",ix,haloBoxIds[ix], sim->nHaloBoxes, sim->nBoxes, nx, ny, nz);
         fflush(stdout);

         hbox = sim->nBoxes + sim->haloSizes[0];

         iface = 0;
         ioff2 = sim->haloOffsets[iface]-nbx;
         ioff = sim->haloSizes[iface];
         for (ix=ioff2; ix<ioff2+sim->haloSizes[iface]; ix++)
         {
            iBox = haloBoxIds[ix];
            fprintf(fp,"to minus face %d %10d %10d dc %f %f %f hbox %d\n",iface,ix,iBox, 
               sim->dCenter[iBox][0], sim->dCenter[iBox][1], sim->dCenter[iBox][2],
               hbox);
            hbox++;
         }

         hbox = sim->nBoxes;

         for (ix=ioff2+sim->haloSizes[iface]; ix<ioff2+2*sim->haloSizes[iface]; ix++) 
         {
            iBox = haloBoxIds[ix];
            fprintf(fp,"to plus face %d %10d %10d dc %f %f %f hbox %d\n",iface,ix,iBox,
               sim->dCenter[iBox][0], sim->dCenter[iBox][1], sim->dCenter[iBox][2],
               hbox);
            hbox++;
         }

         hbox = sim->nBoxes + 2*sim->haloSizes[0] + sim->haloSizes[1];

         iface = 1;
         ioff2 = sim->haloOffsets[iface]-nbx;
         ioff = sim->haloSizes[iface];
         for (ix=ioff2; ix<ioff2+sim->haloSizes[iface]; ix++) 
         {
            iBox = haloBoxIds[ix];
            fprintf(fp,"to minus face %d %10d %10d dc %f %f %f hbox %d\n", iface, ix, iBox, 
               sim->dCenter[iBox][0], sim->dCenter[iBox][1], sim->dCenter[iBox][2],
               hbox);
            hbox++;
         }

         hbox = sim->nBoxes + 2*sim->haloSizes[0];

         for (ix=ioff2+sim->haloSizes[iface]; ix<ioff2+2*sim->haloSizes[iface]; ix++)
         {
            iBox = haloBoxIds[ix];
            fprintf(fp,"to plus face %d %10d %10d dc %f %f %f hbox %d\n", iface, ix, iBox,
               sim->dCenter[iBox][0], sim->dCenter[iBox][1], sim->dCenter[iBox][2],
               hbox);
            hbox++;
         }

         hbox = sim->nBoxes + 2 * sim->haloSizes[0] + 2 * sim->haloSizes[1] + sim->haloSizes[2];
 
         iface = 2;
         ioff2 = sim->haloOffsets[iface]-nbx;
         ioff = sim->haloSizes[iface];
         for (ix=ioff2; ix<ioff2+sim->haloSizes[iface]; ix++) 
         {
            iBox = haloBoxIds[ix];
            fprintf(fp,"to minus face %d %10d %10d dc %f %f %f hbox %d\n", iface, ix,iBox, 
               sim->dCenter[iBox][0], sim->dCenter[iBox][1], sim->dCenter[iBox][2],
               hbox);
            hbox++;
         }

         hbox = sim->nBoxes + 2 * sim->haloSizes[0] + 2 * sim->haloSizes[1];

         for (ix=ioff2+sim->haloSizes[iface]; ix<ioff2+2*sim->haloSizes[iface]; ix++) 
         {
           iBox = haloBoxIds[ix];
           fprintf(fp,"to plus face %d %10d %10d dc %f %f %f hbox %d\n", iface, ix,iBox,
              sim->dCenter[iBox][0], sim->dCenter[iBox][1], sim->dCenter[iBox][2],
              hbox);
           hbox++;
         }

         fclose(fp);
      }
   }

   return haloBoxIds;
}

void loadBuffer(SimFlat* sim, int axis) 
{

   int i, j, k, im;
   int m, n, m3;
   int mj, mk, ij, ik;
   int iOff, i2;
   int* d;
   int fminus, fplus;
   int* haloBoxId;

   fminus = fLookup[axis][0];
   fplus = fLookup[axis][1];

   n = MAXATOMS*sizeof(int);
   m = MAXATOMS*sizeof(real_t);

   haloBoxId = getHaloBoxIds(sim);

   im = 0;
   iOff = sim->haloOffsets[axis] - sim->nBoxes;
   i2 = sim->haloSizes[axis];

   for (i=0; i<sim->haloSizes[axis]; i++, i2++, im+=MAXATOMS) 
   {

      // Work with plus and minus side of axis
      j = haloBoxId[i+iOff];
      k = haloBoxId[i2+iOff];

      sppd[axis][i] = sim->nAtoms[j];
      sppd[axis][i2] = sim->nAtoms[k];

      //printf("%d j %d %d k %d %d\n", getMyParallel(), j, sim->nAtoms[j], k, sim->nAtoms[k]);

      j = MAXATOMS*j;
      k = MAXATOMS*k;

      // Copy id for each atom 
      memcpy(sendIbuffer[fminus][0]+im,sim->id+j,n);
      memcpy(sendIbuffer[fplus][0]+im,sim->id+k,n);

      // Copy type for each atom
      memcpy(sendIbuffer[fminus][1]+im,sim->iType+j,n);
      memcpy(sendIbuffer[fplus][1]+im,sim->iType+k,n);
          
      // Copy positions - x,y,z for each atom
      for (int ll = 0; ll < 3; ll++) 
      {
         for (int l = 0; l < MAXATOMS; l++) 
         {
            sendBuffer[fminus][ll][im+l] = sim->r[j+l][ll];
         }
      }
      //printf("%d For Minus box %d r %f %f %f\n", getMyParallel(), j/MAXATOMS, sim->r[j][0], sim->r[j][1], sim->r[j][2]);

      for (int ll = 0; ll < 3; ll++) 
      {
         for (int l = 0; l < MAXATOMS; l++) 
         {
            sendBuffer[fplus][ll][im+l] = sim->r[k+l][ll];
         }
      }

      //printf("%d For Plus box %d r %f %f %f\n", getMyParallel(), k/MAXATOMS, sim->r[k][0], sim->r[k][1], sim->r[k][2]);

      // Copy momenta x, y, z for each atom
      for (int ll = 0; ll < 3; ll++) 
      {
         for (int l = 0; l < MAXATOMS; l++) 
         {
            sendBuffer[fminus][3+ll][im+l] = sim->p[j+l][ll];
         }
      }

      for (int ll = 0; ll < 3; ll++) 
      {
         for (int l = 0; l < MAXATOMS; l++) 
         {
            sendBuffer[fplus][3+ll][im+l] = sim->p[k+l][ll];
         }
      }

      // Copy force for each atom
      for (int ll = 0; ll < 4; ll++) 
      {
         for (int l = 0; l < MAXATOMS; l++) 
         {
            sendBuffer[fminus][6+ll][im+l] = sim->f[j+l][ll];
         }
      }

      for (int ll = 0; ll < 4; ll++) 
      {
         for (int l = 0; l < MAXATOMS; l++) 
         {
            sendBuffer[fplus][6+ll][im+l] = sim->f[k+l][ll];
         }
      }

      // Copy mass, fi rho for each atom
      memcpy(sendBuffer[fminus][10]+im, sim->mass+j, m);
      memcpy(sendBuffer[fplus][10]+im, sim->mass+k, m);

      memcpy(sendBuffer[fminus][11]+im, sim->fi+j, m);
      memcpy(sendBuffer[fplus][11]+im, sim->fi+k, m);

      memcpy(sendBuffer[fminus][12]+im, sim->rho+j, m);
      memcpy(sendBuffer[fplus][12]+im, sim->rho+k, m);

   }
  
   return;
}

int initHalo(SimFlat* sim) 
{
   if ( ! sendBuffer[0][0] ) 
   {
      real_t* b;
      real_t* rb;
      int* bi;
      int i, n, n3;
      int iFace;
      int iXtra;

      iXtra = sim->nHaloBoxes;
      //printf("Number of halo boxes = %d\n", iXtra);
      rb = suAlignedCalloc(4*MAXATOMS*iXtra*sizeof(real_t));
      b = suAlignedCalloc(MAXVARS*MAXATOMS*iXtra*sizeof(real_t));
      bi = suAlignedCalloc(IMAXVARS*MAXATOMS*iXtra*sizeof(int));
      sppd[0] = suAlignedCalloc(iXtra*sizeof(int));
      rppd[0] = suAlignedCalloc(iXtra*sizeof(int));
      scenter[0] = suAlignedCalloc(3*iXtra*sizeof(real_t));
      rcenter[0] = suAlignedCalloc(3*iXtra*sizeof(real_t));

      if ( ! (rb && b && bi && scenter[0] && rcenter[0] && sppd[0] && rppd[0]) ) return 1;

      sim->haloSizes[0] = sim->lNbx[1]*sim->lNbx[2];
      sim->haloSizes[1] = (sim->lNbx[0]+2)*sim->lNbx[2];
      sim->haloSizes[2] = (sim->lNbx[0]+2)*(sim->lNbx[1]+2);
      //printf("haloSizes = %d %d %d\n", sim->haloSizes[0], sim->haloSizes[1], sim->haloSizes[2]);

      sim->haloOffsets[0] = sim->nBoxes;
      sim->haloOffsets[1] = sim->nBoxes+2*sim->haloSizes[0];
      sim->haloOffsets[2] = sim->nBoxes+2*sim->haloSizes[0]+2*sim->haloSizes[1];
      //printf("haloOffsets = %d %d %d\n", sim->haloOffsets[0], sim->haloOffsets[1], sim->haloOffsets[2]);
    
      n = MAXATOMS*sim->haloSizes[0];
      for (i=0; i<4; i++,rb+=n) 
      {
         receiveBuffer[FACE_X_MINUS][i] = rb;
      }
      for (i=0; i<4; i++,rb+=n) 
      {
         receiveBuffer[FACE_X_PLUS][i] = rb;
      }

      for (i=0; i<MAXVARS; i++,b+=n) 
      {
         sendBuffer[FACE_X_MINUS][i] = b;
      }
      for (i=0; i<MAXVARS; i++,b+=n) 
      {
         sendBuffer[FACE_X_PLUS][i] = b;
      }

      for (i=0; i<IMAXVARS; i++,bi+=n) 
      {
         sendIbuffer[FACE_X_MINUS][i] = bi;
      }
      for (i=0; i<IMAXVARS; i++,bi+=n) 
      {
         sendIbuffer[FACE_X_PLUS][i] = bi;
      }

      n = MAXATOMS*sim->haloSizes[1];
      sppd[1] = sppd[0] + 2*sim->haloSizes[0];
      rppd[1] = rppd[0] + 2*sim->haloSizes[0];

      scenter[1] = scenter[0] + 3*2*sim->haloSizes[0];
      rcenter[1] = rcenter[0] + 3*2*sim->haloSizes[0];

      for (i=0; i<4; i++,rb+=n) 
      {
         receiveBuffer[FACE_Y_MINUS][i] = rb;
      }
      for (i=0; i<4; i++,rb+=n) 
      {
         receiveBuffer[FACE_Y_PLUS][i] = rb;
      }

      for (i=0; i<MAXVARS; i++,b+=n) 
      {
         sendBuffer[FACE_Y_MINUS][i] = b;
      }
      for (i=0; i<MAXVARS; i++,b+=n) 
      {
         sendBuffer[FACE_Y_PLUS][i] = b;
      }
    
      for (i=0; i<IMAXVARS; i++,bi+=n) 
      {
         sendIbuffer[FACE_Y_MINUS][i] = bi;
      }
      for (i=0; i<IMAXVARS; i++,bi+=n) 
      {
         sendIbuffer[FACE_Y_PLUS][i] = bi;
      }

      n = MAXATOMS*sim->haloSizes[2];
      sppd[2] = sppd[1] + 2*sim->haloSizes[1];
      rppd[2] = rppd[1] + 2*sim->haloSizes[1];

      scenter[2] = scenter[1] + 3*2*sim->haloSizes[1];
      rcenter[2] = rcenter[1] + 3*2*sim->haloSizes[1];

      for (i=0; i<4; i++,rb+=n) 
      {
         receiveBuffer[FACE_Z_MINUS][i] = rb;
      }
      for (i=0; i<4; i++,rb+=n) 
      {
         receiveBuffer[FACE_Z_PLUS][i] = rb;
      }

      for (i=0; i<MAXVARS; i++,b+=n) 
      {
         sendBuffer[FACE_Z_MINUS][i] = b;
      }
      for (i=0; i<MAXVARS; i++,b+=n) 
      {
         sendBuffer[FACE_Z_PLUS][i] = b;
      }

      for (i=0; i<IMAXVARS; i++,bi+=n) 
      {
         sendIbuffer[FACE_Z_MINUS][i] = bi;
      }
      for (i=0; i<IMAXVARS; i++,bi+=n) 
      {
         sendIbuffer[FACE_Z_PLUS][i] = bi;
      }

      /*
       * do periodic boundary data **/

      memset(shiftR,0,6*sizeof(double));

      /**
       * The X_Minus face **/
      iFace = FACE_X_MINUS;
      procAddress[iFace] = processorNum(sim, -1,0,0);
      shiftFlags[iFace] = (((sim->minLBound[0] - sim->lSize[0]) < sim->minGBound[0])?1:0);
      shiftR[iFace] = sim->maxGBound[0] - sim->minGBound[0];

      /**
       * The X_Plus face **/
      iFace = FACE_X_PLUS;
      procAddress[iFace] = processorNum(sim, 1,0,0);
      shiftFlags[iFace] = (((sim->maxLBound[0] + sim->lSize[0]) > sim->maxGBound[0])?1:0);
      shiftR[iFace] = -(sim->maxGBound[0] - sim->minGBound[0]);
    
      /**
       * The Y_Minus face **/
      iFace = FACE_Y_MINUS;
      procAddress[iFace] = processorNum(sim, 0,-1,0);
      shiftFlags[iFace] = (((sim->minLBound[1] - sim->lSize[1]) < sim->minGBound[1])?1:0);
      shiftR[iFace] = sim->maxGBound[1] - sim->minGBound[1];

      /**
       * The Y_Plus face **/
      iFace = FACE_Y_PLUS;
      procAddress[iFace] = processorNum(sim, 0,1,0);
      shiftFlags[iFace] = (((sim->maxLBound[1] + sim->lSize[1]) > sim->maxGBound[1])?1:0);
      shiftR[iFace] = -(sim->maxGBound[1] - sim->minGBound[1]);
    
      /**
       * The Z_Minus face **/
      iFace = FACE_Z_MINUS;
      procAddress[iFace] = processorNum(sim, 0,0,-1);
      shiftFlags[iFace] = (((sim->minLBound[2] - sim->lSize[2]) < sim->minGBound[2])?1:0);
      shiftR[iFace] = sim->maxGBound[2] - sim->minGBound[2];

      /**
       * The Z_Plus face **/
      iFace = FACE_Z_PLUS;
      procAddress[iFace] = processorNum(sim, 0,0,1);
      shiftFlags[iFace] = (((sim->maxLBound[2] + sim->lSize[2]) > sim->maxGBound[2])?1:0);
      shiftR[iFace] = -(sim->maxGBound[2] - sim->minGBound[2]);
    
      // Calculate halo centers
      adjustHaloCenters(sim);
   }

   return 0;
}

void sendAndReceiveBuffer(SimFlat* sim, int axis) 
{
   int fminus, minusProc;
   int fplus, plusProc;
   int iszppd, iszppd3;
   int isz, isz3;
   int ioff;
   int* d;
   int* mppd;
   int* pppd;
   real_t* mcenter;
   real_t* pcenter;
   real3* dc;

   isz = sizeof(int)*MAXATOMS*sim->haloSizes[axis];
   isz3 = sizeof(real_t)*MAXATOMS*sim->haloSizes[axis];
   iszppd = sizeof(int)*sim->haloSizes[axis];
   iszppd3 = 3*sizeof(real_t)*sim->haloSizes[axis];

   fminus = fLookup[axis][0];
   fplus = fLookup[axis][1];
   minusProc = procAddress[fminus] & PROC_MASK;
   plusProc = procAddress[fplus] & PROC_MASK;

   //printf("For axis %d proc %d minus proc %d plus proc %d\n", axis, getMyParallel(), minusProc, plusProc);

   /**
    * Send atom counts per halo box
    * Send box centers
    * Send id, type
    * Send the x, y, z
    * Send momenta x, y, z
    * Send force
    * Send mass, fi, rho **/
   sendParallel(minusProc,fminus,sppd[axis],iszppd);

   for (int i = 0; i < IMAXVARS; i++) 
   {
      sendParallel(minusProc,fminus,sendIbuffer[fminus][i],isz);
   }

   for (int i = 0; i < MAXVARS; i++) 
   {
      sendParallel(minusProc,fminus,sendBuffer[fminus][i],isz3);
   }

   // Send to plus proc
   sendParallel(plusProc,fplus,sppd[axis]+sim->haloSizes[axis],iszppd);

   for (int i = 0; i < IMAXVARS; i++) 
   {
      sendParallel(plusProc,fplus,sendIbuffer[fplus][i],isz);
   }

   for (int i = 0; i < MAXVARS; i++) 
   {
      sendParallel(plusProc,fplus,sendBuffer[fplus][i],isz3);
   }

   /**
    * Note the swapped 'm' and 'p'
    * in the indices for the receives
    **/

   /**
    * Receive atom counts per box and centers
    * First insert the ppd into correct locations **/
   mppd = rppd[axis];
   pppd = rppd[axis] + sim->haloSizes[axis];

   receiveBlockingParallel(minusProc,fplus,mppd,iszppd);
   receiveBlockingParallel(plusProc,fminus,pppd,iszppd);

   //printf(" %d %d received ppd %d %d\n", getMyParallel(), axis, minusProc, plusProc);

   // Set box sizes in numbers of particles
   ioff = sim->haloOffsets[axis];
   d = sim->nAtoms + ioff;
   for(int i=0, ir=0; i<sim->haloSizes[axis]; i++, ir=ir+3, d++) 
   {
      d[0] = mppd[i];
      d[sim->haloSizes[axis]] = pppd[i];
      //printf("%d mppd %d pppd %d\n", i, mppd[i], pppd[i]);

   }
    
   /**
    * Receive id, type
    * Receive the x, y, and z 
    * Receive momenta x, y, z
    * Receive force
    * Receive mass, fi, rho
    * 
    * Note the swapped 'm' and 'p'
    **/
   ioff = ioff*MAXATOMS;

   receiveBlockingParallel(minusProc,fplus,sim->id+ioff, isz);
   receiveBlockingParallel(minusProc,fplus,sim->iType+ioff, isz);

   //printf("%d id %d %d %d\n", getMyParallel(), sim->id[ioff], sim->id[ioff+1], sim->id[ioff+2]);

   //printf("%d %d minus received id, type %d\n", getMyParallel(), axis, minusProc);

   for (int i = 0;  i < 3; i++) 
   {
      if (receiveBuffer[fplus][i] != NULL)
         receiveBlockingParallel(minusProc,fplus,receiveBuffer[fplus][i],isz3);
      else
         printf("receive %d from minus %d axis %d r NULL buffer\n", getMyParallel(), minusProc, axis);
   }
   for (int ii = 0; ii < MAXATOMS*sim->haloSizes[axis]; ii++) 
   {
      for (int i = 0; i < 3; i++) 
      {
         sim->r[ioff+ii][i] = receiveBuffer[fplus][i][ii];
      }
      //printf("receive %d %d %d %f %f %f\n", getMyParallel(), ioff, ii, sim->r[ioff+ii][0], sim->r[ioff+ii][1], sim->r[ioff+ii][2]);
   }
   //printf("receive %d %d %f %f %f\n", getMyParallel(), ioff, sim->r[ioff][0], sim->r[ioff][1], sim->r[ioff][2]);  
   //printf("minus received r\n");

   for (int i = 0;  i < 3; i++) 
   {
      if (receiveBuffer[fplus][i] != NULL)
         receiveBlockingParallel(minusProc,fplus,receiveBuffer[fplus][i],isz3);
      else
         printf("%d from minus %d axis %d p NULL buffer\n", getMyParallel(), minusProc, axis);
   } 
   for (int ii = 0; ii < MAXATOMS*sim->haloSizes[axis]; ii++) 
   {
      for (int i = 0; i < 3; i++) 
      {
         sim->p[ioff+ii][i] = receiveBuffer[fplus][i][ii];
      }
   }
   //printf("minus received p\n");

   for (int i = 0;  i < 4; i++) 
   {
      receiveBlockingParallel(minusProc,fplus,receiveBuffer[fplus][i],isz3);
   } 
   for (int ii = 0; ii < MAXATOMS*sim->haloSizes[axis]; ii++) 
   {
      for (int i = 0; i < 4; i++) 
      {
         sim->f[ioff+ii][i] = receiveBuffer[fplus][i][ii];
      }
   }
   //printf("minus received f\n");

   receiveBlockingParallel(minusProc,fplus,sim->mass+ioff, isz3);
   receiveBlockingParallel(minusProc,fplus,sim->fi+ioff, isz3);
   receiveBlockingParallel(minusProc,fplus,sim->rho+ioff, isz3);

   //printf("minus received mass, fi, rho\n");

   /**
    * From plus proc
    * Receive id, type
    * Receive x, y, z
    * Receive momenta x, y, z
    * Receive force
    * Receive mass, fi, rho */
   ioff += MAXATOMS*sim->haloSizes[axis];

   receiveBlockingParallel(plusProc,fminus,sim->id+ioff, isz);
   receiveBlockingParallel(plusProc,fminus,sim->iType+ioff, isz);
   //printf("%d %d plus received id, type %d\n", getMyParallel(), axis, plusProc);

   for (int i = 0;  i < 3; i++) 
   {
      if (receiveBuffer[fminus][i] != NULL)        
         receiveBlockingParallel(plusProc,fminus,receiveBuffer[fminus][i],isz3);
      else
         printf("%d from plus %d axis %d r NULL\n", getMyParallel(), plusProc, axis);
   }
   for (int ii = 0; ii < MAXATOMS*sim->haloSizes[axis]; ii++) 
   {
      for (int i = 0; i < 3; i++) 
      {
         sim->r[ioff+ii][i] = receiveBuffer[fminus][i][ii];
      }
   }
   //printf("plus received r\n");

   for (int i = 0;  i < 3; i++) 
   {
      if (receiveBuffer[fminus][i] != NULL)
         receiveBlockingParallel(plusProc,fminus,receiveBuffer[fminus][i],isz3);
      else
         printf("%d from plus %d axis %d p NULL\n", getMyParallel(), plusProc, axis);
   }
   for (int ii = 0; ii < MAXATOMS*sim->haloSizes[axis]; ii++) 
   {
      for (int i = 0; i < 3; i++) 
      {
         sim->p[ioff+ii][i] = receiveBuffer[fminus][i][ii];
      }
   }
   //printf("plus received p\n");

   for (int i = 0;  i < 4; i++) 
   {
      receiveBlockingParallel(plusProc,fminus,receiveBuffer[fminus][i],isz3);
   }
   for (int ii = 0; ii < MAXATOMS*sim->haloSizes[axis]; ii++) 
   {
      for (int i = 0; i < 4; i++) 
      {
         sim->f[ioff+ii][i] = receiveBuffer[fminus][i][ii];
      }
   }
   //printf("plus received f\n");

   receiveBlockingParallel(plusProc,fminus,sim->mass+ioff, isz3);
   receiveBlockingParallel(plusProc,fminus,sim->fi+ioff, isz3);
   receiveBlockingParallel(plusProc,fminus,sim->rho+ioff, isz3);
   //printf("plus received mass,fi,rho\n");

   barrierParallel();

   return;
}

void adjustHaloCenters(SimFlat* sim) 
{
   real3* mcenter;
   real3* pcenter;
   real3 madjust, padjust;
   int mi, pi;
   real3 bval;
   int ix, iy, iz;

   // Calculate adjustment constants
   bval[1] = 0; bval[2] = 0;
   for (int i = 0; i < 3; i++) 
   {
      madjust[i] = sim->minLBound[i] - sim->boxSize[i];
      padjust[i] = sim->minLBound[i] + sim->lNbx[i] * sim->boxSize[i];
      bval[i] = sim->minLBound[i];
   }

   // Adjust centers for XMinus/XPlus faces
   mcenter = sim->dCenter + sim->nBoxes;
   pcenter = mcenter + sim->haloSizes[0];
   ix = iy = iz = 0;
   for (int iBox = 0; iBox < sim->haloSizes[0]; iBox++, mcenter++, pcenter++) 
   {
      mcenter[0][0] = madjust[0];
      mcenter[0][1] = bval[1] + iy * sim->boxSize[1];
      mcenter[0][2] = bval[2] + iz * sim->boxSize[2];
      pcenter[0][0] = padjust[0];
      pcenter[0][1] = bval[1] + iy * sim->boxSize[1];;
      pcenter[0][2] = bval[2] + iz * sim->boxSize[2];
      iy++;
      if (iy == sim->lNbx[1]) 
      {
         iy=0;
         iz++;
         if (iz == sim->lNbx[2]) iz = 0;
      }
      //printf("X m %f %f %f p %f %f %f\n", mcenter[0][0], mcenter[0][1], mcenter[0][2], pcenter[0][0], pcenter[0][1], pcenter[0][2]);
   }

   // Adjust centers for YMinus/YPlus faces and XMinus/XPlus edges
   mcenter = sim->dCenter + sim->nBoxes + 2*sim->haloSizes[0];
   pcenter = mcenter + sim->haloSizes[1];
   for (int iBox = 0; iBox < sim->lNbx[2]; iBox++) 
   {

      // From XMinus edge
      // Update X per XMinus
      mcenter[0][0] = madjust[0];
      pcenter[0][0] = madjust[0];

      // Update Y
      mcenter[0][1] = madjust[1];
      pcenter[0][1] = padjust[1];

      mcenter[0][2] = bval[2] + iBox * sim->boxSize[2];
      pcenter[0][2] = bval[2] + iBox * sim->boxSize[2];

      //printf("Y XMinus m %f %f %f p %f %f %f\n", mcenter[0][0], mcenter[0][1], mcenter[0][2], pcenter[0][0], pcenter[0][1], pcenter[0][2]);

      mcenter++;
      pcenter++;

      // From YMinus/YPlus faces
      for (int xbox = 0; xbox < sim->lNbx[0]; xbox++, mcenter++, pcenter++)
      {
         mcenter[0][1] = madjust[1];
         mcenter[0][0] = bval[0] + xbox * sim->boxSize[0];
         mcenter[0][2] = bval[2] + iBox * sim->boxSize[2];
         pcenter[0][1] = padjust[1];
         pcenter[0][0] = bval[0] + xbox * sim->boxSize[0];
         pcenter[0][2] = bval[2] + iBox * sim->boxSize[2];
         //printf("Y m %f %f %f p %f %f %f\n", mcenter[0][0], mcenter[0][1], mcenter[0][2], pcenter[0][0], pcenter[0][1], pcenter[0][2]);
      }

      // From XPlus edge
      // Update X
      mcenter[0][0] = padjust[0];
      pcenter[0][0] = padjust[0];

      // Update 
      mcenter[0][1] = madjust[1];
      pcenter[0][1] = padjust[1];

      mcenter[0][2] = bval[2] + iBox * sim->boxSize[2];
      pcenter[0][2] = bval[2] + iBox * sim->boxSize[2];

      //printf("Y XPlus m %f %f %f p %f %f %f\n", mcenter[0][0], mcenter[0][1], mcenter[0][2], pcenter[0][0], pcenter[0][1], pcenter[0][2]);

      mcenter++;
      pcenter++;
   } 

   // Adjust centers for Z face, YMinus/YPlus edges, XMinus/XPlus edges, and corners
   mcenter = sim->dCenter + sim->nBoxes + 2 * sim->haloSizes[0] + 2 * sim->haloSizes[1];
   pcenter = mcenter + sim->haloSizes[2];

   // From YMinus edge
   for (int iBox = 0; iBox < (2 + sim->lNbx[0]); iBox++, mcenter++, pcenter++) 
   {
      // Update Y
      mcenter[0][1] = madjust[1];
      pcenter[0][1] = madjust[1];

      // Update Z
      mcenter[0][2] = madjust[2];
      pcenter[0][2] = padjust[2];

      mcenter[0][0] = bval[0] + (iBox-1) * sim->boxSize[0];
      pcenter[0][0] = bval[0] + (iBox-1) * sim->boxSize[0];

      //printf("Z YMinus m %f %f %f p %f %f %f\n", mcenter[0][0], mcenter[0][1], mcenter[0][2], pcenter[0][0], pcenter[0][1], pcenter[0][2]);
   }

   for (int iBox = 0; iBox < sim->lNbx[1]; iBox++)    
   {

      // From XMinus edge
      // Update X
      mcenter[0][0] = madjust[0];
      pcenter[0][0] = madjust[0];

      // Update Z
      mcenter[0][2] = madjust[2];
      pcenter[0][2] = padjust[2];

      mcenter[0][1] = bval[1] + iBox * sim->boxSize[1];
      pcenter[0][1] = bval[1] + iBox * sim->boxSize[1];
      //printf("Z XMinus m %f %f %f p %f %f %f\n", mcenter[0][0], mcenter[0][1], mcenter[0][2], pcenter[0][0], pcenter[0][1], pcenter[0][2]);

      mcenter++;
      pcenter++;

      // From Z face
      for (int xbox = 0; xbox < sim->lNbx[0]; xbox++, mcenter++, pcenter++) 
      {
         mcenter[0][2] = madjust[2];
         pcenter[0][2] = padjust[2];

         mcenter[0][0] = bval[0] + xbox * sim->boxSize[0];
         pcenter[0][0] = bval[0] + xbox * sim->boxSize[0];

         mcenter[0][1] = bval[1] + iBox * sim->boxSize[1];
         pcenter[0][1] = bval[1] + iBox * sim->boxSize[1];
         //printf("Z m %f %f %f p %f %f %f\n", mcenter[0][0], mcenter[0][1], mcenter[0][2], pcenter[0][0], pcenter[0][1], pcenter[0][2]);
      }

      // From XPlus edge
      // Update X
      mcenter[0][0] = padjust[0];
      pcenter[0][0] = padjust[0];

      // Update Z
      mcenter[0][2] = madjust[2];
      pcenter[0][2] = padjust[2];

      mcenter[0][1] = bval[1] + iBox * sim->boxSize[1];
      pcenter[0][1] = bval[1] + iBox * sim->boxSize[1];
      //printf("Z XPlus m %f %f %f p %f %f %f\n", mcenter[0][0], mcenter[0][1], mcenter[0][2], pcenter[0][0], pcenter[0][1], pcenter[0][2]);

      mcenter++;
      pcenter++;
   }

   // From YPlus edge
   for (int iBox = 0; iBox < (2 + sim->lNbx[0]); iBox++, mcenter++, pcenter++) 
   {
      // Update Y
      mcenter[0][1] = padjust[1];
      pcenter[0][1] = padjust[1];

      // Update Z
      mcenter[0][2] = madjust[2];
      pcenter[0][2] = padjust[2];

      mcenter[0][0] = bval[0] + (iBox-1) * sim->boxSize[0];
      pcenter[0][0] = bval[0] + (iBox-1) * sim->boxSize[0];

      //printf("Z YPlus m %f %f %f p %f %f %f\n", mcenter[0][0], mcenter[0][1], mcenter[0][2], pcenter[0][0], pcenter[0][1], pcenter[0][2]);
   }

}

void checkLocalBoxes(SimFlat* sim) 
{
   real3 pr;
   real_t xtemp;
   int xyz[3];
   int tbox;
   int iOff;

   for (int iBox = 0; iBox < sim->nBoxes; iBox++) 
   {
      iOff = iBox * MAXATOMS;

      for (int i = 0; i < sim->nAtoms[iBox]; i++, iOff++) 
      {
         for (int j = 0; j < 3; j++) 
         {
            pr[j] = sim->r[iOff][j] + sim->dCenter[iBox][j];
         }

         for (int j = 0; j < 3; j++) 
         {
            xtemp = pr[j] - sim->minLBound[j];
            if (xtemp < 0 && xtemp >= -sim->boxSize[j])
              xyz[j] = -1;
            else
              xyz[j] = xtemp / sim->boxSize[j];
         }

         tbox = getIBoxFromIxyz3(sim, xyz);

         if (iBox != tbox)
            printf("%d Local iBox %d %d %d r %f %f %f dCenter %f %f %f xyz %d %d %d tbox %d\n", getMyParallel(), iBox, i, iOff, pr[0], pr[1], pr[2], sim->dCenter[iBox][0], sim->dCenter[iBox][1], sim->dCenter[iBox][2], xyz[0], xyz[1], xyz[2], tbox);

         for (int j = 0; j < 3; j++) 
         {
            if (pr[j] < sim->minLBound[j] || pr[j] >= sim->maxLBound[j]) 
            {
               printf("%d iBox %d atom %d out of box\n", getMyParallel(), iBox, i);
            }
         }
      }
   }
}

void checkHaloBoxes(SimFlat* sim) 
{
   real3 pr;
   real_t xtemp;
   int xyz[3];
   int tbox;
   int iOff;

   for (int axis = 0; axis < 3; axis++) 
   {
/*
      printf("axis %d start %d end %d\n", axis, sim->haloOffsets[axis],
         sim->haloOffsets[axis]+2*sim->haloSizes[axis]);
*/
      for (int iBox = sim->haloOffsets[axis]; iBox < sim->haloOffsets[axis]+2*sim->haloSizes[axis]; iBox++)
      {
         iOff = iBox *MAXATOMS;
/*
         printf("iBox %d natoms %d dCenter %f %f %f\n", iBox, sim->nAtoms[iBox],  sim->dCenter[iBox][0],
            sim->dCenter[iBox][1], sim->dCenter[iBox][2]);
*/
         for (int i = 0; i < sim->nAtoms[iBox]; i++, iOff++)
         {
            for (int j = 0; j < 3; j++)
            {
               pr[j] = sim->r[iOff][j] + sim->dCenter[iBox][j];
            }

            for (int j = 0; j < 3; j++)
            {
               xtemp = pr[j] - sim->minLBound[j];
               if (xtemp < 0 && xtemp >= -sim->boxSize[j]) 
                  xyz[j] = -1;
               else
                  xyz[j] = xtemp / sim->boxSize[j];
            }

               tbox = getIBoxFromIxyz3(sim, xyz);
                
               if (getMyParallel() == 0 && iBox != tbox)
                  printf("%d Halo iBox %d %d r %f %f %f dCenter %f %f %f xyz %d %d %d tbox %d\n", getMyParallel(), iBox, i, pr[0], pr[1], pr[2], sim->dCenter[iBox][0], sim->dCenter[iBox][1], sim->dCenter[iBox][2], xyz[0], xyz[1], xyz[2], tbox);

         }
      }
   }
}

void haloExchange(SimFlat* sim) 
{
   int n;

   if(initHalo(sim)) 
   {
      fprintf(stderr,
	      "\n\n"
	      "Unable to allocate space in haloExchange()\n"
	      "\n\n"
	      );
      exit(1);
   }

   // Face X
   loadBuffer(sim, 0);
   sendAndReceiveBuffer(sim, 0);

   // Face Y
   loadBuffer(sim, 1);
   sendAndReceiveBuffer(sim, 1);

   // Face Z
   loadBuffer(sim, 2);
   sendAndReceiveBuffer(sim, 2);

   barrierParallel();

   return;
}

void printBoxes(SimFlat* sim) 
{
   if (0) 
   {
      FILE* fp;
      int rindex;
      char fname[80];

      sprintf(fname, "local_data%d_%d.txt", step, getMyParallel());
      fp = fopen(fname,"w");

      for (int iBox = 0; iBox < sim->nBoxes; iBox++)
      {
         fprintf(fp, "%f %f %f Box %d atoms %d center %f %f %f :\n", sim->dCenter[iBox][0], sim->dCenter[iBox][1], sim->dCenter[iBox][2], iBox, sim->nAtoms[iBox], sim->dCenter[iBox][0], sim->dCenter[iBox][1], sim->dCenter[iBox][2]);
       
         rindex = iBox * MAXATOMS;

         for (int j = 0; j < sim->nAtoms[iBox]; j++, rindex++) 
         {
            fprintf(fp, "  %f %f %f  atom %d: r %f %f %f f %f %f %f p %f %f %f\n", 
               sim->dCenter[iBox][0], sim->dCenter[iBox][1], sim->dCenter[iBox][2], j,
               sim->r[rindex][0], sim->r[rindex][1], sim->r[rindex][2], 
               sim->f[rindex][0], sim->f[rindex][1], sim->f[rindex][2], 
               sim->p[rindex][0], sim->p[rindex][1], sim->p[rindex][2]); 
         } 
      }

      fclose(fp); 

      sprintf(fname, "halo_data%d_%d.txt", step, getMyParallel());
      fp = fopen(fname,"w");

      for (int iBox = sim->nBoxes; iBox < sim->nTotalBoxes; iBox++) 
      {
         fprintf(fp, "Box %d atoms %d center %f %f %f :\n", iBox, sim->nAtoms[iBox], sim->dCenter[iBox][0], sim->dCenter[iBox][1], sim->dCenter[iBox][2]);

         rindex = iBox * MAXATOMS;

         for (int j = 0; j < sim->nAtoms[iBox]; j++, rindex++) 
         {
            fprintf(fp, "  atom %d: r %f %f %f f %f %f %f p %f %f %f\n", j,
               sim->r[rindex][0], sim->r[rindex][1], sim->r[rindex][2],
               sim->f[rindex][0], sim->f[rindex][1], sim->f[rindex][2],
               sim->p[rindex][0], sim->p[rindex][1], sim->p[rindex][2]);
         }
      }

      fclose(fp);

      step++;
   }

   return;
}

int boxLocation(int obox, int j, int bindex, SimFlat* sim) 
{
   int iBox, ix, iy, iz;
   int xyz[3];
   real3 prAdjusted;
   real_t xtemp;

   // Adjust position values by center
   // Check if particle is on processor
   for (int i = 0; i < 3; i++) 
   {
      prAdjusted[i] = sim->r[bindex+j][i] + sim->dCenter[obox][i];

      // Return if off processor
      if ( (prAdjusted[i] <= sim->minLBound[i]) || 
           (prAdjusted[i] > sim->maxLBound[i]) )
         return -1;

      // Calculate contribution to local box number
      xtemp = prAdjusted[i] - sim->minLBound[i];
      if (xtemp < 0 && xtemp >= -sim->boxSize[i])
         xyz[i] = -1;
      else
         xyz[i] = xtemp / sim->boxSize[i];
   }

   iBox = getIBoxFromIxyz3(sim, xyz);

   if (iBox != obox)
   {
      printf("from box %d %d r %f %f %f dCenter %f %f %f to xyz %d %d %d box %d center %f %f %f\n", obox, j, sim->r[bindex+j][0], sim->r[bindex+j][1], sim->r[bindex+j][2], sim->dCenter[obox][0], sim->dCenter[obox][1], sim->dCenter[obox][2], xyz[0], xyz[1], xyz[2], iBox, sim->dCenter[iBox][0], sim->dCenter[iBox][1], sim->dCenter[iBox][2]);
   }   

   return iBox;
}

void redistRemotePeriodic(SimFlat* sim) 
{
   /**
    * gets the remote skin layer and pulls
    * in particles that belong to us.
      **/
   int j, whichBox;
   int index;
   real3* pr;

   haloExchange(sim);

   //printf("h = %d b = %d tot = %d\n", sim->nHaloBoxes, sim->nBoxes, sim->nTotalBoxes);

   index =  MAXATOMS*(sim->nBoxes);

   for(int iBox = sim->nBoxes; iBox<sim->nTotalBoxes; iBox++, index+=MAXATOMS)
   {

      //    intf("dCenter %f %f %f box %d\n", sim->dCenter[iBox][0], sim->dCenter[iBox][1], sim->dCenter[iBox][2], iBox);

      j = 0;
      while (j < sim->nAtoms[iBox]) 
      {
         whichBox = boxLocation(iBox, j, index, sim);

         //fprintf(stdout, "Add to %d from box %d %d\n", whichBox, iBox, j);
         //printf("from box %d %d pr %f %f %f dCenter %f %f %f to box %d center %f %f %f\n", iBox, j, pr[j][0], pr[j][1], pr[j][2], sim->dCenter[iBox][0], sim->dCenter[iBox][1], sim->dCenter[iBox][2], whichBox, sim->dCenter[whichBox][0], sim->dCenter[whichBox][1], sim->dCenter[whichBox][2]);

         if (iBox != whichBox && whichBox >= 0 && whichBox < sim->nTotalBoxes) 
         {
            for (int i = 0; i < 3; i++) 
            {
               sim->r[index+j][i] += sim->dCenter[iBox][i] - sim->dCenter[whichBox][i];
            }

            printf("%d Atom %d moving from box %d to box %d atoms %d newr %f %f %f\n", getMyParallel(), j, iBox,
 whichBox, sim->nAtoms[whichBox], sim->r[index+j][0], sim->r[index+j][1], sim->r[index+j][2]);

            moveAtom(sim, j, iBox, whichBox);
         }
         else 
         {
            j++;
         }

      }
   }

}

void redistribute(SimFlat* sim) 
{
  printBoxes(sim);

  // redistribute local particles
  reBoxAll(sim);

  /**
   * redistribute remote particles **/
  redistRemotePeriodic(sim);

  /**
   * remove local particles **/
  reBoxAll2(sim);

  return;
}
