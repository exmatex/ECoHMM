#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "eamTypes.h"
#include "mytype.h"
#include "pmd.h"
#include "cheby.h"

#define DIAG_CHEBY 0
#define PI 3.141592653589793

// Based on the OCL version
//different call signature than the regular version
static void eamInterpolateDerivlocal(real_t r,
      real_t* values,
      int nValues,
      real_t *range,
      real_t *value1, 
      real_t *f1)
{
   int i1;
   int i;
   real_t gi, gi1;

   // using different data struct here
   real_t x0 = range[0];
   real_t xn = range[1];
   real_t invDx = range[2];

   // identical to Sriram's loop in eam.c
   if ( r<x0) r = x0;
   else if (r>xn) r = xn;

   r = (r-x0)*(invDx) ;
   i1 = (int)floor(r);

   /* reset r to fractional distance */
   r = r - floor(r);

   gi  = values[i1+1] - values[i1-1];
   gi1 = values[i1+2] - values[i1];


   // values[i1-1] is guaranteed(?) inbounds because 
   // a->x0 = x0 + (xn-x0)/(double)n; 
   // appears in allocPotentialArray
   *value1 = values[i1] + 0.5*r*(
	 r*(values[i1+1]+ values[i1-1] -2.0*values[i1]) +
	 gi
	 );
   if (*value1 > 1.0e10) *value1 = 0.0;
   if(i1<=0) 
      *f1 = 0.0;
   else 
      *f1 = 0.5*(gi + r*(gi1-gi))*invDx;

   return;

}

static inline void eamInterpolateDeriv(struct PotentialArray *a, real_t r, int iType, int jType, real_t *value1, real_t *f1)
{
   /**
    *
    * This routine will not crash if r is out of range.
    *
    * if ( r < a->x0) r = a->x0;
    * if ( r > a->xn)   r = a->xn;
    **/

   int i1;
   real_t gi, gi1;

   if ( r<a->x0) r = a->x0;
   else if (r>a->xn) r = a->xn;

   r = (r-a->x0)*(a->invDx) ;
   i1 = (int)floor(r);

   /* reset r to fractional distance */
   r = r - floor(r);

   gi  = a->values[i1+1] - a->values[i1-1];
   gi1 = a->values[i1+2] - a->values[i1];


   *value1 = a->values[i1] + 0.5*r*(
	 r*(a->values[i1+1]+ a->values[i1-1] -2.0*a->values[i1]) +
	 gi
	 );
   if(i1<=0) *f1 = 0.0;
   else *f1 = 0.5*(gi + r*(gi1-gi))*a->invDx;

   return;
}

/** Given an EAM potential, generate the corresponding Chebychev approximation 
 * with n coefficients **/
struct EamCheby *setChebPot(EamPotential *pot, int n)
{
   //printf("Generating Chebychev coefficients:");

   struct EamCheby *retCheb;

   retCheb = (EamCheby*)malloc(sizeof(EamCheby));

   //printf("phi...");
   retCheb->phi = genCheb(pot->phi, n);
   retCheb->dphi = genDer(retCheb->phi);

   //printf("rho...");
   retCheb->rho = genCheb(pot->rho, n);
   retCheb->drho = genDer(retCheb->rho);

   //printf("f...");
   retCheb->f = genCheb(pot->f, n);
   retCheb->df = genDer(retCheb->f);

   //printf("\n");
   return retCheb;
}

/** Given a tabulated potential array, generate the first n coefficients of 
 * the corresponding Chebychev approximation **/
struct PotentialArray *genCheb(PotentialArray *pot, int n)
{
   PotentialArray *chPot = (struct PotentialArray*)malloc(sizeof(struct PotentialArray));
   chPot->x0 = pot->x0;
   chPot->xn = pot->xn;
   chPot->n = n;
   chPot->values = (real_t*)malloc((n)*sizeof(real_t));

   real_t a = pot->x0; 
   real_t b = pot->xn; 
   real_t *c = chPot->values; 
   real_t *values = pot->values;
   int nValues = pot->n;
   real_t invDx = pot->invDx;

#if (DIAG_CHEBY > 0)
   //printf("range is "EMT1", "EMT1", "EMT1"\n", a, b, invDx);
   //printf("nValues = %d\n", nValues);
#endif

   real_t fac,bpa,bma;
   real_t r_dummy;

   real_t range[3];
   range[0] = a;
   range[1] = b;
   range[2] = invDx;

   real_t f[n];
   bma=0.5*(b-a);
   bpa=0.5*(b+a);
   for (int k=0;k<n;k++)
   {
      real_t y=cos(PI*(k+0.5)/n);
#if (DIAG_CHEBY > 0)
      eamInterpolateDeriv(pot, y*bma+bpa, 0, 0, &f[k], &r_dummy);
      //printf("old %d, "EMT1", "EMT1"\n", k, y*bma+bpa, f[k]);
#endif
      eamInterpolateDerivlocal(y*bma+bpa, values, nValues, range, &f[k], &r_dummy);
#if (DIAG_CHEBY > 0)
      //printf("new %d, "EMT1", "EMT1"\n", k, y*bma+bpa, f[k]);
#endif
   }
   fac=2.0/n;
   for (int j=0;j<n;j++)
   {
      double sum=0.0;
      for (int k=0;k<n;k++)
	 sum += f[k]*cos(PI*j*(k+0.5)/n);
      c[j]=fac*sum;
#if (DIAG_CHEBY > 0)
      //printf("%d, "EMT1"\n",j, c[j]);
#endif
   }
   return chPot;
}

/** given n Chebychev coefficients for some F, compute the corresponding coefficients for dF **/
void chDer(real_t a, real_t b, real_t *c, real_t *cder, int n)
{
#if (DIAG_CHEBY > 0)
   //printf("a = %e,  b = %e\n", a, b);
#endif
   double con;
#if (DIAG_CHEBY > 0)
   //printf("cder:\n");
#endif
   cder[n-1]=0.0;
#if (DIAG_CHEBY > 0)
   //printf("%d, %e\n", n-1, cder[n-1]);
#endif
   cder[n-2]=2.0*(n-1)*c[n-1];
#if (DIAG_CHEBY > 0)
   //printf("%d, %e\n", n-2, cder[n-2]);
#endif
   for (int j=n-2;j>0;j--)
   {
      cder[j-1]=cder[j+1]+2*(j)*c[j];
#if (DIAG_CHEBY > 0)
      //printf("%d, %e\n", j-1, cder[j-1]);
#endif
   }
   con=2.0/(b-a);
#if (DIAG_CHEBY > 0)
   //printf("con: %e\n", con);
#endif
   for (int j=0;j<n;j++)
   {
      cder[j]=cder[j]*con;
#if (DIAG_CHEBY > 0)
      //printf("%d, %e, %e\n", j, c[j], cder[j]);
#endif
   }
}

/** Given a Chebychev potential array, generate the corresponding array for the derivative **/
struct PotentialArray *genDer(PotentialArray *ch)
{
   real_t a = ch->x0;
   real_t b = ch->xn;
   real_t* c = ch->values;
   int n = ch->n;

   PotentialArray *chDer = (struct PotentialArray*)malloc(sizeof(struct PotentialArray));
   chDer->x0 = a;
   chDer->xn = b;
   chDer->invDx = 1.0;
   chDer->n = n;
   chDer->values = (real_t*)malloc(n*sizeof(real_t));

   real_t* cder = chDer->values;

#if (DIAG_CHEBY > 0)
   //printf("a = %e,  b = %e\n", a, b);
   //printf("n = %d\n", n);
#endif
   double con;
#if (DIAG_CHEBY > 0)
   //printf("cder:\n");
#endif
   cder[n-1]=0.0;
#if (DIAG_CHEBY > 0)
   //printf("%d, %e\n", n-1, cder[n-1]);
   fflush(stdout);
#endif
   cder[n-2]=2.0*(n-1)*c[n-1];
#if (DIAG_CHEBY > 0)
   //printf("%d, %e\n", n-2, cder[n-2]);
#endif
   for (int j=n-2;j>0;j--)
   {
      cder[j-1]=cder[j+1]+2*(j)*c[j];
#if (DIAG_CHEBY > 0)
      //printf("%d, %e\n", j-1, cder[j-1]);
#endif
   }
   con=2.0/(b-a);
#if (DIAG_CHEBY > 0)
   //printf("con: %e\n", con);
#endif
   for (int j=0;j<n;j++)
   {
      cder[j]=cder[j]*con;
#if (DIAG_CHEBY > 0)
      //printf("%d, %e, %e\n", j, c[j], cder[j]);
#endif
   }
   return chDer;
}

/** given a list of Chebyshev coefficients c, compute the value at x
 * x must be in the range [a, b]
 * Modified the call signature to use the PotentialArray
 * which incorporates a, b, c as x0, xn, values **/
real_t eamCheb(PotentialArray *cheb, real_t x) 
{
   real_t a = cheb->x0;
   real_t b = cheb->xn;
   real_t *c = cheb->values;
   int m = cheb->n;

#if (DIAG_CHEBY > 1)
   //printf("x0 = "EMT1", xn = "EMT1", n = %d \n", a, b, m);
   for (int i=0;i<m;i++) //printf(EMT1"\n", cheb->values[i]);
#endif

   real_t d, dd, sv, y, y2;
   real_t ch;
   if (x < a) x = a;
   if (x > b) b = b;
   /*
      if ((x-a)*(x-b) > 0.0)
      {
      //printf("x not in range in eamCheb, %f\n", x);
      }
      */
   d=0.0;
   dd=0.0;
   y=(2.0*x-a-b)/(b-a);
   y2=2.0*y;
   for (int j=m-1;j>0;j--)
   {
      sv=d;
      d=y2*d-dd+c[j];
      dd=sv;
   }
   ch=y*d-dd+0.5*c[0];
   return ch;
}


/** given a list of Chebyshev coefficients c, compute the value at x 
 * x must be in the range [a, b] **/
real_t chebev(real_t a, real_t b, real_t *c,int m, real_t x) 
{
   real_t d, dd, sv, y, y2;
   real_t ch;
   if ((x-a)*(x-b) > 0.0)
   {
      //printf("x not in range in chebev, %f\n", x);
   }
   d=0.0;
   dd=0.0;
   y=(2.0*x-a-b)/(b-a);
   y2=2.0*y;
   for (int j=m-1;j>0;j--)
   {
      sv=d;
      d=y2*d-dd+c[j];
      dd=sv;
   }
   ch=y*d-dd+0.5*c[0];
   return ch;
}

/** all these read subroutines are no longer used **/
/*
// read array sizes from the top of the file so I can allocate arrays
void readArraySizes(
EamCheby *chPot,
int* numCh,
int* numRef)
{
int i, d_int;

FILE *ChFile;
ChFile = fopen("mishin_ch.txt", "r");

fscanf(ChFile, "%d %d\n", &d_int, &d_int);
//printf("%d\n", d_int);
int num_vars = d_int/2;
for (int i=0;i<num_vars;i++)
{
fscanf(ChFile, "%d\n", &numRef[i]);
fscanf(ChFile, "%d\n", &numCh[i]);
//printf("%d, %d\n", numRef[i], numCh[i]);
}
chPot->phi->n = numCh[0];
chPot->rho->n = numCh[1];
chPot->f->n = numCh[2];
}

// read in the reference arrays and their derivatives
void readMishinRef(
int* numRef,
real_t* phiRange,
real_t* rhoRange,
real_t* FRange,
real_t* phiRef, 
real_t* rhoRef, 
real_t* FRef, 
real_t* xphiRef, 
real_t* xrhoRef, 
real_t* xFRef) 
{
FILE *RefFile;

int i, i_err, d_int;
int count =0;

RefFile = fopen("data/mishin_cu", "r");
fscanf(RefFile, "%d %d\n", &d_int, &d_int);
//printf("%d\n", d_int);

fscanf(RefFile, "%d %d %d\n", &d_int, &d_int, &d_int);
//printf("%d\n", d_int);
fscanf(RefFile, "%d\n", &d_int);
//printf("%d\n", numRef[count]);
fscanf(RefFile, FMT1" "FMT1"\n", &phiRange[0], &phiRange[1]);
phiRange[2] = (phiRange[1]-phiRange[0])/(numRef[count]);
//printf("%e, %e, %e\n", phiRange[0], phiRange[1], phiRange[2]);

for (int i=0;i<numRef[count]+1;i++)
{
i_err = fscanf(RefFile, EMT1" "EMT1"\n", &xphiRef[i], &phiRef[i]);
#if (DIAG_CHEBY > 0)
if (i_err > 0)
//printf("%d, %e, %e\n", i, xphiRef[i], phiRef[i]);
#endif
}
phiRange[2] = (xphiRef[1]-xphiRef[0]);
count ++;

fscanf(RefFile, "%d %d %d\n", &d_int, &d_int, &d_int);
//printf("%d\n", d_int);
fscanf(RefFile, "%d\n", &numRef[count]);
//printf("%d\n", numRef[count]);
fscanf(RefFile, FMT1" "FMT1"\n", &rhoRange[0], &rhoRange[1]);
rhoRange[2] = (rhoRange[1]-rhoRange[0])/(numRef[count]);
//printf("%e, %e, %e\n", rhoRange[0], rhoRange[1], rhoRange[2]);

for (int i=0;i<numRef[count]+1;i++)
{
   i_err = fscanf(RefFile, EMT1" "EMT1"\n", &xrhoRef[i], &rhoRef[i]);
#if (DIAG_CHEBY > 0)
   if (i_err > 0)
      //printf("%d, %e, %e\n", i, xrhoRef[i], rhoRef[i]);
#endif
}
rhoRange[2] = (xrhoRef[1]-xrhoRef[0]);
count ++;

fscanf(RefFile, "%d %d %d\n", &d_int, &d_int, &d_int);
//printf("%d\n", d_int);
fscanf(RefFile, "%d\n", &numRef[count]);
//printf("%d\n", numRef[count]);
fscanf(RefFile, FMT1" "FMT1"\n", &FRange[0], &FRange[1]);
FRange[2] = (FRange[1]-FRange[0])/(numRef[count]);
//printf("%e, %e, %e\n", FRange[0], FRange[1], FRange[2]);

for (int i=0;i<numRef[count]+1;i++)
{
   i_err = fscanf(RefFile, EMT1" "EMT1"\n", &xFRef[i], &FRef[i]);
#if (DIAG_CHEBY > 0)
   if (i_err > 0)
      //printf("%d, %e, %e\n", i, xFRef[i], FRef[i]);
#endif
}
FRange[2] = (xFRef[1]-xFRef[0]);

fclose(RefFile);
}

// read in the Chebychev coefficients
void readMishinCh(EamCheby *chPot)
{
   //chPot = (EamCheby*)malloc(sizeof(EamCheby));

   FILE *ChFile;

   int d_int;
   int count=0;
   int i;
   int i_err;
   char d_str[128];
   real_t d_real_1, d_real_2;
   int d_int_1;
   int numCh[3];
   real_t x0, xn;

   //printf("Reading Chebychev coefficients.. \n");

   ChFile = fopen("data/mishin_ch.txt", "r");

   fscanf(ChFile, "%d %d\n", &d_int, &d_int);
   //printf("%d\n", d_int);
   int num_vars = d_int/2;
   for (int i=0;i<num_vars;i++)
   {
      fscanf(ChFile, "%d\n", &d_int_1);
      fscanf(ChFile, "%d\n", &numCh[i]);
      //printf("%d, %d\n", d_int_1, numCh[i]);
   }

   // allocate all the array space here
   chPot->phi = (struct PotentialArray*)malloc(sizeof(struct PotentialArray));
   chPot->rho = (struct PotentialArray*)malloc(sizeof(struct PotentialArray));
   chPot->f   = (struct PotentialArray*)malloc(sizeof(struct PotentialArray));

   chPot->dphi = (struct PotentialArray*)malloc(sizeof(struct PotentialArray));
   chPot->drho = (struct PotentialArray*)malloc(sizeof(struct PotentialArray));
   chPot->df   = (struct PotentialArray*)malloc(sizeof(struct PotentialArray));

   chPot->phi->values = (real_t*)malloc((numCh[0]+3)*sizeof(real_t));
   chPot->rho->values = (real_t*)malloc((numCh[1]+3)*sizeof(real_t));
   chPot->f->values   = (real_t*)malloc((numCh[2]+3)*sizeof(real_t));

   chPot->dphi->values = (real_t*)malloc((numCh[0]+3)*sizeof(real_t));
   chPot->drho->values = (real_t*)malloc((numCh[1]+3)*sizeof(real_t));
   chPot->df->values   = (real_t*)malloc((numCh[2]+3)*sizeof(real_t));

   chPot->phi->n = numCh[0];
   chPot->rho->n = numCh[1];
   chPot->f->n   = numCh[2];

   chPot->dphi->n = numCh[0];
   chPot->drho->n = numCh[1];
   chPot->df->n   = numCh[2];

   //printf("Interpolation ranges..\n");
   fscanf(ChFile, FMT1" "FMT1"\n", &chPot->phi->x0, &chPot->phi->xn);
   //printf("%lf, %lf\n", chPot->phi->x0, chPot->phi->xn);

   //fscanf(ChFile, FMT1" "FMT1"\n", &x0, &xn);
   //chPot->phi = allocPotentialArray(numCh[0], x0, xn, 1.0);

   ////printf("Reading phi coefficients...");
   // read in Chebychev coefficients
   for (int i=0;i<numCh[count];i++)
   {
      i_err = fscanf(ChFile, EMT1"\n", &chPot->phi->values[i]);
#if (DIAG_CHEBY > 0)
      if (i_err > 0)
	 //printf("%d, "EMT1"\n", i, chPot->phi->values[i]);
#endif
   }
   ////printf("done\n");

   count ++;

   fscanf(ChFile, FMT1" "FMT1"\n", &chPot->rho->x0, &chPot->rho->xn);
   //printf("%lf, %lf\n", chPot->rho->x0, chPot->rho->xn);


   // read in Chebychev coefficients
   for (int i=0;i<numCh[count];i++)
   {
      i_err = fscanf(ChFile, EMT1"\n", &chPot->rho->values[i]);
#if (DIAG_CHEBY > 0)
      if (i_err > 0)
	 //printf("%d, "EMT1"\n", i, chPot->rho->values[i]);
#endif
   }

   count ++;

   fscanf(ChFile, FMT1" "FMT1"\n", &chPot->f->x0, &chPot->f->xn);
   //printf("%lf, %lf\n", chPot->f->x0, chPot->f->xn);

   // read in Chebychev coefficients
   for (int i=0;i<numCh[count];i++)
   {
      i_err = fscanf(ChFile, EMT1"\n", &chPot->f->values[i]);
#if (DIAG_CHEBY > 0)
      if (i_err > 0)
	 //printf("%d, "EMT1"\n", i, chPot->f->values[i]);
#endif
   }

   //printf("Number of Chebychev coefficients:\n");
   for (int i=0;i<3;i++) //printf("%d, %d\n", i, numCh[i]);

   fclose(ChFile);

#if (DIAG_CHEBY > 0)
   printArray(chPot->phi->values, chPot->phi->n, "phi");
   printArray(chPot->rho->values, chPot->rho->n, "rho");
   printArray(chPot->f->values,   chPot->f->n,   "F");
#endif

}

// call chDer for each set of coefficients
void genDerivs(EamCheby *chPot)
{
   chDer(chPot->phi->x0, chPot->phi->xn, chPot->phi->values, chPot->dphi->values, chPot->phi->n);
   chDer(chPot->rho->x0, chPot->rho->xn, chPot->rho->values, chPot->drho->values, chPot->rho->n);
   chDer(chPot->f->x0, chPot->f->xn, chPot->f->values, chPot->df->values, chPot->f->n);
#if (DIAG_CHEBY > 0)
   printArray(chPot->dphi->values, chPot->phi->n, "dphi");
   printArray(chPot->drho->values, chPot->rho->n, "drho");
   printArray(chPot->df->values, chPot->f->n, "dF");
#endif

}

// compute the function value from the Chebychev coefficients, and compare to the reference values.
void computeCompare(
      real_t* range, 
      real_t* ch, 
      real_t* dch,
      real_t* ref, 
      real_t* x, 
      real_t* ref_val_h,
      real_t* ref_dval_h,
      int n_ref, 
      int n_ch)
{
   int i;
   real_t ch_val, ch_deriv;
   real_t ref_val, ref_deriv;

   real_t l1, l2;
   real_t l1d, l2d;

   l1 = 0.0;
   l1d = 0.0;
   l2 = 0.0;
   l2d = 0.0;

   real_t diff, diffd;

   for (int i=1;i<n_ref;i++)
   {
      eamInterpolateDerivlocal(x[i], ref, n_ref, range, &ref_val, &ref_deriv);
      ref_val_h[i] = ref_val;
      ref_dval_h[i] = ref_deriv;
      ch_val = chebev(range[0], range[1], ch, n_ch, x[i]);
      ch_deriv = chebev(range[0], range[1], dch, n_ch, x[i]);
      diff = fabs(ref_val - ch_val);
      if (diff > l1) l1 = diff;
      l2 += diff*diff;
      diffd = fabs(ref_deriv - ch_deriv);
      if (diffd > l1d) l1d = diffd;
      l2d += diffd*diffd;
      ////printf("%d, %17.9e, %17.9e, %17.9e, %17.9e, %17.9e\n", i, x[i], ch_val, ref_val, ch_deriv, ref_deriv);

   }
   //printf("Values l1, l2 errors: "EMT1", "EMT1"\n", l1, l2/(n_ref - 1));
   //printf("Derivs l1, l2 errors: "EMT1", "EMT1"\n", l1d, l2d/(n_ref - 1));
}

/*
   int main(int argc, char **argv)
   {
   real_t phiRange[3];
   real_t rhoRange[3];
   real_t FRange[3];

   int numCh[3];
   int numRef[3];

   oclInit(0);

// get array sizes from file
readArraySizes(numCh, numRef);

size_t phiRefSize = numRef[0]*sizeof(real_t);
size_t rhoRefSize = numRef[1]*sizeof(real_t);
size_t FRefSize = numRef[2]*sizeof(real_t);

//printf("Allocating host side arrays...");

// allocate arrays for Chebychev coefficients
real_t* phiCh=malloc(numCh[0]*sizeof(real_t));
real_t* dphiCh=malloc(numCh[0]*sizeof(real_t));

real_t* rhoCh=malloc(numCh[1]*sizeof(real_t));
real_t* drhoCh=malloc(numCh[1]*sizeof(real_t));

real_t* FCh=malloc(numCh[2]*sizeof(real_t));
real_t* dFCh=malloc(numCh[2]*sizeof(real_t));

// allocate arrays for EAM arrays and ordinates
real_t* phiRef=malloc(phiRefSize);
real_t* xphiRef=malloc(phiRefSize);

real_t* rhoRef=malloc(numRef[1]*sizeof(real_t));
real_t* xrhoRef=malloc(numRef[1]*sizeof(real_t));

real_t* FRef=malloc(numRef[2]*sizeof(real_t));
real_t* xFRef=malloc(numRef[2]*sizeof(real_t));

real_t* phiChVal_h=malloc(phiRefSize);
real_t* phiRefVal_h=malloc(phiRefSize);

real_t* rhoChVal_h=malloc(numRef[1]*sizeof(real_t));
real_t* rhoRefVal_h=malloc(numRef[1]*sizeof(real_t));

real_t* FChVal_h=malloc(numRef[2]*sizeof(real_t));
real_t* FRefVal_h=malloc(numRef[2]*sizeof(real_t));

real_t* dphiChVal_h=malloc(phiRefSize);
real_t* dphiRefVal_h=malloc(phiRefSize);

real_t* drhoChVal_h=malloc(numRef[1]*sizeof(real_t));
real_t* drhoRefVal_h=malloc(numRef[1]*sizeof(real_t));

real_t* dFChVal_h=malloc(numRef[2]*sizeof(real_t));
real_t* dFRefVal_h=malloc(numRef[2]*sizeof(real_t));

//printf("done\n");

// allocate the corresponding OpenCL arrays
cl_mem phiCh_d;
cl_mem rhoCh_d;
cl_mem FCh_d;

cl_mem dphiCh_d;
cl_mem drhoCh_d;
cl_mem dFCh_d;

cl_mem xphiRef_d;
cl_mem xrhoRef_d;
cl_mem xFRef_d;

cl_mem phiRange_d;
cl_mem rhoRange_d;
cl_mem FRange_d;

cl_mem phiRef_d;
cl_mem rhoRef_d;
cl_mem FRef_d;

cl_mem phiChVal_d;
cl_mem rhoChVal_d;
cl_mem FChVal_d;

cl_mem dphiChVal_d;
cl_mem drhoChVal_d;
cl_mem dFChVal_d;

cl_mem phiRefVal_d;
cl_mem rhoRefVal_d;
cl_mem FRefVal_d;

cl_mem dphiRefVal_d;
cl_mem drhoRefVal_d;
cl_mem dFRefVal_d;

//printf("Allocating device side buffers....");

// Chebychev coefficients
oclCreateReadBuffer(&phiCh_d, numCh[0]*sizeof(real_t));
oclCreateReadBuffer(&rhoCh_d, numCh[1]*sizeof(real_t));
oclCreateReadBuffer(&FCh_d,   numCh[2]*sizeof(real_t));

//Chebychev derivative coefficients
oclCreateReadBuffer(&dphiCh_d, numCh[0]*sizeof(real_t));
oclCreateReadBuffer(&drhoCh_d, numCh[1]*sizeof(real_t));
oclCreateReadBuffer(&dFCh_d,   numCh[2]*sizeof(real_t));

// range tables
oclCreateReadBuffer(&phiRange_d, 3*sizeof(real_t));
oclCreateReadBuffer(&rhoRange_d, 3*sizeof(real_t));
oclCreateReadBuffer(&FRange_d,   3*sizeof(real_t));

oclCreateReadBuffer(&phiRef_d, phiRefSize);
oclCreateReadBuffer(&rhoRef_d, numRef[1]*sizeof(real_t));
oclCreateReadBuffer(&FRef_d,   numRef[2]*sizeof(real_t));

// ordinates to evaluate values
oclCreateReadBuffer(&xphiRef_d, phiRefSize);
oclCreateReadBuffer(&xrhoRef_d, numRef[1]*sizeof(real_t));
oclCreateReadBuffer(&xFRef_d,   numRef[2]*sizeof(real_t));

// Chebychev output values
oclCreateWriteBuffer(&phiChVal_d, phiRefSize);
oclCreateWriteBuffer(&rhoChVal_d, numRef[1]*sizeof(real_t));
oclCreateWriteBuffer(&FChVal_d,   numRef[2]*sizeof(real_t));

oclCreateWriteBuffer(&dphiChVal_d, phiRefSize);
oclCreateWriteBuffer(&drhoChVal_d, numRef[1]*sizeof(real_t));
oclCreateWriteBuffer(&dFChVal_d,   numRef[2]*sizeof(real_t));

oclCreateWriteBuffer(&phiRefVal_d, phiRefSize);
oclCreateWriteBuffer(&rhoRefVal_d, numRef[1]*sizeof(real_t));
oclCreateWriteBuffer(&FRefVal_d,   numRef[2]*sizeof(real_t));

oclCreateWriteBuffer(&dphiRefVal_d, phiRefSize);
oclCreateWriteBuffer(&drhoRefVal_d, numRef[1]*sizeof(real_t));
oclCreateWriteBuffer(&dFRefVal_d,   numRef[2]*sizeof(real_t));

//printf("done\n");

// read in the coefficients from the file
readMishinCh(numCh, numRef, phiRange, rhoRange, FRange, phiCh, rhoCh, FCh);

// generate the derivative coefficient values
genDerivs(phiRange, rhoRange, FRange, phiCh, rhoCh, FCh, dphiCh, drhoCh, dFCh, numCh);

// read the reference coefficients from the file
readMishinRef(numRef, phiRange, rhoRange, FRange, phiRef, rhoRef, FRef, xphiRef, xrhoRef, xFRef);

// do the test with the C code
computeCompare(phiRange, phiCh, dphiCh, phiRef, xphiRef, phiRefVal_h, dphiRefVal_h, numRef[0], numCh[0]);
computeCompare(rhoRange, rhoCh, drhoCh, rhoRef, xrhoRef, rhoRefVal_h, drhoRefVal_h, numRef[1], numCh[1]);
computeCompare(FRange,   FCh,   dFCh,   FRef,   xFRef,   FRefVal_h,   dFRefVal_h, numRef[2], numCh[2]);

// copy Chebychev coefficients to device
oclCopyToDevice(phiRange, phiRange_d, 3*sizeof(real_t));
oclCopyToDevice(rhoRange, rhoRange_d, 3*sizeof(real_t));
oclCopyToDevice(FRange,   FRange_d,   3*sizeof(real_t));

oclCopyToDevice(phiCh, phiCh_d, numCh[0]*sizeof(real_t));
oclCopyToDevice(rhoCh, rhoCh_d, numCh[1]*sizeof(real_t));
oclCopyToDevice(FCh,   FCh_d,   numCh[2]*sizeof(real_t));

oclCopyToDevice(dphiCh, dphiCh_d, numCh[0]*sizeof(real_t));
oclCopyToDevice(drhoCh, drhoCh_d, numCh[1]*sizeof(real_t));
oclCopyToDevice(dFCh,   dFCh_d,   numCh[2]*sizeof(real_t));

oclCopyToDevice(xphiRef, xphiRef_d, phiRefSize);
oclCopyToDevice(xrhoRef, xrhoRef_d, numRef[1]*sizeof(real_t));
oclCopyToDevice(xFRef,   xFRef_d,   numRef[2]*sizeof(real_t));


oclBuildProgramFromFile("ch_kernels.c");

cl_kernel computeCompare;

int err;
cl_event e_1;
real_t t_elapsed, t_enqueued;

// build kernel
computeCompare = clCreateKernel(program, "computeCompare", &err);

// kernel args for phi test
clSetKernelArg(computeCompare, 0, sizeof(cl_mem), &phiRange_d);
clSetKernelArg(computeCompare, 1, sizeof(cl_mem), &phiCh_d);
clSetKernelArg(computeCompare, 2, sizeof(cl_mem), &dphiCh_d);
clSetKernelArg(computeCompare, 3, sizeof(cl_mem), &phiRef_d);
clSetKernelArg(computeCompare, 4, sizeof(cl_mem), &xphiRef_d);
clSetKernelArg(computeCompare, 5, sizeof(cl_mem), &phiChVal_d);
clSetKernelArg(computeCompare, 6, sizeof(cl_mem), &dphiChVal_d);
clSetKernelArg(computeCompare, 7, sizeof(cl_mem), &phiRefVal_d);
clSetKernelArg(computeCompare, 8, sizeof(cl_mem), &dphiRefVal_d);
clSetKernelArg(computeCompare, 9, sizeof(int), &numRef[0]);
clSetKernelArg(computeCompare, 10, sizeof(int), &numCh[0]);

size_t nGlobal = numRef[0];
//printf("Running computeCompare for phi..");
fflush(stdout);
clEnqueueNDRangeKernel(commandq, computeCompare, 1, 0, &nGlobal, NULL, 0, NULL, &e_1);
clWaitForEvents(1, &e_1);
//printf("done\n");
fflush(stdout);

GetElapsedTime(e_1, &t_elapsed, &t_enqueued);
//printf("Elapsed time for phi = "EMT1"\n", t_elapsed);

//printf("Copying results to host...");
fflush(stdout);
oclCopyToHost(phiChVal_d, phiChVal_h, phiRefSize);
oclCopyToHost(dphiChVal_d, dphiChVal_h, phiRefSize);
//printf("done\n");
fflush(stdout);

int i;
#if (DIAG_CHEBY > 0 )
for (int i=0;i<numRef[0];i++)
{
   //printf("%d, "EMT1", "EMT1", "EMT1", "EMT1", "EMT1"\n", i, xphiRef[i], 
	 phiChVal_h[i], phiRef[i],
	 dphiChVal_h[i], dphiRefVal_h[i]);
}
#endif

// kernel args for rho test
clSetKernelArg(computeCompare, 0, sizeof(cl_mem), &rhoRange_d);
clSetKernelArg(computeCompare, 1, sizeof(cl_mem), &rhoCh_d);
clSetKernelArg(computeCompare, 2, sizeof(cl_mem), &drhoCh_d);
clSetKernelArg(computeCompare, 3, sizeof(cl_mem), &rhoRef_d);
clSetKernelArg(computeCompare, 4, sizeof(cl_mem), &xrhoRef_d);
clSetKernelArg(computeCompare, 5, sizeof(cl_mem), &rhoChVal_d);
clSetKernelArg(computeCompare, 6, sizeof(cl_mem), &drhoChVal_d);
clSetKernelArg(computeCompare, 7, sizeof(cl_mem), &rhoRefVal_d);
clSetKernelArg(computeCompare, 8, sizeof(cl_mem), &drhoRefVal_d);
clSetKernelArg(computeCompare, 9, sizeof(int), &numRef[1]);
clSetKernelArg(computeCompare, 10, sizeof(int), &numCh[1]);

nGlobal = numRef[1];
//printf("Running computeCompare for rho..");
fflush(stdout);
clEnqueueNDRangeKernel(commandq, computeCompare, 1, 0, &nGlobal, NULL, 0, NULL, &e_1);
clWaitForEvents(1, &e_1);
//printf("done\n");
fflush(stdout);

GetElapsedTime(e_1, &t_elapsed, &t_enqueued);
//printf("Elapsed time for rho = "EMT1"\n", t_elapsed);

//printf("Copying results to host...");
fflush(stdout);
oclCopyToHost(rhoChVal_d, rhoChVal_h, rhoRefSize);
oclCopyToHost(drhoChVal_d, drhoChVal_h, rhoRefSize);
//printf("done\n");
fflush(stdout);

#if (DIAG_CHEBY > 0 )
for (int i=0;i<numRef[0];i++)
{
   //printf("%d, "EMT1", "EMT1", "EMT1", "EMT1", "EMT1"\n", i, xrhoRef[i], 
	 rhoChVal_h[i], rhoRef[i],
	 drhoChVal_h[i], drhoRefVal_h[i]);
}
#endif

// kernel args for F test
clSetKernelArg(computeCompare, 0, sizeof(cl_mem), &FRange_d);
clSetKernelArg(computeCompare, 1, sizeof(cl_mem), &FCh_d);
clSetKernelArg(computeCompare, 2, sizeof(cl_mem), &dFCh_d);
clSetKernelArg(computeCompare, 3, sizeof(cl_mem), &FRef_d);
clSetKernelArg(computeCompare, 4, sizeof(cl_mem), &xFRef_d);
clSetKernelArg(computeCompare, 5, sizeof(cl_mem), &FChVal_d);
clSetKernelArg(computeCompare, 6, sizeof(cl_mem), &dFChVal_d);
clSetKernelArg(computeCompare, 7, sizeof(cl_mem), &FRefVal_d);
clSetKernelArg(computeCompare, 8, sizeof(cl_mem), &dFRefVal_d);
clSetKernelArg(computeCompare, 9, sizeof(int), &numRef[2]);
clSetKernelArg(computeCompare, 10, sizeof(int), &numCh[2]);

nGlobal = numRef[2];
//printf("Running computeCompare for F..");
fflush(stdout);
clEnqueueNDRangeKernel(commandq, computeCompare, 1, 0, &nGlobal, NULL, 0, NULL, &e_1);
clWaitForEvents(1, &e_1);
//printf("done\n");
fflush(stdout);

GetElapsedTime(e_1, &t_elapsed, &t_enqueued);
//printf("Elapsed time for F = "EMT1"\n", t_elapsed);

//printf("Copying results to host...");
fflush(stdout);
oclCopyToHost(FChVal_d, FChVal_h, FRefSize);
oclCopyToHost(dFChVal_d, dFChVal_h, FRefSize);
//printf("done\n");
fflush(stdout);

#if (DIAG_CHEBY > 0 )
for (int i=0;i<numRef[0];i++)
{
   //printf("%d, "EMT1", "EMT1", "EMT1", "EMT1", "EMT1"\n", i, xFRef[i], 
	 FChVal_h[i], FRef[i],
	 dFChVal_h[i], dFRefVal_h[i]);
}
#endif
}
*/
