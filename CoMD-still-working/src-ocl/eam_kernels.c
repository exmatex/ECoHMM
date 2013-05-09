/** Kernels for computing the EAM potential
  Since we can only block on kernel completion, the 3 sweeps done in the original code 
  need to be implemented as 3 separate kernels. 
  Note also that the potential arrays are large enough to require accessing them from 
  global memory.
  Since OpenCL doesn't pick up #include properly, we need to manually switch real_t from 
  float to double in each kernel file individually.

  Note: More careful analysis shows we can consolidate kernels 1 and 2 into a single pass;
  there is a flag PASS_2 in this file which should be set to match the flag PASS_2 in 
  helpers.c. Switching this flag allows you to test that a) the results match and 2) evaluate
  the overhead of adding the extra kernel which only loops over all particles.
 **/

//Initial implementation of the MD code
#define N_MAX_ATOMS 64
#define N_MAX_NEIGHBORS 27
#define PERIODIC 1

#define KERN_DIAG 0
#define USE_SPLINE 0
#define PASS_2 0 // this should match the setting in helpers.c

#if defined(cl_khr_fp64)  // Khronos extension available?
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)  // AMD extension available?
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

/* CL_REAL_T is set to single or double depending on compile time flags */
typedef CL_REAL_T real_t; 
typedef CL_REAL4_T cl_real4;

#if (USE_CHEBY)
// given a list of Chebyshev coefficients c, compute the value at x
// x must be in the range [a, b]
// call signature might change to put it in line with the regular EAM call
real_t chebev(
        real_t a, 
        real_t b, 
        __global const real_t *c,
        int m, 
        real_t x) 
{
   real_t d, dd, sv, y, y2;
   real_t ch;
#if(KERN_DIAG > 0) 
   if ((x-a)*(x-b) > 0.0)
   {
      printf("x not in range in chebev, %f\n", x);
   }
#endif
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
/* Chebyshev interpolation routine modified to match the table lookup signature.
 * Returns the value and derivative, ch and dch, at the position x
 */
void eamCheby(
	real_t x,
        __global const real_t* c,
        const int m, 
	real_t *ch,
	real_t *dch) 
{
    // coefficients array c has Chebyshev coeffs, derivative coeffs and range limits
   real_t a = c[2*m + 0];
   real_t b = c[2*m + 1];

   real_t d, dd, sv, y, y2;
#if(KERN_DIAG > 0) 
   if ((x-a)*(x-b) > 0.0)
   {
      printf("x not in range in chebev, %f\n", x);
   }
#endif
   // compute the value
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
   *ch=y*d-dd+0.5*c[0];

    // compute the derivative
    // identical to above but c values offset by m
   d=0.0;
   dd=0.0;
   y=(2.0*x-a-b)/(b-a);
   y2=2.0*y;
   for (int j=m-1;j>0;j--)
   {
      sv=d;
      d=y2*d-dd+c[j+m];
      dd=sv;
   }
   *dch=y*d-dd+0.5*c[0+m];
}

#endif

#if (USE_SPLINE)
void PhiSpline(real_t r, real_t *f, real_t *df)
{
    // values for copper
    /* 
    // original values from 87 paper
    real_t a_k[6] = { 29.059214 , -140.05681 , 130.07331 , -17.48135 , 31.82546 , 71.58749};
    real_t r_k[6] = { 1.2247449 , 1.1547054 , 1.1180065 , 1.0000000 , 0.8660254 , 0.7071068};
    */
    // new values for smoother potential
    real_t a_k[6] = {61.73525861, -108.18467800, 57.00053948,-12.88796578, 39.16381901, 0.0};
    real_t r_k[6] = {1.225, 1.202, 1.154, 1.050, 0.866, 0.707};
    real_t az = 1.0*3.615;
    real_t az3 = az*az*az;
    real_t az2 = az*az;
    // set output values to zero
    *f=0.0;
    *df=0.0;
    //r = 1.0; //3.615;
    r = r/3.615;
    //printf("r = %e\n", r);
    // sum over all coefficients
    for (int k=0;k<6;k++)
    {
	r_k[k] = r_k[k]*3.615;
	if (r < r_k[k])
	{
	    *f += (r_k[k]-r)*(r_k[k]-r)*(r_k[k]-r)*a_k[k]/az3;
	    *df -= 3.0*(r_k[k]-r)*(r_k[k]-r)*a_k[k]/az2;
	}
    }
}

void RhoSpline(real_t r, real_t *f, real_t *df)
{
    // values for copper
    /*
    // original values from 87 paper
    real_t R_k[2] = { 1.2247449 , 1.0000000 };
    real_t A_k[2] = { 9.806694 , 16.774638 };
    */
    // new values for smoother potential
    real_t R_k[2] = { 1.225, 0.990 };
    real_t A_k[2] = { 10.03718305, 17.06363299 };
    real_t az = 1.0*3.615;
    real_t az3 = az*az*az;
    real_t az2 = az*az;
    // set output values to zero
    *f=0.0;
    *df=0.0;
    //r = 1.0; //3.615;
    r = r/3.615;
    // sum over all coefficients
    for (int k=0;k<2;k++)
    {
	R_k[k] = R_k[k]*3.615;
	if (r < R_k[k])
	{
	    *f += (R_k[k]-r)*(R_k[k]-r)*(R_k[k]-r)*A_k[k]/az3;
	    *df -= 3.0*(R_k[k]-r)*(R_k[k]-r)*A_k[k]/az2;
	}
    }
}

void FSpline(real_t rho, real_t *f, real_t *df)
{
    *f = -1.0*sqrt(rho);
    *df = -0.5/(*f);
}

#endif

void eamInterpolateDeriv(real_t r,
	__global const real_t* values,
	const int nValues,
	real_t *value1, 
	real_t *f1)
{
    int i1;
    int i;
    real_t gi, gi1;

    // extract values from potential 'struct'
    real_t x0 = values[nValues+3];
    real_t xn = values[nValues+4];
    real_t invDx = values[nValues+5];

    // identical to Sriram's loop in eam.c
    if ( r<x0) r = x0;
    else if (r>xn) r = xn;

    r = (r-x0)*(invDx) ;
    i1 = (int)floor(r);

    /* reset r to fractional distance */
    r = r - floor(r);

    // all indices shifted up by one compared to the original code
    gi  = values[i1+2] - values[i1];
    gi1 = values[i1+3] - values[i1+1];

    // Note the shift removes [i1-1] as a possibility
    // values[i1-1] is guaranteed(?) inbounds because 
    // a->x0 = x0 + (xn-x0)/(double)n; 
    // appears in allocPotentialArray
    *value1 = values[i1+1] + 0.5*r*(
	    r*(values[i1+2]+ values[i1] -2.0*values[i1+1]) +
	    gi
	    );
    if(i1<=1) 
	*f1 = 0.0;
    else 
	*f1 = 0.5*(gi + r*(gi1-gi))*invDx;

    return;

}

// Simple version without local blocking to check for correctness
__kernel void EAM_Force_1(
	__global real_t* xPos,
	__global real_t* yPos,
	__global real_t* zPos,

	__global real_t* fx,
	__global real_t* fy,
	__global real_t* fz,

	__global real_t* energy,
	__global real_t* rho,
	__global real_t* rhobar,

	__global real_t* dcx,
	__global real_t* dcy,
	__global real_t* dcz,

	__global real_t* bounds,
	__global const int* neighborList,
	__global const int* nNeighbors,
	__global const int* nAtoms,

	__global const int* nValues,

	__global const real_t* phi_pot, // the potentials are passed in as real arrays: x0, xn, invDx, values[?]
	__global const real_t* rho_pot,
	__global const real_t* F_pot,

	const real_t cutoff) 
{
#if USE_SPLINE
    // values for copper
    real_t a_k[6] = { 29.059214 , -140.05681 , 130.07331 , -17.48135 , 31.82546 , 71.58749};
    real_t r_k[6] = { 1.2247449 , 1.1547054 , 1.1180065 , 1.0000000 , 0.8660254 , 0.7071068};

    real_t R_k[2] = { 1.2247449 , 1.0000000 };
    real_t A_k[2] = { 9.806694 , 16.774638 };
#endif

    int iAtom = get_global_id(0);
    int iBox = get_global_id(1);

    real_t dx, dy, dz;
    real_t r, r2, r6;
    real_t fr, e_i;
    real_t rho_i;

    real_t dxbox, dybox, dzbox;

    // accumulate local force value
    real_t fx_i, fy_i, fz_i;

    real_t rCut = cutoff;
    real_t rCut2 = rCut*rCut;
    real_t rhoTmp;
    real_t phiTmp;
    real_t dPhi, dRho;
    real_t fi, fiprime;

    int i;
    int j_local;

    int i_offset;
    int iParticle;

    // zero out forces on particles
    i_offset = iBox*N_MAX_ATOMS;
    iParticle = i_offset + iAtom;

    fx_i = 0.0;
    fy_i = 0.0;
    fz_i = 0.0;

    rho_i = 0.0;

    e_i = 0.0;

    if (iAtom < nAtoms[iBox])
    {// each thread executes on a single atom in the box

#if(KERN_DIAG > 0) 
	printf("i = %d, %f, %f, %f\n", iParticle, xPos[iParticle], yPos[iParticle], zPos[iParticle]);

	printf("iBox = %d, nNeighbors = %d\n", iBox, nNeighbors[iBox]);
#endif

	for (int j = 0; j<nNeighbors[iBox]; j++)
	{// loop over neighbor cells
	    int jBox = neighborList[iBox*N_MAX_NEIGHBORS + j];
	    int jOffset = jBox*N_MAX_ATOMS;

	    // compute box center offsets
	    dxbox = dcx[iBox] - dcx[jBox];
	    dybox = dcy[iBox] - dcy[jBox];
	    dzbox = dcz[iBox] - dcz[jBox];

	    // correct for periodic 
	    if(PERIODIC)
	    {
		if (dxbox<-0.5*bounds[0]) dxbox += bounds[0];
		else if (dxbox > 0.5*bounds[0] ) dxbox -= bounds[0];
		if (dybox<-0.5*bounds[1]) dybox += bounds[1];
		else if (dybox > 0.5*bounds[1] ) dybox -= bounds[1];
		if (dzbox<-0.5*bounds[2]) dzbox += bounds[2];
		else if (dzbox > 0.5*bounds[2] ) dzbox -= bounds[2];
	    }

	    //printf("dxbox, dybox, dzbox = %f, %f, %f\n", dxbox, dybox, dzbox);

	    for (int jAtom = 0; jAtom<nAtoms[jBox]; jAtom++)
	    {// loop over all groups in neighbor cell 

		int jParticle = jOffset + jAtom; // global offset of particle

		dx = dxbox + xPos[iParticle] - xPos[jParticle];
		dy = dybox + yPos[iParticle] - yPos[jParticle];
		dz = dzbox + zPos[iParticle] - zPos[jParticle];

#if(KERN_DIAG > 0) 
		printf("dx, dy, dz = %f, %f, %f\n", dx, dy, dz);
		printf("i = %d, j = %d, %f, %f, %f\n", iParticle, jParticle, xPos[jParticle], yPos[jParticle], zPos[jParticle]);
#endif

		r2 = dx*dx + dy*dy + dz*dz;

		if ( r2 <= rCut2 && r2 > 0.0)
		{// no divide by zero

#if(KERN_DIAG > 0) 
		    printf("%d, %d, %f\n", iParticle, jParticle, r2);
		    printf("r2, rCut2 = %f, %f\n", r2, rCut2);
#endif

		    r = sqrt(r2);
		    //r = 1.0;

/*
#if(USE_SPLINE)
		    //r = r*r_k[0]/cutoff;
		    //r = r/3.615;
		    PhiSpline(r, &phiTmp, &dPhi);
		    RhoSpline(r, &rhoTmp, &dRho);
		    */
#if (USE_CHEBY)
		    eamCheby(r, phi_pot, nValues[0], &phiTmp, &dPhi);
		    eamCheby(r, rho_pot, nValues[1], &rhoTmp, &dRho);
#else
		    eamInterpolateDeriv(r, phi_pot, nValues[0], &phiTmp, &dPhi);
		    eamInterpolateDeriv(r, rho_pot, nValues[1], &rhoTmp, &dRho);
#endif

#if(KERN_DIAG > 0) 
		    printf("%d, %d, %f\n", iParticle, jParticle, r2);
		    printf("iParticle = %d, jParticle = %d, phiTmp = %f, dPhi = %f\n", iParticle, jParticle, phiTmp, dPhi);
		    printf("iParticle = %d, jParticle = %d, rhoTmp = %f, dRho = %f\n", iParticle, jParticle, rhoTmp, dRho);
#endif

#if(USE_SPLINE)
		    fx_i += (dRho+dPhi)*dx/r;
		    fy_i += (dRho+dPhi)*dy/r;
		    fz_i += (dRho+dPhi)*dz/r;
#else
		    fx_i += dPhi*dx/r;
		    fy_i += dPhi*dy/r;
		    fz_i += dPhi*dz/r;
#endif

		    e_i += phiTmp;

		    rho_i += rhoTmp;

		} else {
		}


	    } // loop over all atoms
	} // loop over neighbor cells

	fx[iParticle] = fx_i;
	fy[iParticle] = fy_i;
	fz[iParticle] = fz_i;

	// since we loop over all particles, each particle contributes 1/2 the pair energy to the total
	//energy[iParticle] = e_i*0.5;

	rho[iParticle] = rho_i;

	 // we can actually include the Force_2 kernel here, to save some time!
	 // skip the next 4 lines if PASS_2 = 1 in helpers.c
#if(PASS_2 == 0)
#if (USE_CHEBY)
	eamCheby(rho_i,F_pot,nValues[2],&fi,&fiprime);
#else
	eamInterpolateDeriv(rho_i,F_pot,nValues[2],&fi,&fiprime);
#endif
	rhobar[iParticle] = fiprime; // update rhoprime 
	//update energy terms 
	energy[iParticle] = e_i*0.5 + fi;
#else
	energy[iParticle] = e_i*0.5;
#endif
    } else { // zero out the energy of the other particles for safety
	energy[iParticle] = 0.0;
    }
}

__kernel void EAM_Force_2(
	__global real_t* rhobar,
	__global real_t* energy,
	__global real_t* rho,
	__global const int* nAtoms,
	__global const real_t* F_pot,
	__global const int* nValues
	)
{

    int iAtom = get_global_id(0);
    int iBox = get_global_id(1);

    real_t fi, fiprime;

    int i_offset;
    int iParticle;

    i_offset = iBox*N_MAX_ATOMS;
    iParticle = i_offset + iAtom;

    /*
    // local copy of F potential
    real_t F_local[nValues+3];

    // load values into local potentials
    for (int i=0;i<nValues+3;i++)
    {
    F_local[i] = F_pot[i];
    }
     */

    if (iAtom < nAtoms[iBox])
    {// each thread executes on a single atom in the box


	//printf("iBox = %d, iAtom = %d\n", iBox, iAtom);

        /*
#if(USE_SPLINE)
	FSpline(rho[iParticle], &fi, &fiprime);
        */
#if (USE_CHEBY)
	eamCheby(rho[iParticle],F_pot,nValues[2],&fi,&fiprime);
#else
	eamInterpolateDeriv(rho[iParticle],F_pot,nValues[2],&fi,&fiprime);
#endif
	rhobar[iParticle] = fiprime; // update rhoprime 
	//update energy terms 
	energy[iParticle] += fi;

    }

}

__kernel void EAM_Force_3(
	__global real_t* xPos,
	__global real_t* yPos,
	__global real_t* zPos,

	__global real_t* fx,
	__global real_t* fy,
	__global real_t* fz,

	__global real_t* fi,

	__global real_t* dcx,
	__global real_t* dcy,
	__global real_t* dcz,

	__global const real_t* bounds,
	__global const int* neighborList,
	__global const int* nNeighbors,
	__global const int* nAtoms,
	__global const int* nValues,

	__global const real_t* rho_pot,
	const real_t cutoff) 
{

    int iAtom = get_global_id(0);
    int iBox = get_global_id(1);

    real_t dx, dy, dz;
    real_t r, r2;

    real_t dxbox, dybox, dzbox;

    // accumulate local force value
    real_t fx_i, fy_i, fz_i;

    real_t rCut = cutoff;
    real_t rCut2 = rCut*rCut;
    real_t rhoTmp, dRho;

    int i;
    int j_local;

    real_t rTmp, rhoijprime;

    // global offset of local thread
    int i_offset = iBox*N_MAX_ATOMS;
    int iParticle = i_offset + iAtom;

    if (iAtom < nAtoms[iBox])
    {// each thread executes on a single atom in the box

	// zero out forces on particles
	fx_i = fx[iParticle];
	fy_i = fy[iParticle];
	fz_i = fz[iParticle];

#if(KERN_DIAG > 0) 

	printf("i = %d, %f, %f, %f\n", iParticle, xPos[iParticle], yPos[iParticle], zPos[iParticle]);

	printf("iBox = %d, nNeighbors = %d\n", iBox, nNeighbors[iBox]);
#endif

	for (int j = 0; j<nNeighbors[iBox]; j++)
	{// loop over neighbor cells
	    int jBox = neighborList[iBox*N_MAX_NEIGHBORS + j];
	    int jOffset = jBox*N_MAX_ATOMS;

	    // compute box center offsets
	    dxbox = dcx[iBox] - dcx[jBox];
	    dybox = dcy[iBox] - dcy[jBox];
	    dzbox = dcz[iBox] - dcz[jBox];

	    // correct for periodic 
	    if(PERIODIC)
	    {
		if (dxbox<-0.5*bounds[0]) dxbox += bounds[0];
		else if (dxbox > 0.5*bounds[0] ) dxbox -= bounds[0];
		if (dybox<-0.5*bounds[1]) dybox += bounds[1];
		else if (dybox > 0.5*bounds[1] ) dybox -= bounds[1];
		if (dzbox<-0.5*bounds[2]) dzbox += bounds[2];
		else if (dzbox > 0.5*bounds[2] ) dzbox -= bounds[2];
	    }

	    // printf("dxbox, dybox, dzbox = %f, %f, %f\n", dxbox, dybox, dzbox);

	    for (int jAtom = 0; jAtom<nAtoms[jBox]; jAtom++)
	    {// loop over all groups in neighbor cell 

		int jParticle = jOffset + jAtom; // global offset of particle

		/*
		   dx = xPos[iParticle] - xPos[jParticle] + dxbox;
		   dy = yPos[iParticle] - yPos[jParticle] + dybox;
		   dz = zPos[iParticle] - zPos[jParticle] + dzbox;
		 */

		dx = dxbox + xPos[iParticle] - xPos[jParticle];
		dy = dybox + yPos[iParticle] - yPos[jParticle];
		dz = dzbox + zPos[iParticle] - zPos[jParticle];

#if(KERN_DIAG > 0) 
		printf("dx, dy, dz = %f, %f, %f\n", dx, dy, dz);
		printf("i = %d, j = %d, %f, %f, %f\n", iParticle, jParticle, xPos[jParticle], yPos[jParticle], zPos[jParticle]);
#endif

		r2 = dx*dx + dy*dy + dz*dz;

		if ( r2 < rCut2 && r2 > 0.0)
		{// no divide by zero

#if(KERN_DIAG > 0) 
		    printf("r2, rCut2 = %f, %f\n", r2, rCut2);
		    printf("%d, %d, %f\n", iParticle, jParticle, r2);
#endif

		    r = sqrt(r2);

/*
#if (USE_SPLINE)
		    RhoSpline(r, &rhoTmp, &dRho);
#else
		    */
#if (USE_CHEBY)
		    eamCheby(r, rho_pot, nValues[1], &rhoTmp, &dRho);
#else
		    eamInterpolateDeriv(r, rho_pot, nValues[1], &rhoTmp, &dRho);
#endif
		    rhoijprime = dRho;

#if(KERN_DIAG > 0) 
		    printf("%d, %d, %f\n", iParticle, jParticle, r2);
#endif

		    rTmp = (fi[iParticle]+fi[jParticle])*rhoijprime/r;

		    fx_i += (fi[iParticle]+fi[jParticle])*rhoijprime*dx/r;
		    fy_i += (fi[iParticle]+fi[jParticle])*rhoijprime*dy/r;
		    fz_i += (fi[iParticle]+fi[jParticle])*rhoijprime*dz/r;
		    /*
		       fx_i += rTmp*dx;
		       fy_i += rTmp*dy;
		       fz_i += rTmp*dz;
		     */

		} else {
		}


	    } // loop over all atoms in jBox
	} // loop over neighbor cells

#if (USE_SPLINE)
#else
	fx[iParticle] = fx_i;
	fy[iParticle] = fy_i;
	fz[iParticle] = fz_i;
#endif

    } // loop over all atoms in iBox
}

// AoS Versions
// Simple version without local blocking to check for correctness
__kernel void EAM_Force_1_AoS(
	__global cl_real4* pos,

	__global cl_real4* f,

	__global real_t* energy,
	__global real_t* rho,
	__global real_t* rhobar,

	__global cl_real4* dc,

	__global cl_real4* bounds,
	__global const int* neighborList,
	__global const int* nNeighbors,
	__global const int* nAtoms,

	__global const int* nValues,

	__global const real_t* phi_pot, // the potentials are passed in as real arrays: x0, xn, invDx, values[?]
	__global const real_t* rho_pot,
	__global const real_t* F_pot,

	const real_t cutoff)
{

    int iAtom = get_global_id(0);
    int iBox = get_global_id(1);

    real_t r, r2, r6;
    real_t fr, e_i;
    real_t rho_i;

    cl_real4 dr;
    cl_real4 drBox;
    cl_real4 f_i;

    real_t rCut = cutoff;
    real_t rCut2 = rCut*rCut;
    real_t rhoTmp;
    real_t phiTmp;
    real_t dPhi, dRho;
    real_t fi, fiprime;

    int i;
    int j_local;

    int i_offset;
    int iParticle;

    i_offset = iBox*N_MAX_ATOMS;
    iParticle = i_offset + iAtom;

    // zero out forces on particles
    f_i.x = 0.0;
    f_i.y = 0.0;
    f_i.z = 0.0;

    rho_i = 0.0;

    e_i = 0.0;


    if (iAtom < nAtoms[iBox])
    {// each thread executes on a single atom in the box

#if(KERN_DIAG > 0) 
	printf("i = %d, %f, %f, %f\n", iParticle, pos[iParticle].x, pos[iParticle].y, pos[iParticle].z);

	printf("iBox = %d, nNeighbors = %d\n", iBox, nNeighbors[iBox]);
#endif

	for (int j = 0; j<nNeighbors[iBox]; j++)
	{// loop over neighbor cells
	    int jBox = neighborList[iBox*N_MAX_NEIGHBORS + j];
	    int jOffset = jBox*N_MAX_ATOMS;

	    // compute box center offsets
	    drBox = dc[iBox] - dc[jBox];

	    // correct for periodic 
	    if(PERIODIC)
	    {
		if (drBox.x<-0.5*bounds[0].x) drBox.x += bounds[0].x;
		else if (drBox.x > 0.5*bounds[0].x ) drBox.x -= bounds[0].x;
		if (drBox.y<-0.5*bounds[0].y) drBox.y += bounds[0].y;
		else if (drBox.y > 0.5*bounds[0].y ) drBox.y -= bounds[0].y;
		if (drBox.z<-0.5*bounds[0].z) drBox.z += bounds[0].z;
		else if (drBox.z > 0.5*bounds[0].z ) drBox.z -= bounds[0].z;
	    }

#if(KERN_DIAG > 0) 
	    printf("dxbox, dybox, dzbox = %f, %f, %f\n", dxbox, dybox, dzbox);
#endif

	    for (int jAtom = 0; jAtom<nAtoms[jBox]; jAtom++)
	    {// loop over all groups in neighbor cell 

		int jParticle = jOffset + jAtom; // global offset of particle

		dr = pos[iParticle] - pos[jParticle] + drBox;

#if(KERN_DIAG > 0) 
		printf("dx, dy, dz = %f, %f, %f\n", dr.x, dr.y, dr.z);
		printf("i = %d, j = %d, %f, %f, %f\n", 
		iParticle, jParticle, pos[jParticle].x, pos[jParticle].y, pos[jParticle].z);
#endif

		r2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;

		if ( r2 <= rCut2 && r2 > 0.0)
		{// no divide by zero

#if(KERN_DIAG > 0) 
		    printf("%d, %d, %f\n", iParticle, jParticle, r2);
		    printf("r2, rCut2 = %f, %f\n", r2, rCut2);
#endif

		    r = sqrt(r2);

		    eamInterpolateDeriv(r, phi_pot, nValues[0], &phiTmp, &dPhi);
		    eamInterpolateDeriv(r, rho_pot, nValues[1], &rhoTmp, &dRho);

#if(KERN_DIAG > 0) 
		    printf("%d, %d, %f\n", iParticle, jParticle, r2);
		    printf("iParticle = %d, jParticle = %d, phiTmp = %f, dPhi = %f\n", iParticle, jParticle, phiTmp, dPhi);
		    printf("iParticle = %d, jParticle = %d, rhoTmp = %f, dRho = %f\n", iParticle, jParticle, rhoTmp, dRho);
#endif

		    f_i.x += dPhi*dr.x/r;
		    f_i.y += dPhi*dr.y/r;
		    f_i.z += dPhi*dr.z/r;

		    e_i += phiTmp;

		    rho_i += rhoTmp;

		} else {
		}


	    } // loop over all atoms
	} // loop over neighbor cells

	f[iParticle] = f_i;

	// since we loop over all particles, each particle contributes 1/2 the pair energy to the total
	//energy[iParticle] = e_i*0.5;

	rho[iParticle] = rho_i;

#if(PASS_2 == 0)
	eamInterpolateDeriv(rho_i,F_pot,nValues[2],&fi,&fiprime);
	rhobar[iParticle] = fiprime; // update rhoprime 
	//update energy terms 
	energy[iParticle] = e_i*0.5 + fi;
#else
	energy[iParticle] = e_i*0.5;
#endif
    }
}

__kernel void EAM_Force_2_AoS(
	__global real_t* rhobar,
	__global real_t* energy,
	__global real_t* rho,
	__global const int* nAtoms,
	__global const real_t* F_pot,
	__global const int* nValues
	)
{

    int iAtom = get_global_id(0);
    int iBox = get_global_id(1);

    real_t fi, fiprime;


    int i;

    int i_offset;
    int iParticle;

    i_offset = iBox*N_MAX_ATOMS;
    iParticle = i_offset + iAtom;

    if (iAtom < nAtoms[iBox])
    {// each thread executes on a single atom in the box


#if(KERN_DIAG > 0) 
	printf("iBox = %d, iAtom = %d\n", iBox, iAtom);
#endif

	eamInterpolateDeriv(rho[iParticle],F_pot,nValues[2],&fi,&fiprime);
	rhobar[iParticle] = fiprime; // update rhoprime 
	//update energy terms 
	energy[iParticle] += fi;

    }

}
__kernel void EAM_Force_3_AoS(
	__global cl_real4* pos,

	__global cl_real4* f,

	__global real_t* fi,

	__global cl_real4* dc,

	__global const cl_real4* bounds,
	__global const int* neighborList,
	__global const int* nNeighbors,
	__global const int* nAtoms,

	__global const int* nValues,
	__global const real_t* rho_pot,
	const real_t cutoff) 
{

    int iAtom = get_global_id(0);
    int iBox = get_global_id(1);

    cl_real4 dr;
    real_t r, r2;

    cl_real4 drBox;

    // accumulate local force value
    cl_real4 f_i;

    real_t rCut = cutoff;
    real_t rCut2 = rCut*rCut;
    real_t rhoTmp, dRho;

    int i;
    int j_local;

    real_t rTmp, rhoijprime;

    // global offset of local thread
    int i_offset = iBox*N_MAX_ATOMS;
    int iParticle = i_offset + iAtom;

    if (iAtom < nAtoms[iBox])
    {// each thread executes on a single atom in the box

	// zero out forces on particles
	f_i.x = f[iParticle].x;
	f_i.y = f[iParticle].y;
	f_i.z = f[iParticle].z;

#if(KERN_DIAG > 0) 
	printf("i = %d, %f, %f, %f\n", iParticle, pos[iParticle].x, pos[iParticle].y, pos[iParticle].z);
	printf("iBox = %d, nNeighbors = %d\n", iBox, nNeighbors[iBox]);
#endif

	for (int j = 0; j<nNeighbors[iBox]; j++)
	{// loop over neighbor cells
	    int jBox = neighborList[iBox*N_MAX_NEIGHBORS + j];
	    int jOffset = jBox*N_MAX_ATOMS;

	    // compute box center offsets
	    drBox = dc[iBox] - dc[jBox];

	    // correct for periodic 
	    if(PERIODIC)
	    {
		if (drBox.x<-0.5*bounds[0].x) drBox.x += bounds[0].x;
		else if (drBox.x > 0.5*bounds[0].x ) drBox.x -= bounds[0].x;
		if (drBox.y<-0.5*bounds[0].y) drBox.y += bounds[0].y;
		else if (drBox.y > 0.5*bounds[0].y ) drBox.y -= bounds[0].y;
		if (drBox.z<-0.5*bounds[0].z) drBox.z += bounds[0].z;
		else if (drBox.z > 0.5*bounds[0].z ) drBox.z -= bounds[0].z;
	    }

	    // printf("dxbox, dybox, dzbox = %f, %f, %f\n", drBox.x, drBox.y, drBox.z);

	    for (int jAtom = 0; jAtom<nAtoms[jBox]; jAtom++)
	    {// loop over all groups in neighbor cell 

		int jParticle = jOffset + jAtom; // global offset of particle

		dr = pos[iParticle] - pos[jParticle] + drBox;

#if(KERN_DIAG > 0) 
		printf("dx, dy, dz = %f, %f, %f\n", dr.x, dr.y, dr.z);
		printf("i = %d, j = %d, %f, %f, %f\n", iParticle, jParticle, pos[jParticle].x, pos[jParticle].y, pos[jParticle].z);
#endif

		r2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;

		if ( r2 < rCut2 && r2 > 0.0)
		{// no divide by zero

#if(KERN_DIAG > 0) 
		    printf("r2, rCut = %f, %f\n", r2, rCut);
		    printf("%d, %d, %f\n", iParticle, jParticle, r2);
#endif

		    r = sqrt(r2);

		    eamInterpolateDeriv(r, rho_pot, nValues[1], &rhoTmp, &dRho);
		    rhoijprime = dRho;

#if(KERN_DIAG > 0) 
		    printf("%d, %d, %f\n", iParticle, jParticle, r2);
#endif

		    rTmp = (fi[iParticle]+fi[jParticle])*rhoijprime/r;

		    /*
		    f_i.x += (fi[iParticle]+fi[jParticle])*rhoijprime*dr.x/r;
		    f_i.y += (fi[iParticle]+fi[jParticle])*rhoijprime*dr.y/r;
		    f_i.z += (fi[iParticle]+fi[jParticle])*rhoijprime*dr.z/r;
		     */
		       f_i.x += rTmp*dr.x;
		       f_i.y += rTmp*dr.y;
		       f_i.z += rTmp*dr.z;

		} else {
		}


	    } // loop over all atoms in jBox
	} // loop over neighbor cells

	f[iParticle] = f_i;

    } // loop over all atoms in iBox
}

