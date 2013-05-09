#ifndef __CHEBY_H_
#define __CHEBY_H_

#define NUM_CHEBY 16

struct EamCheby *setChebPot(EamPotential *pot, int n);

struct PotentialArray *genCheb(PotentialArray *pot, int n);

struct PotentialArray *genDer(PotentialArray *ch);

real_t eamCheb(PotentialArray *cheb, real_t x);

void chDer(real_t a, real_t b, real_t *c, real_t *cder, int n);

void readArraySizes( EamCheby *chPot, int* numCh, int* numRef);

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
	real_t* xFRef) ;

void readMishinCh(EamCheby *chPot);

void genDerivs(EamCheby *chPot);

void computeCompare(
	real_t* range, 
	real_t* ch, 
	real_t* dch,
	real_t* ref, 
	real_t* x, 
	real_t* ref_val_h,
	real_t* ref_dval_h,
	int n_ref, 
	int n_ch);

static void eamInterpolateDerivlocal(real_t r,
	real_t* values,
	int nValues,
	real_t *range,
	real_t *value1, 
	real_t *f1);

#endif
