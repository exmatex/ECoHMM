
#ifndef __EAM_H
#define __EAM_H

EamPotential *eamReadASCII(char *dir, char *potName);
void eamDestroy(void **inppot);
int eamForce(void *s);

void destroyPotentialArray(struct PotentialArray **a, int doubleFlag);

struct BasePotential *setEamPot(char *dir, char *file);

struct EamPotential *getEamPot();

struct PotentialArray *allocPotentialArray(int n, real_t x0, real_t xn, real_t invDx);

struct PotentialArray *getPotentialArrayFromBinaryFile(char *file);

struct PotentialArray *getPotentialArrayFromFile(char *file);

double getMassFromFile(char *file);

double getLatFromFile(char *file);

EamPotential *eamReadASCII(char *dir, char *potName);

real_t eamCheb(PotentialArray *cheb, real_t x);

inline void eamInterpolateDeriv(struct PotentialArray *a, real_t r, int iType, int jType, real_t *value1, real_t *f1);

void adiffpot(char *name, PotentialArray *a, PotentialArray *b);

void eamComparePots(EamPotential *a, EamPotential *b);


#endif
