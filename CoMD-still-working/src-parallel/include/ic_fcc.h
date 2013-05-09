
#ifndef __IC_FCC_H
#define __IC_FCC_H

real_t ran1(long* idum);

real_t gasdev(long* idum);

void applyKineticEnergy(double targetKineticEnergy, SimFlat* s);

void assignTKE(Command cmd, SimFlat* s);

SimFlat* createFccLattice(Command cmd, struct BasePotential* pot);

#endif
