
#ifndef __READ_H
#define __READ_H

SimFlat *fromFileASCII(Command cmd, struct BasePotential *pot);
SimFlat *fromFileGzip(Command cmd, struct BasePotential *pot);
SimFlat *fromFileTim(Command cmd, struct BasePotential *pot);

#endif
