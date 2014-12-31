#ifndef _SOLVER_H
#define _SOLVER_H

#include "bool.h"
#include "myComplex.h"

typedef struct Solver
{
  void (*update)();
  void (*finish)();
  void (*init)();
  void (*reset)();
  
  dcomplex* (*getDataX)();
  dcomplex* (*getDataY)();
  dcomplex* (*getDataZ)();
  dcomplex* (*getEpsX)();
  dcomplex* (*getEpsY)();
  dcomplex* (*getEpsZ)();  
} Solver;

#endif
