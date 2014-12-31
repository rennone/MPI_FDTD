#ifndef CONCENTRIC_CIRCLE_MODEL_H
#define CONCENTRIC_CIRCLE_MODEL_H

#include "bool.h"

extern double ( *concentricCircleModel_EPS(void) )(double, double, int, int);
extern bool concentricCircleModel_isFinish(void);
extern void concentricCircleModel_moveDirectory(void);
//extern double conc_eps(double x, double y, int col, int row);
#endif
