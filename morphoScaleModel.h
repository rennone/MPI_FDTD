#ifndef MORPHO_SCALE_H
#define MORPHO_SCALE_H

#include "bool.h"
extern double ( *morphoScaleModel_EPS(void))(double, double, int, int);
extern bool morphoScaleModel_isFinish(void);
extern void morphoScaleModel_moveDirectory(void);
extern void morphoScaleModel_needSize(int *, int*);
extern void morphoScaleModel_init();
#endif
