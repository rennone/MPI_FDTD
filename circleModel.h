#ifndef _CIRCLE_H
#define _CIRCLE_H

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include "bool.h"
extern double (*circleModel_EPS(void))(double, double, int,int);
extern void (*circleModel_output(void))(FILE *, double complex*);
extern bool circleModel_isFinish(void);
extern void circleModel_moveDirectory(void);
extern void circleModel_init(void);
extern void circleModel_needSize(int *x_nm,int *y_nm);
#endif
