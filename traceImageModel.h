#ifndef TRACE_IMAGE_H
#define TRACE_IMAGE_H

#include "bool.h"
extern double ( *traceImageModel_EPS(void))(double, double, int, int);
extern bool traceImageModel_isFinish(void);
extern void traceImageModel_moveDirectory(void);
extern void traceImageModel_needSize(int *, int*);
extern void traceImageModel_init();
#endif
