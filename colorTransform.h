#ifndef _COLOR_TRANSFORM_H_
#define _COLOR_TRANSFORM_H_
#include "bool.h"

extern void colorTransform_init(); //初期化
extern void colorTransform_trans(int lambda, double reflec, double *r, double *g, double *b);

extern bool colorTransform_rgbTohsv(double r, double g, double b, double *h, double *s, double *v);
#endif
