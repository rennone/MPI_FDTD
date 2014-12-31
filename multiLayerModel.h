#ifndef MULTI_LAYER_MODEL_H
#define MULTI_LAYER_MODEL_H

#include "bool.h"
extern double ( *multiLayerModel_EPS(void))(double, double, int, int);
extern bool multiLayerModel_isFinish(void);
extern void multiLayerModel_moveDirectory();
extern void multiLayerModel_needSize(int *, int*);
extern void multiLayerModel_init();

//評価
extern void multiLayerModel_evaluate(double **reflec, int stLambda, int enLamba);
//毎フレーム実行される
extern void multiLayerModel_update();
#endif
