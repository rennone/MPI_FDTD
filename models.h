#ifndef _MODELS_H
#define _MODELS_H

#include "bool.h"
enum MODEL
{
  NO_MODEL,
  MIE_CYLINDER,
  LAYER,
  MORPHO_SCALE,
  CONCENTRIC_CIRCLE,
  ZIGZAG,
  TRACE_IMAGE
};

enum MODE
{
  D_X, //x方向に線積分
  D_Y, //y方向に線積分
  D_XY //面積分 
};

//モデルを定義
extern void models_setModel(enum MODEL model);

//(x,y)における誘電率を取得
extern double models_eps(double x, double y, enum MODE mode);
extern bool models_isFinish(void);
extern void models_moveDirectory(void); //所定のディレクトリまで移動する為の関数
extern void models_needSize(int *x_nm,int *y_nm);
extern void models_initModel();
extern void models_evaluate(double **reflec, int stLambda, int enLamba); //評価
extern void models_update(void); //アップデート
//データの吐き出し
//extern void models_output(FILE *fp, double complex *data);
#endif
