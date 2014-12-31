#ifndef _SIMULATOR_H
#define _SIMULATOR_H
#include<complex.h>
#include "bool.h"
#include "models.h"
#include "field.h"

enum SOLVER
{
  MPI_TM_UPML,
  MPI_TE_UPML 
};

extern void simulator_init(FieldInfo field_info);
//extern void simulator_init(int width, int height , double h_u, int pml, double lambda, int step, enum MODEL model, enum SOLVER solver);
extern void simulator_calc(void);
extern bool simulator_isFinish(void);
extern void simulator_finish(void);
extern void simulator_reset(void);      //条件はそのままで, 電磁波の状態だけリセットする.
extern void simulator_solverInit(void); //epsとかも計算し直してから再スタートする用
extern void simulator_changeModelAndRestart(void); //モデル(のパラメータ)を変更するので, カレントディレクトリを一段上に行く.
extern double complex* simulator_getDrawingData();
extern double *simulator_getEps();
extern void simulator_moveDirectory();
extern void simulator_setSolver(enum SOLVER solverType);
#endif
