#ifndef MPI_NTFF_TM_H
#define MPI_NTFF_TM_H

#include <stdio.h>
#include "myComplex.h"

extern void mpiNtffTM_init();

//resEzに周波数領域の遠方解を代入する.
extern void mpiNtffTM_Frequency( dcomplex *Hx, dcomplex *Hy, dcomplex *Ez, dcomplex resEz[360]);

//時間領域の遠方解のupdate処理(Ux,Uy,Wzを更新)
extern void mpiNtffTM_TimeCalc(dcomplex *Hx, dcomplex *Hy, dcomplex *Ez, dcomplex *Ux, dcomplex *Uy, dcomplex *Wz);

//時間領域の遠方解 Ux, Uy, WzからEth, Ephを求める.
extern void mpiNtffTM_TimeTranslate(dcomplex *Ux, dcomplex *Uy, dcomplex *Wz, dcomplex *Eth, dcomplex *Eph);

//時間領域の遠方解 Ux, Uy, WzからEthを求め
//それをfftした結果から波長,反射角ごとのEthのノルムを返す
// res[λ][deg]でアクセスできる. 0 <= deg <= 359, stLambda <= λ <= enLamba
double** mpiNtffTM_TimeCalcNorm(dcomplex *Ux, dcomplex *Uy, dcomplex *Wz, int *stLambda, int *enLamba);

#endif
