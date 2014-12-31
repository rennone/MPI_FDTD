#ifndef _FIELD_H
#define _FIELD_H
#include <stdio.h>
#include "myComplex.h"
#include <math.h>
#include "bool.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

//入射波のモード
// NO USING NOW
enum WAVE_MODE{
  POINT_LIGHT_IN_CENTER,  //中心に点光源
  SCATTER //散乱波
};

//計算領域に関する物理パラメータ
typedef struct FieldInfo
{
  int width_nm, height_nm; // 領域のサイズ
  int h_u_nm;              //1セルの大きさ
  int pml;                 //pmlレイヤの大きさ(セル数)
  int lambda_nm;        //波長
  int angle_deg;       //入射角度
  int stepNum;             //計算ステップ
}FieldInfo;

//プログラムで扱う計算領域のパラメータ
typedef struct FieldInfo_S
{
  int N_X, N_Y;   //領域のサイズ(セル)
  int N_PX, N_PY; //PMLレイヤを含めた領域サイズ(セル)
  int N_CELL;     //全セル数
  int N_PML;      //PMLレイヤの層の数
  int DX, DY;     //x,yのインデックスのオフセット量(まだ使ってない) todo
}FieldInfo_S;

//MPI分割後の小領域のパラメータ
typedef struct SubFieldInfo_S
{
  int OFFSET_X, OFFSET_Y; //左下からのオフセット量(セル)
  int SUB_N_X, SUB_N_Y;
  int SUB_N_PX, SUB_N_PY;
  int SUB_N_CELL;
  int Rank; //自身のランク
  int RtRank, LtRank, TpRank, BmRank; //周りの領域のプロセスランク
  int SUB_DX, SUB_DY;
}SubFieldInfo_S;

//入射波のパラメータ
typedef struct WaveInfo_S
{
  double Lambda_s; //波長
  double T_s;      //周期
  double Omega_s;  //角周波数
  double K_s;      //波数(非推奨)
  double K_0_s;    //真空中の波数(波数は媒質に依存するため, 真空中の波数であることを明記するために名前を変える)
  double Angle_deg;   //入射角
} WaveInfo_S;

typedef struct NFFInfo
{
  int top, bottom, left, right;
  int cx,cy;
  double RFperC; 
  int arraySize; //必要な配列サイズ
} NTFFInfo;

#define C_0_S 0.7071 //下の変数名長いからこっちにする
static const double LIGHT_SPEED_S = 0.7071;

static const double EPSILON_0_S = 1.0;
static const double MU_0_S = 1.0/C_0_S/C_0_S;
static const double Z_0_S  = 1.41422712488; // = √(μ/ε);

extern int N_X;
extern int N_Y;
extern int N_CELL;
extern int N_PML;
extern int N_PX;
extern int N_PY;

extern int FIELD_SUB_NX;
extern int FIELD_SUB_NY;
extern int FIELD_SUB_NPX;
extern int FIELD_SUB_NPY;
extern int FIELD_SUB_CELL;
extern int FIELD_SUB_OFFSET_X;
extern int FIELD_SUB_OFFSET_Y;

extern int field_getOffsetX();
extern int field_getOffsetY();
extern int field_getSubNx();
extern int field_getSubNy();
extern int field_getSubNpx();
extern int field_getSubNpy();
extern int field_getSubNcell();

//インデックスを取ってくる 
extern int ind(const int, const int);

//フィールドの横,縦の大きさ, 1セルのサイズ, pmlレイヤの数, 波長(nm), 計算ステップ
extern void field_init(FieldInfo field_info);
extern void field_reset(void);
//extern void setField(const int wid, const int hei, const double h, const int pml, const double lambda, const double step);

//pml用のσを取得する関数群
extern double field_sigmaX(double x, double y);
extern double field_sigmaY(double x, double y);
extern double field_pmlCoef(double x, double y);
extern double field_pmlCoef_LXY(double x, double y);

//ns pml用のβ,β*を取得する関数群
extern double field_ns_beta(double alpha, double alpha_aster);
extern double field_ns_beta_aster(double alpha, double alpha_aster);

extern double field_toCellUnit(const double);
extern double field_toPhisycalUnit(const double);
extern void field_nextStep(void);
extern bool field_isFinish(void);

extern void field_setWaveAngle(int deg);

//:getter
extern double field_getT(void);
extern double field_getK(void);
extern double field_getRayCoef(void);
extern double field_getOmega(void);
extern double field_getLambda(void);
extern double field_getWaveAngle(void);
extern double field_getTime(void);
extern double field_getMaxTime(void);
extern NTFFInfo field_getNTFFInfo(void);
extern WaveInfo_S field_getWaveInfo_S(void);
extern SubFieldInfo_S field_getSubFieldInfo_S(void);
extern FieldInfo_S field_getFieldInfo_S(void);
extern FieldInfo field_getFieldInfo(void);

// gapX, gapY : Ex-z, Hx-zは格子点からずれた位置に配置され散る為,格子点からのずれを送る必要がある.
//dotは波の方向とx(とy,z軸)との内積

extern void field_nsScatteredWaveNotUPML(dcomplex *p, double *eps, double gapX, double gapY, double dot);
extern void field_scatteredWaveNotUPML(dcomplex *p, double *eps, double gapX, double gapY);
//散乱波

//UPML専用
extern void field_scatteredWave(dcomplex *p, double *eps, double gapX, double gapY);
extern void field_scatteredPulse(dcomplex *p, double *eps, double gapX, double gapY, double dot); //dotは波の方向とx(とy,z軸)との内積

//点光源
extern dcomplex field_pointLight(void);

//座標->インデックス変換
//(乗算命令があるので何度も呼び出すような処理では使わない方がいい)
extern int field_index(int i, int j);
extern int field_subIndex(int i, int j);

//output method
extern void field_outputElliptic(const char *fileName,double complex* data); //
extern void field_outputAllDataComplex(const char *fileName,double complex* data); //
extern void field_outputAllDataDouble(const char *fileName,double* data); //

extern void field_initMPI();

// for(i=1..N_PX-1)
//   for(j=1..N_PY-1) と同じ(はず)
#define FAST_FOR_FOR(k, fInfo_s) \
  for(int k=fInfo_s.N_PY+1, last = fInfo_s.N_CELL-fInfo_s.N_PY; k<last; k+=2) \
    for(int endRow = k+fInfo_s.N_PY-2; k<endRow; k++)

/*
//--------------------for debug--------------------//
static inline void field_debugPrint(double complex *A)
{
  for(int i=1; i<N_PX; i++){
    for(int j=1; j<N_PY; j++){
      double complex a = A[index(i,j)];
      if(creal(a) > 1.0)
	printf("(%d, %d) : %lf \n", i, j, creal(a));
    }
  }
  }*/
#endif
