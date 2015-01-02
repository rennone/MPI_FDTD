#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "simulator.h"
#include "field.h"
#include "function.h"

#include "mpiTM_UPML.h"
#include "mpiTE_UPML.h"
#include "models.h"
#include <sys/time.h>

static void (*update)() = NULL;
static double complex* (*getDataX)() = NULL;
static double complex* (*getDataY)() = NULL;
static double complex* (*getDataZ)() = NULL;
static void (* finishMethod)() = NULL;
static void (* initMethod)() = NULL;
static void (* resetMethod)() = NULL;
static double complex* (*getDrawData)() = NULL;
static double* (* getEpsMethod )() = NULL;

static struct timeval timer1, timer2;

static char *solverDir = "";

static void setMPITMupml(){
  update       = mpi_fdtdTM_upml_getUpdate();
  initMethod   = mpi_fdtdTM_upml_getInit();
  finishMethod = mpi_fdtdTM_upml_getFinish();
  resetMethod  = mpi_fdtdTM_upml_getReset();
  
  getEpsMethod = mpi_fdtdTM_upml_getEps;
  
  getDataX = mpi_fdtdTM_upml_getHx;
  getDataY = mpi_fdtdTM_upml_getHy;
  getDataZ = mpi_fdtdTM_upml_getEz;
  
  getDrawData = getDataZ;
//  printf("MPI TM UPML mode \n");

  solverDir = "MPI_TM_UPML";
}

static void setMPITEupml(){
  update       = mpi_fdtdTE_upml_getUpdate();
  initMethod   = mpi_fdtdTE_upml_getInit();
  finishMethod = mpi_fdtdTE_upml_getFinish();
  resetMethod  = mpi_fdtdTE_upml_getReset();
  
  getEpsMethod = mpi_fdtdTE_upml_getEps;

  getDataX = mpi_fdtdTE_upml_getEx;
  getDataY = mpi_fdtdTE_upml_getEy;
  getDataZ = mpi_fdtdTE_upml_getHz;

  getDrawData = getDataY;
  
//  printf("MPI TE UPML mode \n");

  solverDir = "MPI_TE_UPML";
}

static void setSolver(enum SOLVER solver)
{
  switch(solver){
  case MPI_TM_UPML:
    setMPITMupml();
    break;
  case MPI_TE_UPML:
    setMPITEupml();
    break;
  default:
    printf("error, not implement simulator (simulator.c)\n");
    exit(2);
    break;
  }
}

void simulator_moveDirectory()
{
  makeDirectory(solverDir);
  moveDirectory(solverDir);
}

void simulator_calc(){
  (*update)();
  
  field_nextStep();   //時間を一つ進める  
}

void simulator_setSolver(enum SOLVER solver)
{
  setSolver(solver);    //Solverの設定と初期化
}

void simulator_init(FieldInfo field_info)
{
  //横幅(nm), 縦幅(nm), 1セルのサイズ(nm), pmlレイヤの数, 波長(nm), 計算ステップ
  field_init(field_info); //フィールドの初期化
  models_initModel();   //モデルの初期化
  (*initMethod)(); //Solverの初期化, EPS, Coeffの設定
  
  gettimeofday(&timer1, NULL); //開始時間の取得
}

void simulator_solverInit()
{
  makeDirectory(solverDir);
  moveDirectory(solverDir);
  (*initMethod)();
}

void simulator_reset()
{
  printf("simulator_reset \n");

  field_reset();   //フィールド情報をリセット
  (*resetMethod)(); //シミュレータをリセット
  gettimeofday(&timer1, NULL); //開始時間の取得
}

void simulator_changeModelAndRestart(void) //モデル(のパラメータ)を変更するので, カレントディレクトリを一段上に行く.
{
  moveDirectory(".."); //一段上に行く.  
}

void simulator_finish()
{
  printf("simulator_finish at %d step \n", (int)field_getTime());
  gettimeofday(&timer2,NULL);
  printf("time = %lf \n", timer2.tv_sec-timer1.tv_sec+(timer2.tv_usec-timer1.tv_usec)*1e-6);

  (*finishMethod)();   //メモリの解放等, solverの終了処理
}

double complex* simulator_getDrawingData(void){
  return (* getDrawData)();
}

bool simulator_isFinish(void)
{
  return field_isFinish();
}

//屈折率のマップを取ってくる
double *simulator_getEps()
{
  return (*getEpsMethod)();  
}
