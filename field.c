#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "field.h"
#include "models.h"
#include "function.h"

int N_X;
int N_Y;
int N_CELL;
int N_PML;
int N_PX;
int N_PY;

//_u : 物理量変換単位, _s:シミュレーション単位
static const int H_s = 1;
static double H_u;         //1セルの大きさ(nm)
static double time;     //ステップ数
static double ray_coef; //波をゆっくり入れる為の係数;
static double lambda_s; //波長 
static double k_s;      //波数 
static double w_s;      //角周波数
static double T_s;      //周期 

static double maxTime;

static NTFFInfo ntff_info;
static FieldInfo fieldInfo;
static FieldInfo_S    fieldInfo_s;
static WaveInfo_S     waveInfo_s;
static SubFieldInfo_S subFieldInfo_s;

//:public------------------------------------//
void field_setWaveAngle(int deg) {
  fieldInfo.angle_deg = deg;
  waveInfo_s.Angle_deg = deg;
}

double field_getT() {  return waveInfo_s.T_s; }
double  field_getK(){  return waveInfo_s.K_s;}
double  field_getRayCoef(){  return ray_coef;}
double  field_getOmega(){  return waveInfo_s.Omega_s;}
double  field_getLambda(){  return lambda_s;}
double  field_getWaveAngle(){  return waveInfo_s.Angle_deg;}
double  field_getTime(){  return time;}
double  field_getMaxTime(){  return maxTime;}
NTFFInfo  field_getNTFFInfo()          { return ntff_info;}
WaveInfo_S field_getWaveInfo_S()         { return waveInfo_s;}
SubFieldInfo_S field_getSubFieldInfo_S() { return subFieldInfo_s;}
FieldInfo_S field_getFieldInfo_S()       { return fieldInfo_s;}
FieldInfo field_getFieldInfo()   { return fieldInfo;}

int field_getOffsetX(){  return subFieldInfo_s.OFFSET_X; }
int field_getOffsetY(){  return subFieldInfo_s.OFFSET_Y; }
int field_getSubNx()  {  return subFieldInfo_s.SUB_N_X;  }
int field_getSubNy()  {  return subFieldInfo_s.SUB_N_Y;  }
int field_getSubNpx() {  return subFieldInfo_s.SUB_N_PX; }
int field_getSubNpy() {  return subFieldInfo_s.SUB_N_PY; }
int field_getSubNcell(){ return subFieldInfo_s.SUB_N_CELL;}

int field_subIndex(int i, int j){
    return i*subFieldInfo_s.SUB_N_PY + j;
}

int field_index(int i, int j){
    return i*fieldInfo_s.N_PY + j;
}

int ind(const int i, const int j){  return i*N_PY + j;}//1次元配列に変換

double field_toCellUnit(const double phisycalUnit){
  return phisycalUnit/fieldInfo.h_u_nm; //セル単位に変換 
}
double field_toPhisycalUnit(const double cellUnit){
  return cellUnit*fieldInfo.h_u_nm;//物理単位(nm)に変換
}

void field_reset()
{
  time = 0;
  ray_coef = 0;
}

void field_init(FieldInfo field_info)
{
  //フィールド情報の保存(最初にしないとtoCellUnit, PhisicalUnitが使えない.
  fieldInfo = field_info;
  
  //領域のシミュレータパラメータを計算
  fieldInfo_s.N_PX  = field_toCellUnit(fieldInfo.width_nm);
  fieldInfo_s.N_PY  = field_toCellUnit(fieldInfo.height_nm);
  fieldInfo_s.N_PML = fieldInfo.pml;
  fieldInfo_s.N_X   = fieldInfo_s.N_PX - 2*fieldInfo_s.N_PML;
  fieldInfo_s.N_Y   = fieldInfo_s.N_PY - 2*fieldInfo_s.N_PML;
  fieldInfo_s.N_CELL= fieldInfo_s.N_PY*fieldInfo_s.N_PX;
  fieldInfo_s.DX    = fieldInfo_s.N_PY;
  fieldInfo_s.DY    = 1;
  
  //入射波パラメータの計算
  waveInfo_s.Lambda_s = field_toCellUnit(fieldInfo.lambda_nm);
  waveInfo_s.T_s      = waveInfo_s.Lambda_s/C_0_S;
  waveInfo_s.K_s      = 2*M_PI/waveInfo_s.Lambda_s;
  waveInfo_s.K_0_s    = 2*M_PI/waveInfo_s.Lambda_s;
  waveInfo_s.Omega_s  = C_0_S*waveInfo_s.K_s;
  waveInfo_s.Angle_deg   = fieldInfo.angle_deg;


  // 下位バージョンとの互換性の為
  H_u = fieldInfo.h_u_nm;  
  N_PX = field_toCellUnit(fieldInfo.width_nm);
  N_PY = field_toCellUnit(fieldInfo.height_nm);
  N_PML= fieldInfo.pml;  
  N_X = N_PX - 2*N_PML;
  N_Y = N_PY - 2*N_PML;
  N_CELL = N_PX * N_PY; //全セル数 
  time = 0;
  maxTime = fieldInfo.stepNum;
  
  lambda_s = field_toCellUnit(fieldInfo.lambda_nm);
  k_s = 2*M_PI/lambda_s;
  w_s = C_0_S*k_s;
  T_s = 2*M_PI/w_s;

  ray_coef  = 0;  
  
  /* NTFF設定 */
  ntff_info.cx     = N_PX/2;
  ntff_info.cy     = N_PY/2;
  ntff_info.top    = N_PY - N_PML - 5;
  ntff_info.bottom = N_PML + 5;
  ntff_info.left   = N_PML + 5;
  ntff_info.right  = N_PX - N_PML - 5;

  double len = (ntff_info.top - ntff_info.bottom)/2;
  ntff_info.RFperC = len*2;
  ntff_info.arraySize = maxTime + 2*ntff_info.RFperC;
}

dcomplex field_pointLight()
{
  double time = field_getTime();
  double w_s  = field_getOmega();
  double ray_coef = field_getRayCoef();

  return ray_coef*cexp(I*w_s*time);
}

//散乱波
void field_nsScatteredWaveNotUPML(dcomplex *p, double *eps, double gapX, double gapY, double dot)
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  double time = field_getTime();
  double w_s  = field_getOmega();
  double ray_coef = field_getRayCoef() * dot;
  double k_s = field_getK();
  double rad = field_getWaveAngle()*M_PI/180.0;
  double ks_cos = cos(rad)*k_s, ks_sin = sin(rad)*k_s;//毎回計算すると時間かかりそうだから代入しておく  
  for(int i=1; i<fInfo_s.N_PX-1; i++) {
    for(int j=1; j<fInfo_s.N_PY-1; j++) {
      int k = field_index(i,j);
      double kr = (i+gapX)*ks_cos+(j+gapY)*ks_sin;
      double n  = sqrt( eps[k] / EPSILON_0_S );
     // 光速 C_0_S が変化するので, kをそれに合わせる必要がある
      double u0 = sin(w_s*0.5) / sin(  k_s*0.5);
      double u1 = sin(w_s*0.5) / sin(n*k_s*0.5);
      double _n = u0 / u1;
      p[k] += ray_coef*
        ( 1.0/(_n*n) - 1.0)* ( cexp( I*(kr-w_s*(time+1.0)) ) - cexp( I*(kr-w_s*(time)) ));  //p[k] -= かも(岡田さんのメール参照)
    }
  }
}

void field_scatteredWaveNotUPML(dcomplex *p, double *eps, double gapX, double gapY)
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  double time = field_getTime();
  double w_s  = field_getOmega();
  double ray_coef = field_getRayCoef();
  double k_s = field_getK();
  double rad = field_getWaveAngle()*M_PI/180.0;
  double ks_cos = cos(rad)*k_s, ks_sin = sin(rad)*k_s;//毎回計算すると時間かかりそうだから代入しておく  
  for(int i=1; i<fInfo_s.N_PX-1; i++) {
    for(int j=1; j<fInfo_s.N_PY-1; j++) {
      int k = field_index(i,j);
      double kr = (i+gapX)*ks_cos+(j+gapY)*ks_sin;
      p[k] += ray_coef*
        ( EPSILON_0_S / eps[k] - 1.0)* ( cexp( I*(kr-w_s*(time+0.5)) ) - cexp( I*(kr-w_s*(time-0.5)) ));  //p[k] -= かも(岡田さんのメール参照)
    }
  }
}

//単一波長の散乱波
// gapX, gapY : Ex-z, Hx-zは格子点からずれた位置に配置され散る為,格子点からのずれを送る必要がある.
// 分割領域では使えない->gapをさらに領域のオフセットだけずらせば使えそう
// UPML専用
void field_scatteredWave(dcomplex *p, double *eps, double gapX, double gapY)
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  double time = field_getTime();
  double w_s  = field_getOmega();
  double ray_coef = field_getRayCoef();
  double k_s = field_getK();
  double rad = field_getWaveAngle()*M_PI/180.0;
  double ks_cos = cos(rad)*k_s, ks_sin = sin(rad)*k_s;//毎回計算すると時間かかりそうだから代入しておく  
  for(int i=1; i<fInfo_s.N_PX-1; i++) {
    for(int j=1; j<fInfo_s.N_PY-1; j++) {
      int k = field_index(i,j);
      double kr = (i+gapX)*ks_cos+(j+gapY)*ks_sin;
      p[k] += ray_coef*(EPSILON_0_S/eps[k] - 1.0)*cexp( I*(kr-w_s*time) );  //p[k] -= かも(岡田さんのメール参照)
    }
  }
}

//ガウシアンパルス
// gapX, gapY : Ex-z, Hx-zは格子点からずれた位置に配置され散る為,格子点からのずれを送る必要がある.
// 分割領域では使えない->gapをさらに領域のオフセットだけずらせば使えそう
// UPML専用
void field_scatteredPulse(dcomplex *p, double *eps, double gapX, double gapY, double dot)
{
  double time = field_getTime();
  double w_s  = field_getOmega();
  double rad = field_getWaveAngle()*M_PI/180.0;	//ラジアン変換  

  double cos_per_c = cos(rad)/C_0_S, sin_per_c = sin(rad)/C_0_S;
  const double beam_width = 50; //パルスの幅

  FieldInfo_S fInfo_s      = field_getFieldInfo_S();  
  
  //waveAngleにより, t0の値を変えないとちょうどいいところにピークが来ないため,それを計算.
  const double center_peak = (fInfo_s.N_PX/2.0+gapX)*cos_per_c+(fInfo_s.N_PY/2+gapY)*sin_per_c; //スタートから中心へ進むのにかかる時間

  // TODO : h_uのサイズに応じて変化させないと, パルスが不連続になる可能性が有る.
  // t=500 では h_u = 10nm, 5nmでうまく動いているのでこうしている.
//常に t=500 の時に,領域の中心にピークが来るように初期位相を調整
  const double t0 = -center_peak + 500;
  
  for(int i=1; i<fInfo_s.N_PX-1; i++)
  {
    for(int j=1; j<fInfo_s.N_PY-1; j++)
    {
      int k = field_index(i,j);
      if(EPSILON_0_S == eps[k])
        continue;
      
      const double r = (i+gapX)*cos_per_c+(j+gapY)*sin_per_c-(time-t0); // (x*cos+y*sin)/C - (time-t0)
      const double gaussian_coef = exp( -pow(r/beam_width, 2 ) );
      p[k] += dot*gaussian_coef*(EPSILON_0_S/eps[k] - 1)*cexp(I*r*w_s);     //p[k] -= かも(岡田さんのメール参照)
    }
  } 
}

//----------------PMLに用いるσx,σyの計算--------------------//
double field_sigmaX(const double x, const double y)
{
  const int M = 2;
  if(x<N_PML)
    return pow(1.0*(N_PML-x)/N_PML, M);
  
  else if(N_PML <= x && x < (N_X+N_PML))    
    return 0;
  
  else
    return pow(1.0*(x - (N_PX-N_PML-1))/N_PML, M);
}

double field_sigmaY(const double x, double y)
{
  const int M = 2;
  if(y<N_PML)
    return pow(1.0*(N_PML - y)/N_PML,M);
  
  else if(y>=N_PML && y<(N_Y+N_PML))
    return 0.0;

  else
    return pow(1.0*(y - (N_PY-N_PML-1))/N_PML,M);
}

//pml用の係数のひな形 Δt = 1
//ep_mu εかμ(Eの係数->ε, Hの係数-> μ
//sig  σ
double field_pmlCoef(double ep_mu, double sig)
{
  return (1.0 - sig/ep_mu)/(1.0+sig/ep_mu);
}
double field_pmlCoef_LXY(double ep_mu, double sig)
{
  return 1.0/(ep_mu + sig); // 1.0/{ep_mu(1.0 + sig/ep_mu)}と同じ
}

//NS PML用のβを取得
double field_ns_beta(double alpha, double alpha_aster)
{
  double tanh_a     = tanh(alpha);//(exp(alpha) - exp(-alpha)) / (exp(alpha) + exp(-alpha));
  double tanh_a_ast = tanh(alpha_aster);//(exp(alpha_aster) - exp(-alpha_aster)) / (exp(alpha_aster) + exp(-alpha_aster));
  return tanh_a / (1 + tanh_a*tanh_a_ast);
}
double field_ns_beta_aster(double alpha, double alpha_aster)
{
  double tanh_a     = tanh(alpha);//(exp(alpha) - exp(-alpha)) / (exp(alpha) + exp(-alpha));
  double tanh_a_ast = tanh(alpha_aster);//(exp(alpha_aster) - exp(-alpha_aster)) / (exp(alpha_aster) + exp(-alpha_aster));
  return tanh_a_ast / (1 + tanh_a*tanh_a_ast);
}
    
//------------------getter-------------------------//
 void field_nextStep(void){
  time+=1.0;
  ray_coef = 1.0 - exp(-pow(0.01*time, 2));
}

 bool field_isFinish(void){
  return time >= maxTime;
}

//---------------output method---------------//
void field_outputElliptic(const char *fileName, double complex* data)
{
  printf("output start\n");
  //file open
  FILE *fp;

  if( (fp=fopen(fileName, "w") ) == NULL){
    printf("cannot open file %s \n", fileName);
    exit(2);
  }

  for(int ang=180; ang >=0; ang--){
    double rad = ang*M_PI/180.0;
    double x = 1.2*lambda_s*cos(rad)+N_PX/2.0;
    double y = 1.2*lambda_s*sin(rad)+N_PY/2.0;
    double norm = cnorm(cbilinear(data,x,y,N_PX,N_PY));
    fprintf(fp, "%d %lf \n", 180-ang, norm);    
  }
  
  fclose(fp);
  printf("output to %s end\n", fileName);
}

//呼び出し時の計算領域のすべての値を出力
void field_outputAllDataComplex(const char *fileName, double complex *data)
{
   printf("output all data start\n");
  //file open
  FILE *fp;

  if( (fp=fopen(fileName, "w") ) == NULL){
    printf("cannot open file %s \n", fileName);
    exit(2);
  }

  for(int i=0; i<N_PX; i++){
    for(int j=0; j<N_PY; j++){
      fprintf(fp, "%lf \n",cnorm(data[ind(i,j)]));
    }
  }
  
  fclose(fp);
  printf("output all data to %s end\n", fileName);
}


void field_outputAllDataDouble(const char *fileName, double *data)
{
   printf("output all data double start\n");
  //file open
  FILE *fp;

  if( (fp=fopen(fileName, "w") ) == NULL){
    printf("cannot open file %s \n", fileName);
    exit(2);
  }

  for(int i=0; i<N_PX; i++){
    for(int j=0; j<N_PY; j++){
      if(data[ind(i,j)] != 1.0)
        fprintf(fp, "%lf \n",data[ind(i,j)]);
    }
  }
  
  fclose(fp);
  printf("output all data to %s end\n", fileName);
}

//:private
//MPIにより領域を分割する.
void field_initMPI()
{
  int num_proc;
  int dim       = 2;        //number of dimension is 2
  int procs[2]  = {0,0};    //[0]: x方向の分割数, [1]:y方向の分割数
  int period[2] = {0,0};    //境界条件, 固定境界(1=>周期境界)
  int reorder   = 1;   //re-distribute rank flag
  MPI_Comm grid_comm;

  MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
  MPI_Dims_create(num_proc, dim, procs);
  MPI_Cart_create(MPI_COMM_WORLD, 2, procs, period, reorder, &grid_comm);
  //周りの領域のプロセスランクを取得
  MPI_Cart_shift(grid_comm, 0, 1, &subFieldInfo_s.LtRank, &subFieldInfo_s.RtRank);
  MPI_Cart_shift(grid_comm, 1, 1, &subFieldInfo_s.BmRank, &subFieldInfo_s.TpRank);

  //プロセス座標において自分がどの位置に居るのか求める(何行何列に居るか)
  int coordinates[2];
  MPI_Comm_rank(grid_comm, &subFieldInfo_s.Rank);
  MPI_Cart_coords(grid_comm, subFieldInfo_s.Rank, 2, coordinates);

  //小領域がすべて同じ大きさになるかチェック
  if( fieldInfo_s.N_PX%procs[0] !=0 || fieldInfo_s.N_PY%procs[1] != 0){      
    printf("cannot devide size(%d, %d) by proc(%d,%d) \n",
           fieldInfo_s.N_PX, fieldInfo_s.N_PY, procs[0],procs[1]);
    exit(2);
  }

  //小領域のパラメータ
  subFieldInfo_s.SUB_N_X    = fieldInfo_s.N_PX / procs[0];
  subFieldInfo_s.SUB_N_Y    = fieldInfo_s.N_PY / procs[1];

  subFieldInfo_s.SUB_N_PX   = subFieldInfo_s.SUB_N_X + 2; //のりしろの分2大きい
  subFieldInfo_s.SUB_N_PY   = subFieldInfo_s.SUB_N_Y + 2; //のりしろの分2大きい
  subFieldInfo_s.SUB_N_CELL = subFieldInfo_s.SUB_N_PX*subFieldInfo_s.SUB_N_PY;  
  subFieldInfo_s.OFFSET_X  = coordinates[0] * subFieldInfo_s.SUB_N_X; //ランクのインデックスではなく, セル単位のオフセットなのでSUB_N_Xずれる
  subFieldInfo_s.OFFSET_Y  = coordinates[1] * subFieldInfo_s.SUB_N_Y;

  subFieldInfo_s.SUB_DX    = subFieldInfo_s.SUB_N_PY;
  subFieldInfo_s.SUB_DY    = 1;
}










