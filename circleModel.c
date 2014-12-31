#include <math.h>
#include "circleModel.h"
#include "field.h"
#include "function.h"

#define ST_RADIUS_NM 500
#define EN_RADIUS_NM 500
#define DELTA_RADIUS_NM 100

#define N 1.6
static int radius_nm = ST_RADIUS_NM;

//フィールド領域が決まってから決まる
static double radius_s;
static double epsilon_s;
static double posx_s;
static double posy_s;

static double eps(double, double, int, int);

double (*circleModel_EPS())(double, double, int , int)
{  
  return eps;
}

//col : D_Xモード, row : D_Yモード
static double eps(double x, double y, int col, int row)
{
  if(x < N_PML || y < N_PML || x > N_X+N_PML || y > N_Y + N_PML)
    return EPSILON_0_S;

  double dx = x-posx_s;
  double dy = y-posy_s;
  //2乗距離
  double len = dx*dx+dy*dy;

  //中心との距離がr+1セル以上なら,そのセルは完全に媒質の外 
  if(len >= (radius_s+1)*(radius_s+1))
    return EPSILON_0_S;

  //中心との距離がr-1セル以下なら,そのセルは完全に媒質の外 
  if(len <= (radius_s-1)*(radius_s-1))
    return epsilon_s;

  //さらに32*32分割し媒質内と媒質外の数を求めepsilonを決定する
  double split = 10;
  double half_split = split/2;
  double sum=0;
  for(double i=-half_split+0.5; i<half_split; i+=1){
    for(double j=-half_split+0.5; j<half_split; j+=1){
      if(pow(dx+col*i/split, 2.0) + pow(dy+row*j/split, 2.0) <= radius_s*radius_s)
	sum+=1;
    }
  }
  
  sum /= split*split;
  return epsilon_s*sum + EPSILON_0_S*(1-sum);
}

bool circleModel_isFinish()
{
  radius_nm += DELTA_RADIUS_NM;
  return radius_nm > EN_RADIUS_NM;
}

void circleModel_moveDirectory()
{
  char buf[512];
  sprintf(buf,"radius_%dnm", radius_nm);
  makeDirectory(buf);
  moveDirectory(buf);
}

void circleModel_init(void)
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  radius_s  = field_toCellUnit(radius_nm);  
  posx_s    = fInfo_s.N_PX/2;
  posy_s    = fInfo_s.N_PY/2;
  epsilon_s = N*N*EPSILON_0_S;
}

void circleModel_needSize(int *x_nm,int *y_nm)
{
  //検証時には半径*1.2の位置を使っているので多めに取っておく
  *x_nm = 1.5*radius_nm*2;
  *y_nm = 1.5*radius_nm*2;
}
