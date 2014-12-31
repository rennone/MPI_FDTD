#include "concentricCircleModel.h"
#include "field.h"
#include "function.h"
#include "parser.h"
#include <math.h>
#include <stdlib.h>
#include <mpi.h>

//同心円の二つの半径[0]が内側
static double radius_s[2]; // radius_s[0] < radius[1]
static double n[2];
static double ep[2];
static int num[2]; //行の数と列の数
static double padding_s; //円同士の隙間

//一つの円
//円の中心,　座標, 
static double circle_eps(double cx, double cy, double x, double y, int col, int row)
{  
  double dx = x-cx;
  double dy = y-cy;
  //2乗距離
  double len = dx*dx+dy*dy;

  //空気中
  //中心との距離がradius[1]+1セル以上なら,そのセルは完全に外側の媒質の外 
  if(len >= (radius_s[1]+1)*(radius_s[1]+1))
    return EPSILON_0_S;

  //内側の媒質
  //中心との距離がradius[0]-1セル以下なら,そのセルは完全に内側の媒質の外
  if(len <= (radius_s[0]-1)*(radius_s[0]-1))
    return ep[0];

  //外側の媒質
//中心との距離がradius[0]+1セル以上かつ,中心との距離がradius[1]-1セル以下なら,
  //完全に媒質1の中  
  if(len >= (radius_s[1]+1)*(radius_s[1]+1) &&
     len <= (radius_s[1]-1)*(radius_s[1]-1))
    return ep[1];

  //空気との境界, もしくは媒質1,2の境界上では,分割して考える
  //さらに32*32分割し媒質内と媒質外の数を求めepsilonを決定する
  double sum[2]={0,0};
  for(double i=-16+0.5; i<16; i+=1){
    for(double j=-16+0.5; j<16; j+=1){
      double lenSqr = pow(dx+col*i/32.0, 2.0) + pow(dy+row*j/32.0, 2.0);
      if( lenSqr < radius_s[0]*radius_s[0]){
	sum[0]+=1;
      }
      else if( lenSqr == radius_s[1]*radius_s[1] ){
        sum[1]+=0.5;
        sum[0]+=0.5;
      }      
      else if( lenSqr < radius_s[1]*radius_s[1] )
      {
        sum[1]+=1;
      } else if( lenSqr == radius_s[1]*radius_s[1])
      {
        sum[1]+=0.5; //もう半分は空気なので
      }      
    }
  }

  sum[0] /= 32.0*32.0;
  sum[1] /= 32.0*32.0;
  return sum[0]*ep[0] + sum[1]*ep[1] + EPSILON_0_S*(1-sum[0]-sum[1]);
}

static double eps(double x, double y, int col, int row)
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  
  //PML領域は真空と仮定
  if(x < fInfo_s.N_PML || y < fInfo_s.N_PML ||
     x > fInfo_s.N_X+fInfo_s.N_PML || y > fInfo_s.N_Y + fInfo_s.N_PML)
    return EPSILON_0_S;

  return circle_eps(fInfo_s.N_PX/2, fInfo_s.N_PY/2, x, y, col, row);
}

double ( *concentricCircleModel_EPS(void) )(double, double, int, int)
{
  radius_s[0] = field_toCellUnit(500);
  radius_s[1] = field_toCellUnit(1000);

  n[0] = 3.882; //Si
  n[1] = 1.457; //Si02

  ep[0] = n[0]*n[0]*EPSILON_0_S;
  ep[1] = n[1]*n[1]*EPSILON_0_S;

  return eps;
}

bool concentricCircleModel_isFinish()
{
  return true;
}

void concentricCircleModel_moveDirectory()
{  
  return;
}
