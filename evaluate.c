#include "evaluate.h"
#include "colorTransform.h"
#include <math.h>
#include "bool.h"

// UNUSED
// 補完関数 : 0~1 => 0~1へ写像(例 : f(x) = x^n など)
static double func(double a)
{
  return a;
}

//以下どれか一つのみ定義する
//#define RED_EVAL
//#define GREEN_EVAL
#define BLUE_EVAL

// HSV値から評価値を返す
double evaluate_eval(double h, double s, double v, enum EvalKinds k)
{
  // hsv表色系での重み付け
   //(hの)平均, 分散
  double mean, variance;
  switch(k)
  {
  case EVAL_RED:
    mean     = 0;
    variance = 30;
    break;
  case EVAL_GREEN:
    mean     = 135;
    variance = 40;
    break;
  case EVAL_BLUE:
    mean     = 225;
    variance = 30;
    break;
    //ありえない値
  default:
    mean     = -10000;
    variance = 1;
  }

  double c = 100; //中心の重みがcになるように係数設定
  double d = 0; //中心から離れると,dに漸近するように差分設定

  // 円周上における距離なので
  // 例えば0と350との距離が10になるように,差が180を越えると+-360する.  
  if( h < mean - 180 )
    h += 360;
  if( h > mean + 180 )
    h -= 360;
  
  // hによるガウス分布から評価する
  double value = (c-d) * exp( -pow(h-mean,2)/(2*pow(variance,2)) ) + d;

  //彩度と明度は大きい方が良い(鮮やか)のでそのまま重み付けする.
  return value * s * v;
}

//各波長の角反射率
// reflec[lambda][degree] : 波長λのdegree°方向の反射率
double evaluate_evaluate(double **reflec, int stLambda, int enLambda)
{
  double value = 0;
  
// 90°を中心に左右 deg° の反射角における積分値の全体に対する割合を用いて評価する
  const int eval_deg = 40; //左右eval_degの視野角で評価
  
  for(int i = 0; i <= enLambda - stLambda; i++)
  {
    double red, green, blue;    
    for(int deg=90-eval_deg; deg<=90+eval_deg; deg++)
    {
      double r,g,b;
      colorTransform_trans(i+stLambda, reflec[i][deg], &r, &g, &b);
      red   += r;
      green += g;
      blue  += b;      
    }    
    double h,s,v;
    bool res = colorTransform_rgbTohsv(red,green,blue, &h, &s, &v);      
    if( !res )
      continue;

    value += evaluate_eval(h, s, v, EVAL_BLUE);
  }
  // degの範囲を変えても最大が変わらないように正規化
  value = value / (2.0*eval_deg+1);
  return value;
}
