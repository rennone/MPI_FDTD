#include "colorTransform.h"
#include "function.h"

static int stLambda = -1; //最小λ
static int enLambda = -1; //最大λ

#define LAMBDA_RANGE 500
// XYZSの表
static double X[LAMBDA_RANGE];
static double Y[LAMBDA_RANGE];
static double Z[LAMBDA_RANGE];
static double S[LAMBDA_RANGE];
static double K = 0;

// 引数 : rgb, hsv (0~1)
// 戻り値 : false => hがundefined
bool colorTransform_rgbTohsv(double r, double g, double b, double *h, double *s, double *v)
{
  r = max( 0.0, min( 1.0, r) );
  g = max( 0.0, min( 1.0, g) );
  b = max( 0.0, min( 1.0, b) );
  
  double mx = max(r, max(g, b));
  double mn = min(r, min(g, b));

  if( mn == mx || mx <= 0.00000001 ){
    *h = -1;
    *s = 0;
    *v = 0;
    return false;
  }

  if( mn == b )
    *h = 60 * (g - r) / (mx - mn) + 60;    
  else if( mn == r)
    *h = 60 * (b - g) / (mx - mn) + 180;    
  else
    *h = 60 * (r - b) / (mx - mn) + 300;
    
  if(*h < 0)
    *h += 360;
  if(*h > 360)
    *h -= 360;
    
  *v = mx;
  *s = (mx - mn) / mx;
  return true;
}

void colorTransform_init()
{
  FILE *fp = FileOpen("XYZS.csv", "r");

  double _k = 0;
  for(int i=0; i<LAMBDA_RANGE; i++)
  {
    double lambda=0, x=0, y=0, z=0, s=0;
    if( fscanf(fp, "%lf,%lf,%lf,%lf, %lf", &lambda, &x, &y, &z, &s) == EOF )
      break;

    if( i == 0)
      stLambda = (int)lambda;
    
    X[i] = x;
    Y[i] = y;
    Z[i] = z;
    S[i] = s;
    enLambda = (int)lambda; //enLambdaは常に更新
    _k += y*s;
  }

  K = 200.0 / _k;
}

void colorTransform_trans(int lambda, double reflec, double *r, double *g, double *b)
{  
  if( lambda < stLambda || lambda > enLambda)
  {
    *r = 0;
    *g = 0;
    *b = 0;
  }
  else
  {
    double s = S[lambda - stLambda];
    double x = K*s*reflec*X[lambda - stLambda];
    double y = K*s*reflec*Y[lambda - stLambda];
    double z = K*s*reflec*Z[lambda - stLambda];

    //行列変換
    *r =  2.3655*x - 0.8971*y - 0.4683*z;
    *g = -0.5151*x + 1.4264*y + 0.0887*z;
    *b =  0.0052*x - 0.0144*y + 1.0089*z;  
  }
}
