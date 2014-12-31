#include "traceImageModel.h"
#include "field.h"
#include "function.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define FILE_NAME "traceImage1"
#define ST_SCALE_WIDTH_NM 300
#define EN_SCALE_WIDTH_NM 1000
#define DELTA_SCALE_WIDTH_NM 100
static int scale_width_nm = ST_SCALE_WIDTH_NM;
static int scale_width_px;
static double nmPerPixel, pixelPerCell;
static int width_px, height_px;
static double width_s, height_s;
static double *epsilonMap; //読み込まれた誘電率マップ

static double left_s, top_s;

static int trace_ind(int x, int y)
{
  return x*height_px + y;
}

double getEpsilon(double x, double y)
{
  double x_px = min(  width_px-1, max(0, x*pixelPerCell));
  double y_px = min( height_px-1, max(0, y*pixelPerCell));
  if(x_px == width_px-1 || y_px == height_px-1)
    return EPSILON_0_S;
  
  int i = floor(x_px);
  int j = floor(y_px);
  double dx = x_px - i;
  double dy = y_px - j;
  int index = trace_ind(i,j);
  
  return epsilonMap[index]*(1.0-dx)*(1.0-dy)
       + epsilonMap[index+height_px]*dx*(1.0-dy)
       + epsilonMap[index+1]*(1.0-dx)*dy
       + epsilonMap[index+height_px+1]*dx*dy;
}

static double eps(double _x, double _y, int col, int row)
{
  double x =  _x - left_s; //画像と同じ座標系に合わせる
  double y = -_y + top_s;
  
  if( x<-0.5 || x>=width_s+0.5 || y<-0.5 || y>height_s+0.5)
    return EPSILON_0_S;  
  
  double s=0; //n1の分割セルの数が入る  
  double split = 10;
  double half_split = split/2;
  
  for(double i=-half_split+0.5; i<half_split; i+=1)  {
    for(double j=-half_split+0.5; j<half_split; j+=1){
      double sx = x + col*i/split; //細分化したセルの位置
      double sy = y + row*j/split;

      if( sx<0 || sx>=width_s-1 || sy<0 || sy>=height_s-1)
      {
        s += EPSILON_0_S;
      } else {
        s += getEpsilon(sx, sy);
      }
    }
  }
  s = s/split/split;
  if(s < EPSILON_0_S)
  {
    printf("%lf \n",s);
  }
  return s;
}

static void readImage()
{
  FILE *fp = NULL;
  
  if( !(fp = fopen(FILE_NAME".txt", "r")))
  {
    printf("cannot find traceImage.txt of morphoScaleModel\n");
    exit(2);
  }

  //横幅と縦幅, 1ピクセルの大きさを取得
  // nmPerPixel depend on morpho scale width( in this case it's is asumed 300nm)
  fscanf(fp,"%d %d %d",&width_px, &height_px, &scale_width_px);
  nmPerPixel = 1.0 * scale_width_nm / scale_width_px;
  epsilonMap = newDouble(width_px*height_px);
  for(int y=0; y<height_px; y++){
    for(int x=0; x<width_px; x++){
      double n;
      fscanf(fp, "%lf ", &n);
      if(n<1.0)
        printf("%lf \n",n);
      epsilonMap[trace_ind(x,y)] = n*n*EPSILON_0_S;
    }
  }

  printf("%d %d %.2lf", width_px, height_px, nmPerPixel);
}

static bool nextStructure()
{
  scale_width_nm += DELTA_SCALE_WIDTH_NM;
  if(scale_width_nm > EN_SCALE_WIDTH_NM)
    return true;

  nmPerPixel = scale_width_nm / scale_width_px;
  //別の構造は無い
  return false;
}

double ( *traceImageModel_EPS(void))(double, double, int, int)
{
  readImage();
  return eps;
}

bool traceImageModel_isFinish()
{
  return nextStructure();
}

void traceImageModel_moveDirectory()
{
  char buf[512];

  makeAndMoveDirectory(FILE_NAME);
  
  sprintf(buf,"scale_width%d", scale_width_nm);
  makeDirectory(buf);
  moveDirectory(buf); 
}

void traceImageModel_needSize(int *x, int *y)
{
  *x = ceil(width_px*nmPerPixel);
  *y = ceil(height_px*nmPerPixel);
}

void traceImageModel_init()
{
  pixelPerCell = 1.0/field_toCellUnit(nmPerPixel);  
  width_s  = field_toCellUnit(width_px*nmPerPixel);
  height_s = field_toCellUnit(height_px*nmPerPixel);
  
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  left_s = fInfo_s.N_PX/2 - width_s/2;
  top_s  = fInfo_s.N_PY/2 + height_s/2;

  printf("px=(%d,%d), %lf\n",width_px, height_px, nmPerPixel);
  printf("s=(%lf,%lf)\n",width_s, height_s);
  printf("lp=(%lf,%lf)\n",left_s, top_s);
  printf("%d %d\n", fInfo_s.N_X, fInfo_s.N_Y);
}
