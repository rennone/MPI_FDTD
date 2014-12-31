//#include "field.h"
#include "models.h"
#include "noModel.h"
#include "circleModel.h"
#include "concentricCircleModel.h"
#include "multiLayerModel.h"
#include "morphoScaleModel.h"
#include "zigzagModel.h"
#include "traceImageModel.h"
#include "bool.h"
#include "function.h"

static void (*initModelMethod)();
static void (*needSizeMethod)(int*, int*);
static double (*epsMethod)(double, double, int, int);
static bool (*isFinishMethod)(void);
static void (*moveDirectoryMethod)(void);
static char *dir;

static void noModel(void)
{
  //no material
  dir = "NoModel";
  epsMethod = noModel_EPS();
  isFinishMethod = noModel_isFinish;
  needSizeMethod = noModel_needSize;
  initModelMethod     = noModel_init;
  moveDirectoryMethod = noModel_moveDirectory;
}

static void circleModel(void)
{
  dir = "MieCylinderModel";
  epsMethod           = circleModel_EPS();
  isFinishMethod      = circleModel_isFinish;  
  needSizeMethod      = circleModel_needSize;
  initModelMethod     = circleModel_init;
  moveDirectoryMethod = circleModel_moveDirectory;
}

static void multiLayerModel()
{
  dir = "MultiLayerModel";
  epsMethod       = multiLayerModel_EPS();
  isFinishMethod  = multiLayerModel_isFinish;
  needSizeMethod  = multiLayerModel_needSize;
  initModelMethod = multiLayerModel_init;  
  moveDirectoryMethod = multiLayerModel_moveDirectory;
}

static void morphoScaleModel()
{
  dir = "MorphoScaleModel";
  epsMethod = morphoScaleModel_EPS();
  isFinishMethod = morphoScaleModel_isFinish;
  needSizeMethod = morphoScaleModel_needSize;
  initModelMethod = morphoScaleModel_init;  
  moveDirectoryMethod = morphoScaleModel_moveDirectory;
}

static void zigzagModel()
{
  dir = "ZigZagModel";
  epsMethod = zigzagModel_EPS();
  isFinishMethod = zigzagModel_isFinish;
  needSizeMethod = zigzagModel_needSize;
  initModelMethod = zigzagModel_init;  
  moveDirectoryMethod = zigzagModel_moveDirectory;

}

static void traceImageModel()
{
  dir = "TraceImageModel";
  epsMethod = traceImageModel_EPS();
  isFinishMethod = traceImageModel_isFinish;
  needSizeMethod = traceImageModel_needSize;
  initModelMethod = traceImageModel_init;  
  moveDirectoryMethod = traceImageModel_moveDirectory;
}

static void concentricCircleModel()
{
  dir = "ConcentricCircleModel";
  epsMethod = concentricCircleModel_EPS();
  isFinishMethod      = concentricCircleModel_isFinish;
  moveDirectoryMethod = concentricCircleModel_moveDirectory;
  printf("not implemented concentricCircle Model");
  exit(2);
}

bool models_isFinish()
{
  return (*isFinishMethod)();
}

//モデルを変更したときに,一度rootまで戻るので,再度ネストする用
void models_moveDirectory()
{
  makeDirectory(dir);
  moveDirectory(dir);
  (*moveDirectoryMethod)();
}

void models_setModel(enum MODEL model)
{
  switch(model){
  case NO_MODEL:
    noModel();
    break;
  case MIE_CYLINDER:
    circleModel();
    break;
  case LAYER:
    multiLayerModel();
    break;
  case MORPHO_SCALE:
    morphoScaleModel();
    break;
  case CONCENTRIC_CIRCLE:
    concentricCircleModel();
    break;
  case ZIGZAG:
    zigzagModel();
    break;
  case TRACE_IMAGE:
    traceImageModel();
    break;
  }  
}

double models_eps(double x, double y, enum MODE mode)
{
  double epsilon;
  switch(mode){
  case D_X :
    epsilon = (*epsMethod)(x, y, 1, 0);
    break;
  case D_Y :
    epsilon = (*epsMethod)(x, y, 0, 1);
    break;
  case D_XY :
    epsilon = (*epsMethod)(x, y, 1, 1);
    break;
  default:
    epsilon = (*epsMethod)(x, y, 1, 1);
  }
  
  return epsilon;
}

void models_needSize(int *x_nm,int *y_nm)
{
  (*needSizeMethod)(x_nm, y_nm);
}

void models_initModel()
{
  (*initModelMethod)();
}

void models_evaluate(double **reflec, int stLambda, int enLamba)
{
  multiLayerModel_evaluate(reflec, stLambda, enLamba);
}

void models_update()
{
  multiLayerModel_update();
}
