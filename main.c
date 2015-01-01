#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <complex.h>
#include <unistd.h>
#include "simulator.h"
#include "models.h"
#include "field.h"
#include "parser.h"
#include "function.h"
#include "drawer.h"
#include "colorTransform.h"

typedef struct Config
{
  FieldInfo field_info;
  int startAngle, endAngle, deltaAngle;
  enum MODEL ModelType;
  enum SOLVER SolverType;
}Config;

#define ST_PHI 0
#define EN_PHI 0

#define DELTA_PHI 5

// 以下 OPEN_GLの関数
#ifdef USE_OPENGL

#define WINDOW_WIDTH 300
#define WINDOW_HEIGHT 300
#include <GL/glew.h>

//Macの場合
#ifdef MAC_OS
#include <GLUT/glut.h>
#endif

//Mac以外
#ifndef MAC_OS
#include <GL/glut.h>
#endif

//プロトタイプ宣言
static void drawField();
static void drawSubField();
static void display();
static void idle();

#endif

int rank;
int numProc;
Config config;
char root[512]; //ルートディレクトリ


static void moveDir()
{
  moveDirectory(root);
  models_moveDirectory();

  //h_uの大きさでディレクトリを作る=> TimeNTFFではH_Uが重要なため
  char buf[128];
  sprintf(buf,"hu_%dnm",config.field_info.h_u_nm);
  makeDirectory(buf);
  moveDirectory(buf);

  simulator_moveDirectory();  
}

static void calcFieldSize(FieldInfo *fInfo)
{
  int x_nm, y_nm;
  models_needSize(&x_nm, &y_nm);
  
  //モデルのサイズ + (pml + ntff)*2 + 余白
  fInfo->width_nm  = x_nm + fInfo->h_u_nm*(fInfo->pml + 5)*2 + 200;
  fInfo->height_nm = y_nm + fInfo->h_u_nm*(fInfo->pml + 5)*2 + 200;

  bool square = false;
  if(square)
  {
    if(fInfo->width_nm > fInfo->height_nm)
      fInfo->height_nm = fInfo->width_nm;
    else
      fInfo->width_nm = fInfo->height_nm;
  }
}

static void initParameter()
{
  config.field_info.h_u_nm    = 10;
  config.field_info.pml       = 15; //pml領域のセル数
  config.field_info.lambda_nm = 500;
  //領域の全体サイズが変化し,収束にかかる時間が変わるためh_uによりステップ数を変える必要がある.
  config.field_info.stepNum   = 20000 / config.field_info.h_u_nm;
  config.field_info.angle_deg = ST_PHI;
  config.startAngle = ST_PHI;
  config.endAngle   = EN_PHI;
  config.deltaAngle = DELTA_PHI;
  config.SolverType = -1; //Note : 使用してはいけない
  config.ModelType = -1; //Note : 使用してはいけない
  
  calcFieldSize(&config.field_info);

  printf("field size(%d nm, %d nm)\n", config.field_info.width_nm, config.field_info.height_nm);
}

//次のシミュレーションパラメータを探す
//progress : 現在の状態からどれだけ進むかのステップ数
//最初はrankだけ進んで, 2回目以降はnumProcだけ進む
//changeModelはモデルの構造が変化したか => フィールドを書き換える必要があるか
//戻り値 : シミュレーションが終了かどうか
bool nextSimulation(int progress, bool *changeModel)
{
  (*changeModel) = false;
  //プロセス数の分だけ進む
  for(int i=0; i<progress; i++)
  {    
    config.field_info.angle_deg += config.deltaAngle; //入射角度を増やす

    //入射角度が全部終わったら, 構造を変える
    if(config.field_info.angle_deg > config.endAngle)
    {
      config.field_info.angle_deg = config.startAngle;
      (*changeModel) = true;
      if(models_isFinish()) {
        return true;
      }

#ifdef USE_OPENGL
      drawer_clear(); //画像をクリア
#endif
      
    }
  }
  return false;
}

static void screenshot()
{
  //その構造で,角度がST_PHIのやつだけが画像を保存する.
  if(config.field_info.angle_deg == ST_PHI)
  {
    FieldInfo_S fInfo_s = field_getFieldInfo_S();
    drawer_outputImage("image.bmp", simulator_getDrawingData(), simulator_getEps(), fInfo_s.N_PX, fInfo_s.N_PY);
  }
}

int main( int argc, char *argv[] )
{
  getcwd(root, 512); //カレントディレクトリを保存

  MPI_Init( 0, 0 );
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);

  colorTransform_init();

  models_setModel(MIE_CYLINDER);       // MORPHO_SCALE,TRACE_IMAGE, ZIGZAG,NO_MODEL,
  simulator_setSolver(MPI_TM_UPML); 


//  initConfigFromText();
  initParameter();  //パラメータを設定
  moveDir();

  /*
 //プロセスごとに別のシミュレーションをする.
  bool changeModel;
  if( nextSimulation(rank, &changeModel) == true){
    MPI_Finalize();
    exit(0);
  }
  */
  
  //シミュレーションの初期化.
  simulator_init(config.field_info);
  screenshot();
  
  //exitしたプロセスがあると停止してしまうのでMPI_Finalizeは使えない
//  MPI_Barrier(MPI_COMM_WORLD); //(情報表示がずれないように)全員一緒に始める
  printf("rank=%d, angle=%d\n",rank, config.field_info.angle_deg);
  
#ifndef USE_OPENGL
  //only calculate mode
  while(1)
  {
    //シミュレーションをまわす
    while(!simulator_isFinish()) {
      simulator_calc();    
    }

    //numProc次のシミュレーションにする. 角度が変わるか構造(nm変数のみ)が変わるか
    if( nextSimulation(numProc, &changeModel) == true)
    {
      simulator_finish(); //シミュレーションの終わり
      break;
    }

    //モデルが変化したかどうか
    if(changeModel)
    {
      simulator_finish();     //変化したら, メモリを解放させる為にシミュレーションを終了させる
      calcFieldSize(&config.field_info);        //フィールドサイズの再計算
      moveDir(); //ディレクトリの移動
      simulator_init(config.field_info);
     
      screenshot();
  
    } else {
      simulator_reset(); //変化してなければ,データの書き出しと電磁波の値だけ0に戻す.
      field_setWaveAngle(config.field_info.angle_deg); //角度だけ変える.
    }
  }

  MPI_Finalize(); //プロセスごとにFinalizeしてもok
#endif

  
#ifdef USE_OPENGL
  SubFieldInfo_S subInfo = field_getSubFieldInfo_S();
  
  int windowX = WINDOW_WIDTH*(rank%6);
  int windowY = (WINDOW_HEIGHT+50)*(rank/6);
  enum COLOR_MODE colorMode = CABS;
  glutInit(&argc, argv);
  glutInitWindowPosition(windowX,windowY);
  glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
  glutCreateWindow("FDTD Simulator");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  glewInit();
  drawer_init(colorMode);
  glutMainLoop();
#endif

  return 1;
}



// OpenGL関連の実装部分
#ifdef USE_OPENGL

static void drawField()
{
  FieldInfo_S sInfo = field_getFieldInfo_S();
  drawer_paintImage(0, 0, sInfo.N_PX, sInfo.N_PY, sInfo.N_PX, sInfo.N_PY,
                    simulator_getDrawingData());
  
  drawer_paintModel(0, 0, sInfo.N_PX, sInfo.N_PY, sInfo.N_PX, sInfo.N_PY,
                    simulator_getEps());
}

static void drawSubField()
{
  SubFieldInfo_S subInfo = field_getSubFieldInfo_S();
  drawer_paintImage(1,1, subInfo.SUB_N_X, subInfo.SUB_N_Y, subInfo.SUB_N_PX, subInfo.SUB_N_PY,
                    simulator_getDrawingData());
  drawer_paintModel(1,1, subInfo.SUB_N_X, subInfo.SUB_N_Y, subInfo.SUB_N_PX, subInfo.SUB_N_PY,
                    simulator_getEps());  
}

static void display()
{
  glClear(GL_COLOR_BUFFER_BIT);
  glEnableClientState( GL_VERTEX_ARRAY );
  glEnableClientState( GL_TEXTURE_COORD_ARRAY );

  //drawField(); //分割しない全体領域の描画
  drawSubField();  //分割された領域の描画  
  drawer_draw();

  glDisableClientState( GL_VERTEX_ARRAY );
  glDisableClientState( GL_TEXTURE_COORD_ARRAY );
  glutSwapBuffers();
}

static void idle(void)
{
  simulator_calc();

  //Note:
  //GAでルートプロセスが,他のプロセスの終了を確認するために入れている.
  models_update(); //モデルのアップデート処理
  
  //シミュレーションが続くなら再描画して終わり
  if( !simulator_isFinish() ){
    glutPostRedisplay();  //再描画
    return;
  }
  
  bool changeModel;

  //シミュレーションを進める.
  if( nextSimulation(numProc, &changeModel) == true )
  {
    simulator_finish();
    MPI_Finalize();
    exit(0);
  }

  if( changeModel ){
    simulator_finish(); //シミュレーションを終える.
    
    //モデルを変更して再計算する
    calcFieldSize(&config.field_info);   
    moveDir();
    //シミュレーションの初期化.
    simulator_init(config.field_info);

    screenshot();
  } else {
    simulator_reset();
    field_setWaveAngle(config.field_info.angle_deg);
  }

}
// 以上 OPENGLの関数

/*
  //ファイルからのパラメータ呼び出し
static void readConfig(FieldInfo *field_info)
{
  FILE *fp = NULL;
  if( !(fp = fopen("config.txt", "r")) )
  {
    printf("cannot find config.txt of main\n");
  }

  char buf[1024],**tmp;

  parser_nextLine(fp, buf);
  field_info->width_nm    = atoi(buf);

  parser_nextLine(fp, buf);
  field_info->height_nm   = atoi(buf);

  parser_nextLine(fp, buf);
  field_info->h_u_nm      = atoi(buf);

  parser_nextLine(fp, buf);
  field_info->pml         = atoi(buf);

  parser_nextLine(fp, buf);
  field_info->lambda_nm   = strtod(buf,tmp);

  parser_nextLine(fp, buf);
  field_info->stepNum     = atoi(buf);

//入射角度
  parser_nextLine(fp, buf);
  config.startAngle = atoi(buf);

  parser_nextLine(fp, buf);
  config.endAngle = atoi(buf);

  parser_nextLine(fp, buf);
  config.deltaAngle = atoi(buf);

  //モデルの種類
  parser_nextLine(fp, buf);
  config.ModelType = atoi(buf);

  //Solver情報
  parser_nextLine(fp, buf);
  config.SolverType = atoi(buf);

  fclose(fp);
}

static void initConfigFromText()
{
  //同時にファイルを読み込むのはヤバい気がするので, rank0だけが読み込んで同期する
  if(rank == 0)
  {    
    readConfig(&config.field_info);
    printf("===========FieldSetting=======\n");
    printf("fieldSize(nm) = (%d, %d) \nh_u = %d \npml = %d\n", config.field_info.width_nm, config.field_info.height_nm,
           config.field_info.h_u_nm, config.field_info.pml);
    printf("lambda(nm) = %d  \nstep = %d\n", config.field_info.lambda_nm, config.field_info.stepNum);

    printf("angle = %d .. %d (delta = %d)\n", config.startAngle, config.endAngle, config.deltaAngle);
    printf("==============================\n");
  }
    
  //configを同期
  if(rank == 0)
  {
    for(int i=1; i<numProc; i++)
    {
      MPI_Send((int*)&config, sizeof(Config)/sizeof(int), MPI_INT, i, 1, MPI_COMM_WORLD);
    }
  }else{
    MPI_Status status;
    MPI_Recv((int*)&config, sizeof(Config)/sizeof(int), MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
  }
}
*/
#endif
