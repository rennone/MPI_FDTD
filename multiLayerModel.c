#include "multiLayerModel.h"
#include "field.h"
#include "function.h"
#include <math.h>
#include "evaluate.h"
#include <mpi.h>
//屈折率
#define N_0 1.0
#define N_1 1.56
//serikon
//#define N_1 8.4179

//異方性を入れるかのフラグ
#define UNIAXIAL false
#define N_0_X 1.0
#define N_1_X 1.1

//横幅
#define ST_WIDTH_NM 300
#define EN_WIDTH_NM 300
#define DELTA_WIDTH_NM 10

//ラメラ1の厚さ
#define ST_THICK_NM_1 90
#define EN_THICK_NM_1 150
#define DELTA_THICK_NM_1 30

//ラメラ0の厚さ
#define ST_THICK_NM_0 (ST_THICK_NM_1 +  0)
#define EN_THICK_NM_0 (ST_THICK_NM_1 + 50)
#define DELTA_THICK_NM_0 10

//ラメラの枚数
#define ST_LAYER_NUM 4
#define EN_LAYER_NUM 4
#define DELTA_LAYER_NUM 2

//互い違い => 左右で n_0 n_1を入れ替え
#define ASYMMETRY false

//USE_GAP flag 
//左右でずらす => DELTA_LEFT_GAP_Y ~ thickness_nm まで変化する.
#define USE_GAP false
#define DELTA_LEFT_GAP_Y 0

//中心に以下の幅で軸となる枝を入れる => 軸の屈折率はN_1になる
#define ST_BRANCH_NM 0
#define EN_BRANCH_NM 0
#define DELTA_BRANCH_NM 50

//先端における横幅の割合
#define ST_EDGE_RATE 0.0
#define EN_EDGE_RATE 1.0
#define DELTA_EDGE_RATE 1.0

//ラメラの先端を丸める曲率 (0で四角形のまま, 1.0で最もカーブする)
#define CURVE 0.0

static int width_nm[2]     = {ST_WIDTH_NM, ST_WIDTH_NM};
static double width_s[2];     //幅

static int thickness_nm[2] = {ST_THICK_NM_0, ST_THICK_NM_1};
static double thickness_s[2]; //厚さ

static int branch_width_nm = ST_BRANCH_NM; //枝の幅
static double branch_width_s; //枝の幅

static int layerNum = ST_LAYER_NUM;     //枚数

static double curve_rate = CURVE;

//N_0 N_1から計算
static double ep_s[2];        //誘電率 = n*n*ep0
static double ep_x_s[2];      //異方性用のx方向の誘電率

static double edge_width_rate = ST_EDGE_RATE;

//DELTA_LEFT_GAP_Y
#if USE_GAP
static int left_gap_y_nm = DELTA_LEFT_GAP_Y;
#else
static int left_gap_y_nm = 0;
#endif

static double left_gap_y_s;

//CURVE から計算
static double c0, c1; //2次関数の比例定数

static void GAInitialize();
  
static double calc_width(double sx, double sy, double wid, double hei, double modY, int k)
{
  double p = 1 - sy/hei;
  double new_wid = (wid+branch_width_s)*(p + (1-p)*edge_width_rate);

  //ラメラの下を基準とした位置を求める
  double dh = k==0 ? modY : modY - thickness_s[0];
  double c  = k==0 ? c0 : c1;

//互い違いの場合はdhを再計算
  if(ASYMMETRY && sx < 0){
    dh = (k==1 ? modY : modY - thickness_s[1]);
  }

  //2次関数で横幅を計算
  return c*pow((dh-thickness_s[k]/2),2) + new_wid;
}

//col : D_Xモード row : D_Yモード
//x,yを中心に, 計算領域のセルと同じ大きさの領域を調べる
static double eps(double x, double y, int col, int row)
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();

  double width = max(width_s[0], width_s[1]);
  double thick = thickness_s[0] + thickness_s[1];
  double height = thick*layerNum + left_gap_y_s; //オフセットがあれば行う

  //領域の中心から, 下にheight/2ずれた位置がレイヤの下部
  int oy = fInfo_s.N_PY/2 - height/2;
  int ox = fInfo_s.N_PX/2;
  double _x = x-ox;	//ox,oyを座標の原点に
  double _y = y-oy;

  //上下左右に飛び出ていないか確認(細分化したセルがあるため, 0.5の余白をとっている)
  if( fabs(_x) > (width/2.0+0.5) ||  _y < -0.5 || _y > height+0.5 )  
    return EPSILON_0_S;

  height = thick*layerNum; //元にもどす.
  double s[2]={0,0}; //n1,n2それぞれの分割セルの数が入る
  double split = 10;
  double half_split = split/2;
  for(double i=-half_split+0.5; i<half_split; i+=1){
    for(double j=-half_split+0.5; j<half_split; j+=1){
      double sx = _x + col*i/split; //細分化したセルの位置
      double sy = _y + row*j/split;

      //左側はleft_gapだけずらす.
      if(sx < 0 && USE_GAP)
        sy -= left_gap_y_s;
      
      //上下に飛び出ていないか確認
      if(sy < 0 || sy > height)
        continue;

      double p = 1 - sy/height;
      //枝の部分
      if(fabs(sx) < branch_width_s*(p + (1-p)*edge_width_rate))
      {
        s[1] += 1;
        continue;
      }
      
      //thickで割ったあまり(double型なのでこんなやり方をしている)
      double modY = sy - floor(sy/thick)*thick;

      //境界上のときは両方の平均になる(普通は無い).
      if( modY == thickness_s[0]) {
        s[0] += 0.5*(fabs(sx) < width_s[0]/2);
        s[1] += 0.5*(fabs(sx) < width_s[1]/2);
        continue;
      }

      //どっちの屈折率にいるか調べる
      int k;
      if (sx < 0 && ASYMMETRY) {
        k = (modY < thickness_s[1]); //互い違いかつ左側は, 1が下にある
      } else {
        k = (modY > thickness_s[0]); //それ以外は0が下にある
      }
      
      double wid = calc_width(sx, sy, width_s[k], height, modY, k);
      
      if(fabs(sx) < wid/2)
        s[k] +=1;
    }    
  }

  s[0] /= split*split;
  s[1] /= split*split;
  if(UNIAXIAL && col == 0 && row == 1)
  {
    return EPSILON_0_S*(1-s[0]-s[1]) + ep_x_s[0]*s[0] + ep_x_s[1]*s[1];
  }
  else
  {
    return EPSILON_0_S*(1-s[0]-s[1]) + ep_s[0]*s[0] + ep_s[1]*s[1];
  }
}

double ( *multiLayerModel_EPS(void))(double, double, int, int)
{
  GAInitialize();
  return eps;
}


//構造を一つ進める
static bool nextStructure1()
{
  left_gap_y_nm += DELTA_LEFT_GAP_Y;
  if( left_gap_y_nm >= (thickness_nm[0]+thickness_nm[1]) || !USE_GAP)
  {
    left_gap_y_nm = USE_GAP ? DELTA_LEFT_GAP_Y : 0;
    thickness_nm[0] += DELTA_THICK_NM_0;
    thickness_nm[1] += DELTA_THICK_NM_1;

    if(thickness_nm[0] > EN_THICK_NM_0)
    {
      thickness_nm[0] = ST_THICK_NM_0;
      thickness_nm[1] = ST_THICK_NM_1;

      edge_width_rate += DELTA_EDGE_RATE;
      if(edge_width_rate > EN_EDGE_RATE)
      {
        edge_width_rate = ST_EDGE_RATE;
        branch_width_nm += DELTA_BRANCH_NM;
        if(branch_width_nm > EN_BRANCH_NM)
      	{
	  branch_width_nm = ST_BRANCH_NM;
	  layerNum += DELTA_LAYER_NUM;
	  if( layerNum > EN_LAYER_NUM)
          {
            printf("there are no models which hasn't been simulated yet\n");
            return true;
          }
	}
      }
    }
  }
  return false;  
}

//構造を一つ進める
static bool nextStructure2()
{
  left_gap_y_nm += DELTA_LEFT_GAP_Y;
  if( left_gap_y_nm >= (thickness_nm[0]+thickness_nm[1]) || !USE_GAP){
    left_gap_y_nm = USE_GAP ? DELTA_LEFT_GAP_Y : 0;
    thickness_nm[0] += DELTA_THICK_NM_0;

    if(thickness_nm[0] > /*EN_THICK_NM_0*/ thickness_nm[1]+50){
      thickness_nm[1] += DELTA_THICK_NM_1;
      thickness_nm[0] = thickness_nm[1];////ST_THICK_NM_0;
    
      if(thickness_nm[1] > EN_THICK_NM_1){ 
	thickness_nm[1] = ST_THICK_NM_1;
	edge_width_rate += DELTA_EDGE_RATE;
      
	if(edge_width_rate > EN_EDGE_RATE) {
	  edge_width_rate = ST_EDGE_RATE;
	  branch_width_nm += DELTA_BRANCH_NM;
	
	  if(branch_width_nm > EN_BRANCH_NM) {
	    branch_width_nm = ST_BRANCH_NM;
	    layerNum += DELTA_LAYER_NUM;
	  
	    if( layerNum > EN_LAYER_NUM) {
	      printf("there are no models which hasn't been simulated yet\n");
	      return true;
	    }
	  }
	}
      }
    }
  }
  return false;  
}

/*
bool multiLayerModel_isFinish(void)
{
  return nextStructure2();
}
*/
void multiLayerModel_needSize(int *x_nm, int *y_nm)
{
  (*x_nm) = max( width_nm[0], width_nm[1]) + branch_width_nm;

  //最後の項はgapの分(これは固定にしないと,gapによりフィールドの領域が変わるので図が変に見える).
  (*y_nm) = (thickness_nm[0]+thickness_nm[1])*layerNum + (USE_GAP ? thickness_nm[0]+thickness_nm[0] : 0) ;
}

static char  root[512];
void multiLayerModel_moveDirectory()
{
  getcwd(root, 512);
  makeDirectory("GA_Images");
  sprintf(root, "%s/GA_Images",root);
  printf("root == %s\n",root);
  if(ASYMMETRY){
    makeDirectory("asymmetry");
    moveDirectory("asymmetry");
  } else {
    makeDirectory("symmetry");
    moveDirectory("symmetry");
  }  
  char buf[512];
  // make folder by index of reflaction

  if(UNIAXIAL){
    sprintf(buf,"uniaxial_n0y%.2lf_n0x%.2lf_n1x%.2lf_n1x%.2lf", N_0, N_0_X, N_1, N_1_X);
    makeDirectory(buf);
    moveDirectory(buf);
  }
  else
  {
    sprintf(buf,"n_%.2lf_%.2lf", N_0, N_1);
    makeDirectory(buf);
    moveDirectory(buf);
  }
  
  sprintf(buf,"curve_%.2lf", curve_rate);
  makeAndMoveDirectory(buf);

  sprintf(buf, "width%d_%d", width_nm[0], width_nm[0]);
  makeAndMoveDirectory(buf);

  sprintf(buf, "thick%d_%d", thickness_nm[0], thickness_nm[1]);
  makeAndMoveDirectory(buf);
  
  sprintf(buf, "gap%d", left_gap_y_nm);
  makeAndMoveDirectory(buf);

  sprintf(buf, "layer%d", layerNum);
  makeAndMoveDirectory(buf);

  sprintf(buf, "edge%.2lf", edge_width_rate);
  makeAndMoveDirectory(buf);

  sprintf(buf, "branch%d", branch_width_nm);
  makeAndMoveDirectory(buf);

  getcwd(buf, 512);
  printf("move to %s\n", buf);
}

void multiLayerModel_init()
{
  width_s[0]     = field_toCellUnit(width_nm[0]);
  width_s[1]     = field_toCellUnit(width_nm[1]);
  thickness_s[0] = field_toCellUnit(thickness_nm[0]);
  thickness_s[1] = field_toCellUnit(thickness_nm[1]);
  left_gap_y_s   = field_toCellUnit(left_gap_y_nm);
  
  ep_s[0] = N_0*N_0*EPSILON_0_S;
  ep_s[1] = N_1*N_1*EPSILON_0_S;
  ep_x_s[0] = N_0_X*N_0_X*EPSILON_0_S;
  ep_x_s[1] = N_1_X*N_1_X*EPSILON_0_S;
  
  branch_width_s = field_toCellUnit(branch_width_nm);
  
  c0 = -4*width_s[0]*curve_rate/thickness_s[0]/thickness_s[0];
  c1 = -4*width_s[1]*curve_rate/thickness_s[1]/thickness_s[1];
}

typedef enum Kinds
{
  eTHICK_NM_0,
  eTHICK_NM_1,
  eLAYER_NUM,
  eEDGE,    //先端の縮小率%
  eCURVE,   //ラメラのカーブの割合%
  eBRANCH_NM,
  eKIND_NUM
} Kinds;
  
// 0~4bit  : ラメラ0の幅    ( 10nm ~ 320nm )
// 5~9bit  : ラメラ1の幅    ( 10nm ~ 320nm )
//10~13bit : ラメラの枚数   ( 1 ~ 16毎 )
//14~17bit : 先端の縮小率  edge   ( 0, 1/15, 2/15 ~ 15/15)
//18~21bit : ラメラの丸め率 curve ( 0, 1/15, ~  15/15)
//22~26bit : 幹の太さ( ラメラ1の幅の 0, 1/15, ~ 15/15)
//個体
typedef struct Individual
{
  double evals[EVAL_NUM]; // RGBそれぞれに対する. 評価値(適応度)
  // double eval;       //評価値(適応度)
  int cells[eKIND_NUM]; //個体の状態
} Individual;
static void printIndiv(Individual *p);

//2つが同じ個体か調べる
static bool Equal(Individual *a, Individual *b)
{  
  for(int i=0; i<eKIND_NUM; i++)
    if( a->cells[i] != b->cells[i] )
      return false;

  return true;
}

//現在の設定をIndividual化する
Individual SettingToInidividual(double evals[EVAL_NUM])
{
  Individual p;
  p.cells[eTHICK_NM_0] = thickness_nm[0];
  p.cells[eTHICK_NM_1] = thickness_nm[1];
  p.cells[eLAYER_NUM] = layerNum;
  p.cells[eEDGE]      = (int)(edge_width_rate * 100);
  p.cells[eCURVE]     = (int)(curve_rate * 100);
  p.cells[eBRANCH_NM] = branch_width_nm;

//  p.eval = value;
  for(int i=0; i<EVAL_NUM; i++)
    p.evals[i] = evals[i];
  
  return p;
}

//個体を現在の設定に反映する.
void IndividualToSetting(Individual *p)
{
  thickness_nm[0] = p->cells[eTHICK_NM_0];
  thickness_nm[1] = p->cells[eTHICK_NM_1];
  layerNum        = p->cells[eLAYER_NUM];
  edge_width_rate = 1.0*p->cells[eEDGE] / 100.0;
  curve_rate      = 1.0*p->cells[eCURVE] / 100.0;
  branch_width_nm = p->cells[eBRANCH_NM];
}

static int rank = -1;
static int numProc = -1;
static int N_Param = -1;
static Individual *curGeneration = NULL; //現世代
static Individual *nexGeneration = NULL; //次世代
static int indivNoCur = 0; //現世代で実行済みの個体数
static MPI_Datatype MPI_INDIVIDUAL;

#define MEMO_STRAGE 100000
static Individual *Memo;//[MEMO_STRAGE];
static int numOfMemo = 0;

//すでに存在するか
static bool Exist(Individual *p)
{
  for(int i=0; i<numOfMemo; i++)
    if( Equal(p, &Memo[i]) )
      return true;
  
  return false;
}


 //計算した個体は保存しておく
static void AddMemo(Individual p)
{
   Memo[numOfMemo++] = p;
}

static void BuildDerivedType()
{
  MPI_Datatype typelists[2];
  typelists[0] = MPI_DOUBLE;
  typelists[1] = MPI_INT;

  int block_length[2];
  block_length[0] = EVAL_NUM;  //1番目の型は配列でEVAL_NUM個
  block_length[1] = eKIND_NUM; //2番目の型は配列でeKIND_NUM個

  MPI_Aint displacements[2];
  MPI_Aint start_address;
  MPI_Aint address;

  Individual p; //アドレスの差分取得用
  MPI_Address(&(p.evals[0]), &start_address);
  displacements[0] = 0;
    
  MPI_Address(&(p.cells[0]), &address);
  displacements[1] = address - start_address;

  MPI_Type_struct(2, block_length, displacements, typelists, &MPI_INDIVIDUAL);
  MPI_Type_commit(&MPI_INDIVIDUAL);
}

static void printIndiv(Individual *p)
{
  for(int i=0; i<EVAL_NUM; i++)
    printf("v = %lf\n", p->evals[i]);
  
  for(int i=0; i<eKIND_NUM; i++)
    printf("%d\n", p->cells[i]);
  printf("\n");
}


//
typedef enum Tags
{
  tInitIndiv, //初期の個体通信用
  tEval,      //評価通信用
  tNextIndiv, //次の個体通信用
  tImageNo,  //画像番号
} Tags;

static void BuildTypeTest()
{
  //テスト
  Individual p;
  for(int i=0; i<EVAL_NUM; i++)
    p.evals[i] = rank;
  
  for(int i=0; i<eKIND_NUM; i++)
    p.cells[i] = (i+1)*rank;
  
  if(rank == 0){
    for(int i=1; i<numProc; i++)
    {
      MPI_Status status;
      Individual s;
      MPI_Recv(&s, 1, MPI_INDIVIDUAL, i, 0, MPI_COMM_WORLD, &status);
      printIndiv(&s);
    }
  }
  else{
    MPI_Send(&p, 1, MPI_INDIVIDUAL, 0, 0, MPI_COMM_WORLD);
  }
}

//突然変異
static Individual Mutation()
{
  Individual p;
  p.cells[eTHICK_NM_0] = (rand() % 50) * 10 + 10; // 10 ~ 500nm
  p.cells[eTHICK_NM_1] = (rand() % 50) * 10 + 50; // 50 ~ 550nm
  p.cells[eLAYER_NUM]  = rand() % 10 + 1;     // 1 ~ 10
  p.cells[eEDGE]       = (rand() % 21)*5;     // 0 ~ 100%
  p.cells[eCURVE]      = (rand() % 21)*5;     // 0 ~ 100%

  //ブランチの太さ
  p.cells[eBRANCH_NM]  = 10*(rand() % (ST_WIDTH_NM / 40));

  for(int i=0; i<EVAL_NUM; i++)
    p.evals[i] = -1000;

  if( p.cells[eLAYER_NUM] <= 0)
  {
    printf("LayerNum == 0 at Mutation\n");
    MPI_Finalize();
    exit(2);
  }

  return p;
}

//0~Mod-1の乱数を２つ被らずに取得
static void random2(int Mod, int *a, int *b)
{
  *a = rand() % Mod;
  *b = rand() % Mod;
  while( *b == *a)   *b = rand() % Mod;  
}

// 交叉

static void CrossOver(Individual *p1, Individual *p2)
{
  int p2Cross       = 50; //二点交叉の確率
  int pUniformCorss = 50; //一様交叉の確率
  int p = rand()%100;
  //2点交叉
  if(p < p2Cross){
    int a, b;
    random2(eKIND_NUM, &a, &b);

    int st = min(a,b);
    int en = max(a,b);
    for(int i=st; i<=en; i++){
      int tmp = p1->cells[i];
      p1->cells[i] = p2->cells[i];
      p2->cells[i] = tmp;
    }
  }
  //一様交叉
  else
  {
    for(int i=0; i<eKIND_NUM; i++)
    {
      if(rand()%2){
        int tmp = p1->cells[i];
        p1->cells[i] = p2->cells[i];
        p2->cells[i] = tmp;
      }
    }
  }
  
  if(p1->cells[eLAYER_NUM] == 0 || p2->cells[eLAYER_NUM] == 0)
  {
    printf("LayerNum == 0 at Cross\n");
    MPI_Finalize();
    exit(2);
  }
}

//評価値の中で最大のものを評価とする
static double maxEvals(Individual *p)
{  
  double mx = -1000;
  for(int i=0; i<EVAL_NUM; i++)
  {
    mx = max( mx, p->evals[i] );
  }

  //念のため, pの中身が 0~100の間にあるかチェック      
  return mx;
}

//選択
static Individual Select()
{
  Individual res = Memo[0];
  for(int i=1; i<numOfMemo; i++)
  {
    if( maxEvals(&Memo[i])  > maxEvals(&res) )
      res = Memo[i];
  }

  if(res.cells[eLAYER_NUM] <= 0)
  {
    printf("LayerNum == 0 at Select\n");
    MPI_Finalize();
    exit(2);
  }
  
  return res;
}

// 世代交代
// シミュレーションが終わった最新のnumProc個の個体を用いて世代交代する
// 
static void Heterogenesis()
{
  printf("Heterogenesis\n");
  //現世代を取得
  //シミュレーションが終わった最新のN_Param個を現世代とする
  Individual current[N_Param];
  for(int i=0; i<N_Param; i++)
  {
    if( numOfMemo - i <= 0 )
      current[i] = Mutation();
    else
      current[i] = Memo[numOfMemo - 1 - i];
  }
 
  int selectP   = 20; //選択確率
  int crossP    = 65; //交叉確率
  int mutationP = 15; //突然変異確率
  for(int i=0; i<N_Param;)
  {
    int p = rand()%100;
    if( p < selectP)
    {
      nexGeneration[i++] = Select();      
    }
    else if(p < selectP + crossP)
    {
      int a, b;
      random2(N_Param, &a, &b);
      Individual p1 = current[a];
      Individual p2 = current[b];
      CrossOver(&p1, &p2);
      nexGeneration[i++] = p1;
      nexGeneration[i++] = p2;
    }
    else
    {
      nexGeneration[i++] = Mutation();
    }
  }
}

//次の個体を生成
static Individual NextGeneration()
{
  //現世代すべて計算済みなら世代交代する.
  if( indivNoCur >= N_Param ){
    Heterogenesis();
    indivNoCur = 0;
  }
  
  Individual next = nexGeneration[indivNoCur++];
  return next;
}

//初期の設定はランクによりランダムに決定する.
static Individual InitialGeneration(int rank)
{
  double evals[EVAL_NUM] = {-1000};
  Individual p = SettingToInidividual(evals);
  return p;
}

static MPI_Request *requests, *flgRequests;
static bool *finishedFlags;
#include <time.h>
static void GAInitialize()
{
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0)
  {
    //プロセス数だけ,個体数を持つようにする
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);

    Memo = (Individual*)malloc(sizeof(Individual)*MEMO_STRAGE);
    N_Param = numProc;
    curGeneration = (Individual*)malloc(sizeof(Individual)*N_Param);
    nexGeneration = (Individual*)malloc(sizeof(Individual)*N_Param);
    requests      = (MPI_Request*)malloc(sizeof(MPI_Request) * numProc);
    finishedFlags = (bool*)malloc(sizeof(bool) * numProc);
    flgRequests   = (bool*)malloc(sizeof(bool) * numProc);
  }

  BuildDerivedType(); //構造体を送受信するための方を定義

  srand( (unsigned)time( NULL ) );

  //最初はランダムに生成する.
  if(rank == 0)
  {
    //最初は全部突然変異
    for(int i=0; i<N_Param; i++)
      nexGeneration[i] = Mutation();
    
    //初期の設定を入れる
    for(int i=0; i<numProc; i++)
      curGeneration[i] = NextGeneration();
    
    //自分を設定
    IndividualToSetting(&curGeneration[0]);

    //他のプロセスのものを生成して送る.
    for(int i=1; i<numProc; i++) {
      MPI_Send(&curGeneration[i], 1, MPI_INDIVIDUAL, i, tInitIndiv, MPI_COMM_WORLD);
    }
  }
  else
  {
    MPI_Status status;
    Individual next;
    MPI_Recv(&next, 1, MPI_INDIVIDUAL, 0, tInitIndiv, MPI_COMM_WORLD, &status);  //個体を受け取る.
    IndividualToSetting(&next); //設定に反映
  }
  
  //ルートプロセスは評価を受け入れるために, 受信態勢に入る
  if( rank == 0 )
  {
    for(int i=1; i<numProc; i++)
    {
      MPI_Irecv(&(curGeneration[i]),
                1, MPI_INDIVIDUAL, i, tEval, MPI_COMM_WORLD, &requests[i]); 
    }
  }
}

#include "drawer.h"
#include "simulator.h"
static void save( int imgNo, Individual *p, double r[360], double g[360], double b[360])
{
  double limit = 13; //評価値13以上がないと保存しない(無駄なので)
  bool canSave = false;
  for(int i=0; i<EVAL_NUM; i++)
    canSave |= p->evals[EVAL_RED] > 13;

  if( !canSave )
    return;

  static colorf **img = NULL;
  if( img == NULL ){
    img = (colorf**)malloc(sizeof(colorf*)*180);
    for(int i=0; i<180; i++)
      img[i] = (colorf*)malloc(sizeof(colorf)*64);
  }
  
  for(int i=0; i<180; i++){
    for(int j=0; j<64; j++){
      img[i][j].r = r[i];
      img[i][j].g = g[i];
      img[i][j].b = b[i];
    }
  }    
  
  //画像を保存
  //評価値はファイル名には記述しない.
  //.txtに書いているので, 評価値を使いたいときは読み込むプログラムを書けば良い
  char buf[512];

  //色の画像を保存
  sprintf(buf, "%s/%d.bmp", root, imgNo );
  drawer_saveImage(buf, img, 180, 64);

  //領域の画像を保存
  sprintf(buf, "%s/%d_img.bmp", root, imgNo );
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  drawer_outputImage(buf, simulator_getDrawingData(), simulator_getEps(), fInfo_s.N_PX, fInfo_s.N_PY);

  //領域の設定を保存
  sprintf(buf, "%s/%d.txt", root, imgNo);  
  
  FILE *fp = FileOpen(buf, "w");
  fprintf(fp, "evalValueRed = %.18lf\n", p->evals[EVAL_RED]);
  fprintf(fp, "evalValueGreen = %.18lf\n", p->evals[EVAL_GREEN]);
  fprintf(fp, "evalValueBlue = %.18lf\n", p->evals[EVAL_BLUE]);
  fprintf(fp, "thick_nm_0 = %d\n", p->cells[eTHICK_NM_0]);
  fprintf(fp, "thick_nm_1 = %d\n", p->cells[eTHICK_NM_1]);
  fprintf(fp, "layer_num = %d\n", p->cells[eLAYER_NUM]);
  fprintf(fp, "edge = %d\n", p->cells[eEDGE]);
  fprintf(fp, "curve = %d\n", p->cells[eCURVE]);
  fprintf(fp, "branch_nm = %d\n", p->cells[eBRANCH_NM]);
  fclose(fp);
}

static MPI_Request req, reqImg;
//反射率を用いた評価関数
void multiLayerModel_evaluate(double **reflec, int stLambda, int enLambda)
{
  //色に変換する.
  double red[360]   = {};
  double green[360] = {};
  double blue[360]  = {};

  for(int deg=0; deg<360; deg++)
  {
    for(int i = 0; i <= enLambda - stLambda; i++){

      double r=0,g=0,b=0;
      colorTransform_trans(i+stLambda, reflec[i][deg], &r, &g, &b);
      red[deg]   += r;
      green[deg] += g;
      blue[deg]  += b;
    }
  }

  for(int i=0; i<360; i++){
    red[i] = min(1.0, max(0.0, red[i]));
    green[i] = min(1.0, max(0.0, green[i]));
    blue[i] = min(1.0, max(0.0, blue[i]));
  }
  //評価を求める.
  const int eval_deg = 40; //左右eval_degの視野角で評価
  double evals[EVAL_NUM] = {};
  for(int deg=90-eval_deg; deg<=90+eval_deg; deg++)
  {
    double h=0,s=0,v=0;
    bool res = colorTransform_rgbTohsv(red[deg],green[deg],blue[deg], &h, &s, &v);
    if( !res )
      continue;
    evals[EVAL_BLUE] += evaluate_eval(h, s, v, EVAL_BLUE);
    evals[EVAL_RED]  += evaluate_eval(h, s, v, EVAL_RED);
    evals[EVAL_GREEN]+= evaluate_eval(h, s, v, EVAL_GREEN);
  }

  // eval_degの範囲で最大値が変わらないように正規化
  for(int evalKinds=0; evalKinds<EVAL_NUM; evalKinds++)
    evals[evalKinds] /= (2.0*eval_deg+1.0);

  if(rank == 0)
  {
    for(int i=0; i<EVAL_NUM; i++)
      curGeneration[0].evals[i] = evals[i];
    
    //画像を保存
    int imgNo = numOfMemo;
    save(imgNo, &curGeneration[0], red, green, blue);
    AddMemo(curGeneration[0]);
  }
  else{
    Individual p = SettingToInidividual(evals);
    MPI_Isend(&p, 1, MPI_INDIVIDUAL, 0, tEval, MPI_COMM_WORLD, &req);
    printf("-----------------sended\n");
  }

  // ルートランクは自分で次を計算する.
  if( rank == 0)
  {
    curGeneration[0] = NextGeneration();
    IndividualToSetting(&curGeneration[0]);
  }
  
  // ルートランク以外はルートから受け取るようにする.
  // ルートプロセスは毎フレーム,他のプロセスが終わったか確認しに行き
  // 終わっていたら次を計算し送る(update関数中から送信される.)
  else
  {
    //同期通信で受け取る.
    MPI_Status status;
    Individual next;
    
    printf("rank %d recieving--------\n", rank);    
    int imgNo;
    MPI_Recv(&imgNo, 1, MPI_INT, 0, tImageNo, MPI_COMM_WORLD, &status);          //画像の番号を受け取る.    
    MPI_Recv(&next, 1, MPI_INDIVIDUAL, 0, tNextIndiv, MPI_COMM_WORLD, &status);

    Individual cur = SettingToInidividual(evals);
    save(imgNo, &cur, red, green, blue);

    IndividualToSetting(&next); //設定に反映
    printIndiv(&next);
    printf("-----------------recieved\n");
  }
}

bool multiLayerModel_isFinish(void)
{
  
  return numOfMemo >= MEMO_STRAGE;
}

//polling用
void multiLayerModel_update()
{
  if( rank != 0 )
    return;

  for(int i=1; i<numProc; i++)
  {
    MPI_Status status;
    int flag;
    //通信が終了したか確認する.
    MPI_Test(&requests[i], &flag, &status);
    // 他のプロセスから個体を受け取った場合 => 次の個体を送信する.
    if( flag )
    {
      int No = numOfMemo; //現在の番号
      MPI_Isend(&No,1, MPI_INT, i, tImageNo, MPI_COMM_WORLD, &reqImg);
      AddMemo(curGeneration[i]);
      
      Individual p = NextGeneration(); //次の個体を計算
      curGeneration[i] = p;            //コピーする
      
      printf("root ------------- nextGeneration\n");
      printIndiv(&p);
      printf(" ------------- \n");
      
      //非同期通信で次の個体を送る.　＝＞受け取り側は同期通信
      MPI_Isend(&p, 1, MPI_INDIVIDUAL, i, tNextIndiv, MPI_COMM_WORLD, &req);

      //現世代は受け取ったのでrequestを初期化
      requests[i] = MPI_REQUEST_NULL;
      
      //次の評価を待つ
      MPI_Irecv(&(curGeneration[i]),
                1, MPI_INDIVIDUAL, i, tEval, MPI_COMM_WORLD, &requests[i]);
    }
  }
}
