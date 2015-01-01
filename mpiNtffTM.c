#include <unistd.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "ntff.h"
#include "mpiNtffTM.h"
#include "field.h"
#include "function.h"
#include "bool.h"
#include "cfft.h"
/*
  bottom(k) = k-1
  top(k)    = k+1
  left(k)   = k-SUB_N_PY
  right(k)  = k+SUB_N_PY
*/

static double R0;

static int sub_tp, sub_bm, sub_rt, sub_lt;
static bool IN_TP, IN_BM, IN_LT, IN_RT;
static int sub_ylt, sub_yrt;
static int sub_xtp, sub_xbm;
//static dcomplex *Ux, *Uy, *Wz;

void mpiNtffTM_init()
{
  NTFFInfo nInfo = field_getNTFFInfo();
  
  R0 = 1.0e6 * field_toCellUnit(500);
  
  int tp = nInfo.top;    int bm = nInfo.bottom;  //上下
  int rt = nInfo.right;  int lt = nInfo.left;	 //左右

  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  
  sub_tp = tp - subInfo_s.OFFSET_Y;  sub_bm = bm - subInfo_s.OFFSET_Y;
  sub_rt = rt - subInfo_s.OFFSET_X;  sub_lt = lt - subInfo_s.OFFSET_X;

  //以下どれかでも満たせば積分路上に無い
  bool outX = sub_rt <= 0 || sub_lt >= subInfo_s.SUB_N_PX-1; //rtより右, もしくはltより左の小領域
  bool outY = sub_tp <= 0 || sub_bm >= subInfo_s.SUB_N_PY-1; //tpより上, もしくはbmより下の小領域
  
  // 小領域内にどの積分面が存在するか
  IN_TP = (0 < sub_tp && sub_tp < subInfo_s.SUB_N_PY-1) && !outX;
  IN_BM = (0 < sub_bm && sub_bm < subInfo_s.SUB_N_PY-1) && !outX;
  IN_RT = (0 < sub_rt && sub_rt < subInfo_s.SUB_N_PX-1) && !outY;
  IN_LT = (0 < sub_lt && sub_lt < subInfo_s.SUB_N_PX-1) && !outY;

  sub_ylt=-1, sub_yrt=-2;
  if(IN_TP || IN_BM)
  {
    sub_yrt = min(subInfo_s.SUB_N_PX-2, max( 1, sub_rt) );
    sub_ylt = min(subInfo_s.SUB_N_PX-2, max( 1, sub_lt) );
  }

  sub_xtp=-2, sub_xbm=-1;
  if(IN_RT || IN_LT)
  {
    sub_xbm = min(subInfo_s.SUB_N_PY-2, max( 1, sub_bm+1) );  //bm,tpですでに計算しているため, ひとつずれる
    sub_xtp = min(subInfo_s.SUB_N_PY-2, max( 1, sub_tp-1) );  //
  }
}

void mpiNtffTM_TimeTranslate(dcomplex *Ux, dcomplex *Uy, dcomplex *Wz, dcomplex *Eth, dcomplex *Eph)
{
  const double w_s = field_getOmega();
  // don't divide by R0 -> textbook by uno
  const double complex coef = 1.0/(4*M_PI*C_0_S)*csqrt( 2*M_PI*C_0_S/(I*w_s) );
  const int maxTime = field_getMaxTime();
  NTFFInfo nInfo = field_getNTFFInfo();  
  double theta = 0;
  double ToRad = M_PI/180.0;
  //TM Time Translate
  for(int ang=0, k=0; ang<360; ang++, k+=nInfo.arraySize)
  {
    double phi = ang*ToRad;
    double sx = cos(theta)*cos(phi);
    double sy = cos(theta)*sin(phi);
    double sz = -cos(theta); //宇野先生の本では -sin(theta)になってる(式としては本の方が正しいけどこのプログラムではこうしないと動かない)
    double px = -sin(phi);
    double py = cos(phi);
    
    //TODO maxTime = nInfo.arraySize ??
    for(int i=0; i < maxTime; i++)
    {
      double complex WTH = 0 + 0 + Wz[k+i]*sz;
      double complex WPH = 0 + 0;
      double complex UTH = Ux[k+i]*sx + Uy[k+i]*sy + 0;
      double complex UPH = Ux[k+i]*px + Uy[k+i]*py;
      double complex ETH = coef*(-Z_0_S*WTH-UPH);
      double complex EPH = coef*(-Z_0_S*WPH+UTH);
      
      Eth[k+i] = ETH; //TODO : 物理単位に変換が必要かも
      Eph[k+i] = EPH; //Note : TMモードでは必要無い
    }
  }
}

//時間領域の遠方解 Ux, Uy, WzからEthを求め
//それをfftした結果から波長,反射角ごとのEthのノルムを返す
// res[λ][deg]でアクセスできる. 0 <= deg <= 359, stLambda <= λ <= enLamba
double** mpiNtffTM_TimeCalcNorm(dcomplex *Ux, dcomplex *Uy, dcomplex *Wz, int *stLambda, int *enLamba)
{
  *stLambda = LAMBDA_ST_NM; //戻り値の波長の範囲を指定
  *enLamba  = LAMBDA_EN_NM;
  
  const int maxTime = field_getMaxTime();
  NTFFInfo nInfo = field_getNTFFInfo();
  dcomplex *Eth, *Eph;  
  Eth = newDComplex(360*nInfo.arraySize);
  Eph = newDComplex(360*nInfo.arraySize);

  //Eth, Ephを計算
  mpiNtffTM_TimeTranslate(Ux,Uy,Wz,Eth,Eph);
  
  double **out_ref = (double**)malloc(sizeof(double*) * (LAMBDA_EN_NM-LAMBDA_ST_NM+1));

  double *datas = newDouble(360*(LAMBDA_EN_NM-LAMBDA_ST_NM+1));
  for(int l=0; l<=LAMBDA_EN_NM-LAMBDA_ST_NM; l++)
    out_ref[l] = &datas[360*l];

  //fft用に2の累乗の配列を確保
  dcomplex *eth = newDComplex(NTFF_NUM);
  for(int ang=0; ang<360; ang++)
  {
    int k= ang*nInfo.arraySize;

    memset((void*)eth,0,sizeof(dcomplex)*NTFF_NUM);               //0で初期化
    memcpy((void*)eth, (void*)&Eth[k], sizeof(dcomplex)*maxTime); //コピー
    cfft(eth, NTFF_NUM); //FFT

    FieldInfo fInfo = field_getFieldInfo();
    for(int lambda_nm=LAMBDA_ST_NM; lambda_nm<=LAMBDA_EN_NM; lambda_nm++)
    {
      //線形補完 TODO : index = n-1 となるほどの小さいlambdaを取得しようとするとエラー
      double p = C_0_S * fInfo.h_u_nm * NTFF_NUM / lambda_nm;
      int index = floor(p);
      p = p-index;
      out_ref[lambda_nm-LAMBDA_ST_NM][ang] = ((1-p)*cnorm(eth[index]) + p*cnorm(eth[index+1]))/NTFF_NUM;
    }
  }
  freeDComplex(eth);

  free(Eth);
  free(Eph);

  return out_ref;
}

static void calc(double time_plus_timeShift, dcomplex eh,  dcomplex *UW_ang)
{
  int m = floor(time_plus_timeShift+0.5);
  double a = (0.5 + time_plus_timeShift - m);
  double b = 1.0-a;
  double ab = a-b;
  UW_ang[m-1] += eh*b;
  UW_ang[m]   += eh*ab;
  UW_ang[m+1] -= eh*a;
}

void mpiNtffTM_TimeCalc(dcomplex *Hx, dcomplex *Hy, dcomplex *Ez, dcomplex *Ux, dcomplex *Uy, dcomplex *Wz)
{
  FieldInfo_S fInfo_s      = field_getFieldInfo_S();
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  NTFFInfo nInfo = field_getNTFFInfo();
  
//  const double coef = 1.0/(4*M_PI*C_0_S*R0);

  //分割領域系に置ける, 中心の位置
  double cx = fInfo_s.N_PX/2 - subInfo_s.OFFSET_X;
  double cy = fInfo_s.N_PY/2 - subInfo_s.OFFSET_Y;

  //分割領域のインデックスに変換 (0以下, SUB_N_PX, PY以上だと範囲外)
  int tp = sub_tp; //上面
  int bm = sub_bm; //下面
  int rt = sub_rt; //右
  int lt = sub_lt; //左

  double timeE = field_getTime() - 1;   //t - Δt
  double timeH = field_getTime() - 0.5; //t - Δt/2

  for(int ang=0; ang<360; ang++)
  {
    double rad = ang*M_PI/180.0;
    double r1x = cos(rad), r1y = sin(rad);

    //角度angの0番目のインデックス
    int index_ang = ang * nInfo.arraySize;
    
    dcomplex *Ux_ang = &Ux[index_ang];
    dcomplex *Uy_ang = &Uy[index_ang];
    dcomplex *Wz_ang = &Wz[index_ang];

    //bottom normal(0,-1,0) //W = Js = n × H = ( 0, 0, Hx)  U = Ms = E × n = (Ez, 0,  0)
    if ( IN_BM ) {
      //積分路の端で無ければ, 分割領域の端から端までが積分路である
      // 一番外側はのりしろなので 1~SUB_N_PX-2まで, ただし lt <= x < rt なので 1~SUB_N_PY-1になってる
      const int subLeft  = sub_ylt;
      const int subRight = sub_yrt;
      
      for(int i=subLeft; i<subRight; i++)
      {
        //原点との距離
        const double r2x = i  - cx;
        const double r2y = bm - cy;        
        const double timeShift = -(r1x*r2x + r1y*r2y)/C_0_S + nInfo.RFperC;
        
        int k = field_subIndex(i, bm);
        calc(timeE+timeShift, Ez[k], Ux_ang);
        calc(timeH+timeShift, 0.5*(Hx[k]+Hx[k-1]), Wz_ang);
      }
    }
    
    //right normal(1,0,0)    //Js = n × H = (0, 0,Hy)  Ms = E × n = ( 0,Ez,0)
    if ( IN_RT ) {
      int subTop    = sub_xtp;
      int subBottom = sub_xbm;
      for ( int j=subBottom; j<subTop; j++ ) {
        const double r2x = rt-cx;
        const double r2y =  j-cy;
        const double timeShift = -(r1x*r2x + r1y*r2y)/C_0_S + nInfo.RFperC;        
        
        int k = field_subIndex(rt, j);
        calc(timeE+timeShift, Ez[k], Uy_ang);
        calc(timeH+timeShift, 0.5 * ( Hy[k] + Hy[k-subInfo_s.SUB_N_PY] ), Wz_ang);        
      }
    }

    //top normal(0,1,0)  //Js = n × H = (0, 0,-Hx)  Ms = E × n = (-Ez, 0,  0)
    if ( IN_TP ) {
      int subLeft  = sub_ylt;
      int subRight = sub_yrt;
      for ( int i=subLeft; i<subRight; i++ ) {
        const double r2x  =  i-cx;
        const double r2y  = tp-cy;
        const double timeShift = -(r1x*r2x + r1y*r2y)/C_0_S + nInfo.RFperC; 

        int k = field_subIndex(i, tp);
        calc(timeE+timeShift,               -Ez[k], Ux_ang);
        calc(timeH+timeShift, -0.5*(Hx[k]+Hx[k-1]), Wz_ang);
      }
    }

    // (left,top) -> (left,bottom)  normal(-1,0,0)   //Js = n × H = (0,0,-Hy)    Ms = E × n = (0,-Ez,0)
    if ( IN_LT ) {
      int subTop    = sub_xtp;
      int subBottom = sub_xbm;
      for ( int j=subBottom; j<subTop; j++ ) {
        const double r2x = lt-cx;
        const double r2y =  j-cy;
        const double timeShift = -(r1x*r2x + r1y*r2y)/C_0_S + nInfo.RFperC;

        int k = field_subIndex(lt, j);
        calc(timeE+timeShift, -Ez[k], Uy_ang);
        calc(timeH+timeShift, -0.5*(Hy[k]+Hy[k - subInfo_s.SUB_N_PY] ), Wz_ang);        
      }
    }
  }
}

//周波数領域のNTFF
void mpiNtffTM_Frequency( dcomplex *Hx, dcomplex *Hy, dcomplex *Ez, dcomplex resultEz[360])
{
  FieldInfo fInfo          = field_getFieldInfo();
  FieldInfo_S fInfo_s      = field_getFieldInfo_S();
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  WaveInfo_S waveInfo_s    = field_getWaveInfo_S();
  
  //分割領域系に置ける, 中心の位置
  double cx = fInfo_s.N_PX/2 - subInfo_s.OFFSET_X;
  double cy = fInfo_s.N_PY/2 - subInfo_s.OFFSET_Y;

  dcomplex coef = csqrt( I*waveInfo_s.K_s/(8*M_PI*R0) ) * cexp(I*waveInfo_s.K_s*R0);
  
  NTFFInfo nInfo = field_getNTFFInfo();
  //分割領域のインデックスに変換 (0以下, SUB_N_PX, PY以上だと範囲外)
  int tp =    nInfo.top - subInfo_s.OFFSET_Y; //上面
  int bm = nInfo.bottom - subInfo_s.OFFSET_Y; //下面
  int rt = nInfo.right  - subInfo_s.OFFSET_X; //右
  int lt = nInfo.left   - subInfo_s.OFFSET_X; //左

  for(int ang=0; ang<360; ang++)
  {
    double rad = ang*M_PI/180.0;
    double rx  = cos(rad), ry = sin(rad); //方向ベクトル

    dcomplex Nz = 0;
    dcomplex Lx = 0;
    dcomplex Ly = 0;
		
    // (left,bottom) -> (right,bottom)
    // 法線ベクトルはn=(0, -1)
    //積分路が分割領域内にあるか確認 
    if ( IN_BM )
    {
      int subLeft  = max(1, lt);                    //左右の端で無ければ, 分割領域の端から端までが積分路である
      int subRight = min(subInfo_s.SUB_N_PX-1, rt); //領域の端っこでなれば, 全部が積分路なのでSUB_N_PX-1 でなく SUB_N_PX まで

      for ( int i=subLeft; i<subRight; i++ )
      {
        double r2x  =  i-cx; //cx, cyを分割領域空間に変換してるので, そのまま引き算が出来る
        double r2y  = bm-cy;
        int k = field_subIndex(i, bm);
        dcomplex C_EZ = Ez[k];        
        dcomplex C_HX = 0.5 * ( Hx[k] + Hx[k-1] );
        double innerProd = rx*r2x + ry*r2y;  //内積
        Nz  += C_HX * cexp( I * waveInfo_s.K_s * innerProd );  // J = n × H = (0   ,0,C_HX)
        Lx  += C_EZ * cexp( I * waveInfo_s.K_s * innerProd );  // M = E × n = (C_EZ,0,   0)
      }
    }

    // (right,bottom) -> (right,top) n=(1,0)
    if ( IN_RT )
    {
      int subTop    = min(subInfo_s.SUB_N_PY-1, tp);
      int subBottom = max(1, bm);
      for ( int j=subBottom; j<subTop; j++ )
      {
        double r2x  = rt-cx; //cx, cyを分割領域空間に変換してるので, そのまま引き算が出来る
        double r2y  =  j-cy;

        int k = field_subIndex(rt, j);
        dcomplex C_EZ = Ez[k];
        dcomplex C_HY = 0.5 * ( Hy[k] + Hy[k-subInfo_s.SUB_N_PY] );

        double innerProd = rx*r2x + ry*r2y;  //内積
        Nz  += C_HY * cexp( I * waveInfo_s.K_s * innerProd );  // J = n × H = (0,   0, C_HY)
        Ly  += C_EZ * cexp( I * waveInfo_s.K_s * innerProd );  // M = E × n = (0,C_EZ,    0)
      }
    }


    // (right,top) -> (left,top)  n=(0,1)
    if ( IN_TP )
    {      
      //左右の端で無ければ, 分割領域の端から端までが積分路である
      int subLeft  = max(1, lt);
      int subRight = min(subInfo_s.SUB_N_PX-1, rt); //領域の端っこでなれば, 全部が積分路なのでSUB_N_PX-1 でなく SUB_N_PX まで      
      for ( int i=subLeft; i<subRight; i++ )
      {
        double r2x  =  i-cx;
        double r2y  = tp-cy;

        int k = field_subIndex(i, tp);
        dcomplex C_EZ = Ez[k];
        dcomplex C_HX = 0.5 * ( Hx[k] + Hx[k-1] );

        double innerProd = rx*r2x  + ry*r2y;  //内積
        Nz   -= C_HX * cexp( I * waveInfo_s.K_s * innerProd );  // J = n × H = (0,    0, -C_HX)
        Lx   -= C_EZ * cexp( I * waveInfo_s.K_s * innerProd );  // M = E × n = (-C_EZ,0,     0)
      }
    }

    // (left,top) -> (left,bottom)
    if ( IN_LT )
    {
      int subTop    = min(subInfo_s.SUB_N_PY-1, tp);
      int subBottom = max(1, bm); 
      for ( int j=subBottom; j<subTop; j++ )
      {
        double r2x = lt-cx;
        double r2y =  j-cy;

        int k = field_subIndex(lt, j);
        dcomplex C_EZ  = Ez[k];
        dcomplex C_HY  = 0.5 * ( Hy[k] + Hy[k - subInfo_s.SUB_N_PY] );

        double innerProd = rx*r2x  + ry*r2y;  //内積      
        Nz   -= C_HY * cexp( I * waveInfo_s.K_s * innerProd );  // J = n × H = (0,     0, -C_HY)		
        Ly   -= C_EZ * cexp( I * waveInfo_s.K_s * innerProd );  // M = E × n = (0, -C_EZ,     0)
      }
    }
    dcomplex Lphi = -Lx*sin(rad) + Ly*cos(rad); //極座標変換
    resultEz[ang] = coef * ( Z_0_S*Nz + Lphi ) * fInfo.h_u_nm;
  }
  
  FILE *fp = fopen("ntffStr.txt", "w");
  
  for(int ang = 0; ang<360; ang++)
  {
    fprintf(fp, "%.18lf\n", cnorm(resultEz[ang]));
  }
  printf("ntff tm frequency finished\n");
}
