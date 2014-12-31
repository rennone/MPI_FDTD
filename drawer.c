#include "drawer.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include "function.h"
static int putBmpHeader(FILE *s, int x, int y, int c);
static int fputc4LowHigh(unsigned long d, FILE *s);
static int fputc2LowHigh(unsigned short d, FILE *s);

#ifdef USE_OPENGL

#include <GL/glew.h>

#include <time.h>

#ifdef MAC_OS
#include <GLUT/glut.h>
#endif
#ifndef MAC_OS
#include <GL/glut.h>
#endif

static const int vertexNum = 4; //頂点数
#define TEX_NX 256  //テクスチャの横幅
#define TEX_NY 256  //縦幅
static colorf texColor[TEX_NX][TEX_NY];
static GLuint ver_buf, tex_buf;
static GLuint texId;

static double (*colorMode)(double complex);
static void colorTransform(double p, colorf *c);

static GLfloat vertices[] =
  {-1.0f, -1.0f, 0.0f,
   +1.0f, -1.0f, 0.0f, 
   +1.0f, +1.0f, 0.0f, 
   -1.0f, +1.0f, 0.0f};

static GLfloat texCoords[] =
  { 0.0f, 0.0f,
    0.0f, 1.0f,
    1.0f, 1.0f,
    1.0f, 0.0f };

//--------------------prototype--------------------//
void drawer_paintImage(int l, int b,int r, int t, int wid,int hei, double complex*);
void drawer_paintModel(int l, int b,int r, int t, int wid,int hei, double *);
void drawer_draw();
//--------------------------------------//


//--------------public Method-----------------//
void (*drawer_getDraw(void))(void)
{
  return drawer_draw;
}

//--------------------------------------//
void drawer_init(enum COLOR_MODE cm)
{
  if(cm == CREAL)
    colorMode = creal;
  else
    colorMode = cnorm;//cabs;

  //初期化
  memset(texColor, 0 , sizeof(texColor));
  
  glGenBuffers(1, &ver_buf);

  glBindBuffer(GL_ARRAY_BUFFER, ver_buf);
  glBufferData(GL_ARRAY_BUFFER, 3*vertexNum*sizeof(GLfloat), vertices, GL_STATIC_DRAW);

  glGenBuffers(1, &tex_buf);
  glBindBuffer(GL_ARRAY_BUFFER, tex_buf);
  glBufferData(GL_ARRAY_BUFFER, 2*vertexNum*sizeof(GLfloat), texCoords, GL_STATIC_DRAW);

  //
  glEnable( GL_TEXTURE_2D );
  glGenTextures( 1, &texId );

  glActiveTexture( GL_TEXTURE0 );

  glBindTexture( GL_TEXTURE_2D, texId );
  glTexImage2D( GL_TEXTURE_2D, 0, 3, TEX_NX, TEX_NY, 0, GL_RGB, GL_FLOAT, texColor);

    //min, maxフィルタ
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );    //min, maxフィルター

}

void drawer_draw()
{  
  glBindBuffer( GL_ARRAY_BUFFER, ver_buf);
  glVertexPointer( 3, GL_FLOAT, 0, 0);

  glBindBuffer( GL_ARRAY_BUFFER, tex_buf);
  glTexCoordPointer(2, GL_FLOAT, 0, 0);
  
  glBindTexture( GL_TEXTURE_2D, texId );

  glTexImage2D( GL_TEXTURE_2D, 0, 3, TEX_NX, TEX_NY, 0, GL_RGB, GL_FLOAT, texColor);  

  glDrawArrays( GL_POLYGON, 0, vertexNum);  
}

void drawer_clear()
{  
  //初期化
  memset(texColor, 0 , sizeof(texColor));
}

//電磁場の値(phi)からテクスチャの色を変える
void drawer_paintImage(int left, int bottom, int right, int top, int width, int height, double complex *phis)
{
  colorf c;
  double complex cphi;
  double ux = 1.0*(right-left)/TEX_NX;
  double uy = 1.0*(top-bottom)/TEX_NY;
  double u = max(ux,uy);
  int i,j;
  double x,y;
  
  for(i=0,x=left; i<TEX_NX && x<right-1; i++, x+=u){
    for(j=0,y=bottom; j<TEX_NY && y<top-1; j++, y+=u)
    {
      cphi = cbilinear(phis,x,y,width,height);
      colorTransform(colorMode(cphi), &c);
      texColor[i][j] = c;
    }
  }
}

//モデルの誘電率(phi)からテクスチャの色を変える(暗くする)
void drawer_paintModel(int left, int bottom, int right, int top, int width, int height, double *phis)
{
  double dphi;
  double ux = 1.0*(right-left)/TEX_NX;
  double uy = 1.0*(top-bottom)/TEX_NY;
  double u = max(ux,uy);
  int i,j;
  double x,y;
  
  for(i=0,x=left; i<TEX_NX && x<right-1; i++, x+=u)
  {
    for(j=0,y=bottom; j<TEX_NY && y<top-1; j++, y+=u)
    {
      dphi = dbilinear(phis,x,y,width,height);
      double n = 1-1.0/dphi;
      texColor[i][j].r -= n;
      texColor[i][j].g -= n;
      texColor[i][j].b -= n;
    }
  }
}

void drawer_paintTest(void){
  colorf c;
  for(int i=0; i<TEX_NX; i++){
    for(int j=0; j<TEX_NY; j++){
      colorTransform(1.0*i/TEX_NX, &c);
      texColor[i][j] = c;
    }
  }
}

void drawer_finish()
{
  printf("drawer_finish not implemented\n");
}

void drawer_screenshot(const char *fileName)
{
  const int width  = glutGet(GLUT_WINDOW_WIDTH);
  const int height = glutGet(GLUT_WINDOW_HEIGHT);

  const int bpp = 24; //1ピクセルセル24ビット
  const int datasize = height*((((width*bpp/8) + 3) >> 2) << 2);
  
  int format = GL_RGB;
  GLubyte *buf = (GLubyte*)malloc(sizeof(GLubyte)*datasize);

  glReadPixels(0, 0, width, height, format, GL_UNSIGNED_BYTE, &buf[0]);

  //RGB -> BGR
  for(int i=0; i<datasize; i+=bpp/8)
  {
    GLubyte tmp = buf[i];
    buf[i] = buf[i+2];
    buf[i+2] = tmp;
//    buf[i]   = buf[i]^buf[i+2];
//    buf[i+2] = buf[i]^buf[i+2];
//    buf[i]   = buf[i]^buf[i+2];
  }

  FILE *fp = fopen(fileName, "wb");
  if(fp==NULL)
  {
    printf("can not open file %s",fileName);
    return;
  }

  if( !putBmpHeader(fp, width, height, bpp) ) {
    printf("can not write headers");
    fclose(fp);
    return;
  }

  if( fwrite((unsigned char*)buf, sizeof(unsigned char), datasize, fp) != datasize)
  {
    printf("can not write data");
    fclose(fp);
    return;
  }
  
  fclose(fp);
  return;
}

#endif

//
//--------------------public Method--------------------//

//--------------------Color Trancform---------------------//
static void colorTransform(double phi, colorf *col)
{
  double range = 1.0; //波の振幅  
  double ab_phi = phi < 0 ? -phi : phi;
  double a = ab_phi < range ? (ab_phi <  range/3.0 ? 3.0/range*ab_phi : (-3.0/4.0/range*ab_phi+1.25) ) : 0.5;
  
  col->r = phi > 0 ? a:0;
  col->b = phi < 0 ? a:0;
  col->g = min(1.0, max(0.0, -3*ab_phi+2));
}

static void modelColorTransform(double n, colorf *col)
{
  colorTransform(0, col);
  col->r -= (1-1.0/n/n);
  col->g -= (1-1.0/n/n);
  col->b -= (1-1.0/n/n);  
}

void drawer_saveImage(char *fileName, colorf **cells, int width, int height)
{
  const int bpp        = 24; //1ピクセルセル24ビット
  const int data_width = (width>>2)<<2; //4の倍数に切り下げ
  const int datasize   = height*data_width*(bpp>>3);

  unsigned char *buf = (unsigned char*)malloc(sizeof(unsigned char)*datasize);

  memset(buf, 0, sizeof(unsigned char)*datasize);

  int k=0;
  for(int j=0; j<height; j++)
    for(int i=0; i<data_width; i++)
    {
      buf[k]   = max(0, min(255, cells[i][j].b*255));
      buf[k+1] = max(0, min(255, cells[i][j].g*255));
      buf[k+2] = max(0, min(255, cells[i][j].r*255));
      k+= (bpp>>3);      
    }
  FILE *fp = fopen(fileName, "wb");
  if(fp==NULL)
  {
    printf("can not open file %s",fileName);
    free(buf);
    return;
  }

  if( !putBmpHeader(fp, data_width, height, bpp) ) {
    printf("can not write headers");
    fclose(fp);
    free(buf);
    return;
  }

  if( fwrite((unsigned char*)buf, sizeof(unsigned char), datasize, fp) != datasize)
  {
    printf("can not write data");
    fclose(fp);
    free(buf);
    return;
  }

  free(buf);
  fclose(fp);
  return;
}

void drawer_outputLineImage(char *fileName, double red[360], double green[360], double blue[360])
{
  int height = 64;
  int width  = 180;
  const int bpp = 24; //1ピクセルセル24ビット
  const int data_width = (width>>2)<<2;
  const int datasize = height*data_width*(bpp>>3);
//  const int datasize = height*((((width*bpp/8) + 3) >> 2) << 2); //横のバイト列は4の倍数でなければならないので切り上げ


  unsigned char *buf = (unsigned char*)malloc(sizeof(unsigned char)*datasize);
  memset(buf, 0, sizeof(unsigned char)*datasize);

  colorf c;
  int k=0;
  for(int j=0; j<height; j++)
    for(int i=0; i<data_width; i++)
    {
      buf[k]   = max(0, min(255, blue[i]*255));
      buf[k+1] = max(0, min(255, green[i]*255));
      buf[k+2] = max(0, min(255, red[i]*255));
      k+= (bpp>>3);
    }

  FILE *fp = fopen(fileName, "wb");
  if(fp==NULL)
  {
    printf("can not open file %s",fileName);
    free(buf);
    return;
  }

  if( !putBmpHeader(fp, data_width, height, bpp) ) {
    printf("can not write headers");
    fclose(fp);
    free(buf);
    return;
  }

  if( fwrite((unsigned char*)buf, sizeof(unsigned char), datasize, fp) != datasize)
  {
    printf("can not write data");
    fclose(fp);
    free(buf);
    return;
  }

  free(buf);
  fclose(fp);
  return;
}
 

void drawer_outputImage(char *fileName, dcomplex *data, double *model, int width, int height)
{
  const int bpp = 24; //1ピクセルセル24ビット
  const int data_width = (width>>2)<<2;
  const int datasize = height*data_width*(bpp>>3);
//  const int datasize = height*((((width*bpp/8) + 3) >> 2) << 2); //横のバイト列は4の倍数でなければならないので切り上げ


  unsigned char *buf = (unsigned char*)malloc(sizeof(unsigned char)*datasize);
  memset(buf, 0, sizeof(unsigned char)*datasize);

  colorf c;
  int k=0;
  for(int j=0; j<height; j++)
    for(int i=0; i<data_width; i++)
    {
      modelColorTransform( model[i*height+j] , &c);
      buf[k]   = max(0, min(255, c.b*255));
      buf[k+1] = max(0, min(255, c.g*255));
      buf[k+2] = max(0, min(255, c.r*255));
      k+= (bpp>>3);
    }

  FILE *fp = fopen(fileName, "wb");
  if(fp==NULL)
  {
    printf("can not open file %s",fileName);
    free(buf);
    return;
  }

  if( !putBmpHeader(fp, data_width, height, bpp) ) {
    printf("can not write headers");
    fclose(fp);
    free(buf);
    return;
  }

  if( fwrite((unsigned char*)buf, sizeof(unsigned char), datasize, fp) != datasize)
  {
    printf("can not write data");
    fclose(fp);
    free(buf);
    return;
  }

  free(buf);
  fclose(fp);
  return;
}


/*
  putBmpHeader BMPヘッダ書出
	
  BMPファイルのヘッダを書き出す

  戻り値
  int:0…失敗, 1…成功
  
  引数
  FILE *s:[i] 出力ストリーム
  int x:[i] 画像Xサイズ(dot, 1〜)
  int y:[i] 画像Yサイズ(dot, 1〜)
  intc:[i] 色ビット数(bit/dot, 1 or 4 or 8 or 24)
*/
static int putBmpHeader(FILE *s, int x, int y, int c)
{
  int i;
  int color; /* 色数 */
  unsigned long int bfOffBits; /* ヘッダサイズ(byte) */

  /* 画像サイズが異常の場合,エラーでリターン */
  if (x <= 0 || y <= 0) {
    return 0;
  }

  /* 出力ストリーム異常の場合,エラーでリターン */
  if (s == NULL || ferror(s)) {
    return 0;
  }

  /* 色数を計算 */
  if (c == 24) {
    color = 0;
  } else {
    color = 1;
    for (i=1;i<=c;i++) {
      color *= 2;
    }
  }

  /* ヘッダサイズ(byte)を計算 */
  /* ヘッダサイズはビットマップファイルヘッダ(14) + ビットマップ情報ヘッダ(40) + 色数 */
  bfOffBits = 14 + 40 + 4 * color;

  /* ビットマップファイルヘッダ(計14byte)を書出 */
  /* 識別文字列 */
  fputs("BM", s);

  /* bfSize ファイルサイズ(byte) */
  fputc4LowHigh(bfOffBits + (unsigned long)x * y, s);

  /* bfReserved1 予約領域1(byte) */
  fputc2LowHigh(0, s);

  /* bfReserved2 予約領域2(byte) */
  fputc2LowHigh(0, s);

  /* bfOffBits ヘッダサイズ(byte) */
  fputc4LowHigh(bfOffBits, s);

  /* ビットマップ情報ヘッダ(計40byte) */
  /* biSize 情報サイズ(byte) */
  fputc4LowHigh(40, s);

  /* biWidth 画像Xサイズ(dot) */
  fputc4LowHigh(x, s);

  /* biHeight 画像Yサイズ(dot) */
  fputc4LowHigh(y, s);

  /* biPlanes 面数 */
  fputc2LowHigh(1, s);

  /* biBitCount 色ビット数(bit/dot) */
  fputc2LowHigh(c, s);

  /* biCompression 圧縮方式 */
  fputc4LowHigh(0, s);

  /* biSizeImage 圧縮サイズ(byte) */
  fputc4LowHigh(0, s);

  /* biXPelsPerMeter 水平解像度(dot/m) */
  fputc4LowHigh(0, s);

  /* biYPelsPerMeter 垂直解像度(dot/m) */
  fputc4LowHigh(0, s);

  /* biClrUsed 色数 */
  fputc4LowHigh(0, s);

  /* biClrImportant 重要色数 */
  fputc4LowHigh(0, s);

  /* 書出失敗ならエラーでリターン */
  if (ferror(s)) {
    return 0;
  }

  /* 成功でリターン */
  return 1;
}

/*
  fputc2LowHigh 2バイトデータ書出(下位〜上位)
	
  2バイトのデータを下位〜上位の順で書き出す

  戻り値
  int:EOF…失敗, EOF以外…成功
  引数
  unsigned short d:[i] データ
  FILE *s:[i] 出力ストリーム
*/
static int fputc2LowHigh(unsigned short d, FILE *s)
{
  putc(d & 0xFF, s);
  return putc(d >> CHAR_BIT, s);
}

/*
  fputc4LowHigh 4バイトデータ書出(下位〜上位)
	
  4バイトのデータを下位〜上位の順で書き出す

  ●戻り値
  int:EOF…失敗, EOF以外…成功
  ●引数
  unsigned long d:[i] データ
  FILE *s:[i] 出力ストリーム
*/
static int fputc4LowHigh(unsigned long d, FILE *s)
{
  putc(d & 0xFF, s);
  putc((d >> CHAR_BIT) & 0xFF, s);
  putc((d >> CHAR_BIT * 2) & 0xFF, s);
  return putc((d >> CHAR_BIT * 3) & 0xFF, s);
}



