#ifndef _DRAWER_H
#define _DRAWER_H
#include "myComplex.h"


typedef struct {
 float r,g,b;
}colorf;

enum COLOR_MODE
{
  CREAL,
  CABS
};

#ifdef USE_OPENGL
extern void (*drawer_getDraw(void))(void);
extern void drawer_paintImage(int l, int b, int r, int t,int wid, int hei, double complex*);
extern void drawer_paintModel(int l, int b, int r, int t,int wid, int hei, double *);
extern void drawer_paintTest(void);
extern void drawer_init(enum COLOR_MODE);
extern void drawer_finish(void);
extern void drawer_draw(void);
extern void drawer_clear(void);
#endif //USE_OPENGL

extern void drawer_outputImage(char *fileName, dcomplex* data, double *model, int width, int height);
extern void drawer_screenshot(const char* filename);
extern void drawer_outputLineImage(char *fileName, double red[360], double green[360], double blue[360]);
extern void drawer_saveImage(char *fileName, colorf **cells, int width, int height);
#endif //_DRAWER_H
