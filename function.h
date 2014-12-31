#ifndef _FUNCTION_H
#define _FUNCTION_H
#include "myComplex.h"
#include <stdio.h>
#include <stdlib.h>
#include "bool.h"

// TODO : max(i++, j) のような使い方だと,挙動がうまくいかない
#define max(a,b) (a > b ? a : b)
#define min(a,b) (a < b ? a : b)


extern void delete(void *ptr);
extern double dbilinear(double *p, double x, double y, int width, int height);
extern FILE* openFile(const char* file_name);
extern FILE* FileOpen(const char* file_name, const char* mode);
extern bool makeDirectory(const char*);
extern void moveDirectory(const char*);
extern void makeAndMoveDirectory(const char*);
//まだ実装していない部分をexitさせたいときに使う.
#define UN_DONE(msg) printf("not implement "); printf(msg); printf(" yet\n"); exit(2)

#endif
