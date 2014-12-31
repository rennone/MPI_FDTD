#ifndef ZIGZAG_H
#define ZIGZAG_H

#include "bool.h"
extern double ( *zigzagModel_EPS(void))(double, double, int, int);
extern bool zigzagModel_isFinish(void);
extern void zigzagModel_moveDirectory(void);
extern void zigzagModel_needSize(int *, int*);
extern void zigzagModel_init();
#endif
