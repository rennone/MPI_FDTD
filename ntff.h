#ifndef __NTFF_H
#define __NTFF_H

#define LAMBDA_ST_NM 380
#define LAMBDA_EN_NM 700
#define NTFF_NUM 8192 // 1<<13

extern void ntff_outputEnormTxt(double **e_norm, const char *file_name);
extern void ntff_outputEnormBin(double **e_norm, const char *file_name);

//360°方向に値をもつ配列を正規化
extern void ntff_normalize(double array[360]);
#endif
