#include <stdio.h>
#include <stdlib.h>
#include "ntff.h"
#include "function.h"

void ntff_outputEnormTxt(double **e_norm, const char *file_name)
{
  FILE *fp = FileOpen(file_name, "w");
  for(int lambda_nm=LAMBDA_ST_NM; lambda_nm<=LAMBDA_EN_NM; lambda_nm++)
  {
    fprintf(fp, "%d ", lambda_nm);
    for(int ang=0; ang<360; ang++)
    {
      fprintf(fp, "%lf.20 ", e_norm[lambda_nm-LAMBDA_ST_NM][ang]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}

void ntff_outputEnormBin(double **e_norm, const char *file_name)
{
  FILE *fp_b = FileOpen(file_name, "wb");
  for(int l = 0; l<=LAMBDA_EN_NM-LAMBDA_ST_NM; l++)
  {
    if( fwrite( &e_norm[l][0], sizeof(double), 360, fp_b) < 360 )
    {
      printf("error in write binary %s\n", file_name);
      exit(2);
    }
  }
  fclose(fp_b);
}

void ntff_normalize(double array[360])
{
  double sum = 0;
  for(int i=0; i<360; i++)
    sum += array[i];

  if(sum == 0)
    return;
  
  for(int i=0; i<360; i++)
    array[i] /= sum;  
}
