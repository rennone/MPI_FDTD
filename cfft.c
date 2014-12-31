/*		fft.c		*/
#include <stdio.h>
#include <stdlib.h>
#include "cfft.h"
#define _USE_MATH_DEFINES
#include <math.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

static void swapx(Complex *a, Complex *b)
{
  Complex w;
   w = *a;
  *a = *b;
  *b =  w;
  return;
}

static void fft1(double ar[], double ai[], int n, int iter, int flag)
{
  int i, it, j, j1, j2, k, xp, xp2;
  double arg, dr1, dr2, di1, di2, tr, ti, w, wr, wi;

  if(n < 2)
  {
    fprintf(stderr, "Error : n < 2  in fft1()\n");
    return;
  }
  
  if(iter <= 0)
  {
    iter = 0;
    i = n;
    while((i>>=1))
      iter++;
  }
  printf("iter = %d\n",iter);

  if(n != (1<<iter))
  {
    fprintf(stderr, "Error : n != 2 ^ k  in fft1()\n");
    return;
  }

  w = (flag? M_PI: -M_PI) / (double)n;
  xp2 = n;
  for(it = 0; it < iter; it++)
  {
    xp = xp2;
    xp2 /= 2;
    w *= 2;
//    for(k = 0, i = - xp; k < xp2; i++)
    for(k = 0, i = - xp; k < xp2; i++, k++)
    {
//      wr = cos(arg = w * k++);
//      wi = sin(arg);
      wr = cos(w * k);
      wi = sin(w*k);
      for(j = xp; j <= n; j += xp)
      {
        j2 = (j1 = j + i) + xp2;
        tr = (dr1 = ar[j1]) - (dr2 = ar[j2]);
        ti = (di1 = ai[j1]) - (di2 = ai[j2]);
        ar[j1] = dr1 + dr2;
        ai[j1] = di1 + di2;
        ar[j2] = tr * wr - ti * wi;
        ai[j2] = ti * wr + tr * wi;
      }
    }
  }
  j = j1 = n / 2;
  j2 = n - 1;
  for(i = 1; i < j2; i++)
  {
    if(i < j)
    {
      w = ar[i];
      ar[i] = ar[j];
      ar[j] = w;
      w = ai[i];
      ai[i] = ai[j];
      ai[j] = w; 
    }
    k = j1;
    while(k <= j)
    {
      j -= k;
      k /= 2;
    }
    j += k;
  }
  if(flag == 0)	return;
  w = 1. / (double)n;
  for(i = 0; i < n; i++)
  {
    ar[i] *= w;
    ai[i] *= w;
  }
  return;
}

static void fft1x(Complex a[], int n, int iter, int flag)
{
  int i, j, j1, j2, k, xp, xp2;
  double arg, w;
  Complex d1, d2, *p, *q, t, ww;

  if(n < 2)
  {
    fprintf(stderr, "Error : n < 2  in fft1()\n");
    return;
  }
  
  if(iter <= 0){
    iter = 0;
    i = n;
    while(i>>=1)
      iter++;
  }

  // nが2の累乗かチェック
  if(n != (1<<iter))
  {
    fprintf(stderr, "Error : n != 2 ^ k  in fft1()\n");
    return;
  }

  //符号
  double sign = (flag == 1)? 1.: -1.;
  
  xp2 = n;
  for(int it = 0; it < iter; it++)
  {
    xp = xp2;
    xp2 = xp / 2;
    w = -M_PI / xp2;
    for(int k = 0; k < xp2; k++)
    {
      ww = cexp(I*sign*w*k);
      i = k - xp;
      for(j = xp, p = a + i + xp, q = a + i + xp + xp2;
          j <= n;
          j += xp, p += xp, q += xp)
      {
        t  = (*p) - (*q);
        *p = (*p) + (*q);
        *q =    t * (ww);
      }
    }
  }
  
  j1 = n / 2;
  j2 = n - 1;
  j = 1;
  for(i = 1, p = a; i <= j2; i++, p++)
  {
    if(i < j)
      swapx(p, a + j - 1);
    
    k = j1;
    while(k < j)
    {
      j -= k;
      k /= 2;
    }
    j += k;
  }
  
  if(flag == 0)
    return;
  
  w = n;
  for(i = 0, p = a; i < n; i++, p++)
    *p = (*p)/w;
  
  return;
}

static void fft2(double ar[], double ai[], int n, int nmax, int flag)
{
  int i, iter, j, k;
  double *p, *q, *wr, *wi;

  if(n < 2)
  {
    fprintf(stderr, "Error : Illegal parameter in fft2()\n");
    return;
  }
  iter = 0;
  i = n;
  while((i /= 2) != 0)	iter++;
  j = 1;
  for(i = 1; i <= iter; i++)	j *= 2;
  if(n != j)
  {
    fprintf(stderr, "Error : n != 2 ^ k  in fft2()\n");
    return;
  }
  wr = (double *)malloc(n * sizeof(double));
  if(wr == NULL)
  {
    fprintf(stderr, "Error : Out of memory  in fft2()\n");
    return;
  }
  wi = (double *)malloc(n * sizeof(double));
  if(wi == NULL)
  {
    fprintf(stderr, "Error : Out of memory  in fft2()\n");
    free((char *)wr);
    return;
  }

  for(j = 0; j < n; j++, k += nmax)
  {
    for(i = 0, p = ar + j, q = ai + j; i < n; i++, p += nmax, q += nmax)
    {
      *(wr + i) = *p;
      *(wi + i) = *q;
    }
    fft1(wr, wi, n, iter, flag);
    for(i = 0, p = ar + j, q = ai + j; i < n; i++, p += nmax, q += nmax)
    {
      *p = *(wr + i);
      *q = *(wi + i);
    }
  }
  for(i = k = 0; i < n; i++, k += nmax)
  {
    for(j = 0, p = ar + k, q = ai + k; j < n; j++)
    {
      *(wr + j) = *p++;
      *(wi + j) = *q++;
    }
    fft1(wr, wi, n, iter, flag);
    for(j = 0, p = ar + k, q = ai + k; j < n; j++)
    {
      *p++ = *(wr + j);
      *q++ = *(wi + j);
    }
  }
  free((char *)wr);
  free((char *)wi);
  return;
}

static void fft2x(Complex a[], int n, int nmax, int flag)
{
  int i, iter, j, k;
  Complex *p, *w;

  if(n < 2)
  {
    fprintf(stderr, "Error : Illegal parameter in fft2()\n");
    return;
  }
  
  iter = 0;
  i = n;
  while((i /= 2) != 0)
    iter++;
  
  j = 1;
  for(i = 1; i <= iter; i++)
    j *= 2;

  if(n != j)
  {
    fprintf(stderr, "Error : n != 2 ^ k  in fft2x()\n");
    return;
  }
  
  w = (Complex *)malloc(n * sizeof(Complex));
  
  if(w == NULL)
  {
    fprintf(stderr, "Error : Out of memory  in fft2x()\n");
    return;
  }

  for(j = 0; j < n; j++, k += nmax)
  {
    for(i = 0, p = a + j; i < n; i++, p += nmax)
      *(w + i) = *p;
    
    fft1x(w, n, iter, flag);
    for(i = 0, p = a + j; i < n; i++, p += nmax)
      *p = *(w + i);
  }
  
  for(i = k = 0; i < n; i++, k += nmax)
  {
    for(j = 0, p = a + k; j < n; j++)
      *(w + j) = *p++;
    fft1x(w, n, iter, flag);
    for(j = 0, p = a + k; j < n; j++)
      *p++ = *(w + j);
  }
  free((char *)w);
  return;
}

void cfft(Complex *data, int n)
{
  fft1x(data, n, 0, 0);
}
