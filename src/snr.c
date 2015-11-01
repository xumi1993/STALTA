#include <stdio.h>
#include <math.h>
#include <omp.h>

float rssq(float *x, int nlen){
   float sumx, rsq;
   int i;
   sumx = 0.0;
   for(i=0;i<nlen;i++){
      sumx += x[i]*x[i];
   }
   rsq = sqrt(sumx/nlen);
   return rsq;
}

float snr(float *x, int nlenx, float *y, int nleny){
   float spow, npow, rxy;
   spow = (float)pow((double)rssq(x, nlenx),2);
   npow = (float)pow((double)rssq(y, nleny),2);
   rxy = spow / npow;
   return rxy;
}

