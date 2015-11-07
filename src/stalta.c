#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sac.h>
#include "sac.h"

int main(int argc, char **argv) {
   char *kname, *outname, *proto, *staname;
   char filter[2];
   SACHEAD hd;
   float		*src, *srcenv,  *threshold;
   int i, j, k, nlen, sta_len, min_bound, max_bound, lta_len, error,
       sm_factor, time_begin_len, time_end_len, isout, isdetrend,
       isfilter, order, passes;
   float dt, sta_rang, lta_rang, gate, mid_time, sta, lta,
         buf_aft, buf_sta_bef, buf_lta_bef, frac, f1, f2,
         transition_bandwidth, attenuation, low, high;
   float *smooth(float *x, int nlen, int N);
   void	detrend(float *y, int n);
   error = 0;
   isout = 0;
   isdetrend = 0;
   isfilter = 0;
   transition_bandwidth = 0.0;
   attenuation = 0.0;
   passes = 2;
 /* input parameters */
   for (i=1; !error && i < argc; i++) {
      if (argv[i][0] == '-') {
         switch(argv[i][1]) {
            case 'I':
               kname = strdup(&argv[i][2]);
               break;
            case 'S':
               sscanf(&argv[i][2],"%f",&sta_rang);
               break;
            case 'L':
               sscanf(&argv[i][2],"%f",&lta_rang);
               break;
            case 'O':
               outname = strdup(&argv[i][2]);
               isout = 1;
               break;
            case 'M':
               sscanf(&argv[i][2],"%d",&sm_factor);
               break;
            case 'T':
               sscanf(&argv[i][2],"%f",&gate);
               break;
            case 'D':
               isdetrend = 1;
                break;
            case 'F':
               isfilter = 1;
               j = sscanf(&argv[i][2],"%2s/%d/%f/%f",&filter,&order,&f1,&f2);
               if (j<3) {error=TRUE;}
               if (j==4&&!strcmp(filter,"bp")){
                   proto = strdup("BP");
                   low = f1; high = f2;
               }
               if (j==3&&strcmp(filter,"hp")){
                   proto = strdup("HP");
                   low = f1; high = 0;
               }
               if (j==3&&strcmp(filter,"lp")){
                   proto = strdup("LP");
                   low = 0; high = f1;
               }
               break;
            default:
               error = TRUE;
               break;
         }   
      }
   }
   if (argc == 1 || error) {
      fprintf(stderr,"Usage: stalta -Iinput.sac -Ssta_len -Llta_len -Msmooth_factor -Tthreshold -D -Ftype/order/f1[/f2] [-Ooutput.sac]\n");
      return -1;
   }
//   printf("%s\n", kname);
   src = read_sac(kname, &hd);
   nlen = hd.npts;
   dt = hd.delta;
   staname = hd.kstnm;
   if (isdetrend){detrend(src, nlen);}
   if (isfilter){
       xapiir(src, nlen, SAC_BUTTERWORTH,
        transition_bandwidth, attenuation,
        order, proto,
        low, high, 
	    dt, passes);}
   srcenv = (float *) malloc(sizeof(float) * nlen);
   envelope(nlen, src, srcenv);
   
   sta_len = (int)(sta_rang / dt);
   lta_len = (int)(lta_rang / dt);
   frac = lta_len/sta_len;
   threshold = (float *) malloc(sizeof(float) * nlen);
   min_bound = lta_len-sta_len;
   max_bound = lta_len;
   for (i=0;i<lta_len;i++){threshold[i] = 0;}
   sta = 0;
   for (i=min_bound;i<max_bound;i++){
      sta += srcenv[i]*srcenv[i];
   }
   lta = 0;
   for (i=0;i<lta_len;i++){
      lta += srcenv[i]*srcenv[i];
   }
   for (i=lta_len;i<nlen;i++){ 
      buf_aft = srcenv[i]*srcenv[i];
      buf_sta_bef = srcenv[i-sta_len]*srcenv[i-sta_len];
      buf_lta_bef = srcenv[i-lta_len]*srcenv[i-lta_len];
      sta += buf_aft - buf_sta_bef;
      lta += buf_aft - buf_lta_bef;
      threshold[i] = (sta/lta)*frac;
   }
   for (i=nlen-sta_len/2;i<nlen;i++){threshold[i] = 0;}
   threshold = smooth(threshold,nlen,sm_factor);

   k = 0;
   for (i=1;i<nlen-1;i++){
      if (threshold[i-1]<threshold[i] && threshold[i+1]<threshold[i] && threshold[i]>gate){
         ++k;
         printf("%s %d %d %f %f\n",staname,k,i,i*dt,threshold[i]);
      }
   }

/* Output  */    
   if (isout){
      write_sac(outname,hd,threshold);
   }
    
   free(srcenv);
   free(threshold);
   return 0;
}
