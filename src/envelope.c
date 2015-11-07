/*************************************************************************
    > File Name: envelope.c
    > Author: xu_mijian
*************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sac.h>
#include "sac.h"

int main(int argc, char **argv) {
   char *kname, *outname;
   SACHEAD hd;
   float		*src, *srcenv,  *threshold;// *staenv, *ltaenv,
   int i, j, k, nlen, sta_len, min_bound, max_bound, lta_len, mid_sta_len, lta_min_bound, error;
   float dt, sta_rang, lta_rang, gate, mid_time, rsl;
   error = 0;
 /* input parameters */
   strcpy(kname, argv[1]);
   src = read_sac(kname, &hd);
   nlen = hd.npts;
   dt = hd.delta;
   srcenv = (float *) malloc(sizeof(float) * nlen);
   envelope(nlen,src,srcenv);   
   for (i=0;i<nlen;i++){
	   printf("%f\n",srcenv[i]);
   }
   kname = strdup("enve.sac");    
   write_sac(kname,hd,srcenv);
}
