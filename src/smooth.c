#include <stdio.h>
#include <stdlib.h>

float *smooth(float *x, int nlen, int N){
    int i, j, num;
    float *smooth, sum;
    smooth = (float *) malloc(sizeof(float) * nlen);   
    sum = 0;
    for (i=0;i<N;i++){
        sum += x[i];
        smooth[i] = 0;
    }
    for (i=(N-1)/2;i<nlen-(N-1)/2;i++){
        sum += x[i+(N-1)/2] - x[i-(N-1)/2-1];
        smooth[i] = sum/N;
    }
    for (i=nlen-(N-1)/2;i<nlen;i++){
        smooth[i] = 0;
    }
    return smooth;
}
