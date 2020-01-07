/* 
Author: Dillon Bourne
Contact: dbourne@uoguelph.ca
Filename: forbacksubs.c 
*/
#include "bla.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>

double* backsubs(int m, int n, double **R, double *b){
	if(m < n){
		fprintf(stderr, "The dimension m must be greater than or equal to n.\n");
		exit(1);
	}
	double *x;
	x = malloc(m*sizeof(double));
	
	x[m-1] = (b[m-1] / R[m-1][m-1]);
	
	for(int i=m-2; i >= 0; i--){
		double sum = 0.0;
		for(int j=i+1; j <m; j++){
			sum += R[i][j] * x[j];
		}
		x[i] = (b[i] - sum) / R[i][i];
	}
	return x;
}
double* forsubs(int m, int n, double **L, double *b){
	if(m < n){
		fprintf(stderr, "The dimension m must be greater than or equal to n.\n");
		exit(1);
	}
	double *x;
	x = malloc(n*sizeof(double));
	x[0] = b[0];
	
	for(int i=1; i < n; i++){
		double sum=0.0;
		for(int j=0; j <= (i-1); j++){
			sum += L[i][j] * x[j];
		}
		x[i] = (b[i] - sum)/L[i][i];
	}
	
	return x;
}