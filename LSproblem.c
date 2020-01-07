/* 
Author: Dillon Bourne
Contact: dbourne@uoguelph.ca
Filename: LSproblem.c 
*/
#include "bla.h"
#include "qr_factorize.h"
#include "LSproblem.h"
#include "forbacksubs.h"
#define M_PI 3.14159265358979323846
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#define M 3

struct AB *AB_construct(int n){
	struct AB *ab;
	ab = malloc(sizeof(struct AB));
	double *x = malloc(n*sizeof(double));
	ab->A = (double **) malloc(n*sizeof(double*));
	for(int i=0; i < n; i++){
		ab->A[i] = (double *) malloc(M*sizeof(double));
	}
	ab->b = malloc(n*sizeof(double));
	
	x[0] = 0.0;
	for(int i=1; i < n; i++){
		x[i] = ((double)i/(double)n);
	}
	for(int i=0; i < n; i++){
		for(int j=0; j < M; j++){
			if(j == 0){
				ab->A[i][j] = exp(x[i]);
			}
			else if(j == 1){
				ab->A[i][j] = sqrt(x[i]);
			}
			else{
				ab->A[i][j] = log(x[i] + 1);
			}
		}
	}
	for(int i=0; i < n; i++){
		ab->b[i] = cos(M_PI*x[i]);
	}

	return ab;
}
double LS_solve(int n){
	struct AB *abSolve;
	abSolve = malloc(sizeof(struct AB));
	double *c;
	c = malloc(n*sizeof(double));
	abSolve = AB_construct(n);
	qrhh2(n,M,abSolve->A,abSolve->b);
	c = backsubs(n,M,abSolve->A,abSolve->b);
} 