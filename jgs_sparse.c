/**
* Author: Dillon Bourne,
* Contact: dbourne@uoguelph.ca,
* Filename: jgs_sparse.c
**/
#include "bla.h"
#include "jgs_sparse.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double *jacobi_sparse(struct CSR *A, double *b){
	//variable declarations & initializations
	int nonZeroDiag = 0;
	double **x;
	double *sub;
	x = malloc((A->m+1)*sizeof(double*));
	sub = (double *) malloc(A->n*sizeof(double));
	//set x(row0) to all 0s
	x[0] = calloc(A->n,A->n*sizeof(double));
	
	for(int i=1; i < A->m+1; i++){
		x[i] = calloc(A->n,A->n*sizeof(double));
	}
	//check if all diagonal elements are non-zero
	for(int i=0; i < A->m; i++){
		//if diag element is non-zero 
		if(!(A->V[i*(A->n+1)] == 0)){
			nonZeroDiag++;
		}
	}
	//exit program if any diagonal elements are 0
	if(!(nonZeroDiag == A->m)){
		printf("Error: the diagonal entries of the matrix must be non-zero. \n");
		exit(1);
	}
	int j=0;
	double sum = 0.0;
	int k;
	for(k=1; k <= A->m; k++){
		for(int i=0; i < A->m; i++){
			for(j=A->I[i]; j < A->I[i+1]; j++){
				if(j == (i*(A->n+1))) { continue; }
				sum += A->V[j] * x[k-1][A->J[j]];
			}
			//(b[i] - sum)/A[i][i]
			x[k][i] = (b[i] - sum) / A->V[i*(A->n+1)];
			sum = 0.0;
		}
		for(int p=1; p <= A->n; p++){
			sub[p-1] = x[p] - x[p-1];
		}
		if((norm(A->n,sub))/(norm(A->n,x[k])) < 10e-14){
			printf("Program closed: norm condition triggered.\n");
			exit(1);
		}
	}
	free(sub);
	return x[k-1];
}
double *gs_sparse(struct CSR *A, double *b){
	//variable declarations & initializations
	int nonZeroDiag = 0;
	double **x;
	double *sub;
	x = malloc((A->m+1)*sizeof(double*));
	sub = (double *) malloc(A->n*sizeof(double));
	//set x(row0) to all 0s
	x[0] = calloc(A->n,A->n*sizeof(double));
	
	for(int i=1; i < A->m+1; i++){
		x[i] = calloc(A->n,A->n*sizeof(double));
	}
	//check if all diagonal elements are non-zero
	for(int i=0; i < A->m; i++){
		//if diag element is non-zero 
		if(!(A->V[i*(A->n+1)] == 0)){
			nonZeroDiag++;
		}
	}
	//exit program if any diagonal elements are 0
	if(!(nonZeroDiag == A->m)){
		printf("Error: the diagonal entries of the matrix must be non-zero. \n");
		exit(1);
	}
	int j=0;
	double sum1 = 0.0;
	double sum2 = 0.0;
	int k = 0;
	for(k=1; k <= A->m; k++){
		for(int i=0; i < A->m; i++){
			for(j=A->I[i]; j < A->I[i+1]; j++){
				if(A->J[j] < i) {
				sum1 += A->V[j] * x[k][A->J[j]];
				}
			}
			
			for(j=A->I[i]; j < A->I[i+1]; j++){
				if(A->J[j] > i) {
				sum2 += A->V[j] * x[k-1][A->J[j]];
				}
			}
			
			//(b[i] - sum)/A[i][i]
			x[k][i] = (b[i] - sum1 - sum2) / A->V[i*(A->n+1)];
			sum1 = 0.0;
			sum2 = 0.0;
		}
		for(int p=1; p <= A->n; p++){
			sub[p-1] = x[p] - x[p-1];
		}
		if((norm(A->n,sub))/(norm(A->n,x[k])) < 10e-14){
			printf("Program closed: norm condition triggered.\n");
			exit(1);
		}
	}
	free(sub);
	return x[k-1];
}
