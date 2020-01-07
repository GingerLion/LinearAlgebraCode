/**
* Author: Dillon Bourne,
* Contact: dbourne@uoguelph.ca,
* Filename: palu.c
**/
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "palu.h"

int * palu(int m, double **A){
	//permutation vector
	int *p;
	//scale-factor vector
	int *s;
	
	int biggest = 0;
	double rowMax = 0.0;
	double scaleMax = 0.0;
	//malloc permutation vector
	p = malloc(m*sizeof(int));
	//malloc scaling-factor vector
	s = malloc(m*sizeof(int));
	//initialize permutation vector [0,1,2,...,m]
	for(int i=0; i < m; i++){
		p[i] = i;
	}
	//loop through matrix A to find max-absolute scale-factor
	for(int i=0; i < m; i++){
		for(int j=0; j < m; j++){
			if(A[i][j] > rowMax){
				rowMax = fabs(A[i][j]);
			}
		}
		s[i] = rowMax;
		rowMax = 0.0;
	}
	for(int k=0; k < (m-1); k++){
		
		//for loop variables and initialization
		scaleMax = 0.0;
		double tmp = 0.0;
		double multiplier = 0.0;
		double **temp;
		temp = (double **) malloc(m*sizeof(double*));
		for(int i=0; i < m; i++){
			temp[i] = (double *) malloc(m*sizeof(double));
		}

		//find pivot for column k
		for(int i=k; i < m; i++){
			if((A[i][k]/(double)s[i]) > scaleMax){
				scaleMax = (A[i][k]/(double)s[i]);
				biggest = i;
			}
		}
	
		//use biggest & k as indices to swap elements in permutation vector
		if(!(biggest == k)){
			//swap elements of permutation vector
			tmp = p[k];
			p[k] = p[biggest];
			p[biggest] = tmp;
			//swap rows of A (U)
			temp[0] = A[k];
			A[k] = NULL;
			A[k] = A[biggest];
			A[biggest] = temp[0];

		}

		//perform Gaussian Elimination on A
		for(int j=k+1; j < m; j++){
			
			multiplier = A[j][k] / A[k][k];

			//[k:m]
			for(int n=k; n < m; n++){ 
				A[j][n] -= multiplier * A[k][n];
			}

			A[j][k] = multiplier;
		
		} 
	}

	return p; 

}
