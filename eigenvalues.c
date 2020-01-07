/**
* Author: Dillon Bourne,
* Contact: dbourne@uoguelph.ca,
* Filename: eigenvalues.c
**/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "bla.c"
#include "qrhh3.h"
double sign(double x){
	if(x == 0.0){
		return (1.0);
	}
	else{
		return (x/fabs(x));
	}
}
void upphess(int m, double **A){
	//variables

	double *v = malloc(m*sizeof(double));
	double *x = malloc(m*sizeof(double));
	double *e = malloc(m*sizeof(double));
	double *vTimes2 = malloc(m*sizeof(double));
	double *vStarA = malloc(m*sizeof(double));
	double signNormX = 0.0;
	double **twoVStarA = (double **) malloc(m*sizeof(double *));
	//malloc columns of matrix
	for(int i = 0; i < m; i ++){
		twoVStarA[i] = (double *) malloc(m*sizeof(double));
	}	

	//fill e (1,0,0,0,...,)
	for(int i = 0; i < m-1; i++){
		(i == 0) ? (e[i]=1.0) : (e[i]=0.0);
	} 
	//for k=1 to n
	for(int k = 0; k < (m-1); k++){
		//x = A[k:m][k]
		for(int l = k+1; l < m; l++){
			x[l-(k+1)] = A[l][k];
		}
		//Vk = sign(x1)norm(x)e1 + x
		signNormX = sign(x[0]) * norm((m-(k+1)),x);
		for(int i=0; i < (m-k); i++){
			v[i] = (signNormX * e[i]) + x[i];
		}
		//Vk = Vk/norm(Vk)
		normalize((m-(k+1)),v);
		//2Vk
		for(int i = 0; i < (m-(k+1)); i++){
			vTimes2[i] = 2*v[i];
		}
		//(vk*A[k+1:m][k:m])
		for(int i=k; i < m; i++){
			vStarA[i-k] = inner_col((m-(k+1)),i,&A[k+1],v);
		}
		//2vk(vk*A[k+1:m][k:m])
		for(int i = 0; i < (m-(k+1)); i++){
			for(int j = 0; j < (m-k); j++){
				twoVStarA[i][j] = vTimes2[i] * vStarA[j];
			} 
		}	
		//A[k+1:m][k:m] -= twoVStarA[k+1:m][k:m]
		for(int i=0; i < (m-(k+1)); i++){
			for(int j=0; j < (m-k); j++){
				A[k+1+i][k+j] -= twoVStarA[i][j];
			}
		} 
		double twoATimesV[m];
		//2(A[1:m][k+1:m]*v)
		for(int i=0; i < m; i++){
			twoATimesV[i] = 0.0;
			for(int j=k+1; j < m; j++){
				twoATimesV[i] += A[i][j] * v[j-(k+1)];
			}
			twoATimesV[i] *= 2;
		}
		//2(A[1:m][k+1:m]*v)v*
		for(int i=0; i < m; i++){
			for(int j=k+1; j < m; j++){
				A[i][j] -= twoATimesV[i] * v[j-(k+1)];
			}
			//iteratively clear re-initialize array
			twoATimesV[i] = 0.0;
		}
	}
	//free up everything except A
	free(x);
	free(v);
	free(vStarA);
	free(e);
	free(vTimes2);
	//tried to free them in a loop but got core aborted error.
	free(twoVStarA);
	return;
}
int qralg(int n, double **A){
	//variables
	int k=0;		//number of iterations
	long double mu;		//Wilkinson Shift
	double B[2][2];
	double delta;
	double b2;
	while(fabs(A[n-1][n-2]) > 10e-14){
		B[0][0] = A[n-2][n-2];
		B[0][1] = A[n-2][n-1];
		B[1][0] = A[n-1][n-2];
		B[1][1] = A[n-1][n-1];
		
		delta = (B[0][0] - B[1][1]) / 2;
		b2 = B[0][1] * B[1][0];

		mu = (B[1][1] - (sign(delta) * (b2)) / (fabs(delta) + sqrt((delta*delta)+b2)));
		//mu = 0.0;
		if(!(isnan(mu))){
			//printf("mu is %.2LF..\n",mu);
			//subtract off mu each diagonal of A
			for(int i=0; i < (n); i++){
				A[i][i] -= mu;
			}
			//qr factorize A into product RQ
			qrhh3(n,A);
			//add back on mu to each diagonal
			for(int i=0; i < (n); i++){
				A[i][i] += mu;
			}
		}
		else{
			mu = (B[0][0] + B[1][1]) / 2;
			for(int i=0; i < (n); i++){
				A[i][i] -= mu;
			}
			//qr factorize A into product RQ
			qrhh3(n,A);
			//add back on mu to each diagonal
			for(int i=0; i < (n); i++){
				A[i][i] += mu;
			}
		}
		k++;
	}
	A[n-1][n-2] = 0;
	return k;
}
int *eigval(int m, double **A){
	
	int *itrCounts = malloc(m*sizeof(int));
	int tempCounts[m];
	
	upphess(m,A);
	
	for(int i=0; i < m; i++){
		tempCounts[i] = qralg(m-i,A);		
	}
	for(int i=m; i > 0; i--){
		itrCounts[m-i] = tempCounts[i-1];
	}
	return itrCounts;
	
}