/* 
Author: Dillon Bourne
Contact: dbourne@uoguelph.ca
Filename: qr_factorize.c
*/
#include "bla.h"
#include<stdlib.h>
#include<stdio.h>
#include "qr_factorize.h"
#include "math.h"

double sign(double x){
	if(x == 0.0){
		return (0.0);
	}
	else{
		return (x/fabs(x));
	}
}
void qrhh1(int m,int n, double **A){
	
	if(m < n){
		printf("Error: m must be greater than n.\n");
		exit(1);
	}
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
		twoVStarA[i] = (double *) malloc(n*sizeof(double));
	}	

	//fill e (1,0,0,0,...,)
	for(int i = 0; i < m; i++){
		(i == 0) ? (e[i]=1.0) : (e[i]=0.0);
	} 
	//for k=1 to n
	for(int k = 0; k < (n-1); k++){
		//x = A[k:m][k]
		for(int l = k; l < m; l++){
			x[l-k] = A[l][k];
		}
		//Vk = sign(x1)norm(x)e1 + x
		signNormX = sign(x[0]) * norm((m-k),x);
		for(int i=0; i < (m-k); i++){
			v[i] = (signNormX * e[i]) + x[i];
		}
		//Vk = Vk/norm(Vk)
		normalize((m-k),v);
		//2Vk
		for(int i = 0; i < (m-k); i++){
			vTimes2[i] = 2*v[i];
		}
		//(vk*A[k:m][k:n])
		printf("Here is v*A at k=%d...\n",k);
		for(int i=k; i < n; i++){
			vStarA[i-k] = inner_col((m-k),i,&A[k],v);
			printf("%.2lf ",vStarA[i]);
		}
		//2vk(vk*A[k:m][k:n])
		printf("\nHere is 2v(v*A) at k=%d...\n",k);
		for(int i = 0; i < (m-k); i++){
			for(int j = 0; j < (n-k); j++){
				twoVStarA[i][j] = vTimes2[i] * vStarA[j];
				printf("%.2lf ",twoVStarA[i][j]);
			} 
			printf("\n");
		}	
		//A[k:m][k:n] -= twoVStarA[k:m][k:n]
		for(int i=0; i < (m-k); i++){
			for(int j=0; j < (n-k); j++){
				A[k+i][k+j] -= twoVStarA[i][j];
			}
		} 	
	}
	return;
}
void qrhh2(int m,int n, double **A, double *b){
	if(m < n){
		printf("Error: m must be greater than n.\n");
		exit(1);
	}
	//variables
	double *v = malloc(m*sizeof(double));
	double *x = malloc(m*sizeof(double));
	double *e = malloc(m*sizeof(double));
	double *vStarA = malloc(m*sizeof(double));
	double *vTimes2 = malloc(m*sizeof(double));
	double vStarb = 0.0;
	double signNormX= 0.0;
	double *twoVStarb = malloc(m*sizeof(double));
	double **twoVStarA = (double **) malloc(m*sizeof(double *));
	//malloc columns of matrix
	for(int i = 0; i < m; i ++){
		twoVStarA[i] = (double *) malloc(n*sizeof(double));
	}
	//fill e (1,0,0,0,...,)
	for(int i = 0; i < m; i++){
		(i == 0) ? (e[i]=1.0) : (e[i]=0.0);
	} 
	//for k=1 to n
	for(int k = 0; k < (n-1); k++){
		//x = A[k:m][k]
		for(int l = k; l < m; l++){
			x[l-k] = A[l][k];
		}
		//Vk = sign(x1)norm(x)e1 + x
		signNormX = sign(x[0]) * norm((m-k),x);
		
		for(int i=0; i < (m-k); i++){
			v[i] = (signNormX * e[i]) + x[i];
		}
		//Vk = Vk/norm(Vk)
		normalize((m-k),v);
		//2Vk
		for(int i = 0; i < (m-k); i++){
			vTimes2[i] = 2*v[i];
		}
		//(vk*A[k:m][k:n])
		for(int i=k; i < m; i++){
			vStarA[i-k] = inner_col((m-k),i,&A[k],v);
		}
		//(vk*b[k:m])
		vStarb = inner((m-k),v,b);
		//2vk(vk*b[k:m])
		for(int i=0; i < (m-k); i++){
			twoVStarb[i] = vTimes2[i] * vStarb;
			for(int j = 0; j < (n-k); j++){
				//2vk(vk*A[k:m][k:n])
				twoVStarA[i][j] = vTimes2[i] * vStarA[j];
			} 
		}	
		//A[k:m][k:n] -= twoVStarA[k:m][k:n]
		for(int i=0; i < (m-k); i++){
			b[i] -= twoVStarb[i];
			for(int j=0; j < (n-k); j++){
				A[k+i][k+j] -= twoVStarA[i][j];	
			}
		} 
	}
	return;
}
double** mgs1(int m,int n,double **A){
	
	double **r;
	double **v;
	r = (double **) malloc(n*sizeof(double *));
	v = (double **) malloc(n*sizeof(double *));
	for(int i=0; i < n; i++){
		r[i] = (double *) malloc(n*sizeof(double));
		}
	for(int i=0;i < m; i++){
		v[i] = (double *) malloc(m*sizeof(double));
	}

	for(int i=0; i < n; i++){
		for(int j=0; j < m; j++){
			v[j][i] = A[j][i];
		}
	}

	for(int i=0; i < n; i++){
		double *vColNow;
		double *vColNext;
		double *q;
		vColNow = malloc(m*sizeof(double));
		vColNext = malloc(m*sizeof(double));
		q = malloc(m*sizeof(double));
		
		double vNorm = 0.0;
		for(int j=0; j < m; j++){
			vColNow[j] = v[j][i];
			vNorm += (vColNow[j]*vColNow[j]);
		}
		vNorm = sqrt(vNorm);
		//r[i][i] = ||v||
		r[i][i] = vNorm;
		for(int j=0; j < m; j++){
			//qi = vi / rii
			A[j][i] = vColNow[j] / r[i][i];
			q[j] = A[j][i];
		}
		if((i+1) != n){
			for(int j=0; j < m; j++){
				vColNext[j] = v[j][i+1];
			}
			for(int j=(i+1); j < n; j++){
				//rij = qi*vj
				r[i][j] = inner(n,q,vColNext);
				
				//vj = vj - (rij x qi)
				for(int k = 0; k < m; k++){
					v[k][j] -= (r[i][j] * q[k]);
				}
			}
		}
	}
	return r;
	
}