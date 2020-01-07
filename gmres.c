/**
 * Author: Dillon Bourne,
 * Contact: dbourne@uoguelph.ca,
 * Filename: "gmres.c"
**/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "forbacksubs.h"
#include "CSR_util.h"
#include "arnoldi.h"
#include "qrgruhT.h"
#include "bla.h"
#include "givens.h"
#define N 50
double *GMRES(struct CSR *A, double *b){
	
	double *q;
	double *x;
	double *y;
	double *be;
	double *v;
	double **G;
	double **QT;
	double **HT;
	double **Q;
	int a=0;
	
	Q = (double **) malloc(A->m*sizeof(double*));
	q = (double *) malloc(A->m*sizeof(double));
	x = (double *) malloc(A->m*sizeof(double));
	y = (double *) malloc(A->m*sizeof(double));
	G = (double **) malloc(N*sizeof(double*));
	be = (double *) malloc(A->m*sizeof(double));
	v = (double *) malloc(2*sizeof(double));
	for(int i=0; i < N; i++){
			G[i] = (double *) malloc(2*sizeof(double));
	}
	for(int i=0; i < A->m; i++){
		q[i] = b[i];
		Q[i] = (double *) malloc(A->m*sizeof(double *));
	}
	normalize(A->m,q);
	be[0] = norm(A->m,b);
	while(1){
		for(int n=1; n <= A->m; n++){
			a = arnoldi(A,1,&QT,&HT,q);
			qrgruhT(n,HT,n-1,G);
			//construct ||b||e
			be[n] = 0.0;		
			//perform same givens rotations from qrgruhT on be
			for(int p=n-1; p < n; p++){
				v[0]=be[p];
				v[1]=be[p+1];
				givens(G[p],v);
				be[p] = v[0];
				be[p+1] = v[1];
				v = calloc(2,sizeof(double));
				}
			//HTy = QT||b||e (Q here is not arnoldi Q)
			//find Ry=||b||e
			y = backsubs(n,n,HT,be);
			
		}
		if(fabs(be[A->m]) <= 10e-06){break;}
	}
	//xn = Qny (Arnoldi Q)
	//convert QT into Q
	for(int i=0; i < A->m; i++){
		for(int j=0; j < A->m; j++){
			Q[i][j] = QT[j][i];
		}
	}
	x = mat_vec_mult(A->m,A->m,Q,y);
	return x;
	free(q);free(y);free(be);
	for(int i=0; i < A->m; i++){ free(G[i]); free(Q[i]); free(Q[i]); free(HT[i]); free(QT[i]);}
	free(G);free(Q);free(QT);free(HT);
}