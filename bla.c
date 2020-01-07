#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bla.h"

double inner(int m, double *a, double *b) {
	/* This function computes the inner product of a and b. */
	int i;
	double val = 0.0;

	for (i=0; i<m; i++) val += a[i]*b[i];
	return val;
}

double norm(int m, double *a) {
	/* This function computes the 2-norm of a. */
	
	return sqrt(inner(m,a,a));
}

double normalize(int m, double *a) {
	/* This function normalizes a vector in place.  It returns the norm. */
	int i;
	double len;

	len = norm(m,a);
	if (len != 0.0)
		for (i=0; i<m; i++) a[i] /= len;
	return len;
}

double *normalize2(int m, double *a) {
	/* This function normalizes a vector allocating space. */
	int i;
	double len, *b;

	if ((b = (double *) calloc(m,sizeof(double))) == NULL) {
		fprintf(stderr,"calloc failed\n");
		exit(1);
	}
	len = norm(m,a);
	if (len != 0.0)
		for (i=0; i<m; i++) b[i] = a[i]/len;
	return b;
}

double *project(int m, double *a, double *b) {
	/* This function projects a onto b and returns the result. */
	int i;
	double coeff, *c;

	if ((c = (double *) calloc(m,sizeof(double))) == NULL) {
		fprintf(stderr,"calloc failed\n");
		exit(1);
		}
	coeff = inner(m,b,b);
	if (coeff != 0.0) {
		coeff = inner(m,b,a)/coeff;
		for (i=0; i<m; i++) c[i] = b[i]*coeff;
	}
	return c;
}

double *mat_vec_mult(int m, int n, double **A, double *x) {
	/* This function computes the matrix-vector product Ax. */
	int i;
	double *y;

	if ((y = (double *) malloc(m*sizeof(double))) == NULL) {
		fprintf(stderr,"malloc failed\n");
		exit(1);
	}
	for (i=0; i<m; i++) y[i] = inner(n,A[i],x);
	return y;
}

double inner_col(int m, int col, double **A, double *x) {
	/* This function computes the inner product of column col of A with vector x. */
	/* n is the number of rows of A and number of entries of x. */
	int i;
	double val = 0.0;

	for (i=0; i<m; i++){
	     val += A[i][col]*x[i];
		}
	return val;
}

double inner_col2(int m, int col1, double **A, int col2, double **B) {
	/* This function computes the inner product of column col1 of A with column col2 of B. */
	/* m is the number of rows of both A and B. */
	int i;
	double val = 0.0;

	for (i=0; i<m; i++) val += A[i][col1]*B[i][col2];
	return val;
}
