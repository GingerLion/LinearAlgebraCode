#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bla.h"
#include "arnoldi.h"

#define MAX(a,b)  ((a) >= (b) ? (a) : (b))

int arnoldi(struct CSR *A, int n, double ***QT, double ***HT, double *b) {
	/* Compute n more steps of the Arnoldi iteration.
	 * Input:
	 *	A:  a pointer to a CSR structure (matrix stored in CSR format)
	 *	n:  number of steps to iterate this time
	 *	QT: address of an array of row pointers to the matrix Q^T (Q transpose).  
	 *	HT: address of an array of row pointers to the matrix H^T (H transpose). 
	 *	b:  vector to initialize first col of Q
	 * Output:
	 *  N: number of full steps actually performed (may be less than n if
	 *     the iteration breaks down, that is, if the last subdiagonal 
	 *     entry of H is zero)
	 *
	 *	Each step produces a new column of Q and a new column of H, hence a
	 *	new row of *QT and a new row of *HT.  This function will (re-)allocate 
	 *	space in *HT and *QT to store the additional pointers, and will 
	 *	allocate space for each of these rows of *QT and for the part rows 
	 *	of *HT.  Since H is upper hessenberg, the nth column has zeros below 
	 *	the (n+1)st row; we do not bother to allocate space for or store 
	 *	these zeros. 
	 *
	 *	On the very first call to this function, *QT and *HT should be 
	 *	uninitialized.
	 *	*/

	int i,j,k;   double *v,TOL;
	static int n_rowptr = 0;
	static int n_iter = 0;

	if (n_rowptr == 0) {
		n_rowptr = MAX(10,3*n);  // allocate at least 10 or 3n pointers to start
		if (((*HT = (double **) malloc(n_rowptr*sizeof(double *))) == NULL) ||
			((*QT = (double **) malloc((n_rowptr+1)*sizeof(double *))) == NULL)) {
			fprintf(stderr,"malloc failed\n");
			exit(1);
		}
		(*QT)[0] = normalize2(A->m,b);
	}
	else if (n_rowptr < n_iter + n) {
		n_rowptr += MAX(10,3*n);  // allocate at least 10 or 3n more pointers
		/* reallocate space for row pointers in QT and HT. */
		if (((*HT = (double **) realloc(*HT,n_rowptr*sizeof(double *))) == NULL) ||
			((*QT = (double **) realloc(*QT,(n_rowptr+1)*sizeof(double *))) == NULL)) {
			fprintf(stderr,"realloc failed\n");
			exit(1);
		}
	}
	/* Set the tolerance for when a subdiagonal entry of H is considered 
	 * zero. */
	TOL = A->m*1.0e-15;
	/* Do n more steps of Arnoldi starting at n_iter. */
	for (k=n_iter; k<n_iter+n; k++) {
		v = CSRmult(A, (*QT)[k]);
		/* Allocate space for new columns of Q and H. */
		if ((((*QT)[k+1] = (double *) calloc(A->m,sizeof(double))) == NULL) ||
		    (((*HT)[k] = (double *) malloc((k+2)*sizeof(double))) == NULL)) {
			fprintf(stderr,"malloc failed\n");
			exit(1);
		}
		/* Obtain column k of H */
		for (j=0; j<=k; j++) {
			(*HT)[k][j] = inner(A->m,(*QT)[j],v);
			for (i=0; i<A->m; i++) v[i] -= (*HT)[k][j] * (*QT)[j][i];
		}
		(*HT)[k][k+1] = norm(A->m, v);
		/* Check for breakdown. */
		if (fabs((*HT)[k][k+1]) < TOL) {k++; break;}
		/* Set column k+1 of Q */
		for (i=0; i<A->m; i++) (*QT)[k+1][i] = v[i]/(*HT)[k][k+1];
	}
	free(v);  /* v was allocated by CSRmult */
	/* Number of iterations successfully performed this call: */
	k -= n_iter;
	/* Update n_iter */
	n_iter += k;
	/* return number of successful iterations */
	return k;
}
