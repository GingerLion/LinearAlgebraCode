/**
 * Author: Dillon Bourne,
 * Contact: dbourne@uoguelph.ca,
 * Filename: "qrgruhT.c"
**/
#include "bla.h"
#include "givens.h"
#include "qrgruhT.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void qrgruhT(int n, double **HT, int ndone, double **G){
	
	double *v = malloc(2*sizeof(double));
	double r = 0.0;
	double c = 0.0;
	double s = 0.0;
	
	for(int i=ndone; i < n; i++){ //goes up to n but can still grab [1][n,n+1] cols
		
		v = calloc(2,sizeof(double));
		//find 1x2 elements in HT
		v[0] = HT[i][i];
		v[1] = HT[i][i+1];
		//find givens rotation values
		r = norm(2,v);
		c = v[0] / r;
		s = v[1] / r;
		
		G[i][0] = c;
		G[i][1] = s;
		
		//apply givens rotations to 0-out bottom entry of v
		
		givens(G[i],v);
		
		//update HT entries with v values
		HT[i][i] = v[0];
		HT[i][i+1] = v[1];
		
		//multiply the rest of the 1x2 row subsets of HT by the same givens rotation
		for(int j=i+1; j < n; j++){
			v[0] = HT[j][i];
			v[1] = HT[j][i+1];
			givens(G[i],v);
			HT[j][i] = v[0];
			HT[j][i+1] = v[1];
			v = calloc(2,sizeof(double));
		}
		
	}
	free(v);
}