/**
 * Author: Dillon Bourne,
 * Contact: dbourne@uoguelph.ca,
 * Filename: "givens.c"
**/
#include "givens.h"
#include "bla.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void givens(double *G, double *v){
	double temp = 0.0;
	temp = (G[0] * v[0]) + (G[1] * v[1]);
	v[1] = (-G[1] * v[0]) + (G[0] * v[1]);
	v[0] = temp;
}
