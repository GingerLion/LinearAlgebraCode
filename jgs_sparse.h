/**
* Author: Dillon Bourne,
* Contact: dbourne@uoguelph.ca,
* Filename: jgs_sparse.h
**/
#include "CSR_util.h"
double *jacobi_sparse(struct CSR *A, double *b);
double *gs_sparse(struct CSR *A, double *b);