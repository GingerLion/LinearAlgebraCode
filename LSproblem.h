/* 
Author: Dillon Bourne
Contact: dbourne@uoguelph.ca
Filename: LSproblem.h 
*/
struct AB {
	double **A;
	double *b;
};
struct AB *AB_construct(int n);
double LS_solve(int n);