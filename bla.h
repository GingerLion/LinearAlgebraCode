double inner(int m, double *a, double *b);
double norm(int m, double *a);
double normalize(int m, double *a);
double *normalize2(int m, double *a);
double *project(int m, double *a, double *b);
double *mat_vec_mult(int m, int n, double **A, double *x);
double inner_col(int m, int col, double **A, double *x);
double inner_col2(int m, int col1, double **A, int col2, double **B);
