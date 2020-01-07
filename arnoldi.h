#ifndef ARNOLDI_H
#define ARNOLDI_H

#include "CSR_util.h"

int arnoldi(struct CSR *A, int n, double ***QT, double ***HT, double *b);

#endif
