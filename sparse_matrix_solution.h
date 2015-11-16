#include "lib_matrix.h"

void luDecomp_sparse(cs* matrixA, double* matrixB, double* matrixX, int size);
void choleskyDecomp_spearse(cs* matrixA, double* matrixB, double* matrixX, int size);
void call_bi_cg_sparse(cs* matrixA, double* matrixB,double* matrixX, int size);
void call_cg_sparse(cs* matrixA, double* matrixB,double* matrixX, int size);
