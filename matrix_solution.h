#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_linalg.h>
#include "lib_matrix.h"

void luDecomp(gsl_matrix* matrixA, gsl_vector* matrixB, int size);
void choleskyDecomp(gsl_matrix* matrixA, gsl_vector* matrixB, int size);
void call_bi_cg(gsl_matrix* matrixA, gsl_vector* matrixB, int size);
void call_cg(gsl_matrix* matrixA, gsl_vector* matrixB, int size);
