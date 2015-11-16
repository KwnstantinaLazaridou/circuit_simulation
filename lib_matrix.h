#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_linalg.h>

void fill_matrix_B(gsl_vector* matrix);
void fill_matrix_A(gsl_matrix* matrixA, gsl_vector* matrixB,int n);
void print2DMatrix(gsl_matrix* matrix,int size);
void print1DMatrix(gsl_vector* matrix,int size);
 
