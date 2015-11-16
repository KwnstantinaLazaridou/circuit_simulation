#include "csparse.h"

int sparse_elements;
void fill_sparse_matrix_B(double* matrix);
void fill_sparse_matrix_A(cs* matrixA, double* matrixB, int n);
//css *S;				// Pinakas S gia tin LU kai Cholesky
//csn *N;				// Pinakas N gia tin LU kai Cholesky

//double *B_sparse;			//dianisma B ( double, mege8ous : [(n-1)+m2] x 1 ) 	

//double *x_sparse;		//dianisma x
