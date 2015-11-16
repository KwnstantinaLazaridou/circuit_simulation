#include "lib_matrix.h"
#include "lib.h"

void fill_matrix_B(gsl_vector* matrix)
{
        currI=rootI;
 
	while(currI != NULL){
  		if(currI->node1!=0) 
                    gsl_vector_set(matrix, (currI->node1)-1, gsl_vector_get(matrix, (currI->node1)-1) - currI->value);
		if(currI->node2!=0) 
                    gsl_vector_set(matrix, (currI->node2)-1, gsl_vector_get(matrix, (currI->node2)-1) + currI->value);
		currI= currI->next;
  	}
}

void fill_matrix_A(gsl_matrix* matrixA, gsl_vector* matrixB, int n)
{
    int b=1;
	currR=rootR;
	while (currR!=NULL) {
      	if (currR->node1!=0) 
            gsl_matrix_set(matrixA, (currR->node1)-1, (currR->node1)-1, gsl_matrix_get(matrixA, (currR->node1)-1, (currR->node1)-1) + 1/currR->value);
		if(currR->node2!=0) 
            gsl_matrix_set(matrixA, (currR->node2)-1, (currR->node2)-1, gsl_matrix_get(matrixA, (currR->node2)-1, (currR->node2)-1) + 1/currR->value);
		if(currR->node1!=0 && currR->node2!=0) {
		   	gsl_matrix_set(matrixA, (currR->node1)-1, (currR->node2)-1, gsl_matrix_get(matrixA, (currR->node1)-1, (currR->node2)-1) - 1/currR->value);
		    gsl_matrix_set(matrixA, (currR->node2)-1, (currR->node1)-1, gsl_matrix_get(matrixA, (currR->node2)-1, (currR->node1)-1) - 1/currR->value);
		}
		currR=currR->next;
	}
  
    printf("Number of nodes: %d\n", n);
    currV=rootV;
 	while(currV != NULL){
     
		if(currV->node1!=0) 
                    gsl_matrix_set(matrixA, n-2+b, (currV->node1)-1, 1.00);
		if(currV->node2!=0) 
                    gsl_matrix_set(matrixA, n-2+b, (currV->node2)-1, -1.00);
           
		if(currV->node1!=0) 
                    gsl_matrix_set(matrixA, (currV->node1)-1, n-2+b, 1.00); 
		if(currV->node2!=0) 
                    gsl_matrix_set(matrixA, (currV->node2)-1, n-2+b, -1.00);
		if((currV->node1!=0) || (currV->node2!=0))  
                    gsl_vector_set(matrixB, n-2+b, currV->value);
		b++;
		currV = currV->next;
    }

    currL=rootL;
 	while(currL != NULL){
     
		if(currL->node1!=0) 
                    gsl_matrix_set(matrixA, n-2+b, (currL->node1)-1, 1.00); 
		if(currL->node2!=0) 
                    gsl_matrix_set(matrixA, n-2+b, (currL->node2)-1, -1.00);
           
		if(currL->node1!=0) 
                    gsl_matrix_set(matrixA, (currL->node1)-1, n-2+b, 1.00);  
		if(currL->node2!=0) 
                    gsl_matrix_set(matrixA, (currL->node2)-1, n-2+b, -1.00); 
		if((currL->node1!=0) || (currL->node2!=0)) 
                    gsl_vector_set(matrixB, n-2+b, currL->value); 
 		b++;
		currL = currL->next;
    }

}

void print2DMatrix(gsl_matrix* matrix,int size)
{
    int i,j;

    for(i=0; i<size; i++){
       	for(j=0; j<size; j++){
            printf(" %lf ", gsl_matrix_get(matrix,i,j));
        }
       	printf("\n");
    }
}

void print1DMatrix(gsl_vector* matrix,int size)
{
    int i;

    for(i=0; i<size; i++){
       	printf(" %lf ", gsl_vector_get(matrix,i));
    }
    printf("\n");
}

