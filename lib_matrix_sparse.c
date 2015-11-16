#include "lib_matrix_sparse.h"
#include "lib.h"

void fill_sparse_matrix_B(double* matrix)
{
        currI=rootI;
 
	while(currI != NULL){
  		if(currI->node1!=0) 
                    matrix[currI->node1-1]-=currI->value;
		if(currI->node2!=0) 
                    matrix[currI->node2-1]-=currI->value;
		currI= currI->next;
  	}
}

void fill_sparse_matrix_A(cs* matrixA, double* matrixB, int n)
{
        int b=1;
	int i=0;

	currR=rootR;
     	while(currR!=NULL){
      		if(currR->node1!=0){
		    cs_entry(matrixA, currR->node1-1, currR->node1-1, 1/currR->value); 
		}
		if(currR->node2!=0){ 
		    cs_entry(matrixA, currR->node2-1, currR->node2-1, 1/currR->value); 

                }   
		if(currR->node1!=0 && currR->node2!=0) { 
		    cs_entry(matrixA, currR->node1-1, currR->node2-1, -1/currR->value); 
		    cs_entry(matrixA, currR->node2-1, currR->node1-1, -1/currR->value); 
		}
		currR=currR->next;
	}

    currV=rootV;
 	while (currV != NULL) {
		if (currV->node1!=0)
            cs_entry(matrixA, n-2+b, (currV->node1)-1, 1.00);
        if (currV->node2!=0)
            cs_entry(matrixA, n-2+b, (currV->node2)-1, -1.00);
        if (currV->node1!=0)
            cs_entry(matrixA, (currV->node1)-1, n-2+b, 1.00);
        if (currV->node2!=0)
            cs_entry(matrixA, (currV->node2)-1, n-2+b, -1.00);
		if ((currV->node1!=0) || (currV->node2!=0)){
            matrixB[n-2+b]=currV->value;
		}
		b++;
		currV = currV->next;
    }

    currL=rootL;
 	while(currL != NULL){
     
		if(currL->node1!=0){ 
            matrixA->i[i]=n-2+b;
		    matrixA->p[i]=currL->node1-1;
		    matrixA->x[i]=1.000;
		    i++;
		    matrixA->i[i]=currL->node1-1;
		    matrixA->p[i]=n-2+b;
		    matrixA->x[i]=1.000;
		    i++;
		}
		if(currL->node2!=0){ 
                    matrixA->i[i]=n-2+b;
		    matrixA->p[i]=currL->node2-1;
		    matrixA->x[i]=-1.000;
		    i++;
		    matrixA->i[i]=currL->node2-1;
		    matrixA->p[i]=n-2+b;
		    matrixA->x[i]=-1.000;
		    i++;
           	}
		
 		b++;
		currL = currL->next;
        }

}

