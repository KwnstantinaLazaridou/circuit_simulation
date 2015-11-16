//Lazaridou Kwnstantina 915
//Tsakatika Ilitsa 929
//Xatzhvasiliou Eleni 1050

#include <stdio.h>
#include "lib.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include "lib_matrix.h"
#include "matrix_solution.h"
#include "lib_matrix_sparse.h"
#include "csparse.h"


int main(int argc, char *argv[]){

	char c;
	int n=0, m=0;

	init();

	FILE *fp=NULL;
	fp=fopen(argv[1],"r");			//Anoigma tou arxeiou
	if(fp==NULL){                           //Elegxos an to arxeio anoikse kanonika, alliws termatismos..
            printf("\nProblem opening file. Program terminated...\n");
            return -1;
        }		

//Diavasma xaraktira kai antistoixi leitourgia analoga me auton (Provlepsi gia grammi sxoliwn kai gia telos arxeiou)
	c=fgetc(fp);
	do{
		if(c=='V' || c=='v'){
			m++;
			creatVoltList(fp);
		}
		else if(c=='I' || c=='i'){
			creatAmberList(fp);
		}
		else if(c=='R' || c=='r'){
			creatResistanceList(fp);
		}
		else if(c=='C' || c=='c'){
			creatCapacitorList(fp);
		}
		else if(c=='L' || c=='l'){
                        m++;
			creatInductorList(fp);
		}
		else if(c=='D' || c=='d'){
			creatDiodeList(fp);
		}
		else if(c=='M' || c=='m'){
			creatMOSList(fp);
		}
		else if(c=='B' || c=='b'){
			creatBJTList(fp);
		}
		else if(c=='%'){ 
			c=fgetc(fp);
			while(c!='\n'&&(c!=EOF)){c=fgetc(fp);}/*MOVE TO NEXT LINE*/
		}
                else if(c=='.'){
			analysis(fp);
		}
		else{
			
		}
		if(c!=EOF){c=fgetc(fp);}
	}while(!feof(fp));

	fclose(fp);

        if(ground==0)
        {
	     printf("No ground node! Program terminated!");
             return -1;
        }
        else
        { 
             printVoltList(rootV);
             printAmperList(rootI);
             printResistanceList(rootR); 
             printCapacitorList(rootC);
             printInductorList(rootL);
             printDiodeList(rootD);
             printMosList(rootM);
             printBjttList(rootB);
       }
 
	n=count_nodes();
	printf("\n m=%d , n=%d \n\n", m, n);

	int size=0;
 
	size = (n-1)+m;

    if(sparse_option == 0){     
	gsl_matrix* pinakasA = gsl_matrix_calloc(size,size);
	gsl_vector* pinakasB = gsl_vector_calloc(size);
	
        fill_matrix_A(pinakasA, pinakasB, n);
	fill_matrix_B(pinakasB);
       
   	printf(" O pinakas A einai o: \n");
        print2DMatrix(pinakasA, size); 		
    	printf("\n O pinakas B einai o: \n");
    	print1DMatrix(pinakasB,size);
 
    	printf("\n");
        
        if(use_lu==1){
	    if(found_iter==1){
                call_bi_cg(pinakasA, pinakasB, size);
            } 
            else{
                luDecomp(pinakasA, pinakasB, size);
	    }
        }
        else if (use_cholesky==1){
            if (found_iter==1){
                 call_cg(pinakasA, pinakasB, size);
            }
            else{
                choleskyDecomp(pinakasA, pinakasB, size);
            }
        }
    }
    else {
	int i;

	cs* pinakasA_sparse = cs_spalloc(size,size,4*sparse_elements,1,1);

	double* pinakasB_sparse = (double *)calloc(size,sizeof(double));
	double* pinakasX_sparse = (double *)calloc(size,sizeof(double));	
	
	fill_sparse_matrix_A(pinakasA_sparse, pinakasB_sparse,n);
	cs* pinakasC_sparse = cs_compress(pinakasA_sparse);
	cs_spfree(pinakasA_sparse);
	cs_dupl(pinakasC_sparse);

 	if(use_lu==1){
	    if(found_iter==1){
                call_bi_cg_sparse(pinakasC_sparse, pinakasB_sparse, pinakasX_sparse, size);
            } 
            else{
                luDecomp_sparse(pinakasC_sparse, pinakasB_sparse, pinakasX_sparse, size);
	    }
        }
        else if (use_cholesky==1){
            if (found_iter==1){
                 call_cg_sparse(pinakasC_sparse, pinakasB_sparse, pinakasX_sparse, size);
            }
            else{
                choleskyDecomp_sparse(pinakasC_sparse, pinakasB_sparse, pinakasX_sparse, size);
            }
        }   
    }
    return 0;
}
