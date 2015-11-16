#include <string.h>
#include "csparse.h"
#include "lib_matrix.h"
#include "lib.h"
#include "sparse_matrix_solution.h"
#include <math.h>

void luDecomp_sparse(cs* matrixA, double* matrixB, double* matrixX, int size)
{
    int i;
    css *S;
    csn *N;

    S=cs_sqr(2,matrixA,0);
    N=cs_lu(matrixA,S,1);
    cs_spfree(matrixA);

    if(found_dc_sweep==0){
        cs_ipvec(N->pinv, matrixB, matrixX, size);
	cs_lsolve(N->L, matrixX);
	cs_usolve(N->U, matrixX);
	cs_ipvec(S->q, matrixX, matrixB, size);

	printf("X vector \n");
	for(i=0;i<size;i++){
	    printf(" %.6lf ",matrixX[i]);
	}
	printf("\n");
    }
    else{
	double value;
        double *matrixB_temp = (double *)calloc(size,sizeof(double));
	if (source > -1) {
            for(value=start_value;value<=end_value;value=value+step){
                matrixB[source-1]=value;
		cs_ipvec(N->pinv, matrixB, matrixX, size);
		cs_lsolve(N->L, matrixX);
		cs_usolve(N->U, matrixX);
       	    }
        }
        else { 
	    	if (sweep_node1!=0){
				matrixB[sweep_node1-1]+=sweep_value-start_value;
	    	}
	    	if(sweep_node2!=0){
				matrixB[sweep_node2-1]-=sweep_value+start_value;
	    	}
	    	for(value=start_value;value<=end_value;value=value+step) {
				cs_ipvec(N->pinv, matrixB, matrixX, size);
				cs_lsolve(N->L, matrixX);
				cs_usolve(N->U, matrixX);
				cs_ipvec(S->q, matrixX, matrixB_temp, size);
	   			if (sweep_node1!=0){
			    	matrixB[sweep_node1-1]+=sweep_value-start_value;
	    		}
	    		if(sweep_node2!=0){
		    		matrixB[sweep_node2-1]-=sweep_value+start_value;
	    		}
	    
				printf("value= %lf Matrix X: \n",value);
            	for(i=0;i<size;i++){
	            	printf(" %.6lf ",matrixX[i]);
	    		}
				printf("\n");
	    	} 
		}
    } 
    printf("\n");

}

void choleskyDecomp_sparse(cs* matrixA, double* matrixB, double* matrixX,int size)
{
    int i;
    css *S;
    csn *N;
 	
    S=cs_schol(1, matrixA);
    N=cs_chol(matrixA, S);
    cs_spfree(matrixA);

    if(found_dc_sweep==0){
        cs_ipvec(S->pinv, matrixB, matrixX, size);
        cs_lsolve(N->L, matrixX);
	cs_ltsolve(N->L, matrixX);
	cs_pvec(S->pinv, matrixX, matrixB, size);

	printf("X vector \n");
	for(i=0;i<size;i++){
	    printf(" %.6lf ",matrixX[i]);
	}
	printf("\n");
    }
    else{
        double value;
		double *matrixB_temp = (double *)calloc(size,sizeof(double));
		if (source > -1) {
            for(value=start_value;value<=end_value;value=value+step){
                matrixB[source-1]=value;
				cs_ipvec(S->pinv, matrixB, matrixX, size);
				cs_lsolve(N->L, matrixX);
				cs_ltsolve(N->L, matrixX);
				cs_pvec(S->pinv, matrixX, matrixB_temp, size);
		
				printf("value= %lf Matrix X: \n",value);
                for(i=0;i<size;i++){
	    	    	printf(" %.6lf ",matrixX[i]);
				}
				printf("\n");
       	    }
        }
        else { 
	    	if (sweep_node1!=0){
		 		matrixB[sweep_node1-1]+=sweep_value-start_value;
	    	}
	    	if(sweep_node2!=0){
				matrixB[sweep_node2-1]-=sweep_value-start_value;
	    	}
	    	for(value=start_value;value<=end_value;value=value+step) {
				cs_ipvec(S->pinv, matrixB, matrixX, size);
				cs_lsolve(N->L, matrixX);
				cs_ltsolve(N->L, matrixX);
				cs_pvec(S->pinv, matrixX, matrixB_temp, size);
	        
				if (sweep_node1!=0){
		    		matrixB[sweep_node1-1]-=step;
	       	 	}
	        	if(sweep_node2!=0){
		    		matrixB[sweep_node2-1]+=step;
	        	}
				printf("value= %lf Matrix X: \n",value);
                for(i=0;i<size;i++){
			   	    printf(" %.6lf ",matrixX[i]);
				}
				printf("\n");
		    } 
		}
    }                
    printf("\n");      
}

void conjugate_gradient_sparse(cs *A, double *b, double* x, int n, double itol)
{  
    int i,j;
    int iter;
    double rho,rho1,alpha,beta,omega;
     
    double r[n];
    double z[n];
    double q[n], temp_q[n];
    double p[n], temp_p[n];
    double res[n];
    double precond[n];  //Preconditioner
     
    memset(precond, 0, n*sizeof(double));
    memset(r, 0, n*sizeof(double));
    memset(z, 0, n*sizeof(double));
    memset(q, 0, n*sizeof(double));
    memset(temp_q, 0, n*sizeof(double));
    memset(p, 0, n*sizeof(double));
    memset(temp_p, 0, n*sizeof(double));
 
    /* Preconditioner */
    double max;
    int pp;
    for(j = 0; j < n; ++j){
        for(pp = A->p[j], max = fabs(A->x[pp]); pp < A->p[j+1]; pp++)
            if(fabs(A->x[pp]) > max)                  //vriskei to diagonio stoixeio
                max = fabs(A->x[pp]);
        precond[j] = 1/max;    
    }  
 
    cblas_dcopy (n, x, 1, res, 1);
 
    //r=b-Ax
    cblas_dcopy (n, b, 1, r, 1);
    memset(p, 0, n*sizeof(double));
    cs_gaxpy (A, x, p);
    for(i=0;i<n;i++){
        r[i]=r[i]-p[i];
     
    }
     
    double r_norm = cblas_dnrm2 (n, r, 1);
    double b_norm = cblas_dnrm2 (n, b, 1);
    if(!b_norm)
        b_norm = 1;
 
    iter = 0;  
     
    while( r_norm/b_norm > itol && iter < n )
    {
        iter++;
 
        cblas_dcopy (n, r, 1, z, 1);                //gia na min allaksei o r
         
        for(i=0;i<n;i++){
            z[i]=precond[i]*z[i];
     
        }
 
        rho = cblas_ddot (n, z, 1, r, 1);
        if (fpclassify(fabs(rho)) == FP_ZERO){
            printf("RHO aborting CG due to EPS...\n");
            exit(42);
        }
 
        if (iter == 1){
            cblas_dcopy (n, z, 1, p, 1);
        }
        else{      
            beta = rho/rho1;
     
            //p = z + beta*p;
            cblas_dscal (n, beta, p, 1);    //rescale
            cblas_daxpy (n, 1, z, 1, p, 1); //p = 1*z + p
             
        }      
        rho1 = rho;
         
        //q = Ap
        memset(q, 0, n*sizeof(double));
        cs_gaxpy (A, p, q);
 
        omega = cblas_ddot (n, p, 1, q, 1);
        if (fpclassify(fabs(omega)) == FP_ZERO){
            printf("OMEGA aborting CG due to EPS...\n");
            exit(42);
        }
 
        alpha = rho/omega; 
 
        //x = x + aplha*p;
        cblas_dcopy (n, p, 1, temp_p, 1);
        cblas_dscal (n, alpha, temp_p, 1);//rescale by alpha
        cblas_daxpy (n, 1, temp_p, 1, res, 1);// sum x = 1*x + temp_p
 
        //r = r - aplha*q;
        cblas_dcopy (n, q, 1, temp_q, 1);
        cblas_dscal (n, -alpha, temp_q, 1);//rescale by alpha
        cblas_daxpy (n, 1, temp_q, 1, r, 1);// sum r = 1*r - temp_p
 
        //next step
        r_norm = cblas_dnrm2 (n, r, 1);
    }
    cblas_dcopy (n, res, 1, x, 1);
 
}
 
void bi_conjugate_gradient_sparse(cs *A, double *b, double* x, int n, double itol){
   
    int i,j,iter;
     
    double rho,rho1,alpha,beta,omega;
     
    double r[n], r_t[n];
    double z[n], z_t[n];
    double q[n], q_t[n], temp_q[n];
    double p[n], p_t[n], temp_p[n];
    double res[n];                  //NA VGEI!
    double precond[n];
     
    //Initializations      
    memset(precond, 0, n*sizeof(double));
    memset(r, 0, n*sizeof(double));
    memset(r_t, 0, n*sizeof(double));
    memset(z, 0, n*sizeof(double));
    memset(z_t, 0, n*sizeof(double));
    memset(q, 0, n*sizeof(double));
    memset(q_t, 0, n*sizeof(double));
    memset(temp_q, 0, n*sizeof(double));
    memset(p, 0, n*sizeof(double));
    memset(p_t, 0, n*sizeof(double));
    memset(temp_p, 0, n*sizeof(double));
    memset(res, 0, n*sizeof(double));
     
    /* Preconditioner */
    double max;
    int pp;
    for(j = 0; j < n; ++j){
        for(pp = A->p[j], max = fabs(A->x[pp]); pp < A->p[j+1]; pp++)
            if(fabs(A->x[pp]) > max)                  //vriskei to diagonio stoixeio
                max = fabs(A->x[pp]);
        precond[j] = 1/max;    
    }  
    cs *AT = cs_transpose (A, 1) ;
 
    cblas_dcopy (n, x, 1, res, 1);
 
    //r=b-Ax
    cblas_dcopy (n, b, 1, r, 1);
    memset(p, 0, n*sizeof(double));
    cs_gaxpy (A, x, p);
    for(i=0;i<n;i++){
        r[i]=r[i]-p[i];
     
    }
     
    cblas_dcopy (n, r, 1, r_t, 1);
     
    double r_norm = cblas_dnrm2 (n, r, 1);
    double b_norm = cblas_dnrm2 (n, b, 1);
    if(!b_norm)
        b_norm = 1;
 
    iter = 0;  
   
    while( r_norm/b_norm > itol && iter < n ){
       
        iter++;
 
        cblas_dcopy (n, r, 1, z, 1);            //gia na min allaksei o r
        cblas_dcopy (n, r_t, 1, z_t, 1);        //gia na min allaksei o r_t
        for(i=0;i<n;i++){
            z[i]=precond[i]*z[i];
            z_t[i]=precond[i]*z_t[i];
        }
     
        rho = cblas_ddot (n, z, 1, r_t, 1);    
        if (fpclassify(fabs(rho)) == FP_ZERO){
            printf("RHO aborting Bi-CG due to EPS...\n");
            exit(42);
        }
         
        if (iter == 1){
            cblas_dcopy (n, z, 1, p, 1);
            cblas_dcopy (n, z_t, 1, p_t, 1);
        }
        else{      
            //p = z + beta*p;
            beta = rho/rho1;           
 
            cblas_dscal (n, beta, p, 1);        //rescale p by beta
            cblas_dscal (n, beta, p_t, 1);      //rescale p_t by beta
         
            cblas_daxpy (n, 1, z, 1, p, 1);     //p = 1*z + p
            cblas_daxpy (n, 1, z_t, 1, p_t, 1); //p_t = 1*z_t + p_t
        }
         
        rho1 = rho;
         
        //q = Ap
        //q_t = trans(A)*p_t
        memset(q, 0, n*sizeof(double));
        cs_gaxpy (A, p, q);
        memset(q_t, 0, n*sizeof(double));
        cs_gaxpy(AT, p_t, q_t);        
         
        omega = cblas_ddot (n, p_t, 1, q, 1);
        if (fpclassify(fabs(omega)) == FP_ZERO){
            printf("OMEGA aborting Bi-CG due to EPS...\n");
            exit(42);
        }
 
        alpha = rho/omega;     
 
        //x = x + aplha*p;
        cblas_dcopy (n, p, 1, temp_p, 1);
        cblas_dscal (n, alpha, temp_p, 1);//rescale by aplha
        cblas_daxpy (n, 1, temp_p, 1, res, 1);// sum x = 1*x + temp_p
 
        //R = R - aplha*Q;
        cblas_dcopy (n, q, 1, temp_q, 1);
        cblas_dscal (n, -alpha, temp_q, 1);//rescale by -aplha
        cblas_daxpy (n, 1, temp_q, 1, r, 1);// sum r = 1*r - temp_p    
 
        //~r=~r-alpha*~q
        cblas_dcopy (n, q_t, 1, temp_q, 1);
        cblas_dscal (n, -alpha, temp_q, 1);//rescale by -aplha
        cblas_daxpy (n, 1, temp_q, 1, r_t, 1);// sum r = 1*r - temp_p
 
        r_norm = cblas_dnrm2 (n, r, 1); //next step
    }
    cblas_dcopy (n, res, 1, x, 1);
 
    cs_spfree(AT);
}

void call_cg_sparse(cs* matrixA, double* matrixB, double* matrixX, int size){
     int i=0;

     if (found_dc_sweep==0) {
         conjugate_gradient_sparse(matrixA,matrixB, matrixX,size,itol_value);
         printf("Matrix X:\n");
         for(i=0;i<size;i++){
	     	printf(" %.6lf ",matrixX[i]);
	 	 }
	 	 printf("\n");
     }
     else {
         double value;
	     double *matrixB_temp = (double *)calloc(size,sizeof(double));
		 if (source > -1) {
         	for (value=start_value;value<=end_value;value=value+step){
                matrixB[source-1]=value;
                conjugate_gradient_sparse(matrixA,matrixB,matrixX,size,itol_value); 
				printf("value= %lf Matrix X: \n",value);
                for(i=0;i<size;i++){
	            	printf(" %.6lf ",matrixX[i]);
	 			}
	 			printf("\n"); 
       	    }
        }
        else { 
	    	if (sweep_node1!=0){
				matrixB[sweep_node1-1]+=sweep_value-start_value;
	    	}
	    	if(sweep_node2!=0){
				matrixB[sweep_node2-1]-=sweep_value-start_value;
	    	}
	    	for (value=start_value;value<=end_value;value=value+step) {
				conjugate_gradient_sparse(matrixA,matrixB,matrixX,size,itol_value);
	        	if (sweep_node1!=0){
		    		matrixB[sweep_node1-1]-=step;
	 			}
	        	if(sweep_node2!=0){
		    		matrixB[sweep_node2-1]+=step;
	        	}
				printf("value= %lf Matrix X: \n",value);
            
			    for(i=0;i<size;i++){
		            printf(" %.6lf ",matrixX[i]);
			 	}
			 	printf("\n"); 
		    } 
		}
    }
}

void call_bi_cg_sparse(cs* matrixA, double* matrixB, double* matrixX,int size){
     int i=0;

     if(found_dc_sweep==0){
         bi_conjugate_gradient_sparse(matrixA,matrixB,matrixX,size,itol_value);
         printf("Matrix X:\n");
         for(i=0;i<size;i++){
	     printf(" %.6lf ",matrixX[i]);
	 }
	 printf("\n");
     }
     else{
        double value;
		double *matrixB_temp = (double *)calloc(size,sizeof(double));
		if (source > -1) {
        	for (value=start_value;value<=end_value;value=value+step){
        	    matrixB[source-1]=value;
        	    bi_conjugate_gradient_sparse(matrixA,matrixB,matrixX,size,itol_value); 
				printf("value= %lf Matrix X: \n",value);
        	    
				for (i=0;i<size;i++) {
			         printf(" %.6lf ",matrixX[i]);
				}
				printf("\n"); 
       		}
    	}
	    else { 
		    if (sweep_node1!=0) {
				matrixB[sweep_node1-1]+=sweep_value-start_value;
		    }
		    if (sweep_node2!=0) {
				matrixB[sweep_node2-1]-=sweep_value-start_value;
		    }
		    for (value=start_value;value<=end_value;value=value+step) {
				bi_conjugate_gradient_sparse(matrixA,matrixB,matrixX,size,itol_value);
		        if (sweep_node1!=0) {
			    	matrixB[sweep_node1-1]-=step;
		        }
		        if (sweep_node2!=0) {
			    	matrixB[sweep_node2-1]+=step;
		        }
				printf("value= %lf Matrix X: \n",value);
    	    
			    for(i=0;i<size;i++){
		            printf(" %.6lf ",matrixX[i]);
		 		}
		 		printf("\n"); 
		    } 
		}
    }  
}
