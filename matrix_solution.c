#include "matrix_solution.h"
#include "lib.h"
#include <math.h>

void luDecomp(gsl_matrix* matrixA, gsl_vector* matrixB, int size){
    int s;
    gsl_vector* x = gsl_vector_calloc(size);
     
    gsl_permutation* p = gsl_permutation_alloc (size);
    gsl_permutation_init(p);
    gsl_linalg_LU_decomp(matrixA,p,&s);             
    
    printf("Matrix A: \n");   
    print2DMatrix(matrixA,size); 
 
    if(found_dc_sweep==0){
        gsl_linalg_LU_solve(matrixA,p,matrixB,x);      
	printf("Matrix X: \n");
        print1DMatrix(x,size);   
	printf("\n");     
    }
    else{
        double value;
	if (source > -1) {
            for(value=start_value;value<=end_value;value=value+step){
                gsl_vector_set(matrixB,source-1,value);  //agnoei tin timh pou eixe prin          
                gsl_linalg_LU_solve(matrixA,p,matrixB,x); //prepei na kratame ka8e dianusma 
		printf("value= %lf Matrix X: \n",value);
                print1DMatrix(x,size); 
		printf("\n");
       	    }
        }
        else { 
	    if (sweep_node1!=0){
		gsl_vector_set(matrixB,sweep_node1-1,gsl_vector_get(matrixB,sweep_node2-1)+sweep_value-start_value);
	    }
	    if(sweep_node2!=0){
		gsl_vector_set(matrixB,sweep_node2-1,gsl_vector_get(matrixB,sweep_node2-1)-sweep_value+start_value);
	    }
	    for(value=start_value;value<=end_value;value=value+step) {
		gsl_linalg_LU_solve(matrixA,p,matrixB,x);
	        if (sweep_node1!=0){
		    gsl_vector_set(matrixB,sweep_node1-1,gsl_vector_get(matrixB,sweep_node2-1)+sweep_value-start_value);
	        }
	        if(sweep_node2!=0){
		    gsl_vector_set(matrixB,sweep_node2-1,gsl_vector_get(matrixB,sweep_node2-1)-sweep_value+start_value);
	        }
		printf("value= %lf Matrix X: \n",value);
                print1DMatrix(x,size);
		printf("\n");
	    } 
	}
    } 

    printf("\n");
}

void choleskyDecomp(gsl_matrix* matrixA, gsl_vector* matrixB, int size){

    gsl_vector* x = gsl_vector_calloc(size);
    gsl_linalg_cholesky_decomp(matrixA);
  
    printf("Matrix A: \n");   
    print2DMatrix(matrixA,size); 
 
    if(found_dc_sweep==0){
        gsl_linalg_cholesky_solve(matrixA,matrixB,x);               
        printf("Matrix X:\n");
        print1DMatrix(x,size); 
		printf("\n");
    }
    else{
        double value;
		if (source > -1) {
            for(value=start_value;value<=end_value;value=value+step) {
                gsl_vector_set(matrixB,source-1,value);  //agnoei tin timh pou eixe prin          
                gsl_linalg_cholesky_solve(matrixA,matrixB,x); //prepei na kratame ka8e dianusma 
                print1DMatrix(x,size);
       	    }
        }
        else { 
	    	if (sweep_node1!=0){
				gsl_vector_set(matrixB,sweep_node1-1,gsl_vector_get(matrixB,sweep_node2-1)+sweep_value-start_value);
	    	}
	    	if(sweep_node2!=0){
				gsl_vector_set(matrixB,sweep_node2-1,gsl_vector_get(matrixB,sweep_node2-1)-sweep_value+start_value);
	    	}
	    	for(value=start_value;value<=end_value;value=value+step) {
				gsl_linalg_cholesky_solve(matrixA,matrixB,x);
	        	if (sweep_node1!=0) {
		    		gsl_vector_set(matrixB,sweep_node1-1,gsl_vector_get(matrixB,sweep_node2-1)+sweep_value-start_value);
	        	}
	        	if(sweep_node2!=0){
		    		gsl_vector_set(matrixB,sweep_node2-1,gsl_vector_get(matrixB,sweep_node2-1)-sweep_value+start_value);
	        	}
				printf("value= %lf Matrix X: \n",value);
                print1DMatrix(x,size);
				printf("\n");
	    	} 
		}
    }         
    printf("\n");      
} 

void conjugate_gradient(gsl_matrix *a,gsl_vector *b,gsl_vector *X,int n,double tolerance){ 
 
    double rho,rho1,alpha,beta;
    gsl_vector *r; 
    gsl_vector *z; 
    gsl_vector *p;
    gsl_vector *q;     
    gsl_vector *precond;
    gsl_vector *temp_p;
    gsl_vector *temp_q;
    gsl_vector *res;
     
    r = gsl_vector_calloc(n);
    z = gsl_vector_calloc(n);
    p = gsl_vector_calloc(n);
    q = gsl_vector_calloc(n);  
    precond = gsl_vector_calloc(n);
    temp_p = gsl_vector_calloc(n);
    temp_q = gsl_vector_calloc(n);
    res = gsl_vector_calloc(n);
 
    gsl_vector_view d;                      //gia na parw tin diagwnio
 
    d=gsl_matrix_diagonal(a);                   //d=diagwnios tou A
     
 
    int i;
     
    for(i=0;i<n;i++){
        if(gsl_vector_get(&d.vector,i)==0){
            gsl_vector_set(&d.vector,i,1);          //an kapoio stoixeio ths diagwniou einai 0 tote to 8etoume 1
        }
    }
 
    gsl_blas_dcopy(&d.vector,precond);              //gia na min allaksei h diagwnios tou a to antigrafw allou
     
    printf("\n");
    for(i=0;i<n;i++){
         gsl_vector_set(precond,i,1.0/gsl_vector_get(precond,i));   //precontitioner^-1 (M^-1) = 1/diag(A)
 
    }
 
    gsl_blas_dcopy(X,res);                          //Store X sto temp res
 
    gsl_blas_dcopy(b,r);   
    gsl_blas_dgemv(CblasNoTrans,1,a,res,0.0,p);             //prosorina p=A*x  
    gsl_vector_sub(r,p);   
     
    int iter=0;
     
    double r_norm = gsl_blas_dnrm2(r);
    double b_norm = gsl_blas_dnrm2(b);
    if(!b_norm)
        b_norm = 1;
 
    while( r_norm/b_norm > tolerance && iter < n ){
 
        iter++;
        gsl_blas_dcopy(r,z);                        //gia na min allaksei o r
        gsl_vector_mul(z,precond);                 
         
        gsl_blas_ddot(r,z,&rho);                    //r^T * Z
        if(iter==1){
            gsl_blas_dcopy(z,p);           
        }
     
        else{
            beta=rho/rho1;
            gsl_blas_dscal(beta,p);
            gsl_blas_daxpy(1,z,p); 
        }
        rho1=rho;
        gsl_blas_dgemv(CblasNoTrans,1,a,p,0.0,q);               //q=A*p
     
        gsl_blas_ddot(p,q,&alpha);                  //p^T * q
        alpha=rho/alpha;                        //alpha=rho/p^T*q
                 
                 
        gsl_blas_dcopy(p,temp_p);                   //x=x+alpha*p
        gsl_blas_dscal(alpha,temp_p);
        gsl_blas_daxpy(1,temp_p,res);  
 
        gsl_blas_dcopy(q,temp_q);                   //r=r-alpha*q
        gsl_blas_dscal(-alpha,temp_q);
        gsl_blas_daxpy(1,temp_q,r);
         
        r_norm = gsl_blas_dnrm2(r);                 //new r norm
 
    }
    gsl_blas_dcopy(res,X);                          //Restore res back to X
}
 
void bi_conjugate_gradient(gsl_matrix *a,gsl_vector *b,gsl_vector *X,int n,double tolerance){
 
    int i;
    double EPS = 1e-12;
    double rho,rho1,alpha,beta,omega;
    gsl_vector *r;
    gsl_vector *z;
    gsl_vector *p;
    gsl_vector *q;
    gsl_vector *precond;
    gsl_vector *temp_p;
    gsl_vector *temp_q,*temp_qt;
    gsl_vector *res;
     
    gsl_vector *r_t,*z_t,*q_t,*p_t;
    gsl_matrix *aT;
     
    r = gsl_vector_calloc(n);
    z = gsl_vector_calloc(n);
    p = gsl_vector_calloc(n);
    q = gsl_vector_calloc(n);
    precond = gsl_vector_calloc(n);
    temp_p = gsl_vector_calloc(n);
    temp_q = gsl_vector_calloc(n);
    res = gsl_vector_calloc(n);
    aT = gsl_matrix_calloc(n,n);
     
    r_t = gsl_vector_calloc(n);
    z_t = gsl_vector_calloc(n);
    p_t = gsl_vector_calloc(n);
    q_t = gsl_vector_calloc(n);
    temp_qt = gsl_vector_calloc(n);
    gsl_vector_view d;                      //gia na parw tin diagwnio
 
    d=gsl_matrix_diagonal(a);                   //d=diagwnios tou A
     
 
     
    for(i=0;i<n;i++){
        if(gsl_vector_get(&d.vector,i)==0){
            gsl_vector_set(&d.vector,i,1);          //an kapoio stoixeio ths diagwniou einai 0 tote to 8etoume 1
        }
    }
 
    gsl_blas_dcopy(&d.vector,precond);              //gia na min allaksei h diagwnios tou a to antigrafw allou
     
 
    for(i=0;i<n;i++){
         gsl_vector_set(precond,i,1.0/gsl_vector_get(precond,i));   //precontitioner^-1 (M^-1) = 1/diag(A)
    }
 
     
    gsl_blas_dcopy(X,res);                      //Store X sto temp res
 
    //r=b-Ax
    gsl_blas_dcopy(b,r);   
    gsl_blas_dgemv(CblasNoTrans,1,a,res,0.0,p);         //prosorina p=A*x  
    gsl_vector_sub(r,p);   
     
    //transport r_t = r
    gsl_blas_dcopy(r,r_t);
     
    int iter=0;
     
    double r_norm = gsl_blas_dnrm2(r);
    double b_norm = gsl_blas_dnrm2(b);
    if(!b_norm){b_norm = 1;}
 
 
    while(((r_norm/b_norm) > tolerance) && (iter < n)){
        iter++;
        gsl_blas_dcopy(r,z);                    //gia na min allaksei o r
        gsl_vector_mul(z,precond);             
         
         
        // transport
        // z_t = MNA * z_t
        gsl_blas_dcopy(r_t,z_t);
        gsl_vector_mul(z_t,precond);
 
        gsl_blas_ddot(z,r_t,&rho);              //r^T * Z
 
        if(fabs(rho)<EPS){ printf("---------rho < EPS--------\n"); exit(0);}
        if (iter == 1){
            gsl_blas_dcopy(z,p);
            gsl_blas_dcopy(z_t,p_t);
        }
     
        else{
            beta=rho/rho1;
             
            gsl_blas_dscal(beta,p);
            gsl_blas_daxpy(1,z,p); 
             
            gsl_blas_dscal(beta,p_t);
            gsl_blas_daxpy(1,z_t,p_t);
             
        }
 
        gsl_matrix_memcpy(aT,a);
        gsl_matrix_transpose(aT);
        rho1=rho;
        gsl_blas_dgemv(CblasNoTrans,1,a,p,0.0,q);           //q=A*p
         
        gsl_blas_dgemv(CblasNoTrans,1,aT,p_t,0.0,q_t);      //q_t = trans(A)*p_t
         
         
        gsl_blas_ddot(p_t,q,&omega);                //omega = trasn(p_t)*q
        if(fabs(omega)<EPS){printf("--------- omega < EPS --------\n"); exit(0);}
        alpha = rho/omega;
             
        gsl_blas_dcopy(p,temp_p);               //x=x+alpha*p
        gsl_blas_dscal(alpha,temp_p);
        gsl_blas_daxpy(1,temp_p,res);
 
        gsl_blas_dcopy(q,temp_q);               //r=r-alpha*q
        gsl_blas_dscal(-alpha,temp_q);
        gsl_blas_daxpy(1,temp_q,r);
 
        gsl_blas_dcopy(p,temp_qt);              //r_t = r_t-alpha*q_t
        gsl_blas_dscal(-alpha,temp_qt);
        gsl_blas_daxpy(1,temp_qt,r_t);
         
        r_norm = gsl_blas_dnrm2(r);             //new r norm
 
    }
    gsl_blas_dcopy(res,X);                      //Restore res back to X
}

void call_cg(gsl_matrix* matrixA, gsl_vector* matrixB, int size){
     gsl_vector* x = gsl_vector_calloc(size);
     if(found_dc_sweep==0){
         conjugate_gradient(matrixA,matrixB,x,size,itol_value);
         printf("Matrix X:\n");
         print1DMatrix(x,size);
     }
     else{
        double value;
	if (source > -1) {
            for(value=start_value;value<=end_value;value=value+step){
                gsl_vector_set(matrixB,source-1,value);  //agnoei tin timh pou eixe prin          
                conjugate_gradient(matrixA,matrixB,x,size,itol_value); 
		printf("value= %lf Matrix X: \n",value);
                print1DMatrix(x,size); 
       	    }
        }
        else { 
	    if (sweep_node1!=0){
		gsl_vector_set(matrixB,sweep_node1-1,gsl_vector_get(matrixB,sweep_node2-1)+sweep_value-start_value);
	    }
	    if(sweep_node2!=0){
		gsl_vector_set(matrixB,sweep_node2-1,gsl_vector_get(matrixB,sweep_node2-1)-sweep_value+start_value);
	    }
	    for(value=start_value;value<=end_value;value=value+step) {
		conjugate_gradient(matrixA,matrixB,x,size,itol_value);
	        if (sweep_node1!=0){
		    gsl_vector_set(matrixB,sweep_node1-1,gsl_vector_get(matrixB,sweep_node2-1)+sweep_value-start_value);
	        }
	        if(sweep_node2!=0){
		    gsl_vector_set(matrixB,sweep_node2-1,gsl_vector_get(matrixB,sweep_node2-1)-sweep_value+start_value);
	        }
		printf("value= %lf Matrix X: \n",value);
                print1DMatrix(x,size);
	    } 
	}
    }
}

void call_bi_cg(gsl_matrix* matrixA, gsl_vector* matrixB, int size){
     gsl_vector* x = gsl_vector_calloc(size);
     if(found_dc_sweep==0){
         bi_conjugate_gradient(matrixA,matrixB,x,size,itol_value);
         printf("Matrix X:\n");
         print1DMatrix(x,size);
	 printf("\n");
     }
     else{
        double value;
		if (source > -1) {
            for(value=start_value;value<=end_value;value=value+step){
                gsl_vector_set(matrixB,source-1,value);  //agnoei tin timh pou eixe prin          
                bi_conjugate_gradient(matrixA,matrixB,x,size,itol_value); 
       	    }
        }
        else { 
	    	if (sweep_node1!=0){
				gsl_vector_set(matrixB,sweep_node1-1,gsl_vector_get(matrixB,sweep_node2-1)+sweep_value-start_value);
	    	}
	    	if(sweep_node2!=0){
				gsl_vector_set(matrixB,sweep_node2-1,gsl_vector_get(matrixB,sweep_node2-1)-sweep_value+start_value);
	    	}
	    	for(value=start_value;value<=end_value;value=value+step) {
				bi_conjugate_gradient(matrixA,matrixB,x,size,itol_value);
	        	if (sweep_node1!=0){
		    		gsl_vector_set(matrixB,sweep_node1-1,gsl_vector_get(matrixB,sweep_node2-1)+sweep_value-start_value);
	        	}
	        	if(sweep_node2!=0){
		    		gsl_vector_set(matrixB,sweep_node2-1,gsl_vector_get(matrixB,sweep_node2-1)-sweep_value+start_value);
	        	}
				printf("value= %lf Matrix X: \n",value);
                print1DMatrix(x,size);
				printf("\n");
		    } 
		}
    }  
}
