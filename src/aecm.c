#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<R.h>
#include "functions.h"

void lambda_store(double *lam_vec, double *lambda, int p, int q){
    int i,j,k=0;
    for (i=0;i<p;i++){
        for(j=0;j<q;j++){
            lam_vec[k]=lambda[i*q+j];
            k++;
        }
    }
}

void lambda_storeG(double *lam_vec, double **lambda, int G, int p, int q){
    int i,g,k=0;
    for(g=0;g<G;g++){
        for (i=0;i<p*q;i++){
                lam_vec[k]=lambda[g][i];
                k++;
            }
        }
    }


double aecm(double *z, double *x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){
	
        double bic;
	int it=0,stop=0,paras;
	double psi, log_detpsi, log_c=0.0, log_detsig;

        double *pi = malloc(sizeof(double)*G);
        double *n = malloc(sizeof(double)*G);
        double *max_v = malloc(sizeof(double)*N);
        double *at = malloc(sizeof(double)*150000);
        double *l = malloc(sizeof(double)*150000);
        double *sampcovtilde = malloc(sizeof(double)*p*p);
        double *v = malloc(sizeof(double)*N*G);
        double *lambda = malloc(sizeof(double)*p*q);
        double *beta = malloc(sizeof(double)*q*p);
        double *theta = malloc(sizeof(double)*q*q);
        double *mu = malloc(sizeof(double)*G*p);

	
	psi = *psi_vec;	
	get_data(lam_vec,lambda,p,q);
	
    while(stop==0){            
		update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	  
	    update_mu(mu, n, x, z, G, N, p);

		if(it>0){ 
			update_z(v, x, z, lambda, psi, mu, pi, max_v, log_c, N, G, p, q);
			known_z(cls,z,N,G);
		}

	    update_stilde(sampcovtilde,x, z, mu, G, N, p);

	    update_beta1(beta, psi, lambda, p, q);
        
   	    update_theta(theta, beta, lambda, sampcovtilde, p, q); 
    
        update_lambda(lambda, beta, sampcovtilde, theta, p, q);
 
        psi = update_psi(lambda, beta, sampcovtilde, p, q); 
   
 	    log_detpsi=p*log(psi);
	    
	    log_detsig = update_det_sigma_NEW(lambda, psi, log_detpsi, p, q);
	  
	    log_c = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig;
       
		update_z(v, x, z, lambda, psi, mu, pi, max_v, log_c, N, G, p, q);
		known_z(cls,z,N,G);
 
        stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++;
		/*printf("ll=%f\t",l[it-1]);*/
	}
	/*printf("ll=%f\n",l[it-1]);*/

	paras = G-1 + G*p + p*q - q*(q-1)/2 + 1;
    	
   	bic = 2.0*l[it-1] - paras*log(N);
	
	lambda_store(lam_vec, lambda, p, q);
    
        free(lambda); free(mu); free(n); free(beta), free(theta); 
        free(sampcovtilde); free(l); free(at); free(pi);
       
	return bic;
} 

double aecm2(double *z, double *x, int *cls, int q, int p, int G, int N, double *lam_vec, double *Psi, double tol){

        double bic;
	int i,it=0,stop=0;
	double a, log_detpsi,log_detsig, log_c=0.0;

	/* printf("G=%d \t q=%d\n",G,q);*/

   	/* Allocate memory - matrices and vectors */ 
        double *det = malloc(sizeof(double)*G);
        double *pi = malloc(sizeof(double)*G);
        double *n = malloc(sizeof(double)*G);
        double *at = malloc(sizeof(double)*150000);
        double *l = malloc(sizeof(double)*150000);
        double *sampcovtilde = malloc(sizeof(double)*p*p);
        double *lambda = malloc(sizeof(double)*p*q);
        double *beta = malloc(sizeof(double)*q*p);
        double *theta = malloc(sizeof(double)*q*q);
        double *mu = malloc(sizeof(double)*G*p);
        double *w = malloc(sizeof(double)*G*N);
        double *max_v = malloc(sizeof(double)*N);
        double *v = malloc(sizeof(double)*N*G);
		
	get_data(lam_vec,lambda,p,q);
		
   	while(stop==0 && it<25 ){            
	    
	    update_n(n, z, G, N);
	    
	    update_pi(pi, n, G, N);
	    
	    update_mu(mu, n, x, z, G, N, p); 

        if (it>0){ 
			update_z2(v, x, z, lambda, Psi, mu, pi, max_v, log_c, N, G, p, q);
			known_z(cls,z,N,G);
		}
        
		update_stilde(sampcovtilde, x, z, mu, G, N, p);        
        
	    update_beta2(beta, Psi, lambda, p, q);

	    update_theta(theta, beta, lambda, sampcovtilde, p, q);

	    update_lambda(lambda, beta, sampcovtilde, theta, p, q);
			    
	    update_psi2(Psi, lambda, beta, sampcovtilde, p, q);

		log_detpsi = 0.0; for(i=0;i<p;i++) log_detpsi += log(Psi[i]);

	    log_detsig = update_det_sigma_NEW2(lambda, Psi, log_detpsi, p, q);
	    
		log_c = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig; 
   
		update_z2(v, x, z, lambda, Psi, mu, pi, max_v, log_c, N, G, p, q);
		known_z(cls,z,N,G);
 
        stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	
	    it++; 
    }
	/*printf("ll=%f\n",l[it-1]);*/
    
	/* no of parameters {pq - q(q - 1)/2} + p */
    a = G-1 + G*p + p*q - q*(q-1)/2 + p;
    
    /* BIC = 2log_likelihood - mlog n; */
    bic = 2.0*l[it-1] - a*log(N);
    
	lambda_store(lam_vec, lambda, p, q);
    
    /* Deallocate memory */
    free(lambda); free(mu); free(w); free(n); free(det);
    free(beta); free(theta); free(sampcovtilde); free(l);
    free(at); free(pi);

	return bic;
}

double aecm3(double *z, double *x, int *cls, int q, int p, int G, int N, double *lam_vec, double *Psi, double tol){

        double bic, a;
	int g,it=0,stop=0;

	/* printf("G=%d \t q=%d\n",G,q);*/
   double *log_det = malloc(sizeof(double)*G);
   double *log_detpsi = malloc(sizeof(double)*G);
   double *log_detsig = malloc(sizeof(double)*G);
   double *pi = malloc(sizeof(double)*G);
   double *n = malloc(sizeof(double)*G);
   double *at = malloc(sizeof(double)*150000);
   double *l = malloc(sizeof(double)*150000);
   double *lambda = malloc(sizeof(double)*p*q);
   double **sampcov    = malloc(sizeof(double*)*G);
   double **beta    = malloc(sizeof(double*)*G);
   double **theta    = malloc(sizeof(double*)*G);
   for(g=0; g < G; g++) {
      sampcov[g]      = malloc(sizeof(double)*p*p); 
      beta[g]      = malloc(sizeof(double)*q*p); 
      theta[g]      = malloc(sizeof(double)*q*q);
    } 
   double *mu = malloc(sizeof(double)*G*p);
   double *max_v = malloc(sizeof(double)*N);
   double *v = malloc(sizeof(double)*N*G);
    
	get_data(lam_vec,lambda,p,q);
		
    while(stop==0){            
           
        update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	
	    update_mu(mu, n, x, z, G, N, p); 
                          
		if(it>0){ 
			update_z3(v, x, z, lambda, Psi, mu, pi, max_v, log_det, N, G, p, q);
			known_z(cls,z,N,G);
        }

	    update_sg(sampcov, x, z, mu, n, p, G, N);
	    
	    for(g=0;g<G;g++) update_beta1(beta[g], Psi[g], lambda, p, q);    
	
    	for (g=0;g<G;g++) update_theta(theta[g], beta[g], lambda, sampcov[g], p, q);
        
	    update_lambda2(lambda, beta, sampcov, theta, n, Psi, p, q, G);
	    
        for(g=0;g<G;g++) Psi[g] = update_psi3(lambda, beta[g], sampcov[g], theta[g], p, q);
    
        for(g=0;g<G;g++){
	    	log_detpsi[g] = p*log(Psi[g]);
			log_detsig[g] = update_det_sigma_NEW(lambda, Psi[g], log_detpsi[g], p, q);
			log_det[g] = (p/2.0)*log(2*M_PI) + 0.5*log_detsig[g];
		} 
	    
		update_z3(v, x, z, lambda, Psi, mu, pi, max_v, log_det, N, G, p, q);
		known_z(cls,z,N,G);

	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
		it++;
    }
	
	/*printf("ll=%f\n",l[it-1]);*/
	
   	a = G-1 + G*p+ p*q - q*(q-1)/2 + G;
    
   	bic = 2.0*l[it-1] - a*log(N);
    
	lambda_store(lam_vec, lambda, p, q);
    
   	/* Deallocate memory */
        free(lambda); free(mu); free(v); free(n); free(log_det); free(max_v);
        free(l); free(at); free(pi); free(log_detpsi);
        for(g=0; g < G; g++) {
           free(beta[g]);
           free(theta[g]);
           free(sampcov[g]);
        }
        free(beta); free(theta); free(sampcov);    
	return bic;
}

double aecm4(double *z, double *x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){

        double bic;
	int it=0,stop=0,paras, g, j;
	
	/* printf("G=%d \t q=%d\n",G,q);*/

   	/* Allocate memory - matrices and vectors */ 
        double *max_v = malloc(sizeof(double)*N);
        double *v = malloc(sizeof(double)*N*G);
        double *pi = malloc(sizeof(double)*G);
        double *n = malloc(sizeof(double)*G);
        double *at = malloc(sizeof(double)*150000);
        double *l = malloc(sizeof(double)*150000);
        double *lambda = malloc(sizeof(double)*p*q);
        double **sampcov    = malloc(sizeof(double*)*G);
        double **beta    = malloc(sizeof(double*)*G);
        double **theta    = malloc(sizeof(double*)*G);
        for(g=0; g < G; g++) {
           sampcov[g]      = malloc(sizeof(double)*p*p);
           beta[g]      = malloc(sizeof(double)*q*p);
           theta[g]      = malloc(sizeof(double)*q*q);
        }
        double *mu = malloc(sizeof(double)*G*p);
        double *w = malloc(sizeof(double)*G*N);
        double *Psi = malloc(sizeof(double)*G*p);
        double *log_detpsi = malloc(sizeof(double)*G);
        double *log_detsig = malloc(sizeof(double)*G);
        double *log_c = malloc(sizeof(double)*G);
        double *Psi0 = malloc(sizeof(double)*p);

	get_data(psi_vec,Psi,G,p);
	get_data(lam_vec,lambda,p,q);
        
   	while(stop==0  ){            
   
	    update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	
	    update_mu(mu, n, x, z, G, N, p); 

        if(it>0){ 
			update_z4(v, x, z, lambda, Psi, mu, pi, max_v, log_c, N, G, p, q);
			known_z(cls,z,N,G);
		} 

 
	    update_sg(sampcov,x, z, mu, n, p, G, N);
        
	    for(g=0;g<G;g++){
                for(j=0; j<p; j++){
                   Psi0[j] = Psi[g*p+j];
                }
                 update_beta2(beta[g], Psi0, lambda, p, q);
             } 
   	    for(g=0;g<G;g++) update_theta(theta[g], beta[g], lambda, sampcov[g], p, q); 
    
	    update_lambda_cuu(lambda, beta, sampcov, theta, n, Psi, p, q, G);
	 
        update_psi_cuu(Psi, lambda, beta, sampcov, theta, p, q, G);

	    for (g=0;g<G;g++){
			log_detpsi[g]=0.0;
			for(j=0;j<p;j++) log_detpsi[g] += log(Psi[g*p+j]);
	    }
	
	    for(g=0;g<G;g++){
               for(j=0; j<p; j++){
                   Psi0[j] = Psi[g*p+j];
                }
                log_detsig[g] = update_det_sigma_NEW2(lambda, Psi0, log_detpsi[g], p, q);
             }
		for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];

        update_z4(v, x, z, lambda, Psi, mu, pi, max_v, log_c, N, G, p, q);
		known_z(cls,z,N,G);
            
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++;
	}
	
	/* no of parameters {pq - q(q - 1)/2} + Gp */
	paras = G-1 + G*p + p*q - q*(q-1)/2 + G*p;
    
   	/* BIC = 2log_likelihood - mlog n; */
   	bic = 2.0*l[it-1] - paras*log(N);
    
	lambda_store(lam_vec, lambda, p, q);
    
    /* Deallocate memory */

        free(lambda); free(mu); free(w); free(n); free(l); free(at); free(pi);
        free(log_detsig); free(log_c); free(log_detpsi); free(Psi); free(Psi0);

        for(g=0; g < G; g++) {
           free(beta[g]);
           free(theta[g]);
           free(sampcov[g]);
        }
        free(beta); free(theta); free(sampcov);
	return bic;
} 

double aecm5(double *z, double *x, int*cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){

        double bic, log_detpsi, psi;
	int it=0,stop=0,paras, g,i;

	/* printf("G=%d \t q=%d\n",G,q);*/
        double *max_v = malloc(sizeof(double)*N);
        double *v = malloc(sizeof(double)*N*G);
        double *log_c = malloc(sizeof(double)*G);
        double *log_detsig = malloc(sizeof(double)*G);
        double *pi = malloc(sizeof(double)*G);
        double *n = malloc(sizeof(double)*G);
        double *at = malloc(sizeof(double)*150000);
        double *l = malloc(sizeof(double)*150000);
        double **sampcov    = malloc(sizeof(double*)*G);
        double **lambda    = malloc(sizeof(double*)*G);
        double **beta    = malloc(sizeof(double*)*G);
        double **theta    = malloc(sizeof(double*)*G);
        for(g=0; g < G; g++) {
           sampcov[g]      = malloc(sizeof(double)*p*p);
           lambda[g]      = malloc(sizeof(double)*p*q);
           beta[g]      = malloc(sizeof(double)*q*p);
           theta[g]      = malloc(sizeof(double)*q*q);
        }
        double *mu = malloc(sizeof(double)*G*p);
        double *w = malloc(sizeof(double)*G*N);

	psi=*psi_vec;
	get_data2(lam_vec,lambda,G,p,q);
       
   	while(stop==0){            
    
	    update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	
	    update_mu(mu, n, x, z, G, N, p); 

		if(it>0){ 
			update_z5(v, x, z, lambda, psi, mu, pi, max_v, log_c, N, G, p, q);
			known_z(cls,z,N,G);
        }
 
	    update_sg(sampcov,x, z, mu, n, p, G, N);
        
	    for(g=0;g<G;g++) update_beta1(beta[g], psi, lambda[g], p, q);
        
   	    for(g=0;g<G;g++) update_theta(theta[g], beta[g], lambda[g], sampcov[g], p, q); 
    
        for (g=0;g<G;g++) update_lambda(lambda[g], beta[g], sampcov[g], theta[g], p, q);
	 
	    psi = update_psi_ucc(lambda, beta, sampcov, p, q, pi, G);
   
		log_detpsi = 0.0;
		for(i=0;i<p;i++) log_detpsi += log(psi);
	    
	    for(g=0;g<G;g++) log_detsig[g] = update_det_sigma_NEW(lambda[g], psi, log_detpsi, p, q);
	    
		for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];
       
		update_z5(v, x, z, lambda, psi, mu, pi, max_v, log_c, N, G, p, q);
		known_z(cls,z,N,G);
            
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++; 
	}
	/*printf("ll=%f\n",l[it-1]);*/

	paras = G-1 + G*p + G*(p*q - q*(q-1)/2) + 1;
    
   	bic = 2.0*l[it-1] - paras*log(N);
    
    lambda_storeG(lam_vec, lambda, G, p, q);

    free(mu); free(w); free(n); free(l); free(at); free(pi); free(log_detsig); free(log_c);
        for(g=0; g < G; g++) {
           free(beta[g]);
           free(lambda[g]);
           free(theta[g]);
           free(sampcov[g]);
        }
        free(beta); free(lambda); free(theta); free(sampcov);

	return bic;
}

double aecm6(double *z, double *x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi, double tol){

        
        double bic, log_detpsi;
	int it=0,stop=0,paras, g, j;

	/* printf("G=%d \t q=%d\n",G,q);*/
        double *max_v = malloc(sizeof(double)*N);
        double *v = malloc(sizeof(double)*N*G);
        double *log_detsig = malloc(sizeof(double)*G);
        double *log_c = malloc(sizeof(double)*G);
        double *pi = malloc(sizeof(double)*G);
        double *n = malloc(sizeof(double)*G);
        double *at = malloc(sizeof(double)*150000);
        double *l = malloc(sizeof(double)*150000);
        double **sampcov    = malloc(sizeof(double*)*G);
        double **lambda    = malloc(sizeof(double*)*G);
        double **beta    = malloc(sizeof(double*)*G);
        double **theta    = malloc(sizeof(double*)*G);
        for(g=0; g < G; g++) {
           sampcov[g]      = malloc(sizeof(double)*p*p);
           lambda[g]      = malloc(sizeof(double)*p*q);
           beta[g]      = malloc(sizeof(double)*q*p);
           theta[g]      = malloc(sizeof(double)*q*q);
        }
        double *mu = malloc(sizeof(double)*G*p);

	
	get_data2(lam_vec,lambda,G,p,q);
      
   	while(stop==0){            
    
	    update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	
	    update_mu(mu, n, x, z, G, N, p); 

        if(it>0){ 
			update_z6(v, x, z, lambda, psi, mu, pi, max_v, log_c, N, G, p, q);
			known_z(cls,z,N,G);
		} 

	    update_sg(sampcov,x, z, mu, n, p, G, N);
        
	    for(g=0;g<G;g++) update_beta2(beta[g], psi, lambda[g], p, q);
        
   	    for(g=0;g<G;g++) update_theta(theta[g], beta[g], lambda[g], sampcov[g], p, q); 
    
        for (g=0;g<G;g++) update_lambda(lambda[g], beta[g], sampcov[g], theta[g], p, q);
	 
	    update_psi_ucu(psi, lambda, beta, sampcov, p, q, pi, G);
                    
		log_detpsi = 0.0;
		for(j=0;j<p;j++) log_detpsi += log(psi[j]);
	    
	    for(g=0;g<G;g++) log_detsig[g] = update_det_sigma_NEW2(lambda[g], psi, log_detpsi, p, q);
	    
		for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];

        update_z6(v, x, z, lambda, psi, mu, pi, max_v, log_c, N, G, p, q);
		known_z(cls,z,N,G);
            
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++; 
	}
	/*printf("ll=%f\n",l[it-1]);*/
	
	paras = G-1 + G*p + G*(p*q - q*(q-1)/2) + p;
    
   	bic = 2.0*l[it-1] - paras*log(N);
    
    lambda_storeG(lam_vec, lambda, G, p, q);
    
        free(mu); free(v); free(n); free(max_v); free(l); free(at); free(pi); free(log_detsig); free(log_c);
        for(g=0; g < G; g++) {
           free(beta[g]);
           free(lambda[g]);
           free(theta[g]);
           free(sampcov[g]);
        }
        free(beta); free(lambda); free(theta); free(sampcov);

	return bic;
}

double aecm7(double *z, double *x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi, double tol){

        double bic;
	int it=0,stop=0,paras, g;

	/* printf("G=%d \t q=%d\n",G,q); */
        double *max_v = malloc(sizeof(double)*N);
        double *v = malloc(sizeof(double)*N*G);
        double *log_detpsi = malloc(sizeof(double)*G);
        double *log_detsig = malloc(sizeof(double)*G);
        double *log_c = malloc(sizeof(double)*G);
        double *pi = malloc(sizeof(double)*G);
        double *n = malloc(sizeof(double)*G);
        double *at = malloc(sizeof(double)*150000);
        double *l = malloc(sizeof(double)*150000);
        double **sampcov    = malloc(sizeof(double*)*G);
        double **lambda    = malloc(sizeof(double*)*G);
        double **beta    = malloc(sizeof(double*)*G);
        double **theta    = malloc(sizeof(double*)*G);
        for(g=0; g < G; g++) {
           sampcov[g]      = malloc(sizeof(double)*p*p);
           lambda[g]      = malloc(sizeof(double)*p*q);
           beta[g]      = malloc(sizeof(double)*q*p);
           theta[g]      = malloc(sizeof(double)*q*q);
        }
        double *mu = malloc(sizeof(double)*G*p);

	
	get_data2(lam_vec,lambda,G,p,q);

   	while(stop==0){            
    
	    update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	
	    update_mu(mu, n, x, z, G, N, p); 

        if(it>0){ 
			update_z7(v, x, z, lambda, psi, mu, pi, max_v, log_c, N, G, p, q);
			known_z(cls,z,N,G);
		}

	    update_sg(sampcov,x, z, mu, n, p, G, N);
        
	    for(g=0;g<G;g++) update_beta1(beta[g], psi[g], lambda[g], p, q);
        
   	    for(g=0;g<G;g++) update_theta(theta[g], beta[g], lambda[g], sampcov[g], p, q); 
    
        for (g=0;g<G;g++) update_lambda(lambda[g], beta[g], sampcov[g], theta[g], p, q);
	 
        for (g=0;g<G;g++) psi[g]=update_psi(lambda[g], beta[g], sampcov[g], p, q); 
           
	    for(g=0;g<G;g++) log_detpsi[g] = p*log(psi[g]);
	    
	    for(g=0;g<G;g++) log_detsig[g] = update_det_sigma_NEW(lambda[g], psi[g], log_detpsi[g], p, q);
	   
		for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];

		update_z7(v, x, z, lambda, psi, mu, pi, max_v, log_c, N, G, p, q);
		known_z(cls,z,N,G);
            
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++;
	}
	/*printf("ll=%f\n",l[it-1]);*/

	paras = G-1 + G*p + G*(p*q - q*(q-1)/2) + G;
    
   	bic = 2.0*l[it-1] - paras*log(N);
    
    lambda_storeG(lam_vec, lambda, G, p, q);

    free(mu); free(v); free(n); free(l); free(at); free(pi);  free(log_detpsi); 
     free(log_detsig); free(log_c);
        for(g=0; g < G; g++) {
           free(beta[g]);
           free(lambda[g]);
           free(theta[g]);
           free(sampcov[g]);
        }
        free(beta); free(lambda); free(theta); free(sampcov);
	return bic;
}

double aecm8(double *z, double *x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){

        double bic;
	int it=0,stop=0,paras, g, j;

	/* printf("G=%d \t q=%d\n",G,q); */

        double *max_v = malloc(sizeof(double)*N);
        double *v = malloc(sizeof(double)*N*G);
        double *pi = malloc(sizeof(double)*G);
        double *n = malloc(sizeof(double)*G);
        double *at = malloc(sizeof(double)*150000);
        double *l = malloc(sizeof(double)*150000);
        double **sampcov    = malloc(sizeof(double*)*G);
        double **lambda    = malloc(sizeof(double*)*G);
        double **beta    = malloc(sizeof(double*)*G);
        double **theta    = malloc(sizeof(double*)*G);
        for(g=0; g < G; g++) {
           sampcov[g]      = malloc(sizeof(double)*p*p);
           lambda[g]      = malloc(sizeof(double)*p*q);
           beta[g]      = malloc(sizeof(double)*q*p);
           theta[g]      = malloc(sizeof(double)*q*q);
        }
        double *mu = malloc(sizeof(double)*G*p);
        double *Psi = malloc(sizeof(double)*G*p);
        double *log_detpsi = malloc(sizeof(double)*G);
        double *log_detsig = malloc(sizeof(double)*G);
        double *log_c = malloc(sizeof(double)*G);
        double *Psi0 = malloc(sizeof(double)*p);
	
	get_data(psi_vec,Psi,G,p);
	get_data2(lam_vec,lambda,G,p,q);
	       
   	while(stop==0){            
    
	    update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	
	    update_mu(mu, n, x, z, G, N, p);

        if(it>0){ 
			update_z8(v, x, z, lambda, Psi, mu, pi, max_v, log_c, N, G, p, q);
			known_z(cls,z,N,G);
		}	
       
	    update_sg(sampcov,x, z, mu, n, p, G, N);
	    
	    for(g=0;g<G;g++){
               for(j=0; j<p; j++){
                   Psi0[j] = Psi[g*p+j];
                }

                update_beta2(beta[g], Psi0, lambda[g], p, q);
             }
        
   	    for(g=0;g<G;g++) update_theta(theta[g], beta[g], lambda[g], sampcov[g], p, q); 
    
        for (g=0;g<G;g++) update_lambda(lambda[g], beta[g], sampcov[g], theta[g], p, q);
	 
        for (g=0;g<G;g++){
         
             update_psi2(Psi0, lambda[g], beta[g], sampcov[g], p, q);
             for(j=0; j<p; j++){
                   Psi[g*p+j] = Psi0[j];
                }
        } 
           
	    for(g=0;g<G;g++){
		    log_detpsi[g]=0.0;
		    for(j=0;j<p;j++) log_detpsi[g] += log(Psi[g*p+j]);
	    }
	    
		for(g=0;g<G;g++){
                   for(j=0; j<p; j++){
                   Psi0[j] = Psi[g*p+j];
                }

                    log_detsig[g] = update_det_sigma_NEW2(lambda[g], Psi0, log_detpsi[g], p, q);
	    }	    
		for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];

		update_z8(v, x, z, lambda, Psi, mu, pi, max_v, log_c, N, G, p, q);
		known_z(cls,z,N,G);
            
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++; 
	}
	/*printf("ll=%f\n",l[it-1]);*/
	paras = G-1 + G*p + G*(p*q - q*(q-1)/2) + G*p;
    
   	bic = 2.0*l[it-1] - paras*log(N);
    
    lambda_storeG(lam_vec, lambda, G, p, q);
    free(mu); free(v); free(n); free(log_detpsi); free(l); free(at); free(pi);
    free(log_detsig); free(log_c); free(Psi); free(max_v); free(Psi0); 

        for(g=0; g < G; g++) {
           free(beta[g]);
           free(lambda[g]);
           free(theta[g]);
           free(sampcov[g]);
        }
        free(beta); free(lambda); free(theta); free(sampcov);
	return bic;
}

double aecm9(double *z, double *x, int *cls, int q, int p, int G, int N, double *lam_vec, double *omega, double tol){

        double bic, a;
	int g,i,it=0,stop=0;
	       
	/* printf("G=%d \t q=%d\n",G,q);*/

   	/* Allocate memory - matrices and vectors */ 
        double *max_v = malloc(sizeof(double)*N);
        double *v = malloc(sizeof(double)*N*G);
        double *log_detpsi = malloc(sizeof(double)*G);
        double *log_detsig = malloc(sizeof(double)*G);
        double *log_c = malloc(sizeof(double)*G);
        double *pi = malloc(sizeof(double)*G);
        double *n = malloc(sizeof(double)*G);
        double *at = malloc(sizeof(double)*150000);
        double *l = malloc(sizeof(double)*150000);
        double *lambda = malloc(sizeof(double)*p*q);
        double **sampcov    = malloc(sizeof(double*)*G);
        double **beta    = malloc(sizeof(double*)*G);
        double **theta    = malloc(sizeof(double*)*G);
        for(g=0; g < G; g++) {
           sampcov[g]      = malloc(sizeof(double)*p*p);
           beta[g]      = malloc(sizeof(double)*q*p);
           theta[g]      = malloc(sizeof(double)*q*q);
        }
        double *mu = malloc(sizeof(double)*G*p);
        double *delta = malloc(sizeof(double)*p);
        double *psi = malloc(sizeof(double)*p);

	get_data(lam_vec,lambda,p,q);
	for(i=0;i<p;i++) delta[i] = 1.0; 
	
    while(stop==0){            
           
        update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	
	    update_mu(mu, n, x, z, G, N, p); 
                          
        if(it>0){ 
			update_z9(v, x, z, lambda, omega, delta, mu, pi, max_v, log_c, N, G,p,q);
			known_z(cls,z,N,G);
		}	
		
	    update_sg(sampcov, x, z, mu, n, p, G, N);

	    for(g=0;g<G;g++){ 
			for (i=0;i<p;i++) psi[i] = omega[g]*delta[i];
			update_beta2(beta[g], psi, lambda, p, q);    
		}	

    	for (g=0;g<G;g++) update_theta(theta[g], beta[g], lambda, sampcov[g], p, q);
        
	    update_lambda2(lambda, beta, sampcov, theta, n, omega, p, q, G);
	    
        for(g=0;g<G;g++) omega[g] = update_omega(lambda, delta, beta[g], sampcov[g], theta[g], p, q);

	    update_delta(delta, lambda, omega, beta, sampcov, theta, n, p, q, N, G); 
      
        for(g=0;g<G;g++){
			for (i=0;i<p;i++) psi[i] = omega[g]*delta[i];
	    	log_detpsi[g] = p*log(omega[g]);
			log_detsig[g] = update_det_sigma_NEW2(lambda, psi, log_detpsi[g], p, q);
			log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];
		}

	    update_z9(v, x, z, lambda, omega, delta, mu, pi, max_v, log_c, N, G,p,q);
        known_z(cls,z,N,G);
		
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
        it++;  
    }
	/*printf("ll=%f\n",l[it-1]);*/

    /* no of parameters */
    a = G-1 + G*p+ p*q - q*(q-1)/2 + G + (p-1);
    
    /* BIC = 2log_likelihood - mlog n; */
    bic = 2.0*l[it-1] - a*log(N);
    
	lambda_store(lam_vec, lambda, p, q);
    
    /* Deallocate memory */
     free(lambda); free(mu); free(v);  free(n); free(log_c);
      free(l); free(at); free(pi); free(log_detpsi); free(delta); free(log_detsig);
        for(g=0; g < G; g++) {
           free(beta[g]);
           free(theta[g]);
           free(sampcov[g]);
        }

       free(beta); free(theta); free(sampcov);

	
	return bic;
}

double aecm10(double *z, double *x, int *cls, int q, int p, int G, int N, double *lam_vec, double *omega, double tol){

        double bic,a;
	int g,i,it=0,stop=0;
	       
	/* printf("G=%d \t q=%d\n",G,q);*/

        double *max_v = malloc(sizeof(double)*N);
        double *v = malloc(sizeof(double)*N*G);
        double *log_detpsi = malloc(sizeof(double)*G);
        double *log_detsig = malloc(sizeof(double)*G);
        double *log_c = malloc(sizeof(double)*G);
        double *pi = malloc(sizeof(double)*G);
        double *n = malloc(sizeof(double)*G);
        double *at = malloc(sizeof(double)*150000);
        double *l = malloc(sizeof(double)*150000);
        double **sampcov    = malloc(sizeof(double*)*G);
        double **lambda    = malloc(sizeof(double*)*G);
        double **beta    = malloc(sizeof(double*)*G);
        double **theta    = malloc(sizeof(double*)*G);
        for(g=0; g < G; g++) {
           sampcov[g]      = malloc(sizeof(double)*p*p);
           lambda[g]      = malloc(sizeof(double)*p*q);
           beta[g]      = malloc(sizeof(double)*q*p);
           theta[g]      = malloc(sizeof(double)*q*q);
        }
        double *mu = malloc(sizeof(double)*G*p);
        double *delta = malloc(sizeof(double)*p);
        double *psi = malloc(sizeof(double)*p);

	
	get_data2(lam_vec,lambda,G,p,q);
	for(i=0;i<p;i++) delta[i] = 1.0; 
	
   	while(stop==0){            
           
        update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	
	    update_mu(mu, n, x, z, G, N, p); 
                          
        if(it>0){ 
			update_z10(v, x, z, lambda, omega, delta, mu, pi, max_v, log_c, N, G,p,q);
            known_z(cls,z,N,G);
		}	
		
	    update_sg(sampcov, x, z, mu, n, p, G, N);

	    for(g=0;g<G;g++){ 
			for (i=0;i<p;i++) psi[i] = omega[g]*delta[i];
			update_beta2(beta[g], psi, lambda[g], p, q);    
		}
	
    	for (g=0;g<G;g++) update_theta(theta[g], beta[g], lambda[g], sampcov[g], p, q);
        
	    for(g=0;g<G;g++) update_lambda(lambda[g], beta[g], sampcov[g], theta[g], p, q);
	    
        for(g=0;g<G;g++) omega[g] = update_omega2(lambda[g], delta, beta[g], sampcov[g], p, q);
	   	
	    update_delta2(delta, lambda, omega, beta, sampcov, theta, n, p, q, N, G); 
	    
        for(g=0;g<G;g++){
			for (i=0;i<p;i++) psi[i] = omega[g]*delta[i];
			log_detpsi[g] = p*log(omega[g]);
	    	log_detsig[g] = update_det_sigma_NEW2(lambda[g], psi, log_detpsi[g], p, q);
			log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];
		}

	    update_z10(v, x, z, lambda, omega, delta, mu, pi, max_v, log_c, N, G,p,q);
        known_z(cls,z,N,G);
		
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
        it++;
   	}
	/*printf("ll=%f\n",l[it-1]);*/
   
   	a = G-1 + G*p+ G*(p*q - q*(q-1)/2) + G + (p-1);
    
   	bic = 2.0*l[it-1] - a*log(N);
    
    lambda_storeG(lam_vec, lambda, G, p, q);
    free(mu); free(v);  free(n); free(log_c); free(l); free(at); free(pi); free(log_detpsi); free(delta); free(log_detsig);    

        for(g=0; g < G; g++) {
           free(beta[g]);
           free(lambda[g]);
           free(theta[g]);
           free(sampcov[g]);
        }

      free(beta); free(lambda); free(theta); free(sampcov);

	return bic;
}

double aecm11(double *z, double *x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){

        double bic, log_detpsi, a, omega;
	int g,i,j,it=0,stop=0;

	/* printf("G=%d \t q=%d\n",G,q);*/
        double *max_v = malloc(sizeof(double)*N);
        double *v = malloc(sizeof(double)*N*G);
        double *pi = malloc(sizeof(double)*G);
        double *n = malloc(sizeof(double)*G);
        double *at = malloc(sizeof(double)*150000);
        double *l = malloc(sizeof(double)*150000);
        double *lambda = malloc(sizeof(double)*p*q);
        double **sampcov    = malloc(sizeof(double*)*G);
        double **beta    = malloc(sizeof(double*)*G);
        double **theta    = malloc(sizeof(double*)*G);
        for(g=0; g < G; g++) {
           sampcov[g]      = malloc(sizeof(double)*p*p);
           beta[g]      = malloc(sizeof(double)*q*p);
           theta[g]      = malloc(sizeof(double)*q*q);
        }
        double *mu = malloc(sizeof(double)*G*p);
        double *delta = malloc(sizeof(double)*G*p);
        double *log_detsig = malloc(sizeof(double)*G);
        double *log_c = malloc(sizeof(double)*G);
        double *psi = malloc(sizeof(double)*p);
        double *delta0 = malloc(sizeof(double)*p);

	
	omega = *psi_vec;
	get_data(lam_vec,lambda,p,q);
	
	for(g=0;g<G;g++)
		for(i=0;i<p;i++) 
			delta[g*p+i] = 1.0; 
	
   	while(stop==0){            
           
        update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	
	    update_mu(mu, n, x, z, G, N, p); 
                          
        if(it>0){ 
			update_z11(v, x, z, lambda, omega, delta, mu, pi, max_v, log_c, N, G,p,q);
            known_z(cls,z,N,G);
		}	
		
	    update_sg(sampcov, x, z, mu, n, p, G, N);

	    for(g=0;g<G;g++){ 
			for (i=0;i<p;i++) psi[i] = omega*delta[g*p+i];
			update_beta2(beta[g], psi, lambda, p, q);    
		}	
	
    	for (g=0;g<G;g++) update_theta(theta[g], beta[g], lambda, sampcov[g], p, q);
        
	    update_lambda_cuu(lambda, beta, sampcov, theta, n, delta, p, q, G);
	    
        omega =0.0;
	    for(g=0;g<G;g++){
                for(j=0; j<p; j++){
                   delta0[j] = delta[g*p+j];
                }

               omega += pi[g]*update_omega(lambda, delta0, beta[g], sampcov[g], theta[g], p, q);
 	    }	
	    for(g=0;g<G;g++){
               for(j=0; j<p; j++){
                   delta0[j] = delta[g*p+j];
                }

                update_delta3(delta0, lambda, omega, beta[g], sampcov[g], theta[g], n[g], p, q);
               for(j=0; j<p; j++){
                  delta[g*p+j] = delta0[j];
               }

	     }	
		log_detpsi = p*log(omega);

        for(g=0;g<G;g++){
			for (i=0;i<p;i++) psi[i] = omega*delta[g*p+i];
	    	log_detsig[g] = update_det_sigma_NEW2(lambda, psi, log_detpsi, p, q);
			log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];
		}

	    update_z11(v, x, z, lambda, omega, delta, mu, pi, max_v, log_c, N, G,p,q);
        known_z(cls,z,N,G);

	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
        it++;
    }
	/*printf("ll=%f\n",l[it-1]);*/

   	a = G-1 + G*p+ p*q - q*(q-1)/2 + 1 + G*(p-1);
    
   	bic = 2.0*l[it-1] - a*log(N);
    
	lambda_store(lam_vec, lambda, p, q);
    

        free(lambda); free(mu); free(v);  free(n); free(log_c); free(l); free(at); free(pi); free(delta);  free(log_detsig);
        free(delta0);
        for(g=0; g < G; g++) {
           free(beta[g]);
           free(theta[g]);
           free(sampcov[g]);
        }
        free(beta); free(theta); free(sampcov);

	
	return bic;
} 

double aecm12(double *z, double *x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){

        double bic, omega,  log_detpsi, a;
	int g,i,j,it=0,stop=0;
	       
	/* printf("G=%d \t q=%d\n",G,q);*/
        double *max_v = malloc(sizeof(double)*N);
        double *v = malloc(sizeof(double)*N*G);
        double *log_detsig = malloc(sizeof(double)*G);
        double *log_c = malloc(sizeof(double)*G);
        double *pi = malloc(sizeof(double)*G);
        double *n = malloc(sizeof(double)*G);
        double *at = malloc(sizeof(double)*150000);
        double *l = malloc(sizeof(double)*150000);
        double **sampcov    = malloc(sizeof(double*)*G);
        double **lambda    = malloc(sizeof(double*)*G);
        double **beta    = malloc(sizeof(double*)*G);
        double **theta    = malloc(sizeof(double*)*G);
        for(g=0; g < G; g++) {
           sampcov[g]      = malloc(sizeof(double)*p*p);
           lambda[g]      = malloc(sizeof(double)*p*q);
           beta[g]      = malloc(sizeof(double)*q*p);
           theta[g]      = malloc(sizeof(double)*q*q);
        }
        double *mu = malloc(sizeof(double)*G*p);
        double *delta = malloc(sizeof(double)*G*p);
        double *psi = malloc(sizeof(double)*p);
        double *delta0 = malloc(sizeof(double)*p);


	
	omega = *psi_vec;
	get_data2(lam_vec,lambda,G,p,q);
	
	for(g=0;g<G;g++)
		for(i=0;i<p;i++) 
			delta[g*p+i] = 1.0; 
	
   	while(stop==0){            
           
		update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	
	    update_mu(mu, n, x, z, G, N, p); 
                          
	    if(it>0){
			update_z12(v, x, z, lambda, omega, delta, mu, pi, max_v, log_c, N, G,p,q);
            known_z(cls,z,N,G);
		}	
		
	    update_sg(sampcov, x, z, mu, n, p, G, N);
	    
	    for(g=0;g<G;g++){ 
			for (i=0;i<p;i++) psi[i] = omega*delta[g*p+i];
			update_beta2(beta[g], psi, lambda[g], p, q);    
		}	
	
    	for (g=0;g<G;g++) update_theta(theta[g], beta[g], lambda[g], sampcov[g], p, q);
        
    	for (g=0;g<G;g++) update_lambda(lambda[g], beta[g], sampcov[g], theta[g], p, q);
		    
        omega =0.0;
	    for(g=0;g<G;g++){
               for(j=0; j<p; j++){
                   delta0[j] = delta[g*p+j];
                }

               omega += pi[g]*update_omega2(lambda[g], delta0, beta[g], sampcov[g], p, q);
	    }
	    for(g=0;g<G;g++){
               for(j=0; j<p; j++){
                   delta0[j] = delta[g*p+j];
                }

                update_delta3(delta0, lambda[g], omega, beta[g], sampcov[g], theta[g], n[g], p, q); 
              for(j=0; j<p; j++){
                 delta[g*p+j] = delta0[j];
              }

	     }
		log_detpsi = p*log(omega);        

		for(g=0;g<G;g++){
			for (i=0;i<p;i++) psi[i] = omega*delta[g*p+i];
			log_detsig[g] = update_det_sigma_NEW2(lambda[g], psi, log_detpsi, p, q);
			log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];
		}

	    update_z12(v, x, z, lambda, omega, delta, mu, pi, max_v, log_c, N, G,p,q);
        known_z(cls,z,N,G);

	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
        it++;
    }
	/*printf("ll=%f\n",l[it-1]);*/
   
   	a = G-1 + G*p+ G*(p*q - q*(q-1)/2) + 1+ G*(p-1);
    
   	bic = 2.0*l[it-1] - a*log(N);
    
    lambda_storeG(lam_vec, lambda, G, p, q);
    
        free(mu); free(v); free(n); free(l); free(at); free(pi); free(delta);
        free(log_c); free(log_detsig);


        for(g=0; g < G; g++) {
           free(beta[g]);
           free(theta[g]);
           free(lambda[g]);
           free(sampcov[g]);
        }
       free(beta); free(theta); free(lambda); free(sampcov);


	return bic;
}

double claecm(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){
    double bic;
	int it=0,stop=0,paras;
	double psi, log_detpsi, log_c=0.0, log_detsig;

        double *pi = malloc(sizeof(double)*G);
        double *n = malloc(sizeof(double)*G);
        double *at = malloc(sizeof(double)*150000);
        double *l = malloc(sizeof(double)*150000);
        double *sampcovtilde = malloc(sizeof(double)*p*p);
        double *max_v = malloc(sizeof(double)*N);
        double *v = malloc(sizeof(double)*N*G);	
        double *lambda = malloc(sizeof(double)*p*q);	
        double *beta = malloc(sizeof(double)*q*p);	
        double *theta = malloc(sizeof(double)*q*q);	
        double *mu = malloc(sizeof(double)*G*p);	

	psi = *psi_vec;	
	get_data(lam_vec,lambda,p,q);
	
    while(stop==0){            
    
		update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	  
	    update_mu(mu, n, x, z, G, N, p);

		if(it>0){ 
			update_z(v, x, z, lambda, psi, mu, pi, max_v, log_c, N, G, p, q);
		}

	    update_stilde(sampcovtilde,x, z, mu, G, N, p);

	    update_beta1(beta, psi, lambda, p, q);
        
   	    update_theta(theta, beta, lambda, sampcovtilde, p, q); 
    
        update_lambda(lambda, beta, sampcovtilde, theta, p, q);
 
        psi = update_psi(lambda, beta, sampcovtilde, p, q); 
   
 	    log_detpsi=p*log(psi);
	    
	    log_detsig = update_det_sigma_NEW(lambda, psi, log_detpsi, p, q);
	  
	    log_c = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig;
       
		update_z(v, x, z, lambda, psi, mu, pi, max_v, log_c, N, G, p, q);
        stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++;
		/*printf("ll=%f\t",l[it-1]);*/
	}
	/*printf("ll=%f\n",l[it-1]);*/

	paras = G-1 + G*p + p*q - q*(q-1)/2 + 1;
    	
   	bic = 2.0*l[it-1] - paras*log(N);

    lambda_store(lam_vec, lambda, p, q);
    psi_vec[0]=psi;
    
    free(lambda); free(mu); free(n); free(beta); free(theta); free(sampcovtilde);
    free(l); free(at); free(pi);

       
	return bic;
} 

double claecm2(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *Psi, double tol){

        double bic;
	int i,it=0,stop=0;
	double a, log_detpsi,log_detsig, log_c=0.0;

	/* printf("G=%d \t q=%d\n",G,q);*/

   	/* Allocate memory - matrices and vectors */ 
        double *det = malloc(sizeof(double)*G);
        double *pi = malloc(sizeof(double)*G);
        double *n = malloc(sizeof(double)*G);
        double *at = malloc(sizeof(double)*150000);
        double *l = malloc(sizeof(double)*150000); 
        double *sampcovtilde = malloc(sizeof(double)*p*p);
        double *lambda = malloc(sizeof(double)*p*q);
        double *beta = malloc(sizeof(double)*q*p); 
        double *theta = malloc(sizeof(double)*q*q);
        double *mu = malloc(sizeof(double)*G*p);
        double *w = malloc(sizeof(double)*G*N);
        double *max_v = malloc(sizeof(double)*N);
        double *v = malloc(sizeof(double)*N*G);
                

	get_data(lam_vec,lambda,p,q);
	
   	while(stop==0){            
	    
	    update_n(n, z, G, N);
	    
	    update_pi(pi, n, G, N);
	    
	    update_mu(mu, n, x, z, G, N, p); 

        if (it>0){ 
			update_z2(v, x, z, lambda, Psi, mu, pi, max_v, log_c, N, G, p, q);
		}
        
		update_stilde(sampcovtilde, x, z, mu, G, N, p);        
        
	    update_beta2(beta, Psi, lambda, p, q);

	    update_theta(theta, beta, lambda, sampcovtilde, p, q);

	    update_lambda(lambda, beta, sampcovtilde, theta, p, q);
			    
	    update_psi2(Psi, lambda, beta, sampcovtilde, p, q);

		log_detpsi = 0.0; for(i=0;i<p;i++) log_detpsi += log(Psi[i]);

	    log_detsig = update_det_sigma_NEW2(lambda, Psi, log_detpsi, p, q);
	    
		log_c = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig; 
   
		update_z2(v, x, z, lambda, Psi, mu, pi, max_v, log_c, N, G, p, q);
 
        stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	
	    it++; 
    }
	/*printf("ll=%f\n",l[it-1]);*/
    
	/* no of parameters {pq - q(q - 1)/2} + p */
    a = G-1 + G*p + p*q - q*(q-1)/2 + p;
    
    /* BIC = 2log_likelihood - mlog n; */
    bic = 2.0*l[it-1] - a*log(N);
    
	lambda_store(lam_vec, lambda, p, q);
    
    /* Deallocate memory */
    free(lambda); free(mu); free(w);  free(n); free(det); free(beta);
    free(theta); free(sampcovtilde); free(l); free(at); free(pi);

	return bic;
}

double claecm3(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *Psi, double tol){

        double bic, a;
	int g,it=0,stop=0;

	/* printf("G=%d \t q=%d\n",G,q);*/
   double *log_det = malloc(sizeof(double)*G);
   double *log_detpsi = malloc(sizeof(double)*G);
   double *log_detsig = malloc(sizeof(double)*G);
   double *pi = malloc(sizeof(double)*G);
   double *n = malloc(sizeof(double)*G);
   double *at = malloc(sizeof(double)*150000);
   double *l = malloc(sizeof(double)*150000);
   double *lambda = malloc(sizeof(double)*p*q);
   double **sampcov    = malloc(sizeof(double*)*G);
   double **beta    = malloc(sizeof(double*)*G);
   double **theta    = malloc(sizeof(double*)*G);
   for(g=0; g < G; g++) {
      sampcov[g]      = malloc(sizeof(double)*p*p);
      beta[g]      = malloc(sizeof(double)*q*p);
      theta[g]      = malloc(sizeof(double)*q*q);
    } 
   double *mu = malloc(sizeof(double)*G*p);
   double *max_v = malloc(sizeof(double)*N);
   double *v = malloc(sizeof(double)*N*G);


	get_data(lam_vec,lambda,p,q);
	
    while(stop==0){            
           
        update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	
	    update_mu(mu, n, x, z, G, N, p); 
                          
		if(it>0){ 
			update_z3(v, x, z, lambda, Psi, mu, pi, max_v, log_det, N, G, p, q);
        }

	    update_sg(sampcov, x, z, mu, n, p, G, N);
	    
	    for(g=0;g<G;g++) update_beta1(beta[g], Psi[g], lambda, p, q);    
	
    	for (g=0;g<G;g++) update_theta(theta[g], beta[g], lambda, sampcov[g], p, q);
        
	    update_lambda2(lambda, beta, sampcov, theta, n, Psi, p, q, G);
	    
        for(g=0;g<G;g++) Psi[g] = update_psi3(lambda, beta[g], sampcov[g], theta[g], p, q);
    
        for(g=0;g<G;g++){
	    	log_detpsi[g] = p*log(Psi[g]);
			log_detsig[g] = update_det_sigma_NEW(lambda, Psi[g], log_detpsi[g], p, q);
			log_det[g] = (p/2.0)*log(2*M_PI) + 0.5*log_detsig[g];
		} 
	    
		update_z3(v, x, z, lambda, Psi, mu, pi, max_v, log_det, N, G, p, q);

	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
		it++;
    }
	
	/*printf("ll=%f\n",l[it-1]);*/
	
   	a = G-1 + G*p+ p*q - q*(q-1)/2 + G;
    
   	bic = 2.0*l[it-1] - a*log(N);
    
	lambda_store(lam_vec, lambda, p, q);
    
   	/* Deallocate memory */

        free(lambda); free(mu); free(v); free(n); free(log_det); free(max_v);
        free(l); free(at); free(pi); free(log_detpsi);
        for(g=0; g < G; g++) {
           free(beta[g]);
           free(theta[g]);
           free(sampcov[g]);
        }

        free(beta); free(theta); free(sampcov);
    
	return bic;
}

double claecm4(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){

        double bic;
	int it=0,stop=0,paras, g, j;
	
   	/* Allocate memory - matrices and vectors */ 
        double *max_v = malloc(sizeof(double)*N);
        double *v = malloc(sizeof(double)*N*G);
        double *pi = malloc(sizeof(double)*G);
        double *n = malloc(sizeof(double)*G);
        double *at = malloc(sizeof(double)*150000);
        double *l = malloc(sizeof(double)*150000);
        double *lambda = malloc(sizeof(double)*p*q);
        double **sampcov    = malloc(sizeof(double*)*G);
        double **beta    = malloc(sizeof(double*)*G);
        double **theta    = malloc(sizeof(double*)*G);
        for(g=0; g < G; g++) {
           sampcov[g]      = malloc(sizeof(double)*p*p);
           beta[g]      = malloc(sizeof(double)*q*p);
           theta[g]      = malloc(sizeof(double)*q*q);
        }
        double *mu = malloc(sizeof(double)*G*p);
        double *w = malloc(sizeof(double)*G*N);
        double *Psi = malloc(sizeof(double)*G*p);
        double *log_detpsi = malloc(sizeof(double)*G);
        double *log_detsig = malloc(sizeof(double)*G);
        double *log_c = malloc(sizeof(double)*G);
        double *Psi0 = malloc(sizeof(double)*p);

	
	get_data(psi_vec,Psi,G,p);
	get_data(lam_vec,lambda,p,q);
        
   	while(stop==0){            
   
	    update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	
	    update_mu(mu, n, x, z, G, N, p); 

        if(it>0){ 
			update_z4(v, x, z, lambda, Psi, mu, pi, max_v, log_c, N, G, p, q);
		} 
 
	    update_sg(sampcov,x, z, mu, n, p, G, N);
        
               
	    for(g=0;g<G;g++){
               for(j=0; j<p; j++){
                   Psi0[j] = Psi[g*p+j];
               }

                update_beta2(beta[g], Psi0, lambda, p, q);
            }
   	    for(g=0;g<G;g++) update_theta(theta[g], beta[g], lambda, sampcov[g], p, q); 
    
	    update_lambda_cuu(lambda, beta, sampcov, theta, n, Psi, p, q, G);
	 
        update_psi_cuu(Psi, lambda, beta, sampcov, theta, p, q, G);

	    for (g=0;g<G;g++){
			log_detpsi[g]=0.0;
			for(j=0;j<p;j++) log_detpsi[g] += log(Psi[g*p+j]);
	    }
	    for(g=0;g<G;g++){
                for(j=0; j<p; j++){
                   Psi0[j] = Psi[g*p+j];
                }
                 log_detsig[g] = update_det_sigma_NEW2(lambda, Psi0, log_detpsi[g], p, q);
            }

		for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];

        update_z4(v, x, z, lambda, Psi, mu, pi, max_v, log_c, N, G, p, q);
            
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++;
	}
	
	/* no of parameters {pq - q(q - 1)/2} + Gp */
	paras = G-1 + G*p + p*q - q*(q-1)/2 + G*p;
/*new section*/	
    
   	/* BIC = 2log_likelihood - mlog n; */
   	bic = 2.0*l[it-1] - paras*log(N);
    
    lambda_store(lam_vec, lambda, p, q);
    /*{int i,k=0;
    for(g=0;g<G;g++){
        for (i=0;i<p;i++){
            psi_vec[k]=Psi[g][i];
            k++;
            Rprintf("%f ",Psi[g][i]);
        }Rprintf("\n");
    }Rprintf("\n");
    }*/
   lambda_store(psi_vec,Psi,G,p);
     
    /* Deallocate memory */
        free(lambda); free(mu); free(w); free(n); free(l); free(at); free(pi);
        free(log_detsig); free(log_c); free(log_detpsi); free(Psi0);
        free(max_v); free(v); free(Psi);

        for(g=0; g < G; g++) {
           free(beta[g]);
           free(theta[g]);
           free(sampcov[g]);
        }
        free(beta); free(theta); free(sampcov);

	return bic;
} 

double claecm5(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){

        double bic, log_detpsi, psi;
	int it=0,stop=0,paras, g,i;

	/*Rprintf("G=%d \t q=%d\n",G,q);*/
        double *max_v = malloc(sizeof(double)*N);
        double *v = malloc(sizeof(double)*N*G);
        double *log_c = malloc(sizeof(double)*G);
        double *log_detsig = malloc(sizeof(double)*G);
        double *pi = malloc(sizeof(double)*G);
        double *n = malloc(sizeof(double)*G);
        double *at = malloc(sizeof(double)*150000);
        double *l = malloc(sizeof(double)*150000);
        double **sampcov    = malloc(sizeof(double*)*G);
        double **lambda    = malloc(sizeof(double*)*G);
        double **beta    = malloc(sizeof(double*)*G);
        double **theta    = malloc(sizeof(double*)*G);
        for(g=0; g < G; g++) {
           sampcov[g]      = malloc(sizeof(double)*p*p);
           lambda[g]      = malloc(sizeof(double)*p*q);
           beta[g]      = malloc(sizeof(double)*q*p);
           theta[g]      = malloc(sizeof(double)*q*q);
        }
        double *mu = malloc(sizeof(double)*G*p);
        double *w = malloc(sizeof(double)*G*N);


	psi=*psi_vec;
	get_data2(lam_vec,lambda,G,p,q);
       
   	while(stop==0){            
    
	    update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	
	    update_mu(mu, n, x, z, G, N, p); 

		if(it>0){ 
			update_z5(v, x, z, lambda, psi, mu, pi, max_v, log_c, N, G, p, q);
        }
 
	    update_sg(sampcov,x, z, mu, n, p, G, N);
        
	    for(g=0;g<G;g++) update_beta1(beta[g], psi, lambda[g], p, q);
        
   	    for(g=0;g<G;g++) update_theta(theta[g], beta[g], lambda[g], sampcov[g], p, q); 
    
        for (g=0;g<G;g++) update_lambda(lambda[g], beta[g], sampcov[g], theta[g], p, q);
	 
	    psi = update_psi_ucc(lambda, beta, sampcov, p, q, pi, G);
   
		log_detpsi = 0.0;
		for(i=0;i<p;i++) log_detpsi += log(psi);
	    
	    for(g=0;g<G;g++) log_detsig[g] = update_det_sigma_NEW(lambda[g], psi, log_detpsi, p, q);
	    
		for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];
       
		update_z5(v, x, z, lambda, psi, mu, pi, max_v, log_c, N, G, p, q);
            
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++; 
	}
	/*printf("ll=%f\n",l[it-1]);*/

	paras = G-1 + G*p + G*(p*q - q*(q-1)/2) + 1;
    
   	bic = 2.0*l[it-1] - paras*log(N);
    
    lambda_storeG(lam_vec, lambda, G, p, q);
    /*Rprintf("psi=%f\n",psi);*/
    psi_vec[0]=psi;
    free(mu); free(w); free(n); free(l); free(at); free(pi); free(log_detsig); free(log_c);
        for(g=0; g < G; g++) {
           free(beta[g]);
           free(lambda[g]);
           free(theta[g]);
           free(sampcov[g]);
        }

        free(beta); free(lambda); free(theta); free(sampcov);
    
	return bic;
}

double claecm6(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *psi, double tol){

 
        double bic, log_detpsi;
	int it=0,stop=0,paras, g, j;

	/* printf("G=%d \t q=%d\n",G,q);*/
        double *max_v = malloc(sizeof(double)*N);
        double *v = malloc(sizeof(double)*N*G);
        double *log_detsig = malloc(sizeof(double)*G);
        double *log_c = malloc(sizeof(double)*G);
        double *pi = malloc(sizeof(double)*G);
        double *n = malloc(sizeof(double)*G);
        double *at = malloc(sizeof(double)*150000);
        double *l = malloc(sizeof(double)*150000);
        double **sampcov    = malloc(sizeof(double*)*G);
        double **lambda    = malloc(sizeof(double*)*G);
        double **beta    = malloc(sizeof(double*)*G);
        double **theta    = malloc(sizeof(double*)*G);
        for(g=0; g < G; g++) {
           sampcov[g]      = malloc(sizeof(double)*p*p);
           lambda[g]      = malloc(sizeof(double)*p*q);
           beta[g]      = malloc(sizeof(double)*q*p);
           theta[g]      = malloc(sizeof(double)*q*q);
        }
        double *mu = malloc(sizeof(double)*G*p);

	
	get_data2(lam_vec,lambda,G,p,q);
  
   	while(stop==0){            
    
	    update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	
	    update_mu(mu, n, x, z, G, N, p); 

        if(it>0){ 
			update_z6(v, x, z, lambda, psi, mu, pi, max_v, log_c, N, G, p, q);
		} 

	    update_sg(sampcov,x, z, mu, n, p, G, N);
        
	    for(g=0;g<G;g++) update_beta2(beta[g], psi, lambda[g], p, q);
        
   	    for(g=0;g<G;g++) update_theta(theta[g], beta[g], lambda[g], sampcov[g], p, q); 
    
        for (g=0;g<G;g++) update_lambda(lambda[g], beta[g], sampcov[g], theta[g], p, q);
	 
	    update_psi_ucu(psi, lambda, beta, sampcov, p, q, pi, G);
                    
		log_detpsi = 0.0;
		for(j=0;j<p;j++) log_detpsi += log(psi[j]);
	    
	    for(g=0;g<G;g++) log_detsig[g] = update_det_sigma_NEW2(lambda[g], psi, log_detpsi, p, q);
	    
		for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];

        update_z6(v, x, z, lambda, psi, mu, pi, max_v, log_c, N, G, p, q);
            
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++; 
	}
	/*printf("ll=%f\n",l[it-1]);*/
	
	paras = G-1 + G*p + G*(p*q - q*(q-1)/2) + p;
    
   	bic = 2.0*l[it-1] - paras*log(N);
    
    lambda_storeG(lam_vec, lambda, G, p, q);
    
        free(mu); free(v); free(n); free(max_v); free(l); free(at); free(pi); free(log_detsig); free(log_c);
        for(g=0; g < G; g++) {
           free(beta[g]);
           free(lambda[g]);
           free(theta[g]);
           free(sampcov[g]);
        }
        free(beta); free(lambda); free(theta); free(sampcov);
	return bic;
}

double claecm7(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *psi, double tol){

        double bic;
	int it=0,stop=0,paras, g;

	/* printf("G=%d \t q=%d\n",G,q); */
        double *max_v = malloc(sizeof(double)*N);
        double *v = malloc(sizeof(double)*N*G);
        double *log_detpsi = malloc(sizeof(double)*G);
        double *log_detsig = malloc(sizeof(double)*G);
        double *log_c = malloc(sizeof(double)*G);
        double *pi = malloc(sizeof(double)*G);
        double *n = malloc(sizeof(double)*G);
        double *at = malloc(sizeof(double)*150000);
        double *l = malloc(sizeof(double)*150000);
        double **sampcov    = malloc(sizeof(double*)*G);
        double **lambda    = malloc(sizeof(double*)*G);
        double **beta    = malloc(sizeof(double*)*G);
        double **theta    = malloc(sizeof(double*)*G);
        for(g=0; g < G; g++) {
           sampcov[g]      = malloc(sizeof(double)*p*p);
           lambda[g]      = malloc(sizeof(double)*p*q);
           beta[g]      = malloc(sizeof(double)*q*p);
           theta[g]      = malloc(sizeof(double)*q*q);
        }
        double *mu = malloc(sizeof(double)*G*p);
	
	get_data2(lam_vec,lambda,G,p,q);

   	while(stop==0){            
    
	    update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	
	    update_mu(mu, n, x, z, G, N, p); 

        if(it>0){ 
			update_z7(v, x, z, lambda, psi, mu, pi, max_v, log_c, N, G, p, q);
		}

	    update_sg(sampcov,x, z, mu, n, p, G, N);
        
	    for(g=0;g<G;g++) update_beta1(beta[g], psi[g], lambda[g], p, q);
        
   	    for(g=0;g<G;g++) update_theta(theta[g], beta[g], lambda[g], sampcov[g], p, q); 
    
        for (g=0;g<G;g++) update_lambda(lambda[g], beta[g], sampcov[g], theta[g], p, q);
	 
        for (g=0;g<G;g++) psi[g]=update_psi(lambda[g], beta[g], sampcov[g], p, q); 
           
	    for(g=0;g<G;g++) log_detpsi[g] = p*log(psi[g]);
	    
	    for(g=0;g<G;g++) log_detsig[g] = update_det_sigma_NEW(lambda[g], psi[g], log_detpsi[g], p, q);
	   
		for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];

		update_z7(v, x, z, lambda, psi, mu, pi, max_v, log_c, N, G, p, q);
            
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++;
	}
	/*printf("ll=%f\n",l[it-1]);*/

	paras = G-1 + G*p + G*(p*q - q*(q-1)/2) + G;
    
   	bic = 2.0*l[it-1] - paras*log(N);
    
    lambda_storeG(lam_vec, lambda, G, p, q);
    free(mu); free(v); free(n); free(l); free(at); free(pi);  free(log_detpsi);
     free(log_detsig); free(log_c);
        for(g=0; g < G; g++) {
           free(beta[g]);
           free(lambda[g]);
           free(theta[g]);
           free(sampcov[g]);
        }
    
        free(beta); free(lambda); free(theta); free(sampcov);	
	return bic;
}

double claecm8(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){

        double bic;
	int it=0,stop=0,paras, g, j;

	/* printf("G=%d \t q=%d\n",G,q); */
        double *max_v = malloc(sizeof(double)*N);
        double *v = malloc(sizeof(double)*N*G);
        double *pi = malloc(sizeof(double)*G);
        double *n = malloc(sizeof(double)*G);
        double *at = malloc(sizeof(double)*150000);
        double *l = malloc(sizeof(double)*150000);
        double **sampcov    = malloc(sizeof(double*)*G);
        double **lambda    = malloc(sizeof(double*)*G);
        double **beta    = malloc(sizeof(double*)*G);
        double **theta    = malloc(sizeof(double*)*G);
        for(g=0; g < G; g++) {
           sampcov[g]      = malloc(sizeof(double)*p*p);
           lambda[g]      = malloc(sizeof(double)*p*q);
           beta[g]      = malloc(sizeof(double)*q*p);
           theta[g]      = malloc(sizeof(double)*q*q);
        }
        double *mu = malloc(sizeof(double)*G*p);
        double *Psi = malloc(sizeof(double)*G*p);
        double *log_detpsi = malloc(sizeof(double)*G);
        double *log_detsig = malloc(sizeof(double)*G);
        double *log_c = malloc(sizeof(double)*G);
        double *Psi0 = malloc(sizeof(double)*p);


	
	get_data(psi_vec,Psi,G,p);
	get_data2(lam_vec,lambda,G,p,q);
	
   	while(stop==0){            

    
	    update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	
	    update_mu(mu, n, x, z, G, N, p);

        if(it>0){ 
			update_z8(v, x, z, lambda, Psi, mu, pi, max_v, log_c, N, G, p, q);
		}	
       
	    update_sg(sampcov,x, z, mu, n, p, G, N);
	    
	    for(g=0;g<G;g++){
                for(j=0; j<p; j++){
                   Psi0[j] = Psi[g*p+j];
                }

              update_beta2(beta[g], Psi0, lambda[g], p, q);
           }
        
   	    for(g=0;g<G;g++) update_theta(theta[g], beta[g], lambda[g], sampcov[g], p, q); 
    
        for (g=0;g<G;g++) update_lambda(lambda[g], beta[g], sampcov[g], theta[g], p, q);
	 
        for (g=0;g<G;g++){

              update_psi2(Psi0, lambda[g], beta[g], sampcov[g], p, q); 
             for(j=0; j<p; j++){
                   Psi[g*p+j] = Psi0[j];
                }
           }
	    for(g=0;g<G;g++){
		    log_detpsi[g]=0.0;
		    for(j=0;j<p;j++) log_detpsi[g] += log(Psi[g*p+j]);
	    }
	    
		for(g=0;g<G;g++){
                   for(j=0; j<p; j++){
                   Psi0[j] = Psi[g*p+j];
                }
                    log_detsig[g] = update_det_sigma_NEW2(lambda[g], Psi0, log_detpsi[g], p, q);
                }
	    	    
		for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];

		update_z8(v, x, z, lambda, Psi, mu, pi, max_v, log_c, N, G, p, q);
            
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++; 
	}
	/*printf("ll=%f\n",l[it-1]);*/
	paras = G-1 + G*p + G*(p*q - q*(q-1)/2) + G*p;
   	
	bic = 2.0*l[it-1] - paras*log(N);
    
    lambda_storeG(lam_vec, lambda, G, p, q);
    
    free(mu); free(v); free(n); free(log_detpsi); free(l); free(at); free(pi);
    free(log_detsig); free(log_c); free(Psi); free(max_v); free(Psi0);

        for(g=0; g < G; g++) {
           free(beta[g]);
           free(lambda[g]);
           free(theta[g]);
           free(sampcov[g]);
        }
        free(beta); free(lambda); free(theta); free(sampcov);
	return bic;
}

double claecm9(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *omega, double tol){

        double bic,a;
	int g,i,it=0,stop=0;
	       
	/* printf("G=%d \t q=%d\n",G,q);*/

   	/* Allocate memory - matrices and vectors */ 
        double *max_v = malloc(sizeof(double)*N);
        double *v = malloc(sizeof(double)*N*G);
        double *log_detpsi = malloc(sizeof(double)*G);
        double *log_detsig = malloc(sizeof(double)*G);
        double *log_c = malloc(sizeof(double)*G);
        double *pi = malloc(sizeof(double)*G);
        double *n = malloc(sizeof(double)*G);
        double *at = malloc(sizeof(double)*150000);
        double *l = malloc(sizeof(double)*150000);
        double *lambda = malloc(sizeof(double)*p*q);
        double **sampcov    = malloc(sizeof(double*)*G);
        double **beta    = malloc(sizeof(double*)*G);
        double **theta    = malloc(sizeof(double*)*G);
        for(g=0; g < G; g++) {
           sampcov[g]      = malloc(sizeof(double)*p*p);
           beta[g]      = malloc(sizeof(double)*q*p);
           theta[g]      = malloc(sizeof(double)*q*q);
        }
        double *mu = malloc(sizeof(double)*G*p);
        double *delta = malloc(sizeof(double)*p);
        double *psi = malloc(sizeof(double)*p);

	
	get_data(lam_vec,lambda,p,q);
	for(i=0;i<p;i++) delta[i] = 1.0; 
	
    while(stop==0){            
           
        update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	
	    update_mu(mu, n, x, z, G, N, p); 
                          
        if(it>0){ 
			update_z9(v, x, z, lambda, omega, delta, mu, pi, max_v, log_c, N, G,p,q);
		}	
		
	    update_sg(sampcov, x, z, mu, n, p, G, N);

	    for(g=0;g<G;g++){ 
			for (i=0;i<p;i++) psi[i] = omega[g]*delta[i];
			update_beta2(beta[g], psi, lambda, p, q);    
		}	
 

    	for (g=0;g<G;g++) update_theta(theta[g], beta[g], lambda, sampcov[g], p, q);
        
	    update_lambda2(lambda, beta, sampcov, theta, n, omega, p, q, G);
	    
        for(g=0;g<G;g++) omega[g] = update_omega(lambda, delta, beta[g], sampcov[g], theta[g], p, q);

	    update_delta(delta, lambda, omega, beta, sampcov, theta, n, p, q, N, G); 
      
        for(g=0;g<G;g++){
			for (i=0;i<p;i++) psi[i] = omega[g]*delta[i];
	    	log_detpsi[g] = p*log(omega[g]);
			log_detsig[g] = update_det_sigma_NEW2(lambda, psi, log_detpsi[g], p, q);
			log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];
		}

	    update_z9(v, x, z, lambda, omega, delta, mu, pi, max_v, log_c, N, G,p,q);
		
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
        it++;  
    }
	/*printf("ll=%f\n",l[it-1]);*/

    /* no of parameters */
    a = G-1 + G*p+ p*q - q*(q-1)/2 + G + (p-1);
    
    /* BIC = 2log_likelihood - mlog n; */
    bic = 2.0*l[it-1] - a*log(N);
    
	lambda_store(lam_vec, lambda, p, q);
    for(i=0;i<p;i++){omega[G+i]=delta[i];
                     /*Rprintf("omega is %f\n", omega[G+i]);*/
	}
    
    /* Deallocate memory: psi not freed? */
     free(lambda); free(mu); free(v);  free(n); free(log_c);
      free(l); free(at); free(pi); free(log_detpsi); free(delta); free(log_detsig);
        for(g=0; g < G; g++) {
           free(beta[g]);
           free(theta[g]);
           free(sampcov[g]);
        }

        free(beta); free(theta); free(sampcov);
	return bic;
}

double claecm10(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *omega, double tol){

        double bic, a;
	int g,i,it=0,stop=0;
	       
	/* printf("G=%d \t q=%d\n",G,q);*/
        double *max_v = malloc(sizeof(double)*N);
        double *v = malloc(sizeof(double)*N*G);
        double *log_detpsi = malloc(sizeof(double)*G);
        double *log_detsig = malloc(sizeof(double)*G);
        double *log_c = malloc(sizeof(double)*G);
        double *pi = malloc(sizeof(double)*G);
        double *n = malloc(sizeof(double)*G);
        double *at = malloc(sizeof(double)*150000);
        double *l = malloc(sizeof(double)*150000);
        double **sampcov    = malloc(sizeof(double*)*G);
        double **lambda    = malloc(sizeof(double*)*G);
        double **beta    = malloc(sizeof(double*)*G);
        double **theta    = malloc(sizeof(double*)*G);
        for(g=0; g < G; g++) {
           sampcov[g]      = malloc(sizeof(double)*p*p);
           lambda[g]      = malloc(sizeof(double)*p*q);
           beta[g]      = malloc(sizeof(double)*q*p);
           theta[g]      = malloc(sizeof(double)*q*q);
        }
        double *mu = malloc(sizeof(double)*G*p);
        double *delta = malloc(sizeof(double)*p);
        double *psi = malloc(sizeof(double)*p);

	
	get_data2(lam_vec,lambda,G,p,q);
	for(i=0;i<p;i++) delta[i] = 1.0; 
	
   	while(stop==0){            
           
        update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	
	    update_mu(mu, n, x, z, G, N, p); 
                          
        if(it>0){ 
			update_z10(v, x, z, lambda, omega, delta, mu, pi, max_v, log_c, N, G,p,q);
		}	
		
	    update_sg(sampcov, x, z, mu, n, p, G, N);

	    for(g=0;g<G;g++){ 
			for (i=0;i<p;i++) psi[i] = omega[g]*delta[i];
			update_beta2(beta[g], psi, lambda[g], p, q);    
		}
	
    	for (g=0;g<G;g++) update_theta(theta[g], beta[g], lambda[g], sampcov[g], p, q);
        
	    for(g=0;g<G;g++) update_lambda(lambda[g], beta[g], sampcov[g], theta[g], p, q);
	    
        for(g=0;g<G;g++) omega[g] = update_omega2(lambda[g], delta, beta[g], sampcov[g], p, q);
	   	
	    update_delta2(delta, lambda, omega, beta, sampcov, theta, n, p, q, N, G); 
	    
        for(g=0;g<G;g++){
			for (i=0;i<p;i++) psi[i] = omega[g]*delta[i];
			log_detpsi[g] = p*log(omega[g]);
	    	log_detsig[g] = update_det_sigma_NEW2(lambda[g], psi, log_detpsi[g], p, q);
			log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];
		}

	    update_z10(v, x, z, lambda, omega, delta, mu, pi, max_v, log_c, N, G,p,q);
		
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
        it++;
   	}
   
   	a = G-1 + G*p+ G*(p*q - q*(q-1)/2) + G + (p-1);
    
   	bic = 2.0*l[it-1] - a*log(N);
    
    lambda_storeG(lam_vec, lambda, G, p, q);
    for(i=0;i<p;i++){omega[G+i]=delta[i];} 
    
        free(mu); free(v);  free(n); free(log_c); free(l); free(at); free(pi); free(log_detpsi); free(max_v); free(psi);
        free(delta); free(log_detsig);
        for(g=0; g < G; g++) {
           free(lambda[g]);
           free(beta[g]);
           free(theta[g]);
           free(sampcov[g]);
        }
       free(lambda); free(beta); free(theta); free(sampcov);

	return bic;
}

double claecm11(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){

        double bic,  log_detpsi, a, omega;
	int g,i,j,k=1,it=0,stop=0;

	/* printf("G=%d \t q=%d\n",G,q);*/
        double *max_v = malloc(sizeof(double)*N);
        double *v = malloc(sizeof(double)*N*G);
        double *pi = malloc(sizeof(double)*G);
        double *n = malloc(sizeof(double)*G);
        double *at = malloc(sizeof(double)*150000);
        double *l = malloc(sizeof(double)*150000);
        double *lambda = malloc(sizeof(double)*p*q);
        double **sampcov    = malloc(sizeof(double*)*G);
        double **beta    = malloc(sizeof(double*)*G);
        double **theta    = malloc(sizeof(double*)*G);
        for(g=0; g < G; g++) {
           sampcov[g]      = malloc(sizeof(double)*p*p);
           beta[g]      = malloc(sizeof(double)*q*p);
           theta[g]      = malloc(sizeof(double)*q*q);
        }
        double *mu = malloc(sizeof(double)*G*p);
        double *delta = malloc(sizeof(double)*G*p);
        double *log_detsig = malloc(sizeof(double)*G);
        double *log_c = malloc(sizeof(double)*G);
        double *psi = malloc(sizeof(double)*p);
         double *delta0 = malloc(sizeof(double)*p);

	
	omega = *psi_vec;
	get_data(lam_vec,lambda,p,q);
	
	for(g=0;g<G;g++)
		for(i=0;i<p;i++) 
			delta[g*p+i] = 1.0; 
	
   	while(stop==0){            
           
        update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	
	    update_mu(mu, n, x, z, G, N, p); 
                          
        if(it>0){ 
			update_z11(v, x, z, lambda, omega, delta, mu, pi, max_v, log_c, N, G,p,q);
		}	
		
	    update_sg(sampcov, x, z, mu, n, p, G, N);

	    for(g=0;g<G;g++){ 
			for (i=0;i<p;i++) psi[i] = omega*delta[g*p+i];
			update_beta2(beta[g], psi, lambda, p, q);    
		}	
	
    	for (g=0;g<G;g++) update_theta(theta[g], beta[g], lambda, sampcov[g], p, q);
        
	    update_lambda_cuu(lambda, beta, sampcov, theta, n, delta, p, q, G);
	    
        omega =0.0;
	    for(g=0;g<G;g++){
               for(j=0; j<p; j++){
                   delta0[j] = delta[g*p+j];
                }

               omega += pi[g]*update_omega(lambda, delta0, beta[g], sampcov[g], theta[g], p, q);
            }	   	
	    for(g=0;g<G;g++){
               for(j=0; j<p; j++){
                   delta0[j] = delta[g*p+j];
                }

                update_delta3(delta0, lambda, omega, beta[g], sampcov[g], theta[g], n[g], p, q); 
               for(j=0; j<p; j++){
                  delta[g*p+j] = delta0[j];
               }
	    }	
		log_detpsi = p*log(omega);

        for(g=0;g<G;g++){
			for (i=0;i<p;i++) psi[i] = omega*delta[g*p+i];
	    	log_detsig[g] = update_det_sigma_NEW2(lambda, psi, log_detpsi, p, q);
			log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];
		}

	    update_z11(v, x, z, lambda, omega, delta, mu, pi, max_v, log_c, N, G,p,q);

	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
        it++;
    }
	/*printf("ll=%f\n",l[it-1]);*/

   	a = G-1 + G*p+ p*q - q*(q-1)/2 + 1 + G*(p-1);
    
   	bic = 2.0*l[it-1] - a*log(N);
    
	lambda_store(lam_vec, lambda, p, q);
    psi_vec[0]=omega;
    for(g=0;g<G;g++){
        for(i=0;i<p;i++){
            psi_vec[k]=delta[g*p+i];
            k++;
        }
    }
    
	
        free(lambda); free(mu); free(v);  free(n); free(log_c); free(l); free(at); free(pi); free(delta);  free(log_detsig);
        free(delta0);
        for(g=0; g < G; g++) {
           free(beta[g]);
           free(theta[g]);
           free(sampcov[g]);
        }
      free(beta); free(theta); free(sampcov);
	return bic;
} 

double claecm12(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){

        double bic, log_detpsi, a, omega;
	int g,i,k=1,j,it=0,stop=0;
	       
	/* printf("G=%d \t q=%d\n",G,q);*/
        double *max_v = malloc(sizeof(double)*N);
        double *v = malloc(sizeof(double)*N*G);
        double *log_detsig = malloc(sizeof(double)*G);
        double *log_c = malloc(sizeof(double)*G);
        double *pi = malloc(sizeof(double)*G);
        double *n = malloc(sizeof(double)*G);
        double *at = malloc(sizeof(double)*150000);
        double *l = malloc(sizeof(double)*150000);
        double **sampcov    = malloc(sizeof(double*)*G);
        double **lambda    = malloc(sizeof(double*)*G);
        double **beta    = malloc(sizeof(double*)*G);
        double **theta    = malloc(sizeof(double*)*G);
        for(g=0; g < G; g++) {
           sampcov[g]      = malloc(sizeof(double)*p*p);
           lambda[g]      = malloc(sizeof(double)*p*q);
           beta[g]      = malloc(sizeof(double)*q*p);
           theta[g]      = malloc(sizeof(double)*q*q);
        }
        double *mu = malloc(sizeof(double)*G*p);
        double *delta = malloc(sizeof(double)*G*p);
        double *psi = malloc(sizeof(double)*p);
        double *delta0 = malloc(sizeof(double)*p);

	
	omega = *psi_vec;
	get_data2(lam_vec,lambda,G,p,q);
		
	for(g=0;g<G;g++)
		for(i=0;i<p;i++) 
			delta[g*p+i] = 1.0; 
	
   	while(stop==0){            
           
		update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	
	    update_mu(mu, n, x, z, G, N, p); 
                          
	    if(it>0){
			update_z12(v, x, z, lambda, omega, delta, mu, pi, max_v, log_c, N, G,p,q);
		}	
		
	    update_sg(sampcov, x, z, mu, n, p, G, N);
	    
	    for(g=0;g<G;g++){ 
			for (i=0;i<p;i++) psi[i] = omega*delta[g*p+i];
			update_beta2(beta[g], psi, lambda[g], p, q);    
		}	
	
    	for (g=0;g<G;g++) update_theta(theta[g], beta[g], lambda[g], sampcov[g], p, q);
        
    	for (g=0;g<G;g++) update_lambda(lambda[g], beta[g], sampcov[g], theta[g], p, q);
		    
        omega =0.0;
	    for(g=0;g<G;g++){
               for(j=0; j<p; j++){
                   delta0[j] = delta[g*p+j];
                }
             omega += pi[g]*update_omega2(lambda[g], delta0, beta[g], sampcov[g], p, q);
	   } 
	    for(g=0;g<G;g++){
               for(j=0; j<p; j++){
                   delta0[j] = delta[g*p+j];
                }
            update_delta3(delta0, lambda[g], omega, beta[g], sampcov[g], theta[g], n[g], p, q); 
              for(j=0; j<p; j++){
                 delta[g*p+j] = delta0[j];
              }
            }
	    
		log_detpsi = p*log(omega);        

		for(g=0;g<G;g++){
			for (i=0;i<p;i++) psi[i] = omega*delta[g*p+i];
			log_detsig[g] = update_det_sigma_NEW2(lambda[g], psi, log_detpsi, p, q);
			log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];
		}

	    update_z12(v, x, z, lambda, omega, delta, mu, pi, max_v, log_c, N, G,p,q);

	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
        it++;
    }
	/*printf("ll=%f\n",l[it-1]);*/
   
   	a = G-1 + G*p+ G*(p*q - q*(q-1)/2) + 1+ G*(p-1);
    
   	bic = 2.0*l[it-1] - a*log(N);
    
    lambda_storeG(lam_vec, lambda, G, p, q);
    psi_vec[0]=omega;
    for(g=0;g<G;g++){
        for(i=0;i<p;i++){
            psi_vec[k]=delta[g*p+i];
            k++;
        }
    }
        free(mu); free(v); free(n); free(l); free(at); free(pi); free(delta);
        free(log_c); free(log_detsig); free(delta0);


        for(g=0; g < G; g++) {
           free(beta[g]);
           free(theta[g]);
           free(lambda[g]);
           free(sampcov[g]);
        }
       free(beta); free(theta); free(lambda); free(sampcov);    


	return bic;
}
