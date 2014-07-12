#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<R.h>
#include "functions.h"

void lambda_store(double *lam_vec, double **lambda, int p, int q){
    int i,j,k=0;
    for (i=0;i<p;i++){
        for(j=0;j<q;j++){
            lam_vec[k]=lambda[i][j];
            k++;
            /*Rprintf("%f ",lambda[i][j]);*/
        }/*Rprintf("\n");*/
    }
}

void lambda_storeG(double *lam_vec, double ***lambda, int G, int p, int q){
    int i,j,g,k=0;
    for(g=0;g<G;g++){
        for (i=0;i<p;i++){
            for(j=0;j<q;j++){
                lam_vec[k]=lambda[g][i][j];
                k++;
                /*Rprintf("%f ",lambda[g][i][j]);*/
            }/*Rprintf("\n");*/
        }/*Rprintf("\n");*/
    }
}


double aecm(double **z, double **x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){
	
	double **lambda, **mu, bic, **beta, **theta, **sampcovtilde, **v, *max_v, *l, *at, *pi, *n; 
	int it=0,stop=0,paras;
	double psi, log_detpsi, log_c=0.0, log_detsig;

   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc(&sampcovtilde,p,p);		
	my_alloc_vec(&max_v,N);
   	my_alloc(&v,N,G);		
   	my_alloc(&lambda,p,q);
   	my_alloc(&beta,q,p);
   	my_alloc(&theta,q,q);
   	my_alloc(&mu,G,p);
	
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
    
   	my_free(lambda,p); my_free(mu,G); free(n); my_free(beta,q); 
	my_free(theta,q); my_free(sampcovtilde,p); free(l); free(at); free(pi);
       
	return bic;
} 

double aecm2(double **z, double **x, int *cls, int q, int p, int G, int N, double *lam_vec, double *Psi, double tol){

   	double **lambda, **mu, bic, **w, **beta, **theta, 
	       **sampcovtilde, **v, *max_v, *l, *at, *pi, *det, *n; 
	int i,it=0,stop=0;
	double a, log_detpsi,log_detsig, log_c=0.0;

	/* printf("G=%d \t q=%d\n",G,q);*/

   	/* Allocate memory - matrices and vectors */ 
   	my_alloc_vec(&det,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc(&sampcovtilde,p,p);		
   	my_alloc(&lambda,p,q);			
   	my_alloc(&beta,q,p);
   	my_alloc(&theta,q,q);
   	my_alloc(&mu,G,p);
   	my_alloc(&w,G,N);
	my_alloc_vec(&max_v,N);
   	my_alloc(&v,N,G);	
		
	get_data(lam_vec,lambda,p,q);
		
   	while(stop==0){            
	    
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
    my_free(lambda,p); my_free(mu,G); my_free(w,G); free(n); free(det); 
    my_free(beta,q); my_free(theta,q); my_free(sampcovtilde,p); free(l); 
    free(at); free(pi);   
	return bic;
}

double aecm3(double **z, double **x, int *cls, int q, int p, int G, int N, double *lam_vec, double *Psi, double tol){

	double **lambda, **mu, bic, ***beta, ***theta, ***sampcov, 
	       *l, *at, *pi, *n, *log_detpsi, *log_detsig, *log_det,
			**v, *max_v, a; 
	int g,it=0,stop=0;

	/* printf("G=%d \t q=%d\n",G,q);*/

    my_alloc_vec(&log_det,G);
    my_alloc_vec(&log_detpsi,G);
    my_alloc_vec(&log_detsig,G);
   	my_alloc_vec(&pi,G);
    my_alloc_vec(&n,G);
    my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
    my_alloc_3d(&sampcov,G,p,p);		
    my_alloc(&lambda,p,q);		
    my_alloc_3d(&beta,G,q,p);		
    my_alloc_3d(&theta,G,q,q);
    my_alloc(&mu,G,p);
	my_alloc_vec(&max_v,N);
   	my_alloc(&v,N,G);		
    
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
   	my_free(lambda,p); my_free(mu,G); my_free(v,G); free(n); free(log_det); free(max_v); 
	my_free_3d(beta,G,q); my_free_3d(theta,G,q); my_free_3d(sampcov,G,p); 
	free(l); free(at); free(pi); free(log_detpsi);
    
	return bic;
}

double aecm4(double **z, double **x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){

	double **lambda, **mu, bic, **w, ***beta, ***theta, ***sampcov, *l, *at, *pi, *n,
	       **Psi, *log_detpsi, *log_detsig, *log_c, **v, *max_v; 
	int it=0,stop=0,paras, g, j;
	
	/* printf("G=%d \t q=%d\n",G,q);*/

   	/* Allocate memory - matrices and vectors */ 
   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&sampcov,G,p,p);		
   	my_alloc(&lambda,p,q);
   	my_alloc_3d(&beta,G,q,p);
   	my_alloc_3d(&theta,G,q,q);
   	my_alloc(&mu,G,p);
   	my_alloc(&w,G,N);
	my_alloc(&Psi,G,p);
	my_alloc_vec(&log_detpsi,G);
	my_alloc_vec(&log_detsig,G);
	my_alloc_vec(&log_c,G);
	
	get_data(psi_vec,Psi,G,p);
	get_data(lam_vec,lambda,p,q);
        
   	while(stop==0){            
   
	    update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	
	    update_mu(mu, n, x, z, G, N, p); 

        if(it>0){ 
			update_z4(v, x, z, lambda, Psi, mu, pi, max_v, log_c, N, G, p, q);
			known_z(cls,z,N,G);
		} 
 
	    update_sg(sampcov,x, z, mu, n, p, G, N);
        
	    for(g=0;g<G;g++) update_beta2(beta[g], Psi[g], lambda, p, q);
        
   	    for(g=0;g<G;g++) update_theta(theta[g], beta[g], lambda, sampcov[g], p, q); 
    
	    update_lambda_cuu(lambda, beta, sampcov, theta, n, Psi, p, q, G);
	 
        update_psi_cuu(Psi, lambda, beta, sampcov, theta, p, q, G);

	    for (g=0;g<G;g++){
			log_detpsi[g]=0.0;
			for(j=0;j<p;j++) log_detpsi[g] += log(Psi[g][j]);
	    }
	
	    for(g=0;g<G;g++) log_detsig[g] = update_det_sigma_NEW2(lambda, Psi[g], log_detpsi[g], p, q);

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
    my_free(lambda,p); my_free(mu,G); my_free(w,G); free(n); my_free_3d(beta,G,q); 
	my_free_3d(theta,G,q); my_free_3d(sampcov,G,p); free(l); free(at); free(pi);
	free(log_detsig); free(log_c); free(log_detpsi);
	return bic;
} 

double aecm5(double **z, double **x, int*cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){

	double ***lambda, **mu, bic, **w, ***beta, ***theta, **v, *max_v, ***sampcov, *l, *at, *pi, 
	*n, *log_detsig, *log_c, log_detpsi, psi; 
	int it=0,stop=0,paras, g,i;

	/* printf("G=%d \t q=%d\n",G,q);*/

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
   	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&log_detsig,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&sampcov,G,p,p);		
   	my_alloc_3d(&lambda,G,p,q);
   	my_alloc_3d(&beta,G,q,p);
   	my_alloc_3d(&theta,G,q,q);
   	my_alloc(&mu,G,p);
   	my_alloc(&w,G,N);
	
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
    
   	my_free_3d(lambda,G,p); my_free(mu,G); my_free(w,G); free(n); 
	my_free_3d(beta,G,q); my_free_3d(theta,G,q); my_free_3d(sampcov,G,p); 
	free(l); free(at); free(pi);free(log_detsig); free(log_c);
	return bic;
}

double aecm6(double **z, double **x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi, double tol){

	double ***lambda, **mu, bic, ***beta, ***theta,**v,*max_v, ***sampcov, *l, *at, *pi, 
	*n, *log_detsig,*log_c, log_detpsi; 
	int it=0,stop=0,paras, g, j;

	/* printf("G=%d \t q=%d\n",G,q);*/

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
	my_alloc_vec(&log_detsig,G);
	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&sampcov,G,p,p);		
   	my_alloc_3d(&lambda,G,p,q);
   	my_alloc_3d(&beta,G,q,p);
   	my_alloc_3d(&theta,G,q,q);
   	my_alloc(&mu,G,p);
	
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
    
   	my_free_3d(lambda,G,p); my_free(mu,G); my_free(v,N); free(n);free(max_v); 
	my_free_3d(beta,G,q); my_free_3d(theta,G,q); my_free_3d(sampcov,G,p); 
	free(l); free(at); free(pi); free(log_detsig); free(log_c);
	return bic;
}

double aecm7(double **z, double **x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi, double tol){

	double ***lambda, **mu, bic, ***beta, ***theta, **v, *max_v, ***sampcov, *l, *at, *pi, 
		*n, *log_detpsi, *log_detsig, *log_c; 
	int it=0,stop=0,paras, g;

	/* printf("G=%d \t q=%d\n",G,q); */

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
	my_alloc_vec(&log_detpsi,G);
	my_alloc_vec(&log_detsig,G);
	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&sampcov,G,p,p);		
   	my_alloc_3d(&lambda,G,p,q);
   	my_alloc_3d(&beta,G,q,p);
   	my_alloc_3d(&theta,G,q,q);
   	my_alloc(&mu,G,p);
	
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
    
    my_free_3d(lambda,G,p); my_free(mu,G); my_free(v,N); free(n); 
	my_free_3d(beta,G,q); free(log_detpsi); my_free_3d(theta,G,q); 
	my_free_3d(sampcov,G,p); free(l); free(at); free(pi);
	free(log_detsig); free(log_c);
	
	return bic;
}

double aecm8(double **z, double **x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){

	double ***lambda, **mu, bic, ***beta, ***theta, **v, *max_v, ***sampcov, *l, *at, *pi, 
	*n, **Psi, *log_c,*log_detsig,*log_detpsi; 
	int it=0,stop=0,paras, g, j;

	/* printf("G=%d \t q=%d\n",G,q); */

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&sampcov,G,p,p);		
   	my_alloc_3d(&lambda,G,p,q);
   	my_alloc_3d(&beta,G,q,p);
   	my_alloc_3d(&theta,G,q,q);
   	my_alloc(&mu,G,p);
	my_alloc(&Psi,G,p);
	my_alloc_vec(&log_detpsi,G);
	my_alloc_vec(&log_detsig,G);
	my_alloc_vec(&log_c,G);
	
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
	    
	    for(g=0;g<G;g++) update_beta2(beta[g], Psi[g], lambda[g], p, q);
        
   	    for(g=0;g<G;g++) update_theta(theta[g], beta[g], lambda[g], sampcov[g], p, q); 
    
        for (g=0;g<G;g++) update_lambda(lambda[g], beta[g], sampcov[g], theta[g], p, q);
	 
        for (g=0;g<G;g++) update_psi2(Psi[g], lambda[g], beta[g], sampcov[g], p, q); 
           
	    for(g=0;g<G;g++){
		    log_detpsi[g]=0.0;
		    for(j=0;j<p;j++) log_detpsi[g] += log(Psi[g][j]);
	    }
	    
		for(g=0;g<G;g++) log_detsig[g] = update_det_sigma_NEW2(lambda[g], Psi[g], log_detpsi[g], p, q);
	    	    
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
    
   	my_free_3d(lambda,G,p); my_free(mu,G); my_free(v,N); free(n); 
	my_free_3d(beta,G,q); free(log_detpsi); my_free_3d(theta,G,q); 
	my_free_3d(sampcov,G,p); free(l); free(at); free(pi);
	free(log_detsig); free(log_c); my_free(Psi,G); free(max_v);
	return bic;
}

double aecm9(double **z, double **x, int *cls, int q, int p, int G, int N, double *lam_vec, double *omega, double tol){

	double **lambda, **mu, bic, ***beta, ***theta, ***sampcov, **v, *max_v, *l, *at, *pi, 
	*n, *delta, *log_c, *log_detpsi, *log_detsig, a, *psi; 
	int g,i,it=0,stop=0;
	       
	/* printf("G=%d \t q=%d\n",G,q);*/

   	/* Allocate memory - matrices and vectors */ 
   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
	my_alloc_vec(&log_detpsi,G);
	my_alloc_vec(&log_detsig,G);
	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&sampcov,G,p,p);		
   	my_alloc(&lambda,p,q);		
   	my_alloc_3d(&beta,G,q,p);		
   	my_alloc_3d(&theta,G,q,q);
   	my_alloc(&mu,G,p);
	my_alloc_vec(&delta,p);
	my_alloc_vec(&psi,p);
	
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
    my_free(lambda,p); my_free(mu,G); my_free(v,N); free(n); free(log_c); 
	my_free_3d(beta,G,q); my_free_3d(theta,G,q); my_free_3d(sampcov,G,p); 
	free(l); free(at); free(pi); free(log_detpsi); free(delta); free(log_detsig);
	
	return bic;
}

double aecm10(double **z, double **x, int *cls, int q, int p, int G, int N, double *lam_vec, double *omega, double tol){

	double ***lambda, **mu, bic, ***beta, ***theta, ***sampcov,**v, *max_v, *l, *at, *pi, 
	*n, *delta, *log_detsig, *log_c, *log_detpsi, a, *psi; 
	int g,i,it=0,stop=0;
	       
	/* printf("G=%d \t q=%d\n",G,q);*/

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
	my_alloc_vec(&log_detpsi,G);
	my_alloc_vec(&log_detsig,G);
	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&sampcov,G,p,p);		
   	my_alloc_3d(&lambda,G,p,q);		
   	my_alloc_3d(&beta,G,q,p);		
   	my_alloc_3d(&theta,G,q,q);
   	my_alloc(&mu,G,p);
	my_alloc_vec(&delta,p);
	my_alloc_vec(&psi,p);
	
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
    
   	my_free_3d(lambda,G,p); my_free(mu,G); my_free(v,N); free(n); free(log_c); 
	my_free_3d(beta,G,q); my_free_3d(theta,G,q); my_free_3d(sampcov,G,p); 
	free(l); free(at); free(pi); free(log_detpsi); free(delta); free(log_detsig);

	return bic;
}

double aecm11(double **z, double **x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){

	double **lambda, **mu, bic, ***beta, ***theta, ***sampcov, **v, *max_v, *l, *at, *pi, 
	*n, **delta, *log_c, *log_detsig, log_detpsi, a, omega, *psi; 
	int g,i,it=0,stop=0;

	/* printf("G=%d \t q=%d\n",G,q);*/

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&sampcov,G,p,p);		
   	my_alloc(&lambda,p,q);		
   	my_alloc_3d(&beta,G,q,p);		
	my_alloc_3d(&theta,G,q,q);
   	my_alloc(&mu,G,p);
	my_alloc(&delta,G,p);
	my_alloc_vec(&log_detsig,G);
	my_alloc_vec(&log_c,G);
	my_alloc_vec(&psi,p);
	
	omega = *psi_vec;
	get_data(lam_vec,lambda,p,q);
	
	for(g=0;g<G;g++)
		for(i=0;i<p;i++) 
			delta[g][i] = 1.0; 
	
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
			for (i=0;i<p;i++) psi[i] = omega*delta[g][i];
			update_beta2(beta[g], psi, lambda, p, q);    
		}	
	
    	for (g=0;g<G;g++) update_theta(theta[g], beta[g], lambda, sampcov[g], p, q);
        
	    update_lambda_cuu(lambda, beta, sampcov, theta, n, delta, p, q, G);
	    
        omega =0.0;
	    for(g=0;g<G;g++) omega += pi[g]*update_omega(lambda, delta[g], beta[g], sampcov[g], theta[g], p, q);
 	   	
	    for(g=0;g<G;g++) update_delta3(delta[g], lambda, omega, beta[g], sampcov[g], theta[g], n[g], p, q);
		
		log_detpsi = p*log(omega);

        for(g=0;g<G;g++){
			for (i=0;i<p;i++) psi[i] = omega*delta[g][i];
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
    
   	my_free(lambda,p); my_free(mu,G); my_free(v,N); free(n); free(log_c); my_free_3d(beta,G,q); my_free_3d(theta,G,q); 
	my_free_3d(sampcov,G,p); free(l); free(at); free(pi); my_free(delta,G); free(log_detsig); 
	
	return bic;
} 

double aecm12(double **z, double **x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){

	double ***lambda, **mu, bic, ***beta, ***theta, ***sampcov, **v, *max_v, *l, *at, *pi, 
	*n, **delta, *log_detsig, *log_c, log_detpsi, a, omega, *psi; 
	int g,i,it=0,stop=0;
	       
	/* printf("G=%d \t q=%d\n",G,q);*/

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
	my_alloc_vec(&log_detsig,G);
	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&sampcov,G,p,p);		
   	my_alloc_3d(&lambda,G,p,q);		
   	my_alloc_3d(&beta,G,q,p);		
	my_alloc_3d(&theta,G,q,q);
   	my_alloc(&mu,G,p);
	my_alloc(&delta,G,p);
	my_alloc_vec(&psi,p);
	
	omega = *psi_vec;
	get_data2(lam_vec,lambda,G,p,q);
	
	for(g=0;g<G;g++)
		for(i=0;i<p;i++) 
			delta[g][i] = 1.0; 
	
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
			for (i=0;i<p;i++) psi[i] = omega*delta[g][i];
			update_beta2(beta[g], psi, lambda[g], p, q);    
		}	
	
    	for (g=0;g<G;g++) update_theta(theta[g], beta[g], lambda[g], sampcov[g], p, q);
        
    	for (g=0;g<G;g++) update_lambda(lambda[g], beta[g], sampcov[g], theta[g], p, q);
		    
        omega =0.0;
	    for(g=0;g<G;g++) omega += pi[g]*update_omega2(lambda[g], delta[g], beta[g], sampcov[g], p, q);
	    
	    for(g=0;g<G;g++) update_delta3(delta[g], lambda[g], omega, beta[g], sampcov[g], theta[g], n[g], p, q); 
	    
		log_detpsi = p*log(omega);        

		for(g=0;g<G;g++){
			for (i=0;i<p;i++) psi[i] = omega*delta[g][i];
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
    
   	my_free_3d(lambda,G,p); my_free(mu,G); my_free(v,N); free(n);  
	my_free_3d(beta,G,q); my_free_3d(theta,G,q); my_free_3d(sampcov,G,p); 
	free(l); free(at); free(pi); my_free(delta,G); free(log_c);
	free(log_detsig);

	return bic;
}

double claecm(double **z, double **x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){
	
	double **lambda, **mu, bic, **beta, **theta, **sampcovtilde, **v, *max_v, *l, *at, *pi, *n; 
	int it=0,stop=0,paras;
	double psi, log_detpsi, log_c=0.0, log_detsig;
	
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc(&sampcovtilde,p,p);		
	my_alloc_vec(&max_v,N);
   	my_alloc(&v,N,G);		
   	my_alloc(&lambda,p,q);
   	my_alloc(&beta,q,p);
   	my_alloc(&theta,q,q);
   	my_alloc(&mu,G,p);
	
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
    
   	my_free(lambda,p); my_free(mu,G); free(n); my_free(beta,q); 
	my_free(theta,q); my_free(sampcovtilde,p); free(l); free(at); free(pi);
       
	return bic;
} 

double claecm2(double **z, double **x, int q, int p, int G, int N, double *lam_vec, double *Psi, double tol){

   	double **lambda, **mu, bic, **w, **beta, **theta, 
	       **sampcovtilde, **v, *max_v, *l, *at, *pi, *det, *n; 
	int i,it=0,stop=0;
	double a, log_detpsi,log_detsig, log_c=0.0;

	/* printf("G=%d \t q=%d\n",G,q);*/

   	/* Allocate memory - matrices and vectors */ 
   	my_alloc_vec(&det,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc(&sampcovtilde,p,p);		
   	my_alloc(&lambda,p,q);			
   	my_alloc(&beta,q,p);
   	my_alloc(&theta,q,q);
   	my_alloc(&mu,G,p);
   	my_alloc(&w,G,N);
	my_alloc_vec(&max_v,N);
   	my_alloc(&v,N,G);		

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
    my_free(lambda,p); my_free(mu,G); my_free(w,G); free(n); free(det); 
    my_free(beta,q); my_free(theta,q); my_free(sampcovtilde,p); free(l); 
    free(at); free(pi);   
	return bic;
}

double claecm3(double **z, double **x, int q, int p, int G, int N, double *lam_vec, double *Psi, double tol){

	double **lambda, **mu, bic, ***beta, ***theta, ***sampcov, 
	       *l, *at, *pi, *n, *log_detpsi, *log_detsig, *log_det,
			**v, *max_v, a; 
	int g,it=0,stop=0;

	/* printf("G=%d \t q=%d\n",G,q);*/

    my_alloc_vec(&log_det,G);
    my_alloc_vec(&log_detpsi,G);
    my_alloc_vec(&log_detsig,G);
   	my_alloc_vec(&pi,G);
    my_alloc_vec(&n,G);
    my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
    my_alloc_3d(&sampcov,G,p,p);		
    my_alloc(&lambda,p,q);		
    my_alloc_3d(&beta,G,q,p);		
    my_alloc_3d(&theta,G,q,q);
    my_alloc(&mu,G,p);
	my_alloc_vec(&max_v,N);
   	my_alloc(&v,N,G);		
        
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
   	my_free(lambda,p); my_free(mu,G); my_free(v,G); free(n); free(log_det); free(max_v); 
	my_free_3d(beta,G,q); my_free_3d(theta,G,q); my_free_3d(sampcov,G,p); 
	free(l); free(at); free(pi); free(log_detpsi);
    
	return bic;
}

double claecm4(double **z, double **x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){

	double **lambda, **mu, bic, **w, ***beta, ***theta, ***sampcov, *l, *at, *pi, *n,
	       **Psi, *log_detpsi, *log_detsig, *log_c, **v, *max_v; 
	int it=0,stop=0,paras, g, j;
	
   	/* Allocate memory - matrices and vectors */ 
   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&sampcov,G,p,p);		
   	my_alloc(&lambda,p,q);
   	my_alloc_3d(&beta,G,q,p);
   	my_alloc_3d(&theta,G,q,q);
   	my_alloc(&mu,G,p);
   	my_alloc(&w,G,N);
	my_alloc(&Psi,G,p);
	my_alloc_vec(&log_detpsi,G);
	my_alloc_vec(&log_detsig,G);
	my_alloc_vec(&log_c,G);
	
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
        
	    for(g=0;g<G;g++) update_beta2(beta[g], Psi[g], lambda, p, q);
        
   	    for(g=0;g<G;g++) update_theta(theta[g], beta[g], lambda, sampcov[g], p, q); 
    
	    update_lambda_cuu(lambda, beta, sampcov, theta, n, Psi, p, q, G);
	 
        update_psi_cuu(Psi, lambda, beta, sampcov, theta, p, q, G);

	    for (g=0;g<G;g++){
			log_detpsi[g]=0.0;
			for(j=0;j<p;j++) log_detpsi[g] += log(Psi[g][j]);
	    }
	
	    for(g=0;g<G;g++) log_detsig[g] = update_det_sigma_NEW2(lambda, Psi[g], log_detpsi[g], p, q);

		for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];

        update_z4(v, x, z, lambda, Psi, mu, pi, max_v, log_c, N, G, p, q);
            
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++;
	}
	
	/* no of parameters {pq - q(q - 1)/2} + Gp */
	paras = G-1 + G*p + p*q - q*(q-1)/2 + G*p;
    
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
    my_free(lambda,p); my_free(mu,G); my_free(w,G); free(n); my_free_3d(beta,G,q); 
	my_free_3d(theta,G,q); my_free_3d(sampcov,G,p); free(l); free(at); free(pi);
	free(log_detsig); free(log_c); free(log_detpsi);
	return bic;
} 

double claecm5(double **z, double **x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){

	double ***lambda, **mu, bic, **w, ***beta, ***theta, **v, *max_v, ***sampcov, *l, *at, *pi, 
	*n, *log_detsig, *log_c, log_detpsi, psi; 
	int it=0,stop=0,paras, g,i;

	/*Rprintf("G=%d \t q=%d\n",G,q);*/

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
   	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&log_detsig,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&sampcov,G,p,p);		
   	my_alloc_3d(&lambda,G,p,q);
   	my_alloc_3d(&beta,G,q,p);
   	my_alloc_3d(&theta,G,q,q);
   	my_alloc(&mu,G,p);
   	my_alloc(&w,G,N);
	
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
    
   	my_free_3d(lambda,G,p); my_free(mu,G); my_free(w,G); free(n); 
	my_free_3d(beta,G,q); my_free_3d(theta,G,q); my_free_3d(sampcov,G,p); 
	free(l); free(at); free(pi);free(log_detsig); free(log_c);
	return bic;
}

double claecm6(double **z, double **x, int q, int p, int G, int N, double *lam_vec, double *psi, double tol){

	double ***lambda, **mu, bic, ***beta, ***theta,**v,*max_v, ***sampcov, *l, *at, *pi, 
	*n, *log_detsig,*log_c, log_detpsi; 
	int it=0,stop=0,paras, g, j;

	/* printf("G=%d \t q=%d\n",G,q);*/

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
	my_alloc_vec(&log_detsig,G);
	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&sampcov,G,p,p);		
   	my_alloc_3d(&lambda,G,p,q);
   	my_alloc_3d(&beta,G,q,p);
   	my_alloc_3d(&theta,G,q,q);
   	my_alloc(&mu,G,p);
	
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
    
   	my_free_3d(lambda,G,p); my_free(mu,G); my_free(v,N); free(n);free(max_v); 
	my_free_3d(beta,G,q); my_free_3d(theta,G,q); my_free_3d(sampcov,G,p); 
	free(l); free(at); free(pi); free(log_detsig); free(log_c);
	return bic;
}

double claecm7(double **z, double **x, int q, int p, int G, int N, double *lam_vec, double *psi, double tol){

	double ***lambda, **mu, bic, ***beta, ***theta, **v, *max_v, ***sampcov, *l, *at, *pi, 
		*n, *log_detpsi, *log_detsig, *log_c; 
	int it=0,stop=0,paras, g;

	/* printf("G=%d \t q=%d\n",G,q); */

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
	my_alloc_vec(&log_detpsi,G);
	my_alloc_vec(&log_detsig,G);
	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&sampcov,G,p,p);		
   	my_alloc_3d(&lambda,G,p,q);
   	my_alloc_3d(&beta,G,q,p);
   	my_alloc_3d(&theta,G,q,q);
   	my_alloc(&mu,G,p);
	
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
    /*for(g=0;g<G;g++){
        Rprintf("%f ",psi[g]);
    }*/
    
    my_free_3d(lambda,G,p); my_free(mu,G); my_free(v,N); free(n); 
	my_free_3d(beta,G,q); free(log_detpsi); my_free_3d(theta,G,q); 
	my_free_3d(sampcov,G,p); free(l); free(at); free(pi);
	free(log_detsig); free(log_c);
	
	return bic;
}

double claecm8(double **z, double **x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){

	double ***lambda, **mu, bic, ***beta, ***theta, **v, *max_v, ***sampcov, *l, *at, *pi, 
	*n, **Psi, *log_c,*log_detsig,*log_detpsi; 
	int it=0,stop=0,paras, g, j;

	/* printf("G=%d \t q=%d\n",G,q); */

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&sampcov,G,p,p);		
   	my_alloc_3d(&lambda,G,p,q);
   	my_alloc_3d(&beta,G,q,p);
   	my_alloc_3d(&theta,G,q,q);
   	my_alloc(&mu,G,p);
	my_alloc(&Psi,G,p);
	my_alloc_vec(&log_detpsi,G);
	my_alloc_vec(&log_detsig,G);
	my_alloc_vec(&log_c,G);
	
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
	    
	    for(g=0;g<G;g++) update_beta2(beta[g], Psi[g], lambda[g], p, q);
        
   	    for(g=0;g<G;g++) update_theta(theta[g], beta[g], lambda[g], sampcov[g], p, q); 
    
        for (g=0;g<G;g++) update_lambda(lambda[g], beta[g], sampcov[g], theta[g], p, q);
	 
        for (g=0;g<G;g++) update_psi2(Psi[g], lambda[g], beta[g], sampcov[g], p, q); 
           
	    for(g=0;g<G;g++){
		    log_detpsi[g]=0.0;
		    for(j=0;j<p;j++) log_detpsi[g] += log(Psi[g][j]);
	    }
	    
		for(g=0;g<G;g++) log_detsig[g] = update_det_sigma_NEW2(lambda[g], Psi[g], log_detpsi[g], p, q);
	    	    
		for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];

		update_z8(v, x, z, lambda, Psi, mu, pi, max_v, log_c, N, G, p, q);
            
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++; 
	}
	/*printf("ll=%f\n",l[it-1]);*/
	paras = G-1 + G*p + G*(p*q - q*(q-1)/2) + G*p;
   	
	bic = 2.0*l[it-1] - paras*log(N);
    
    lambda_storeG(lam_vec, lambda, G, p, q);
    
   	my_free_3d(lambda,G,p); my_free(mu,G); my_free(v,N); free(n); 
	my_free_3d(beta,G,q); free(log_detpsi); my_free_3d(theta,G,q); 
	my_free_3d(sampcov,G,p); free(l); free(at); free(pi);
	free(log_detsig); free(log_c); my_free(Psi,G); free(max_v);
	return bic;
}

double claecm9(double **z, double **x, int q, int p, int G, int N, double *lam_vec, double *omega, double tol){

	double **lambda, **mu, bic, ***beta, ***theta, ***sampcov, **v, *max_v, *l, *at, *pi, 
	*n, *delta, *log_c, *log_detpsi, *log_detsig, a, *psi; 
	int g,i,it=0,stop=0;
	       
	/* printf("G=%d \t q=%d\n",G,q);*/

   	/* Allocate memory - matrices and vectors */ 
   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
	my_alloc_vec(&log_detpsi,G);
	my_alloc_vec(&log_detsig,G);
	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&sampcov,G,p,p);		
   	my_alloc(&lambda,p,q);		
   	my_alloc_3d(&beta,G,q,p);		
   	my_alloc_3d(&theta,G,q,q);
   	my_alloc(&mu,G,p);
	my_alloc_vec(&delta,p);
	my_alloc_vec(&psi,p);
	
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
    for(i=0;i<p;i++){omega[G+i]=delta[i];}
    
    /* Deallocate memory: psi not freed? */
    my_free(lambda,p); my_free(mu,G); my_free(v,N); free(n); free(log_c); 
	my_free_3d(beta,G,q); my_free_3d(theta,G,q); my_free_3d(sampcov,G,p); 
	free(l); free(at); free(pi); free(log_detpsi); free(delta); free(log_detsig);
	
	return bic;
}

double claecm10(double **z, double **x, int q, int p, int G, int N, double *lam_vec, double *omega, double tol){

	double ***lambda, **mu, bic, ***beta, ***theta, ***sampcov,**v, *max_v, *l, *at, *pi, 
	*n, *delta, *log_detsig, *log_c, *log_detpsi, a, *psi; 
	int g,i,it=0,stop=0;
	       
	/* printf("G=%d \t q=%d\n",G,q);*/

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
	my_alloc_vec(&log_detpsi,G);
	my_alloc_vec(&log_detsig,G);
	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&sampcov,G,p,p);		
   	my_alloc_3d(&lambda,G,p,q);		
   	my_alloc_3d(&beta,G,q,p);		
   	my_alloc_3d(&theta,G,q,q);
   	my_alloc(&mu,G,p);
	my_alloc_vec(&delta,p);
	my_alloc_vec(&psi,p);
	
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
	/*printf("ll=%f\n",l[it-1]);*/
   
   	a = G-1 + G*p+ G*(p*q - q*(q-1)/2) + G + (p-1);
    
   	bic = 2.0*l[it-1] - a*log(N);
    
    lambda_storeG(lam_vec, lambda, G, p, q);
    for(i=0;i<p;i++){omega[G+i]=delta[i];}
    
   	my_free_3d(lambda,G,p); my_free(mu,G); my_free(v,N); free(n); free(log_c); 
	my_free_3d(beta,G,q); my_free_3d(theta,G,q); my_free_3d(sampcov,G,p); 
	free(l); free(at); free(pi); free(log_detpsi); free(delta); free(log_detsig);

	return bic;
}

double claecm11(double **z, double **x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){

	double **lambda, **mu, bic, ***beta, ***theta, ***sampcov, **v, *max_v, *l, *at, *pi, 
	*n, **delta, *log_c, *log_detsig, log_detpsi, a, omega, *psi; 
	int g,i,k=1,it=0,stop=0;

	/* printf("G=%d \t q=%d\n",G,q);*/

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&sampcov,G,p,p);		
   	my_alloc(&lambda,p,q);		
   	my_alloc_3d(&beta,G,q,p);		
	my_alloc_3d(&theta,G,q,q);
   	my_alloc(&mu,G,p);
	my_alloc(&delta,G,p);
	my_alloc_vec(&log_detsig,G);
	my_alloc_vec(&log_c,G);
	my_alloc_vec(&psi,p);
	
	omega = *psi_vec;
	get_data(lam_vec,lambda,p,q);
	
	for(g=0;g<G;g++)
		for(i=0;i<p;i++) 
			delta[g][i] = 1.0; 
	
   	while(stop==0){            
           
        update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	
	    update_mu(mu, n, x, z, G, N, p); 
                          
        if(it>0){ 
			update_z11(v, x, z, lambda, omega, delta, mu, pi, max_v, log_c, N, G,p,q);
		}	
		
	    update_sg(sampcov, x, z, mu, n, p, G, N);

	    for(g=0;g<G;g++){ 
			for (i=0;i<p;i++) psi[i] = omega*delta[g][i];
			update_beta2(beta[g], psi, lambda, p, q);    
		}	
	
    	for (g=0;g<G;g++) update_theta(theta[g], beta[g], lambda, sampcov[g], p, q);
        
	    update_lambda_cuu(lambda, beta, sampcov, theta, n, delta, p, q, G);
	    
        omega =0.0;
	    for(g=0;g<G;g++) omega += pi[g]*update_omega(lambda, delta[g], beta[g], sampcov[g], theta[g], p, q);
 	   	
	    for(g=0;g<G;g++) update_delta3(delta[g], lambda, omega, beta[g], sampcov[g], theta[g], n[g], p, q); 
		
		log_detpsi = p*log(omega);

        for(g=0;g<G;g++){
			for (i=0;i<p;i++) psi[i] = omega*delta[g][i];
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
            psi_vec[k]=delta[g][i];
            k++;
        }
    }
    
   	my_free(lambda,p); my_free(mu,G); my_free(v,N); free(n); free(log_c); my_free_3d(beta,G,q); my_free_3d(theta,G,q); 
	my_free_3d(sampcov,G,p); free(l); free(at); free(pi); my_free(delta,G); free(log_detsig); 
	
	return bic;
} 

double claecm12(double **z, double **x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol){

	double ***lambda, **mu, bic, ***beta, ***theta, ***sampcov, **v, *max_v, *l, *at, *pi, 
	*n, **delta, *log_detsig, *log_c, log_detpsi, a, omega, *psi; 
	int g,i,k=1,it=0,stop=0;
	       
	/* printf("G=%d \t q=%d\n",G,q);*/

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
	my_alloc_vec(&log_detsig,G);
	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&sampcov,G,p,p);		
   	my_alloc_3d(&lambda,G,p,q);		
   	my_alloc_3d(&beta,G,q,p);		
	my_alloc_3d(&theta,G,q,q);
   	my_alloc(&mu,G,p);
	my_alloc(&delta,G,p);
	my_alloc_vec(&psi,p);
	
	omega = *psi_vec;
	get_data2(lam_vec,lambda,G,p,q);
		
	for(g=0;g<G;g++)
		for(i=0;i<p;i++) 
			delta[g][i] = 1.0; 
	
   	while(stop==0){            
           
		update_n(n, z, G, N);
	
	    update_pi(pi, n, G, N);
	
	    update_mu(mu, n, x, z, G, N, p); 
                          
	    if(it>0){
			update_z12(v, x, z, lambda, omega, delta, mu, pi, max_v, log_c, N, G,p,q);
		}	
		
	    update_sg(sampcov, x, z, mu, n, p, G, N);
	    
	    for(g=0;g<G;g++){ 
			for (i=0;i<p;i++) psi[i] = omega*delta[g][i];
			update_beta2(beta[g], psi, lambda[g], p, q);    
		}	
	
    	for (g=0;g<G;g++) update_theta(theta[g], beta[g], lambda[g], sampcov[g], p, q);
        
    	for (g=0;g<G;g++) update_lambda(lambda[g], beta[g], sampcov[g], theta[g], p, q);
		    
        omega =0.0;
	    for(g=0;g<G;g++) omega += pi[g]*update_omega2(lambda[g], delta[g], beta[g], sampcov[g], p, q);
	    
	    for(g=0;g<G;g++) update_delta3(delta[g], lambda[g], omega, beta[g], sampcov[g], theta[g], n[g], p, q); 
	    
		log_detpsi = p*log(omega);        

		for(g=0;g<G;g++){
			for (i=0;i<p;i++) psi[i] = omega*delta[g][i];
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
            psi_vec[k]=delta[g][i];
            k++;
        }
    }
    
   	my_free_3d(lambda,G,p); my_free(mu,G); my_free(v,N); free(n);  
	my_free_3d(beta,G,q); my_free_3d(theta,G,q); my_free_3d(sampcov,G,p); 
	free(l); free(at); free(pi); my_free(delta,G); free(log_c);
	free(log_detsig);

	return bic;
}
