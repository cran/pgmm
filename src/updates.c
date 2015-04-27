#include <R.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "functions.h"

void update_n(double *n, double *z, int G, int N){    
	int g,i;	    
	for(g=0; g<G; g++){
            n[g]=0.0;
            for(i=0; i<N; i++){
                n[g] += z[g+i*G];
            }
    }
}
    
void update_pi(double *pi, double *n, int G, int N){    
	int g;	    
	for(g=0; g<G; g++) pi[g] = n[g]/N; 
}
    
void update_mu(double *mu, double *n, double *x, double *z,int G, 
	    int N, int p){    
	int i,j,g;
	for(g=0; g<G; g++){
            for(j=0; j<p; j++){
                mu[g*p+j]=0.0;
                for(i=0; i<N; i++) {
                    mu[g*p+j] += z[g+ i*G]*x[j + i*p]; 
				}
                mu[g*p+j] /= n[g];
            }
    }

}


void update_stilde(double *sampcovtilde, double *x, double *z, double *mu,
	int G, int N, int p){
        int i,j,k,g;
        
	for(j=0; j<p; j++){
            for(k=0; k<p; k++){
                sampcovtilde[j*p+k]=0.0;
                for(g=0;g<G;g++)    
                    for(i=0; i<N; i++){
                        
sampcovtilde[j*p+k]+=z[g+i*G]*(x[j+i*p] -mu[j+g*p])*(x[k+i*p]-mu[k+g*p]);
                     }
				sampcovtilde[j*p+k]/=N;
            }
    }
}

void update_sg(double **sg, double *x, double *z, double *mu,
	double *n, int p, int G, int N){
        int i,j,k,g;
	
	for(g=0;g<G;g++){    
	   for(j=0; j<p; j++){
               for(k=0; k<p; k++){
                    sg[g][j*p+k]=0.0;
                    for(i=0; i<N; i++){
                        sg[g][j*p+k]+=z[g+i*G]*(x[j+i*p]-mu[g*p+j])*(x[k+i*p]-mu[g*p+k])/n[g];
                    }
               }
           }
	}

}

void update_beta1(double *beta, double psi, double *lambda, int p, int q){
	int i,j;
	double det[1];

        double *lhs = malloc(sizeof(double)*q*p);
        double *rhs = malloc(sizeof(double)*p*p);
        double *cp = malloc(sizeof(double)*q*q);
        double *result = malloc(sizeof(double)*p*p);


 	/* Lambda'/psi */
    mx_trans(p, q, lambda, lhs);
    for(i=0;i<q;i++)
		for(j=0;j<p;j++)
            lhs[i*p+j]/=psi; /*psi[j]*/
 	/* Lambda'/psi* Lambda */
    mx_mult(q, p, q, lhs, lambda, cp);

    /* (I + Lambda'/psi*Lambda)^{-1} */
    for(i=0; i<q;i++){
		for(j=0;j<q;j++){
			result[i*q+j] = cp[i*q+j];
			if (i==j) result[i*q+i] += 1.0;
		}
	}
    GaussJordan(q, result, rhs, det);
    /* work out rhs*/
	mx_mult(q,q,q,cp,rhs,result);
    mx_mult(q,q,p,result,lhs,rhs);
	
	for(i=0;i<q;i++){
		for(j=0;j<p;j++){
			beta[i*p+j] = lhs[i*p+j] - rhs[i*p+j]; 
                }
        }      
        free(lhs); free(result); free(rhs); free(cp); 

}

void update_beta2(double *beta, double *Psi, double *lambda, int p, int q){
	int i,j;
	double det[1];
        double *lhs = malloc(sizeof(double)*p*p);
        double *rhs = malloc(sizeof(double)*p*p);
        double *cp = malloc(sizeof(double)*p*p);
        double *result = malloc(sizeof(double)*p*p);

 	/* Lambda'/psi */
    mx_trans(p, q, lambda, lhs);
    for(i=0;i<q;i++)
		for(j=0;j<p;j++)
            lhs[i*p+j]/=Psi[j];
 	/* Lambda'/psi* Lambda */
    mx_mult(q, p, q, lhs, lambda, cp);
        
    /* (I + Lambda'/psi*Lambda)^{-1} */
    for(i=0; i<q;i++){
		for(j=0;j<q;j++){
			result[i*q+j] = cp[i*q+j];
			if (i==j) result[i*q+i] += 1.0;
		}
	}
    GaussJordan(q, result, rhs, det);
    /* work out rhs*/
	mx_mult(q,q,q,cp,rhs,result);
    mx_mult(q,q,p,result,lhs,rhs);
	
	for(i=0;i<q;i++){
		for(j=0;j<p;j++){
			beta[i*p+j] = lhs[i*p+j] - rhs[i*p+j];
                }
        }       
        free(lhs); free(result); free(rhs); free(cp); 

}

void update_theta(double *theta, double *beta, double *lambda, double *sampcovtilde, int p, int q){ 
	int i,j;
	
        double *tmp = malloc(sizeof(double)*p*p);
        double *r_1 = malloc(sizeof(double)*q*q);
        double *r_2 = malloc(sizeof(double)*q*p);
        double *r_3 = malloc(sizeof(double)*q*q);
        double *id_q = malloc(sizeof(double)*q*q);
	generate_identity(q, id_q);	
	
	/* Work out beta*lambda */
    mx_mult(q, p, q, beta, lambda, r_1);
    
    /* Now, work out beta*sampcovtilde*beta' */
    mx_mult(q, p, p, beta, sampcovtilde, r_2);
    mx_trans(q, p, beta, tmp);    
    mx_mult(q, p, q, r_2, tmp, r_3);    
                 
    for(i=0;i<q;i++){
		for(j=0;j<q;j++){
			theta[i*q+j] = id_q[i*q+j]-r_1[i*q+j]+r_3[i*q+j];
                }
     }
        free(id_q); free(tmp); free(r_1); free(r_2); free(r_3);
       
}

void update_lambda(double *lambda, double *beta, double *s, double *theta, int p, int q){
	
	int i,j;
        double  det[1];
	
        double *tran = malloc(sizeof(double)*p*q);
        double *res1 = malloc(sizeof(double)*p*q);
        double *res2 = malloc(sizeof(double)*q*q);
        double *res3 = malloc(sizeof(double)*q*q);
	
	/* sampcovtilde*beta'*/
    mx_trans(q,p,beta,tran);
    mx_mult(p,p,q,s,tran,res1);
    	
	/* Make of copy of theta ahead of Gauss-Jordan */
	for(i=0;i<q;i++)
        for(j=0;j<q;j++)
            res3[i*q+j] = theta[i*q+j];   
        
	/* lambda=sampcovtilde*beta'*theta^{-1} */
    GaussJordan(q, res3, res2,det);
    mx_mult(p,q,q,res1,res2,lambda);
    free(tran); free(res1); free(res2); free(res3);
    
}

/* For CUU */
void update_lambda_cuu(double *lambda, double **beta, double **s, 
	double **theta, double *n, double *Psi, int p, int q, int G){
	
	int i,j,k,g;
        double det[1];
        double *tran = malloc(sizeof(double)*p*q);
        double *res1 = malloc(sizeof(double)*p*q);
        double *res2 = malloc(sizeof(double)*p*q);
        double *res3 = malloc(sizeof(double)*q*q);
        double *result = malloc(sizeof(double)*q);
        double *lambda0 = malloc(sizeof(double)*q);
	
	
    /* Compute the RHS --- only needs to happen once */
	/* sum_g[sampcov_g*beta_g'*n_g/psi_g] */
    for(g=0;g<G;g++){
		mx_trans(q,p,beta[g],tran);
		mx_mult(p,p,q,s[g],tran,res1);
		if(g==0){ 
			for(i=0;i<p;i++){
				for(j=0;j<q;j++){
					res2[i*q+j] = res1[i*q+j]*n[g]/Psi[g*p+i];
                                }
                       }
		}else{
			for(i=0;i<p;i++)
				for(j=0;j<q;j++)
					res2[i*q+j] += res1[i*q+j]*n[g]/Psi[g*p+i];    
		}
    }
                         

    /* Now solve for lambda row-by-row */
	for (i=0;i<p;i++){
		/* First, compute the theta sum*/
		for(g=0;g<G;g++){
			if(g==0){ 
				for(k=0;k<q;k++)
					for(j=0;j<q;j++)
						res3[k*q+j] = theta[g][k*q+j]*n[g]/Psi[i+g*p];
			}else{
				for(k=0;k<q;k++)
					for(j=0;j<q;j++)
						res3[k*q+j] += theta[g][k*q+j]*n[g]/Psi[i+g*p];    
			}
		}
		/* Invert theta sum */
		GaussJordan(q, res3, res1, det);
		
		/* Now solve for row i of lambda */
           for(j=0;j<q;j++)
               {result[j] = res2[i*q+j];}
		vec_mx_mult(q, q, result, res1, lambda0);
           for(j=0;j<q;j++)
               {lambda[i*q+j] = lambda0[j];}
	}     	

        free(tran); free(res1); free(res2); free(res3); free(result); free(lambda0);
}

/* For CUC */
void update_lambda2(double *lambda, double **beta, double **s, double **theta, double *n, double *Psi, int p, 
		int q, int G){
	
	int i,j,g;
        double det[1];
        double *tran = malloc(sizeof(double)*p*q);
        double *res1 = malloc(sizeof(double)*p*q);
        double *res2 = malloc(sizeof(double)*p*q);
        double *res3 = malloc(sizeof(double)*q*q);
	
	
    /* sum_g[sampcov_g*beta_g'*n_g/psi_g] */
    for(g=0;g<G;g++){
		mx_trans(q,p,beta[g],tran);
		mx_mult(p,p,q,s[g],tran,res1);
		if(g==0){ 
			for(i=0;i<p;i++)
				for(j=0;j<q;j++)
					res2[i*q+j] = res1[i*q+j]*n[g]/Psi[g];
		}else{
			for(i=0;i<p;i++)
				for(j=0;j<q;j++)
					res2[i*q+j] += res1[i*q+j]*n[g]/Psi[g];    
		}
    }

   	/* sum_g[theta_g'*n_g/psi_g]^{-1} */
    for(g=0;g<G;g++){
		if(g==0){ 
			for(i=0;i<q;i++)
				for(j=0;j<q;j++)
					res3[i*q+j] = theta[g][i*q+j]*n[g]/Psi[g];
		}else{
			for(i=0;i<q;i++)
				for(j=0;j<q;j++)
					res3[i*q+j] += theta[g][i*q+j]*n[g]/Psi[g];    
		}
    }
        
    GaussJordan(q, res3, res1,det); /* inverting the theta sum */
    mx_mult(p,q,q,res2,res1,lambda); /* this gives lambda */ 
	
    free(tran); free(res1); free(res2); free(res3);
}

double update_psi(double *lambda, double *beta, double *sampcovtilde, int p, int q){ 
	
	int i;
	double psi;
        double *result_1 = malloc(sizeof(double)*p*p);
        double *result_2 = malloc(sizeof(double)*p);


	/* lambda*beta*sampcovtilde */
    mx_mult(p,q,p,lambda, beta, result_1);
    mx_mult_diag1(p,p,result_1, sampcovtilde, result_2);
       
    psi = 0.0;
    for(i=0;i<p;i++)
        psi += sampcovtilde[i*p+i]-result_2[i];
    psi /= p;
	
	free(result_1); free(result_2); 
	
	return psi;
}

/* Account for Psi_p*/
void update_psi2(double *psi, double *lambda, double *beta, double *sampcovtilde, int p, int q){ 
	
	int i;
        double *result_1 = malloc(sizeof(double)*p*p);
        double *result_2 = malloc(sizeof(double)*p);

	/* lambda*beta*sampcovtilde */
	mx_mult(p,q,p,lambda, beta, result_1);
    mx_mult_diag1(p,p,result_1, sampcovtilde, result_2);
        
    for(i=0;i<p;i++){
		psi[i] = sampcovtilde[i*p+i]-result_2[i];
    }
	free(result_1); free(result_2); 
}

/* CUC case */
double update_psi3(double *lambda, double *beta, double *sampcovg, double *theta, int p, int q){ 
	
	int i;
	double psi;
        double *temp = malloc(sizeof(double)*q*p);
        double *result_1 = malloc(sizeof(double)*p*p);
        double *result_2 = malloc(sizeof(double)*p);
        double *result_3 = malloc(sizeof(double)*p);


	/* lambda*beta*sampcov */
    mx_mult(p,q,p,lambda, beta, result_1);
    mx_mult_diag1(p,p,result_1, sampcovg, result_2);
       
	/* lambda*theta*lambda' */
    mx_trans(p,q,lambda,temp);
    mx_mult(p,q,q,lambda,theta,result_1); 
    mx_mult_diag1(p,q,result_1,temp,result_3);
        
    psi = 0.0;
    for(i=0;i<p;i++)
        psi += sampcovg[i*p+i]-2*result_2[i]+result_3[i];
    psi /= p;
	
	free(temp); free(result_1); free(result_2); 
	free(result_3);
	
	return psi;
}

/* CUU case */
void update_psi_cuu(double *psi, double *lambda, double **beta, double **sampcovg, double 
	**theta, int p, int q, int G){ 
	
	int i, g, j;
        double *temp = malloc(sizeof(double)*q*p);
        double *result_1 = malloc(sizeof(double)*p*p);
        double *result_2 = malloc(sizeof(double)*G*p);
        double *result_3 = malloc(sizeof(double)*G*p);
        double *result = malloc(sizeof(double)*p);
        

	/* lambda*beta*sampcov */
    for (g=0;g<G;g++){
		mx_mult(p,q,p,lambda, beta[g], result_1);
    	mx_mult_diag1(p,p,result_1, sampcovg[g], result);
        for(j=0; j<p; j++) {result_2[g*p+j] = result[j];}
	}
       
	/* lambda*theta*lambda' */
    for(g=0;g<G;g++){
		mx_trans(p,q,lambda,temp);
        mx_mult(p,q,q,lambda,theta[g],result_1); 
        mx_mult_diag1(p,q,result_1,temp,result);
        for(j=0; j<p; j++) {result_3[g*p+j] = result[j];}
	}
        
    for(g=0;g<G;g++)
	    for(i=0;i<p;i++)
             	psi[g*p+i] = sampcovg[g][i*p+i]-2*result_2[g*p+i]+result_3[g*p+i];
	
	free(temp); free(result_1); free(result_2); 
	free(result_3); free(result);
}

/* UCU*/
void update_psi_ucu(double *psi, double **lambda, double **beta, double **sampcov, int p, int q, double *pi, int G){ 
	
	int i,g,j;
        double *result_1 = malloc(sizeof(double)*p*p);
        double *result_2 = malloc(sizeof(double)*G*p);
        double *result = malloc(sizeof(double)*p);


	/* lambda*beta*sampcovtilde */
    for (g=0; g<G; g++){
	    mx_mult(p,q,p,lambda[g], beta[g], result_1);
        mx_mult_diag1(p,p,result_1, sampcov[g], result);
            for(j=0; j<p; j++) {result_2[g*p+j] = result[j];}
	}
	
    for(i=0;i<p;i++){
		psi[i]=0.0;
		for(g=0;g<G;g++) psi[i] += pi[g]*(sampcov[g][i*p+i]-result_2[i+g*p]); 
	}
	
	free(result_1); free(result_2); free(result);
}

/* UCC*/
double update_psi_ucc(double **lambda, double **beta, double **sampcov, int p, int q, double *pi, int G){ 
	
	int i,g,j;
        double psi;
        double *result_1 = malloc(sizeof(double)*p*p);
        double *result_2 = malloc(sizeof(double)*G*p);
        double *result = malloc(sizeof(double)*p);

	/* lambda*beta*sampcovtilde */
	for (g=0; g<G; g++){
		mx_mult(p,q,p,lambda[g], beta[g], result_1);
        mx_mult_diag1(p,p,result_1, sampcov[g], result);
            for(j=0; j<p; j++) {result_2[g*p+j] = result[j];}
	}
	
    psi=0.0;
	for(g=0;g<G;g++)
	    for(i=0;i<p;i++)
			psi += pi[g]*(sampcov[g][i*p+i]-result_2[g*p+i]);
	psi /= p;
	
	free(result_1); free(result_2); free(result);

	return psi;
}

double update_det_sigma_NEW(double *lambda, double psi, double log_detpsi, int p, int q){

	int i, j;
        double det[1];

        double *tmp = malloc(sizeof(double)*p*p);
        double *tmp2 = malloc(sizeof(double)*p*p);

	update_beta1(tmp2, psi, lambda, p, q);
	mx_mult(q,p,q,tmp2, lambda, tmp);
    for(i=0;i<q;i++){
		for(j=0;j<q;j++){
			tmp2[i*q+j] = tmp[i*q+j]*(-1.0);
			if(i==j) tmp2[i*q+i] += 1.0;
		}
	}
	GaussJordan(q,tmp2,tmp,det);

	free(tmp); free(tmp2);

	return (log_detpsi - log(det[0]));
}

double update_det_sigma_NEW2(double *lambda, double *psi, double log_detpsi, int p, int q){

	int i, j;
        double det[1];

        double *tmp = malloc(sizeof(double)*p*p);
        double *tmp2 = malloc(sizeof(double)*p*p);

	update_beta2(tmp2, psi, lambda, p, q);
	mx_mult(q,p,q,tmp2, lambda, tmp);
    for(i=0;i<q;i++){
		for(j=0;j<q;j++){
			tmp2[i*q+j] = tmp[i*q+j]*(-1.0);
			if(i==j) tmp2[i*q+i] += 1.0;
		}
	}
	GaussJordan(q,tmp2,tmp,det);

	free(tmp); free(tmp2);

	return (log_detpsi - log(det[0]));
}

double update_det_sigma(double *lambda, double *sigma, double log_detpsi, int p, int q){

	int i,j;
	double siginv, det[1];
	
        double *tmp = malloc(sizeof(double)*p*p);
        double *r_1 = malloc(sizeof(double)*p*p);
        double *r_2 = malloc(sizeof(double)*p*p);
        double *r_3 = malloc(sizeof(double)*p*p);

	
	/* |I-lambda'*sigma*lambda| */
    mx_trans(p, q, lambda, tmp);
    mx_mult(q, p, p, tmp, sigma, r_1);
    mx_mult(q,p,q,r_1,lambda,r_3);
    for(i=0;i<q;i++){
       	for(j=0;j<q;j++){
        	r_3[i*q+j] = 0-r_3[i*q+j];
            if(i==j) r_3[i*q+i] += 1.0;
        }
    }
    GaussJordan(q,r_3,r_2,det);
    siginv = log_detpsi - log(det[0]);
	
	free(tmp); free(r_1); free(r_2); free (r_3);
	/*if (siginv<0) printf("-");*/
	return siginv;
}

/***************** Update Z values WITHOUT Sigma ************************/
int update_z(double *v, double *x, double *z, double *lambda, double psi, double *mu, 
	      double *pi, double *max_v, double log_c, int N, int G, int p, int q){
    
    int i,g,j;
    double a,e,d_alt;        
    double *x0 = malloc(sizeof(double)*p);
    double *mu0 = malloc(sizeof(double)*p);
    double *v0 = malloc(sizeof(double)*G);
 

    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
           for(j=0; j<p; j++){
               x0[j] = x[i*p+j];
               mu0[j]= mu[g*p+j];
           }
			a = woodbury(x0, lambda, psi, mu0, p, q);
			e = a/2.0*(-1.0);
			v[i*G+g] = e + log(pi[g]) - log_c;
        }/* End g loop */
                for(j=0;j<G;j++) {v0[j]=v[i*G+j];}
		max_v[i] = maximum_array(v0,G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i*G+g]-max_v[i]); 
		for(g=0;g<G;g++) z[i*G+g]=exp(v[i*G+g]-max_v[i])/d_alt;
                                 
    }/* End i loop */  
     free(x0); free(mu0); free(v0); 

    return 0;
}

/*************************** Woodbury trick **********************************/
double woodbury(double *x, double *lambda, double psi, double *mu, int p, int q){
	int i,j;
        double det[1], lhs, rhs;
        double *lvec = malloc(sizeof(double)*p);
        double *tvec = malloc(sizeof(double)*p);
        double *cvec = malloc(sizeof(double)*p);
        double *temp = malloc(sizeof(double)*q*p);
        double *cp = malloc(sizeof(double)*q*p);
        double *result = malloc(sizeof(double)*q*p);


	/* Compute LHS */
	lhs = 0.0;
	for (j=0;j<p;j++)
		lhs += (x[j]-mu[j])*(x[j]-mu[j]);
	lhs /= psi; 

    /* Compute RHS */
	for (j=0;j<p;j++)
		lvec[j] = (x[j]-mu[j])/psi;
	vec_mx_mult(p,q,lvec,lambda,tvec);
 	/* Lambda'/psi */
    mx_trans(p, q, lambda, temp);
    for(i=0;i<q;i++)
		for(j=0;j<p;j++)
            temp[i*p+j]/=psi;
 	/* Lambda'/psi* Lambda */
    mx_mult(q, p, q, temp, lambda, result);
    /* (I + Lambda'/psi*Lambda)^{-1} */
    for(i=0; i<q;i++) result[i*q+i] += 1.0;
    GaussJordan(q, result, cp, det);
    mx_trans(p, q, lambda, temp);
    mx_mult(q, q, p, cp, temp, result);
	vec_mx_mult(q,p,tvec,result,cvec);
	rhs = 0.0;
	for (j=0;j<p;j++)
		rhs += cvec[j]*(x[j]-mu[j]);
	rhs /= psi; 
	
	free(lvec); free(cvec); free(tvec); free(result); free(temp); free(cp);

	return (lhs-rhs);
}

/***************** Update Z values WITHOUT Sigma ************************/
int update_z2(double *v, double *x, double *z, double *lambda, double *psi, double *mu, 
	      double *pi, double *max_v, double log_c, int N, int G, int p, int q){
    
    int i,j,g;
    double a,e,d_alt;        
    double *x0 = malloc(sizeof(double)*p);    
    double *mu0 = malloc(sizeof(double)*p);    
    double *v0 = malloc(sizeof(double)*G);    
    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
           for(j=0; j<p; j++){
               x0[j] = x[i*p+j];
               mu0[j]= mu[g*p+j];
           }
			a = woodbury2(x0, lambda, psi, mu0, p, q);
			e = a/2.0*(-1.0);
			v[i*G+g] = e + log(pi[g]) - log_c;
        }/* End g loop */
                for(j=0;j<G;j++) {v0[j]=v[i*G+j];}
		max_v[i] = maximum_array(v0,G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i*G+g]-max_v[i]); 
		for(g=0;g<G;g++) z[i*G+g]=exp(v[i*G+g]-max_v[i])/d_alt;
    }/* End i loop */
    free(x0); free(mu0); free(v0);  
    return 0;
}

/*************************** Woodbury trick **********************************/
double woodbury2(double *x, double *lambda, double *psi, double *mu, int p, int q){
	int i,j;
        double det[1], lhs, rhs;
        double *lvec = malloc(sizeof(double)*p);
        double *tvec = malloc(sizeof(double)*p);
        double *cvec = malloc(sizeof(double)*p);
        double *temp = malloc(sizeof(double)*q*p);
        double *cp = malloc(sizeof(double)*q*p);
        double *result = malloc(sizeof(double)*q*p);

	/* Compute LHS */
	lhs = 0.0;
	for (j=0;j<p;j++)
		lhs += (x[j]-mu[j])*(x[j]-mu[j])/psi[j];

    /* Compute RHS */
	for (j=0;j<p;j++)
		lvec[j] = (x[j]-mu[j])/psi[j];
	vec_mx_mult(p,q,lvec,lambda,tvec);
 	/* Lambda'/psi */
    mx_trans(p, q, lambda, temp);
    for(i=0;i<q;i++)
		for(j=0;j<p;j++)
            temp[i*p+j]/=psi[j];
 	/* Lambda'/psi* Lambda */
    mx_mult(q, p, q, temp, lambda, result);
    /* (I + Lambda'/psi*Lambda)^{-1} */
    for(i=0; i<q;i++) result[i*q+i] += 1.0;
    GaussJordan(q, result, cp, det);
    mx_trans(p, q, lambda, temp);
    mx_mult(q, q, p, cp, temp, result);
	vec_mx_mult(q,p,tvec,result,cvec);
	rhs = 0.0;
	for (j=0;j<p;j++)
		rhs += cvec[j]*(x[j]-mu[j])/psi[j];
	
	free(lvec); free(cvec); free(tvec); free(result); free(temp); free(cp);

	return (lhs-rhs);
}

/***************** Update Z values WITHOUT Sigma ************************/
int update_z3(double *v, double *x, double *z, double *lambda, double *psi, double *mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,j,g;
    double a,e,d_alt;        
    double *x0 = malloc(sizeof(double)*p);
    double *mu0 = malloc(sizeof(double)*p);
    double *v0 = malloc(sizeof(double)*G);
    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
           for(j=0; j<p; j++){
               x0[j] = x[i*p+j];
               mu0[j]= mu[g*p+j];
           } 
			a = woodbury(x0, lambda, psi[g], mu0, p, q);
			e = a/2.0*(-1.0);
			v[i*G+g] = e + log(pi[g]) - log_c[g];
        }/* End g loop */
                for(j=0;j<G;j++) {v0[j]=v[i*G+j];}
		max_v[i] = maximum_array(v0,G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i*G+g]-max_v[i]); 
		for(g=0;g<G;g++) z[i*G+g]=exp(v[i*G+g]-max_v[i])/d_alt;
    }/* End i loop */  
    free(x0); free(mu0); free(v0);
    return 0;
}

/***************** Update Z values WITHOUT Sigma ************************/
int update_z4(double *v, double *x, double *z, double *lambda, double *psi, double *mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,j,g;
    double a,e,d_alt;        
    double *x0 = malloc(sizeof(double)*p);
    double *mu0 = malloc(sizeof(double)*p);
    double *v0 = malloc(sizeof(double)*G);
    double *psi0 = malloc(sizeof(double)*p);
    
    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
           for(j=0; j<p; j++){
               x0[j] = x[i*p+j];
               mu0[j]= mu[g*p+j];
               psi0[j]= psi[g*p+j];
           } 

			a = woodbury2(x0, lambda, psi0, mu0, p, q);
			e = a/2.0*(-1.0);
			v[i*G+g] = e + log(pi[g]) - log_c[g];
        }/* End g loop */
                for(j=0;j<G;j++) {v0[j]=v[i*G+j];}
		max_v[i] = maximum_array(v0,G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i*G+g]-max_v[i]); 
		for(g=0;g<G;g++) z[i*G+g]=exp(v[i*G+g]-max_v[i])/d_alt;
    }/* End i loop */  
    free(x0); free(mu0); free(v0); free(psi0);

    return 0;
}

int update_z5(double *v, double *x, double *z, double **lambda, double psi, double *mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,j,g;
    double a,e,d_alt;        
    double *x0 = malloc(sizeof(double)*p);
    double *mu0 = malloc(sizeof(double)*p);
    double *v0 = malloc(sizeof(double)*G);
   
    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
           for(j=0; j<p; j++){
               x0[j] = x[i*p+j];
               mu0[j]= mu[g*p+j];
           } 

			a = woodbury(x0, lambda[g], psi, mu0, p, q);
			e = a/2.0*(-1.0);
			v[i*G+g] = e + log(pi[g]) - log_c[g];
        }/* End g loop */
                for(j=0;j<G;j++) {v0[j]=v[i*G+j];}
		max_v[i] = maximum_array(v0,G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i*G+g]-max_v[i]); 
		for(g=0;g<G;g++) z[i*G+g]=exp(v[i*G+g]-max_v[i])/d_alt;
    }/* End i loop */  
    free(x0); free(mu0); free(v0);
    return 0;
}

int update_z6(double *v, double *x, double *z, double **lambda, double *psi, double *mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,j,g;
    double a,e,d_alt;        
    double *x0 = malloc(sizeof(double)*p);
    double *mu0 = malloc(sizeof(double)*p);
    double *v0 = malloc(sizeof(double)*G);

    
    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
           for(j=0; j<p; j++){
               x0[j] = x[i*p+j];
               mu0[j]= mu[g*p+j];
           } 

			a = woodbury2(x0, lambda[g], psi, mu0, p, q);
			e = a/2.0*(-1.0);
			v[i*G+g] = e + log(pi[g]) - log_c[g];
        }/* End g loop */
                for(j=0;j<G;j++) {v0[j]=v[i*G+j];}
		max_v[i] = maximum_array(v0,G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i*G+g]-max_v[i]); 
		for(g=0;g<G;g++) z[i*G+g]=exp(v[i*G+g]-max_v[i])/d_alt;
    }/* End i loop */  
    free(x0); free(mu0); free(v0);
    return 0;
}

int update_z7(double *v, double *x, double *z, double **lambda, double *psi, double *mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,j,g;
    double a,e,d_alt;        
    double *x0 = malloc(sizeof(double)*p);
    double *mu0 = malloc(sizeof(double)*p);
    double *v0 = malloc(sizeof(double)*G);
 
    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
           for(j=0; j<p; j++){
               x0[j] = x[i*p+j];
               mu0[j]= mu[g*p+j];
           } 

			a = woodbury(x0, lambda[g], psi[g], mu0, p, q);
			e = a/2.0*(-1.0);
			v[i*G+g] = e + log(pi[g]) - log_c[g];
        }/* End g loop */
                for(j=0;j<G;j++) {v0[j]=v[i*G+j];}
		max_v[i] = maximum_array(v0,G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i*G+g]-max_v[i]); 
		for(g=0;g<G;g++) z[i*G+g]=exp(v[i*G+g]-max_v[i])/d_alt;
    }/* End i loop */  
    free(x0); free(mu0); free(v0);
    return 0;
}

int update_z8(double *v, double *x, double *z, double **lambda, double *psi, double *mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,j,g;
    double a,e,d_alt;        
    double *x0 = malloc(sizeof(double)*p);
    double *mu0 = malloc(sizeof(double)*p);
    double *v0 = malloc(sizeof(double)*G);
    double *psi0 = malloc(sizeof(double)*p);
 
    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
           for(j=0; j<p; j++){
               x0[j] = x[i*p+j];
               mu0[j]= mu[g*p+j];
               psi0[j]= psi[g*p+j];
           } 

			a = woodbury2(x0, lambda[g], psi0, mu0, p, q);
			e = a/2.0*(-1.0);
			v[i*G+g] = e + log(pi[g]) - log_c[g];
        }/* End g loop */
                for(j=0;j<G;j++) {v0[j]=v[i*G+j];}
		max_v[i] = maximum_array(v0,G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i*G+g]-max_v[i]); 
		for(g=0;g<G;g++) {z[i*G+g]=exp(v[i*G+g]-max_v[i])/d_alt;
//                                  Rprintf(" z is %f\n", z[i*G+g]);
                                 }
    }/* End i loop */  
    free(x0); free(mu0); free(v0); free(psi0);
    return 0;
}

int update_z9(double *v, double *x, double *z, double *lambda, double *omega, double *delta, double *mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g,j,k;
    double a,e,d_alt;
    double *psi = malloc(sizeof(double)*p);
    double *x0 = malloc(sizeof(double)*p);
    double *mu0 = malloc(sizeof(double)*p);
    double *v0 = malloc(sizeof(double)*G);

    
    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
			for(j=0;j<p;j++) psi[j]=omega[g]*delta[j];
           for(k=0; k<p; k++){
               x0[k] = x[i*p+k];
               mu0[k]= mu[g*p+k];
           } 

			a = woodbury2(x0, lambda, psi, mu0, p, q);
			e = a/2.0*(-1.0);
			v[i*G+g] = e + log(pi[g]) - log_c[g];
        }/* End g loop */
                for(k=0;k<G;k++) {v0[k]=v[i*G+k];}
		max_v[i] = maximum_array(v0,G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i*G+g]-max_v[i]); 
		for(g=0;g<G;g++) {z[i*G+g]=exp(v[i*G+g]-max_v[i])/d_alt;
                                 }
    }/* End i loop */  
     free(psi); free(x0); free(mu0); free(v0);
    return 0;
}

int update_z10(double *v, double *x, double *z, double **lambda, double *omega, double *delta, double *mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g,j,k;
    double a,e,d_alt;

    double *psi = malloc(sizeof(double)*p);
    double *x0 = malloc(sizeof(double)*p);
    double *mu0 = malloc(sizeof(double)*p);
    double *v0 = malloc(sizeof(double)*G);

    
    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
			for(j=0;j<p;j++) psi[j]=omega[g]*delta[j];
           for(k=0; k<p; k++){
               x0[k] = x[i*p+k];
               mu0[k]= mu[g*p+k];
           }

			a = woodbury2(x0, lambda[g], psi, mu0, p, q);
			e = a/2.0*(-1.0);
			v[i*G+g] = e + log(pi[g]) - log_c[g];
        }/* End g loop */
                for(k=0;k<G;k++) {v0[k]=v[i*G+k];}
		max_v[i] = maximum_array(v0,G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i*G+g]-max_v[i]); 
		for(g=0;g<G;g++) z[i*G+g]=exp(v[i*G+g]-max_v[i])/d_alt;
    }/* End i loop */  
     free(psi); free(x0); free(mu0); free(v0);
	return 0;
}

int update_z11(double *v, double *x, double *z, double *lambda, double omega, double *delta, double *mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g,j,k;
    double a,e,d_alt;

    double *psi = malloc(sizeof(double)*p);
    double *x0 = malloc(sizeof(double)*p);
    double *mu0 = malloc(sizeof(double)*p);
    double *v0 = malloc(sizeof(double)*G);

    
    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
			for(j=0;j<p;j++) psi[j]=omega*delta[g*p+j];
           for(k=0; k<p; k++){
               x0[k] = x[i*p+k];
               mu0[k]= mu[g*p+k];
           }

			a = woodbury2(x0, lambda, psi, mu0, p, q);
			e = a/2.0*(-1.0);
			v[i*G+g] = e + log(pi[g]) - log_c[g];
        }/* End g loop */
                for(k=0;k<G;k++) {v0[k]=v[i*G+k];}
		max_v[i] = maximum_array(v0,G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i*G+g]-max_v[i]); 
		for(g=0;g<G;g++) z[i*G+g]=exp(v[i*G+g]-max_v[i])/d_alt;
    }/* End i loop */  
    free(psi); free(x0); free(mu0); free(v0);
	return 0;
}

int update_z12(double *v, double *x, double *z, double **lambda, double omega, double *delta, double *mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g,j,k;
    double a,e,d_alt;

    double *psi = malloc(sizeof(double)*p);
    double *x0 = malloc(sizeof(double)*p);
    double *mu0 = malloc(sizeof(double)*p);
    double *v0 = malloc(sizeof(double)*G);

    
    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
			for(j=0;j<p;j++) psi[j]=omega*delta[g*p+j];
           for(k=0; k<p; k++){
               x0[k] = x[i*p+k];
               mu0[k]= mu[g*p+k];
           }

			a = woodbury2(x0, lambda[g], psi, mu0, p, q);
			e = a/2.0*(-1.0);
			v[i*G+g] = e + log(pi[g]) - log_c[g];
        }/* End g loop */
                for(k=0;k<G;k++) {v0[k]=v[i*G+k];}
		max_v[i] = maximum_array(v0,G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i*G+g]-max_v[i]); 
		for(g=0;g<G;g++) z[i*G+g]=exp(v[i*G+g]-max_v[i])/d_alt;
    }/* End i loop */  
      free(psi);  free(x0); free(mu0); free(v0);
	return 0;
}

double update_omega(double *lambda, double *delta, double *beta, double *sampcovg, double *theta, int p, int q){ 
	
	int i;
	double omega;

        double *temp = malloc(sizeof(double)*q*p);
        double *result_1 = malloc(sizeof(double)*p*p);
        double *result_2 = malloc(sizeof(double)*p);
        double *result_3 = malloc(sizeof(double)*p);


	/* lambda*beta*sampcov */
    mx_mult(p,q,p,lambda, beta, result_1);
    mx_mult_diag1(p,p,result_1, sampcovg, result_2);
       
	/* lambda*theta*lambda' */
    mx_trans(p,q,lambda,temp);
    mx_mult(p,q,q,lambda,theta,result_1); 
    mx_mult_diag1(p,q,result_1,temp,result_3);
        
    omega = 0.0;
    for(i=0;i<p;i++)
        omega += (sampcovg[i*p+i]-2*result_2[i]+result_3[i])/delta[i];
    omega /= p;
	
	free(temp); free(result_1); free(result_2);free(result_3);
	
	return omega;
}

double update_omega2(double *lambda, double *delta, double *beta, double *sampcovg, int p, int q){ 
	
	int i;
	double omega;

        double *result_1 = malloc(sizeof(double)*p*p);
        double *result_2 = malloc(sizeof(double)*p);

	/* lambda*beta*sampcov */
    mx_mult(p,q,p,lambda, beta, result_1);
    mx_mult_diag1(p,p,result_1, sampcovg, result_2);
       
    omega = 0.0;
    for(i=0;i<p;i++)
        omega += (sampcovg[i*p+i]-result_2[i])/delta[i];
    omega /= p;
	
	free(result_1); free(result_2); 
	
	return omega;
}


void update_delta(double *delta, double *lambda, double *omega, double **beta, double **sampcov, double **theta, 
		double *n, int p, int q, int N, int G){

	int i,g,j;
	double lagrange = 0.0; /* Lagrange Multiplier */
        double *temp = malloc(sizeof(double)*q*p);
        double *result_1 = malloc(sizeof(double)*p*p);
        double *result_2 = malloc(sizeof(double)*G*p);
        double *result_3 = malloc(sizeof(double)*G*p);
        double *temp1 = malloc(sizeof(double)*p);
        double *result = malloc(sizeof(double)*p);
	/* lambda*beta_g*sampcov_g */
    for(g=0;g<G;g++){
		mx_mult(p,q,p,lambda, beta[g], result_1);
       	mx_mult_diag1(p,p,result_1, sampcov[g], result);
        for(j=0;j<p;j++)
           result_2[g*p+j]=result[j]; 
    }
	
	/* lambda*theta_g*lambda' */
    for(g=0;g<G;g++){
		mx_trans(p,q,lambda,temp);
       	mx_mult(p,q,q,lambda,theta[g],result_1); 
       	mx_mult_diag1(p,q,result_1,temp,result);
        for(j=0;j<p;j++)
           result_3[g*p+j]=result[j]; 
	}
	
    for(i=0;i<p;i++){
        temp1[i] = 0.0;
	    for(g=0;g<G;g++)
		    temp1[i] += (sampcov[g][i*p+i]-2.0*result_2[i+g*p]+result_3[i+g*p])*n[g]/omega[g];
	    lagrange += log(temp1[i]);
	}	
	
	/* Compute Lagrange Multiplier */
	lagrange /= (double)p;
	lagrange = exp(lagrange);
	lagrange -= (double)N; 
	lagrange /= (double)2; 

	for(i=0;i<p;i++){ delta[i] = temp1[i]/((double)N+(2.0*lagrange));
                        }
		
	free(temp); free(result_1); free(result_2); 
	free(result_3); free(temp1); free(result);
}

void update_delta2(double *delta, double **lambda, double *omega, double **beta, double **sampcov, double **theta, 
		double *n, int p, int q, int N, int G){

	int i,g,j;
        double toprint, lagrange_alt=0.0;
        double *temp = malloc(sizeof(double)*q*p);
        double *result_1 = malloc(sizeof(double)*p*p);
        double *result_2 = malloc(sizeof(double)*G*p);
        double *result_3 = malloc(sizeof(double)*G*p);
        double *temp1 = malloc(sizeof(double)*p);
        double *result1 = malloc(sizeof(double)*p);
        double *result2 = malloc(sizeof(double)*p);


	/* lambda*beta_g*sampcov_g */
    for(g=0;g<G;g++){
		mx_mult(p,q,p,lambda[g], beta[g], result_1);
       	mx_mult_diag1(p,p,result_1, sampcov[g], result1);
        for(j=0;j<p;j++)
           result_2[g*p+j]=result1[j]; 
    }
	
	/* lambda*theta_g*lambda' */
    for(g=0;g<G;g++){
		mx_trans(p,q,lambda[g],temp);
        mx_mult(p,q,q,lambda[g],theta[g],result_1); 
        mx_mult_diag1(p,q,result_1,temp,result2);
        for(j=0;j<p;j++)
           result_3[g*p+j]=result2[j]; 
	}
	
    for(i=0;i<p;i++){
        temp1[i] = 0.0;
	    for(g=0;g<G;g++){
		temp1[i] += ((sampcov[g][i*p+i]-2*result_2[i+g*p]+result_3[i+g*p])*n[g]/omega[g]);
	    lagrange_alt += log(temp1[i]);
	}	
	}
	lagrange_alt /= p;
	toprint = exp(lagrange_alt);
	lagrange_alt = toprint;
	lagrange_alt -= N; 
	lagrange_alt /= 2; 
	
	for(i=0;i<p;i++){
	    delta[i] = temp1[i]/(N+(2*lagrange_alt));
        } 
	
	free(temp); free(result_1); free(result_2); 
	free(result_3); free(temp1); free(result1); free(result2);
}

/* Omega not a vector */
void update_delta3(double *delta, double *lambda, double omega, double *beta, double *sampcov, double *theta, 
		double n, int p, int q){

	int i;
	double lagrange = 0.0; /* Lagrange Multiplier THIS WAS CHANGED FROM 1.0 NOV 2008*/
        double *temp = malloc(sizeof(double)*p*p);
        double *result_1 = malloc(sizeof(double)*p*p);
        double *result_2 = malloc(sizeof(double)*p);
        double *result_3 = malloc(sizeof(double)*p);
        double *temp1 = malloc(sizeof(double)*p);


	/* lambda*beta_g*sampcov_g */
    mx_mult(p,q,p,lambda, beta, result_1);
    mx_mult_diag1(p,p,result_1, sampcov, result_2);
       		
	/* lambda*theta_g*lambda' */
    mx_trans(p,q,lambda,temp);
    mx_mult(p,q,q,lambda,theta,result_1); 
    mx_mult_diag1(p,q,result_1,temp,result_3);
		
    for(i=0;i<p;i++) temp1[i] = sampcov[i*p+i]-2*result_2[i]+result_3[i];
	
	for(i=0;i<p;i++) lagrange += log(temp1[i]);
		
	/* Compute Lagrange Multiplier */
	lagrange /= p; 
	lagrange = exp(lagrange);
	lagrange /= omega;
	lagrange -= 1; 
	lagrange *= (n/2); 

	for(i=0;i<p;i++) {delta[i] = temp1[i]/((1+(2*lagrange/n))*omega);
                         }
		
	free(temp); free(result_1); free(result_2); 
	free(result_3); free(temp1);
}


