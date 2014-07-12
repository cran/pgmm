#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "functions.h"

void update_n(double *n, double **z, int G, int N){    
	int g,i;	    
	for(g=0; g<G; g++){
            n[g]=0.0;
            for(i=0; i<N; i++)
                n[g] += z[i][g];
    }
}
    
void update_pi(double *pi, double *n, int G, int N){    
	int g;	    
	for(g=0; g<G; g++) pi[g] = n[g]/N; 
}
    
void update_mu(double **mu, double *n, double **x, double **z,int G, 
	    int N, int p){    
	int i,j,g;
	for(g=0; g<G; g++){
            for(j=0; j<p; j++){
                mu[g][j]=0.0;
                for(i=0; i<N; i++) {
                    mu[g][j] += z[i][g]*x[i][j]; 
				}
                mu[g][j] /= n[g];
            }
    }
}

void update_stilde(double **sampcovtilde, double **x, double **z, double **mu,
	int G, int N, int p){
        int i,j,k,g;
	for(j=0; j<p; j++){
            for(k=0; k<=j; k++){
                sampcovtilde[j][k]=0.0;
                for(g=0;g<G;g++)    
                    for(i=0; i<N; i++)
                        sampcovtilde[j][k]+=z[i][g]*(x[i][j]-mu[g][j])*(x[i][k]-mu[g][k]);
				sampcovtilde[j][k]/=N;
                sampcovtilde[k][j] = sampcovtilde[j][k];		
            }
    }
}

void update_sg(double ***sg, double **x, double **z, double **mu,
	double *n, int p, int G, int N){
        int i,j,k,g;
	
	for(g=0;g<G;g++){    
	   for(j=0; j<p; j++){
               for(k=0; k<=j; k++){
                    sg[g][j][k]=0.0;
                    for(i=0; i<N; i++){
                        sg[g][j][k]+=z[i][g]*(x[i][j]-mu[g][j])*(x[i][k]-mu[g][k])/n[g];
                    }
                    sg[g][k][j]=sg[g][j][k];
               }
        }
	}
}

void update_beta1(double **beta, double psi, double **lambda, int p, int q){
	int i,j;
	double **lhs, **rhs, **cp, **result;
	double det[1];
	my_alloc(&lhs,p,p);
	my_alloc(&rhs,p,p);
	my_alloc(&cp,p,p);
	my_alloc(&result,p,p);

 	/* Lambda'/psi */
    mx_trans(p, q, lambda, lhs);
    for(i=0;i<q;i++)
		for(j=0;j<p;j++)
            lhs[i][j]/=psi; /*psi[j]*/
 	/* Lambda'/psi* Lambda */
    mx_mult(q, p, q, lhs, lambda, cp);
        
    /* (I + Lambda'/psi*Lambda)^{-1} */
    for(i=0; i<q;i++){
		for(j=0;j<q;j++){
			result[i][j] = cp[i][j];
			if (i==j) result[i][i] += 1.0;
		}
	}
    GaussJordan(q, result, rhs, det);
    /* work out rhs*/
	mx_mult(q,q,q,cp,rhs,result);
    mx_mult(q,q,p,result,lhs,rhs);
	
	for(i=0;i<q;i++)
		for(j=0;j<p;j++)
			beta[i][j] = lhs[i][j] - rhs[i][j];       
 
	my_free(lhs,p); my_free(result,p); my_free(rhs,p); my_free(cp,p);
}

void update_beta2(double **beta, double *Psi, double **lambda, int p, int q){
	int i,j;
	double **lhs, **rhs, **cp, **result;
	double det[1];
	my_alloc(&lhs,p,p);
	my_alloc(&rhs,p,p);
	my_alloc(&cp,p,p);
	my_alloc(&result,p,p);

 	/* Lambda'/psi */
    mx_trans(p, q, lambda, lhs);
    for(i=0;i<q;i++)
		for(j=0;j<p;j++)
            lhs[i][j]/=Psi[j];
 	/* Lambda'/psi* Lambda */
    mx_mult(q, p, q, lhs, lambda, cp);
        
    /* (I + Lambda'/psi*Lambda)^{-1} */
    for(i=0; i<q;i++){
		for(j=0;j<q;j++){
			result[i][j] = cp[i][j];
			if (i==j) result[i][i] += 1.0;
		}
	}
    GaussJordan(q, result, rhs, det);
    /* work out rhs*/
	mx_mult(q,q,q,cp,rhs,result);
    mx_mult(q,q,p,result,lhs,rhs);
	
	for(i=0;i<q;i++)
		for(j=0;j<p;j++)
			beta[i][j] = lhs[i][j] - rhs[i][j];       
 
	my_free(lhs,p); my_free(result,p); my_free(rhs,p); my_free(cp,p);
}

void update_theta(double **theta, double **beta, double **lambda, double **sampcovtilde, int p, int q){ 
	int i,j;
	double **id_q, **r_1, **r_2, **r_3, **tmp;
	
	my_alloc(&tmp,p,p);
	my_alloc(&r_1,q,q);
	my_alloc(&r_2,q,p);
	my_alloc(&r_3,q,q);
	my_alloc(&id_q,q,q);
	generate_identity(q, id_q);	
	
	/* Work out beta*lambda */
    mx_mult(q, p, q, beta, lambda, r_1);
    
    /* Now, work out beta*sampcovtilde*beta' */
    mx_mult(q, p, p, beta, sampcovtilde, r_2);
    mx_trans(q, p, beta, tmp);    
    mx_mult(q, p, q, r_2, tmp, r_3);    
                 
    for(i=0;i<q;i++)
		for(j=0;j<q;j++)
			theta[i][j] = id_q[i][j]-r_1[i][j]+r_3[i][j];
       
	my_free(id_q,q); my_free(tmp,p); my_free(r_1,q); my_free(r_2,q); 
	my_free(r_3,q);
}

void update_lambda(double **lambda, double **beta, double **s, double **theta, int p, int q){
	
	int i,j;
	double **tran, **res1, **res2, **res3, det[1];
	
	my_alloc(&tran,p,q);
	my_alloc(&res1,p,q);
	my_alloc(&res2,q,q);
	my_alloc(&res3,q,q);
	
	/* sampcovtilde*beta'*/
    mx_trans(q,p,beta,tran);
    mx_mult(p,p,q,s,tran,res1);
    	
	/* Make of copy of theta ahead of Gauss-Jordan */
	for(i=0;i<q;i++)
        for(j=0;j<q;j++)
            res3[i][j] = theta[i][j];   
        
	/* lambda=sampcovtilde*beta'*theta^{-1} */
    GaussJordan(q, res3, res2,det);
    mx_mult(p,q,q,res1,res2,lambda);
    
	my_free(tran,p); my_free(res1,p); my_free(res2,q); my_free(res3,q);
}

/* For CUU */
void update_lambda_cuu(double **lambda, double ***beta, double ***s, 
	double ***theta, double *n, double **Psi, int p, int q, int G){
	
	int i,j,k,g;
	double **tran, **res1, **res2, **res3, det[1];
	
	my_alloc(&tran,p,q);
	my_alloc(&res1,p,q);
	my_alloc(&res2,p,q);
	my_alloc(&res3,q,q);	
	
    /* Compute the RHS --- only needs to happen once */
	/* sum_g[sampcov_g*beta_g'*n_g/psi_g] */
    for(g=0;g<G;g++){
		mx_trans(q,p,beta[g],tran);
		mx_mult(p,p,q,s[g],tran,res1);
		if(g==0){ 
			for(i=0;i<p;i++)
				for(j=0;j<q;j++)
					res2[i][j] = res1[i][j]*n[g]/Psi[g][i];
		}else{
			for(i=0;i<p;i++)
				for(j=0;j<q;j++)
					res2[i][j] += res1[i][j]*n[g]/Psi[g][i];    
		}
    }

    /* Now solve for lambda row-by-row */
	for (i=0;i<p;i++){
		/* First, compute the theta sum*/
		for(g=0;g<G;g++){
			if(g==0){ 
				for(k=0;k<q;k++)
					for(j=0;j<q;j++)
						res3[k][j] = theta[g][k][j]*n[g]/Psi[g][i];
			}else{
				for(k=0;k<q;k++)
					for(j=0;j<q;j++)
						res3[k][j] += theta[g][k][j]*n[g]/Psi[g][i];    
			}
		}
		/* Invert theta sum */
		GaussJordan(q, res3, res1, det);
		
		/* Now solve for row i of lambda */
		vec_mx_mult(q, q, res2[i], res1, lambda[i]);
	}     	
	my_free(tran,p); my_free(res1,p); my_free(res2,p); my_free(res3,q);
}

/* For CUC */
void update_lambda2(double **lambda, double ***beta, double ***s, double ***theta, double *n, double *Psi, int p, 
		int q, int G){
	
	int i,j,g;
	double **tran, **res1, **res2, **res3, det[1];
	
	my_alloc(&tran,p,q);
	my_alloc(&res1,p,q);
	my_alloc(&res2,p,q);
	my_alloc(&res3,q,q);
	
    /* sum_g[sampcov_g*beta_g'*n_g/psi_g] */
    for(g=0;g<G;g++){
		mx_trans(q,p,beta[g],tran);
		mx_mult(p,p,q,s[g],tran,res1);
		if(g==0){ 
			for(i=0;i<p;i++)
				for(j=0;j<q;j++)
					res2[i][j] = res1[i][j]*n[g]/Psi[g];
		}else{
			for(i=0;i<p;i++)
				for(j=0;j<q;j++)
					res2[i][j] += res1[i][j]*n[g]/Psi[g];    
		}
    }

   	/* sum_g[theta_g'*n_g/psi_g]^{-1} */
    for(g=0;g<G;g++){
		if(g==0){ 
			for(i=0;i<q;i++)
				for(j=0;j<q;j++)
					res3[i][j] = theta[g][i][j]*n[g]/Psi[g];
		}else{
			for(i=0;i<q;i++)
				for(j=0;j<q;j++)
					res3[i][j] += theta[g][i][j]*n[g]/Psi[g];    
		}
    }
        
    GaussJordan(q, res3, res1,det); /* inverting the theta sum */
    mx_mult(p,q,q,res2,res1,lambda); /* this gives lambda */ 
	
	my_free(tran,p); my_free(res1,p); my_free(res2,p); my_free(res3,q);
}

double update_psi(double **lambda, double **beta, double **sampcovtilde, int p, int q){ 
	
	int i;
	double **result_1, *result_2; 
	double psi;
	my_alloc(&result_1,p,p);
	my_alloc_vec(&result_2,p);

	/* lambda*beta*sampcovtilde */
    mx_mult(p,q,p,lambda, beta, result_1);
    mx_mult_diag1(p,p,result_1, sampcovtilde, result_2);
       
    psi = 0.0;
    for(i=0;i<p;i++)
        psi += sampcovtilde[i][i]-result_2[i];
    psi /= p;
	
	my_free(result_1,p); free(result_2); 
	
	return psi;
}

/* Account for Psi_p*/
void update_psi2(double *psi, double **lambda, double **beta, double **sampcovtilde, int p, int q){ 
	
	int i;
	double **result_1, *result_2;
	my_alloc(&result_1,p,p);
	my_alloc_vec(&result_2,p);

	/* lambda*beta*sampcovtilde */
	mx_mult(p,q,p,lambda, beta, result_1);
    mx_mult_diag1(p,p,result_1, sampcovtilde, result_2);
        
    for(i=0;i<p;i++)
		psi[i] = sampcovtilde[i][i]-result_2[i];
	
	my_free(result_1,p); free(result_2); 
}

/* CUC case */
double update_psi3(double **lambda, double **beta, double **sampcovg, double **theta, int p, int q){ 
	
	int i;
	double **temp, **result_1, *result_2, *result_3; 
	double psi;
	my_alloc(&temp,q,p);
	my_alloc(&result_1,p,p);
	my_alloc_vec(&result_2,p);
	my_alloc_vec(&result_3,p);

	/* lambda*beta*sampcov */
    mx_mult(p,q,p,lambda, beta, result_1);
    mx_mult_diag1(p,p,result_1, sampcovg, result_2);
       
	/* lambda*theta*lambda' */
    mx_trans(p,q,lambda,temp);
    mx_mult(p,q,q,lambda,theta,result_1); 
    mx_mult_diag1(p,q,result_1,temp,result_3);
        
    psi = 0.0;
    for(i=0;i<p;i++)
        psi += sampcovg[i][i]-2*result_2[i]+result_3[i];
    psi /= p;
	
	my_free(temp,q); my_free(result_1,p); free(result_2); 
	free(result_3);
	
	return psi;
}

/* CUU case */
void update_psi_cuu(double **psi, double **lambda, double ***beta, double ***sampcovg, double 
	***theta, int p, int q, int G){ 
	
	int i, g;
	double **temp, **result_1, **result_2, **result_3;
	my_alloc(&temp,q,p);
	my_alloc(&result_1,p,p);
	my_alloc(&result_2,G,p);
	my_alloc(&result_3,G,p);

	/* lambda*beta*sampcov */
    for (g=0;g<G;g++){
		mx_mult(p,q,p,lambda, beta[g], result_1);
    	mx_mult_diag1(p,p,result_1, sampcovg[g], result_2[g]);
	}
       
	/* lambda*theta*lambda' */
    for(g=0;g<G;g++){
		mx_trans(p,q,lambda,temp);
        mx_mult(p,q,q,lambda,theta[g],result_1); 
        mx_mult_diag1(p,q,result_1,temp,result_3[g]);
	}
        
    for(g=0;g<G;g++)
	    for(i=0;i<p;i++)
             	psi[g][i] = sampcovg[g][i][i]-2*result_2[g][i]+result_3[g][i];
	
	my_free(temp,q); my_free(result_1,p); my_free(result_2,G); 
	my_free(result_3,G);
}

/* UCU*/
void update_psi_ucu(double *psi, double ***lambda, double ***beta, double ***sampcov, int p, int q, double *pi, int G){ 
	
	int i,g;
	double **result_1, **result_2;
	my_alloc(&result_1,p,p);
	my_alloc(&result_2,G,p);

	/* lambda*beta*sampcovtilde */
    for (g=0; g<G; g++){
	    mx_mult(p,q,p,lambda[g], beta[g], result_1);
        mx_mult_diag1(p,p,result_1, sampcov[g], result_2[g]);
	}
	
    for(i=0;i<p;i++){
		psi[i]=0.0;
		for(g=0;g<G;g++) psi[i] += pi[g]*(sampcov[g][i][i]-result_2[g][i]);
	}
	
	my_free(result_1,p); my_free(result_2,G); 
}

/* UCC*/
double update_psi_ucc(double ***lambda, double ***beta, double ***sampcov, int p, int q, double *pi, int G){ 
	
	int i,g;
	double **result_1, **result_2, psi;
	my_alloc(&result_1,p,p);
	my_alloc(&result_2,G,p);

	/* lambda*beta*sampcovtilde */
	for (g=0; g<G; g++){
		mx_mult(p,q,p,lambda[g], beta[g], result_1);
        mx_mult_diag1(p,p,result_1, sampcov[g], result_2[g]);
	}
	
    psi=0.0;
	for(g=0;g<G;g++)
	    for(i=0;i<p;i++)
			psi += pi[g]*(sampcov[g][i][i]-result_2[g][i]);
	psi /= p;
	
	my_free(result_1,p); my_free(result_2,G);

	return psi;
}

double update_det_sigma_NEW(double **lambda, double psi, double log_detpsi, int p, int q){

	int i, j;
	double **tmp, **tmp2, det[1];

	my_alloc(&tmp,p,p);
	my_alloc(&tmp2,p,p);

	update_beta1(tmp2, psi, lambda, p, q);
	mx_mult(q,p,q,tmp2, lambda, tmp);
    for(i=0;i<q;i++){
		for(j=0;j<q;j++){
			tmp2[i][j] = tmp[i][j]*(-1.0);
			if(i==j) tmp2[i][i] += 1.0;
		}
	}
	GaussJordan(q,tmp2,tmp,det);

	my_free(tmp,p); my_free(tmp2,p);

	return (log_detpsi - log(det[0]));
}

double update_det_sigma_NEW2(double **lambda, double *psi, double log_detpsi, int p, int q){

	int i, j;
	double **tmp, **tmp2, det[1];

	my_alloc(&tmp,p,p);
	my_alloc(&tmp2,p,p);

	update_beta2(tmp2, psi, lambda, p, q);
	mx_mult(q,p,q,tmp2, lambda, tmp);
    for(i=0;i<q;i++){
		for(j=0;j<q;j++){
			tmp2[i][j] = tmp[i][j]*(-1.0);
			if(i==j) tmp2[i][i] += 1.0;
		}
	}
	GaussJordan(q,tmp2,tmp,det);

	my_free(tmp,p); my_free(tmp2,p);

	return (log_detpsi - log(det[0]));
}

double update_det_sigma(double **lambda, double **sigma, double log_detpsi, int p, int q){

	int i,j;
	double siginv, det[1];
	double **r_1, **r_2, **r_3, **tmp;
	
	my_alloc(&tmp,p,p);
	my_alloc(&r_1,p,p);
	my_alloc(&r_2,p,p);
	my_alloc(&r_3,p,p);
	
	/* |I-lambda'*sigma*lambda| */
    mx_trans(p, q, lambda, tmp);
    mx_mult(q, p, p, tmp, sigma, r_1);
    mx_mult(q,p,q,r_1,lambda,r_3);
    for(i=0;i<q;i++){
       	for(j=0;j<q;j++){
        	r_3[i][j] = 0-r_3[i][j];
            if(i==j) r_3[i][i] += 1.0;
        }
    }
    GaussJordan(q,r_3,r_2,det);
    siginv = log_detpsi - log(det[0]);
	
	my_free(tmp,p); my_free(r_1,p); my_free(r_2,p); my_free (r_3,p);
	/*if (siginv<0) printf("-");*/
	return siginv;
}

/***************** Update Z values WITHOUT Sigma ************************/
int update_z(double **v, double **x, double **z, double **lambda, double psi, double **mu, 
	      double *pi, double *max_v, double log_c, int N, int G, int p, int q){
    
    int i,g;
    double a,e,d_alt;        
   
    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
			a = woodbury(x[i], lambda, psi, mu[g], p, q);
			e = a/2.0*(-1.0);
			v[i][g] = e + log(pi[g]) - log_c;
        }/* End g loop */
		max_v[i] = maximum_array(v[i],G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i][g]-max_v[i]); 
		for(g=0;g<G;g++) z[i][g]=exp(v[i][g]-max_v[i])/d_alt;
    }/* End i loop */  
    return 0;
}

/*************************** Woodbury trick **********************************/
double woodbury(double *x, double **lambda, double psi, double *mu, int p, int q){
	int i,j;
	double *lvec, *tvec, *cvec, **temp, **cp, **result, det[1], lhs, rhs;
	my_alloc_vec(&lvec,p);
	my_alloc_vec(&tvec,p);
	my_alloc_vec(&cvec,p);
	my_alloc(&temp,q,p);
	my_alloc(&cp,q,p);
	my_alloc(&result,q,p);

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
            temp[i][j]/=psi;
 	/* Lambda'/psi* Lambda */
    mx_mult(q, p, q, temp, lambda, result);
    /* (I + Lambda'/psi*Lambda)^{-1} */
    for(i=0; i<q;i++) result[i][i] += 1.0;
    GaussJordan(q, result, cp, det);
    mx_trans(p, q, lambda, temp);
    mx_mult(q, q, p, cp, temp, result);
	vec_mx_mult(q,p,tvec,result,cvec);
	rhs = 0.0;
	for (j=0;j<p;j++)
		rhs += cvec[j]*(x[j]-mu[j]);
	rhs /= psi; 
	
	free(lvec); free(cvec); free(tvec); my_free(result,q); my_free(temp,q); my_free(cp,q);

	return (lhs-rhs);
}

/***************** Update Z values WITHOUT Sigma ************************/
int update_z2(double **v, double **x, double **z, double **lambda, double *psi, double **mu, 
	      double *pi, double *max_v, double log_c, int N, int G, int p, int q){
    
    int i,g;
    double a,e,d_alt;        
    
    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
			a = woodbury2(x[i], lambda, psi, mu[g], p, q);
			e = a/2.0*(-1.0);
			v[i][g] = e + log(pi[g]) - log_c;
        }/* End g loop */
		max_v[i] = maximum_array(v[i],G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i][g]-max_v[i]); 
		for(g=0;g<G;g++) z[i][g]=exp(v[i][g]-max_v[i])/d_alt;
    }/* End i loop */  
    return 0;
}

/*************************** Woodbury trick **********************************/
double woodbury2(double *x, double **lambda, double *psi, double *mu, int p, int q){
	int i,j;
	double *lvec, *tvec, *cvec, **temp, **cp, **result, det[1], lhs, rhs;
	my_alloc_vec(&lvec,p);
	my_alloc_vec(&tvec,p);
	my_alloc_vec(&cvec,p);
	my_alloc(&temp,q,p);
	my_alloc(&cp,q,p);
	my_alloc(&result,q,p);

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
            temp[i][j]/=psi[j];
 	/* Lambda'/psi* Lambda */
    mx_mult(q, p, q, temp, lambda, result);
    /* (I + Lambda'/psi*Lambda)^{-1} */
    for(i=0; i<q;i++) result[i][i] += 1.0;
    GaussJordan(q, result, cp, det);
    mx_trans(p, q, lambda, temp);
    mx_mult(q, q, p, cp, temp, result);
	vec_mx_mult(q,p,tvec,result,cvec);
	rhs = 0.0;
	for (j=0;j<p;j++)
		rhs += cvec[j]*(x[j]-mu[j])/psi[j];
	
	free(lvec); free(cvec); free(tvec); my_free(result,q); my_free(temp,q); my_free(cp,q);

	return (lhs-rhs);
}

/***************** Update Z values WITHOUT Sigma ************************/
int update_z3(double **v, double **x, double **z, double **lambda, double *psi, double **mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g;
    double a,e,d_alt;        
    
    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
			a = woodbury(x[i], lambda, psi[g], mu[g], p, q);
			e = a/2.0*(-1.0);
			v[i][g] = e + log(pi[g]) - log_c[g];
        }/* End g loop */
		max_v[i] = maximum_array(v[i],G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i][g]-max_v[i]); 
		for(g=0;g<G;g++) z[i][g]=exp(v[i][g]-max_v[i])/d_alt;
    }/* End i loop */  
    return 0;
}

/***************** Update Z values WITHOUT Sigma ************************/
int update_z4(double **v, double **x, double **z, double **lambda, double **psi, double **mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g;
    double a,e,d_alt;        
    
    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
			a = woodbury2(x[i], lambda, psi[g], mu[g], p, q);
			e = a/2.0*(-1.0);
			v[i][g] = e + log(pi[g]) - log_c[g];
        }/* End g loop */
		max_v[i] = maximum_array(v[i],G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i][g]-max_v[i]); 
		for(g=0;g<G;g++) z[i][g]=exp(v[i][g]-max_v[i])/d_alt;
    }/* End i loop */  
    return 0;
}

int update_z5(double **v, double **x, double **z, double ***lambda, double psi, double **mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g;
    double a,e,d_alt;        
   
    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
			a = woodbury(x[i], lambda[g], psi, mu[g], p, q);
			e = a/2.0*(-1.0);
			v[i][g] = e + log(pi[g]) - log_c[g];
        }/* End g loop */
		max_v[i] = maximum_array(v[i],G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i][g]-max_v[i]); 
		for(g=0;g<G;g++) z[i][g]=exp(v[i][g]-max_v[i])/d_alt;
    }/* End i loop */  
    return 0;
}

int update_z6(double **v, double **x, double **z, double ***lambda, double *psi, double **mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g;
    double a,e,d_alt;        
    
    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
			a = woodbury2(x[i], lambda[g], psi, mu[g], p, q);
			e = a/2.0*(-1.0);
			v[i][g] = e + log(pi[g]) - log_c[g];
        }/* End g loop */
		max_v[i] = maximum_array(v[i],G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i][g]-max_v[i]); 
		for(g=0;g<G;g++) z[i][g]=exp(v[i][g]-max_v[i])/d_alt;
    }/* End i loop */  
    return 0;
}

int update_z7(double **v, double **x, double **z, double ***lambda, double *psi, double **mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g;
    double a,e,d_alt;        
    
    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
			a = woodbury(x[i], lambda[g], psi[g], mu[g], p, q);
			e = a/2.0*(-1.0);
			v[i][g] = e + log(pi[g]) - log_c[g];
        }/* End g loop */
		max_v[i] = maximum_array(v[i],G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i][g]-max_v[i]); 
		for(g=0;g<G;g++) z[i][g]=exp(v[i][g]-max_v[i])/d_alt;
    }/* End i loop */  
    return 0;
}

int update_z8(double **v, double **x, double **z, double ***lambda, double **psi, double **mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g;
    double a,e,d_alt;        
    
    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
			a = woodbury2(x[i], lambda[g], psi[g], mu[g], p, q);
			e = a/2.0*(-1.0);
			v[i][g] = e + log(pi[g]) - log_c[g];
        }/* End g loop */
		max_v[i] = maximum_array(v[i],G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i][g]-max_v[i]); 
		for(g=0;g<G;g++) z[i][g]=exp(v[i][g]-max_v[i])/d_alt;
    }/* End i loop */  
    return 0;
}

int update_z9(double **v, double **x, double **z, double **lambda, double *omega, double *delta, double **mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g,j;
    double a,e,d_alt, *psi;

	my_alloc_vec(&psi,p);
    
    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
			for(j=0;j<p;j++) psi[j]=omega[g]*delta[j];
			a = woodbury2(x[i], lambda, psi, mu[g], p, q);
			e = a/2.0*(-1.0);
			v[i][g] = e + log(pi[g]) - log_c[g];
        }/* End g loop */
		max_v[i] = maximum_array(v[i],G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i][g]-max_v[i]); 
		for(g=0;g<G;g++) z[i][g]=exp(v[i][g]-max_v[i])/d_alt;
    }/* End i loop */  
    return 0;
}

int update_z10(double **v, double **x, double **z, double ***lambda, double *omega, double *delta, double **mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g,j;
    double a,e,d_alt, *psi;

	my_alloc_vec(&psi,p);
    
    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
			for(j=0;j<p;j++) psi[j]=omega[g]*delta[j];
			a = woodbury2(x[i], lambda[g], psi, mu[g], p, q);
			e = a/2.0*(-1.0);
			v[i][g] = e + log(pi[g]) - log_c[g];
        }/* End g loop */
		max_v[i] = maximum_array(v[i],G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i][g]-max_v[i]); 
		for(g=0;g<G;g++) z[i][g]=exp(v[i][g]-max_v[i])/d_alt;
    }/* End i loop */  
	return 0;
}

int update_z11(double **v, double **x, double **z, double **lambda, double omega, double **delta, double **mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g,j;
    double a,e,d_alt, *psi;

	my_alloc_vec(&psi,p);
    
    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
			for(j=0;j<p;j++) psi[j]=omega*delta[g][j];
			a = woodbury2(x[i], lambda, psi, mu[g], p, q);
			e = a/2.0*(-1.0);
			v[i][g] = e + log(pi[g]) - log_c[g];
        }/* End g loop */
		max_v[i] = maximum_array(v[i],G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i][g]-max_v[i]); 
		for(g=0;g<G;g++) z[i][g]=exp(v[i][g]-max_v[i])/d_alt;
    }/* End i loop */  
	return 0;
}

int update_z12(double **v, double **x, double **z, double ***lambda, double omega, double **delta, double **mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g,j;
    double a,e,d_alt, *psi;

	my_alloc_vec(&psi,p);
    
    for(i=0; i<N; i++){
        for(g=0; g<G; g++){
			for(j=0;j<p;j++) psi[j]=omega*delta[g][j];
			a = woodbury2(x[i], lambda[g], psi, mu[g], p, q);
			e = a/2.0*(-1.0);
			v[i][g] = e + log(pi[g]) - log_c[g];
        }/* End g loop */
		max_v[i] = maximum_array(v[i],G);
		d_alt=0.0;
		for(g=0;g<G;g++) d_alt += exp(v[i][g]-max_v[i]); 
		for(g=0;g<G;g++) z[i][g]=exp(v[i][g]-max_v[i])/d_alt;
    }/* End i loop */  
	return 0;
}

double update_omega(double **lambda, double *delta, double **beta, double **sampcovg, double **theta, int p, int q){ 
	
	int i;
	double **temp, **result_1, *result_2, *result_3; 
	double omega;

	my_alloc(&temp,q,p);
	my_alloc(&result_1,p,p);
	my_alloc_vec(&result_2,p);
	my_alloc_vec(&result_3,p);

	/* lambda*beta*sampcov */
    mx_mult(p,q,p,lambda, beta, result_1);
    mx_mult_diag1(p,p,result_1, sampcovg, result_2);
       
	/* lambda*theta*lambda' */
    mx_trans(p,q,lambda,temp);
    mx_mult(p,q,q,lambda,theta,result_1); 
    mx_mult_diag1(p,q,result_1,temp,result_3);
        
    omega = 0.0;
    for(i=0;i<p;i++)
        omega += (sampcovg[i][i]-2*result_2[i]+result_3[i])/delta[i];
    omega /= p;
	
	my_free(temp,q); my_free(result_1,p); free(result_2);free(result_3);
	
	return omega;
}

double update_omega2(double **lambda, double *delta, double **beta, double **sampcovg, int p, int q){ 
	
	int i;
	double **result_1, *result_2; 
	double omega;

	my_alloc(&result_1,p,p);
	my_alloc_vec(&result_2,p);

	/* lambda*beta*sampcov */
    mx_mult(p,q,p,lambda, beta, result_1);
    mx_mult_diag1(p,p,result_1, sampcovg, result_2);
       
    omega = 0.0;
    for(i=0;i<p;i++)
        omega += (sampcovg[i][i]-result_2[i])/delta[i];
    omega /= p;
	
	my_free(result_1,p); free(result_2); 
	
	return omega;
}


void update_delta(double *delta, double **lambda, double *omega, double ***beta, double ***sampcov, double ***theta, 
		double *n, int p, int q, int N, int G){

	int i,g;
	double **temp, **result_1, **result_2, **result_3, *temp1;
	double lagrange = 0.0; /* Lagrange Multiplier */
	my_alloc(&temp,q,p);
	my_alloc(&result_1,p,p);
	my_alloc(&result_2,G,p);
	my_alloc(&result_3,G,p);
	my_alloc_vec(&temp1,p);

	/* lambda*beta_g*sampcov_g */
    for(g=0;g<G;g++){
		mx_mult(p,q,p,lambda, beta[g], result_1);
       	mx_mult_diag1(p,p,result_1, sampcov[g], result_2[g]);
    }
	
	/* lambda*theta_g*lambda' */
    for(g=0;g<G;g++){
		mx_trans(p,q,lambda,temp);
       	mx_mult(p,q,q,lambda,theta[g],result_1); 
       	mx_mult_diag1(p,q,result_1,temp,result_3[g]);
	}
	
    for(i=0;i<p;i++){
        temp1[i] = 0.0;
	    for(g=0;g<G;g++)
		    temp1[i] += (sampcov[g][i][i]-2.0*result_2[g][i]+result_3[g][i])*n[g]/omega[g];
	    lagrange += log(temp1[i]);
	}	
	
	/* Compute Lagrange Multiplier */
	lagrange /= (double)p;
	lagrange = exp(lagrange);
	lagrange -= (double)N; 
	lagrange /= (double)2; 

	for(i=0;i<p;i++) delta[i] = temp1[i]/((double)N+(2.0*lagrange));
		
	my_free(temp,q); my_free(result_1,p); my_free(result_2,G); 
	my_free(result_3,G); free(temp1);
}

void update_delta2(double *delta, double ***lambda, double *omega, double ***beta, double ***sampcov, double ***theta, 
		double *n, int p, int q, int N, int G){

	int i,g;
	double **temp, **result_1, **result_2, **result_3, *temp1, toprint, lagrange_alt=0.0;
	my_alloc(&temp,q,p);
	my_alloc(&result_1,p,p);
	my_alloc(&result_2,G,p);
	my_alloc(&result_3,G,p);
	my_alloc_vec(&temp1,p);

	/* lambda*beta_g*sampcov_g */
    for(g=0;g<G;g++){
		mx_mult(p,q,p,lambda[g], beta[g], result_1);
       	mx_mult_diag1(p,p,result_1, sampcov[g], result_2[g]);
    }
	
	/* lambda*theta_g*lambda' */
    for(g=0;g<G;g++){
		mx_trans(p,q,lambda[g],temp);
        mx_mult(p,q,q,lambda[g],theta[g],result_1); 
        mx_mult_diag1(p,q,result_1,temp,result_3[g]);
	}
	
    for(i=0;i<p;i++){
        temp1[i] = 0.0;
	    for(g=0;g<G;g++)
		temp1[i] += ((sampcov[g][i][i]-2*result_2[g][i]+result_3[g][i])*n[g]/omega[g]);
	    lagrange_alt += log(temp1[i]);
	}	
	
	lagrange_alt /= p;
	toprint = exp(lagrange_alt);
	lagrange_alt = toprint;
	lagrange_alt -= N; 
	lagrange_alt /= 2; 
	
	for(i=0;i<p;i++)
	    delta[i] = temp1[i]/(N+(2*lagrange_alt)); 
	
	my_free(temp,q); my_free(result_1,p); my_free(result_2,G); 
	my_free(result_3,G); free(temp1);
}

/* Omega not a vector */
void update_delta3(double *delta, double **lambda, double omega, double **beta, double **sampcov, double **theta, 
		double n, int p, int q){

	int i;
	double **temp, **result_1, *result_2, *result_3, *temp1;
	double lagrange = 0.0; /* Lagrange Multiplier THIS WAS CHANGED FROM 1.0 NOV 2008*/
	my_alloc(&temp,p,p);
	my_alloc(&result_1,p,p);
	my_alloc_vec(&result_2,p);
	my_alloc_vec(&result_3,p);
	my_alloc_vec(&temp1,p);

	/* lambda*beta_g*sampcov_g */
    mx_mult(p,q,p,lambda, beta, result_1);
    mx_mult_diag1(p,p,result_1, sampcov, result_2);
       		
	/* lambda*theta_g*lambda' */
    mx_trans(p,q,lambda,temp);
    mx_mult(p,q,q,lambda,theta,result_1); 
    mx_mult_diag1(p,q,result_1,temp,result_3);
		
    for(i=0;i<p;i++) temp1[i] = sampcov[i][i]-2*result_2[i]+result_3[i];
	
	for(i=0;i<p;i++) lagrange += log(temp1[i]);
		
	/* Compute Lagrange Multiplier */
	lagrange /= p; 
	lagrange = exp(lagrange);
	lagrange /= omega;
	lagrange -= 1; 
	lagrange *= (n/2); 

	for(i=0;i<p;i++) delta[i] = temp1[i]/((1+(2*lagrange/n))*omega);
		
	my_free(temp,p); my_free(result_1,p); free(result_2); 
	free(result_3); free(temp1);
}
