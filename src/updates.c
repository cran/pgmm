#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "functions.h"

void p5M_n(double *n, double **z, int G, int N){    
	int g,i;	    
	for(g=0; g<G; g++){
            n[g]=0.0;
            for(i=0; i<N; i++)
                n[g] += z[i][g];
    }
}
    
void p5M_pi(double *pi, double *n, int G, int N){    
	int g;	    
	for(g=0; g<G; g++) pi[g] = n[g]/N; 
}
    
void updau8m(double **mu, double *n, double **x, double **z,int G, 
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

void p5M_stilde(double **y2etilde, double **x, double **z, double **mu,
	int G, int N, int p){
        int i,j,k,g;
	for(j=0; j<p; j++){
            for(k=0; k<=j; k++){
                y2etilde[j][k]=0.0;
                for(g=0;g<G;g++)    
                    for(i=0; i<N; i++)
                        y2etilde[j][k]+=z[i][g]*(x[i][j]-mu[g][j])*(x[i][k]-mu[g][k]);
				y2etilde[j][k]/=N;
                y2etilde[k][j] = y2etilde[j][k];		
            }
    }
}

void p5M_sg(double ***sg, double **x, double **z, double **mu,
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

void p5M_kf41(double **kf4, double bx1, double **zR6, int p, int q){
	int i,j;
	double **lhs, **rhs, **cp, **result;
	double det[1];
	my_alloc(&lhs,p,p);
	my_alloc(&rhs,p,p);
	my_alloc(&cp,p,p);
	my_alloc(&result,p,p);

 	/* zR6'/bx1 */
    mx_trans(p, q, zR6, lhs);
    for(i=0;i<q;i++)
		for(j=0;j<p;j++)
            lhs[i][j]/=bx1; /*bx1[j]*/
 	/* zR6'/bx1* zR6 */
    mx_mult(q, p, q, lhs, zR6, cp);
        
    /* (I + zR6'/bx1*zR6)^{-1} */
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
			kf4[i][j] = lhs[i][j] - rhs[i][j];       
 
	my_free(lhs,p); my_free(result,p); my_free(rhs,p); my_free(cp,p);
}

void p5M_kf42(double **kf4, double *bx1, double **zR6, int p, int q){
	int i,j;
	double **lhs, **rhs, **cp, **result;
	double det[1];
	my_alloc(&lhs,p,p);
	my_alloc(&rhs,p,p);
	my_alloc(&cp,p,p);
	my_alloc(&result,p,p);

 	/* zR6'/bx1 */
    mx_trans(p, q, zR6, lhs);
    for(i=0;i<q;i++)
		for(j=0;j<p;j++)
            lhs[i][j]/=bx1[j];
 	/* zR6'/bx1* zR6 */
    mx_mult(q, p, q, lhs, zR6, cp);
        
    /* (I + zR6'/bx1*zR6)^{-1} */
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
			kf4[i][j] = lhs[i][j] - rhs[i][j];       
 
	my_free(lhs,p); my_free(result,p); my_free(rhs,p); my_free(cp,p);
}

void p5M_po1(double **po1, double **kf4, double **zR6, double **y2etilde, int p, int q){ 
	int i,j;
	double **id_q, **r_1, **r_2, **r_3, **tmp;
	
	my_alloc(&tmp,p,p);
	my_alloc(&r_1,q,q);
	my_alloc(&r_2,q,p);
	my_alloc(&r_3,q,q);
	my_alloc(&id_q,q,q);
	generate_identity(q, id_q);	
	
	/* Work out kf4*zR6 */
    mx_mult(q, p, q, kf4, zR6, r_1);
    
    /* Now, work out kf4*y2etilde*kf4' */
    mx_mult(q, p, p, kf4, y2etilde, r_2);
    mx_trans(q, p, kf4, tmp);    
    mx_mult(q, p, q, r_2, tmp, r_3);    
                 
    for(i=0;i<q;i++)
		for(j=0;j<q;j++)
			po1[i][j] = id_q[i][j]-r_1[i][j]+r_3[i][j];
       
	my_free(id_q,q); my_free(tmp,p); my_free(r_1,q); my_free(r_2,q); 
	my_free(r_3,q);
}

void p5M_zR6(double **zR6, double **kf4, double **s, double **po1, int p, int q){
	
	int i,j;
	double **tran, **res1, **res2, **res3, det[1];
	
	my_alloc(&tran,p,q);
	my_alloc(&res1,p,q);
	my_alloc(&res2,q,q);
	my_alloc(&res3,q,q);
	
	/* y2etilde*kf4'*/
    mx_trans(q,p,kf4,tran);
    mx_mult(p,p,q,s,tran,res1);
    	
	/* Make of copy of po1 ahead of Gauss-Jordan */
	for(i=0;i<q;i++)
        for(j=0;j<q;j++)
            res3[i][j] = po1[i][j];   
        
	/* zR6=y2etilde*kf4'*po1^{-1} */
    GaussJordan(q, res3, res2,det);
    mx_mult(p,q,q,res1,res2,zR6);
    
	my_free(tran,p); my_free(res1,p); my_free(res2,q); my_free(res3,q);
}

/* For CUU */
void p5M_zR6_cuu(double **zR6, double ***kf4, double ***s, 
	double ***po1, double *n, double **bx1, int p, int q, int G){
	
	int i,j,k,g;
	double **tran, **res1, **res2, **res3, det[1];
	
	my_alloc(&tran,p,q);
	my_alloc(&res1,p,q);
	my_alloc(&res2,p,q);
	my_alloc(&res3,q,q);	
	
    /* Compute the RHS --- only needs to happen once */
	/* sum_g[y2e_g*kf4_g'*n_g/bx1_g] */
    for(g=0;g<G;g++){
		mx_trans(q,p,kf4[g],tran);
		mx_mult(p,p,q,s[g],tran,res1);
		if(g==0){ 
			for(i=0;i<p;i++)
				for(j=0;j<q;j++)
					res2[i][j] = res1[i][j]*n[g]/bx1[g][i];
		}else{
			for(i=0;i<p;i++)
				for(j=0;j<q;j++)
					res2[i][j] += res1[i][j]*n[g]/bx1[g][i];    
		}
    }

    /* Now solve for zR6 row-by-row */
	for (i=0;i<p;i++){
		/* First, compute the po1 sum*/
		for(g=0;g<G;g++){
			if(g==0){ 
				for(k=0;k<q;k++)
					for(j=0;j<q;j++)
						res3[k][j] = po1[g][k][j]*n[g]/bx1[g][i];
			}else{
				for(k=0;k<q;k++)
					for(j=0;j<q;j++)
						res3[k][j] += po1[g][k][j]*n[g]/bx1[g][i];    
			}
		}
		/* Invert po1 sum */
		GaussJordan(q, res3, res1, det);
		
		/* Now solve for row i of zR6 */
		vec_mx_mult(q, q, res2[i], res1, zR6[i]);
	}     	
	my_free(tran,p); my_free(res1,p); my_free(res2,p); my_free(res3,q);
}

/* For CUC */
void p5M_zR62(double **zR6, double ***kf4, double ***s, double ***po1, double *n, double *bx1, int p, 
		int q, int G){
	
	int i,j,g;
	double **tran, **res1, **res2, **res3, det[1];
	
	my_alloc(&tran,p,q);
	my_alloc(&res1,p,q);
	my_alloc(&res2,p,q);
	my_alloc(&res3,q,q);
	
    /* sum_g[y2e_g*kf4_g'*n_g/bx1_g] */
    for(g=0;g<G;g++){
		mx_trans(q,p,kf4[g],tran);
		mx_mult(p,p,q,s[g],tran,res1);
		if(g==0){ 
			for(i=0;i<p;i++)
				for(j=0;j<q;j++)
					res2[i][j] = res1[i][j]*n[g]/bx1[g];
		}else{
			for(i=0;i<p;i++)
				for(j=0;j<q;j++)
					res2[i][j] += res1[i][j]*n[g]/bx1[g];    
		}
    }

   	/* sum_g[po1_g'*n_g/bx1_g]^{-1} */
    for(g=0;g<G;g++){
		if(g==0){ 
			for(i=0;i<q;i++)
				for(j=0;j<q;j++)
					res3[i][j] = po1[g][i][j]*n[g]/bx1[g];
		}else{
			for(i=0;i<q;i++)
				for(j=0;j<q;j++)
					res3[i][j] += po1[g][i][j]*n[g]/bx1[g];    
		}
    }
        
    GaussJordan(q, res3, res1,det); /* inverting the po1 sum */
    mx_mult(p,q,q,res2,res1,zR6); /* this gives zR6 */ 
	
	my_free(tran,p); my_free(res1,p); my_free(res2,p); my_free(res3,q);
}

double p5M_bx1(double **zR6, double **kf4, double **y2etilde, int p, int q){ 
	
	int i;
	double **result_1, *result_2; 
	double bx1;
	my_alloc(&result_1,p,p);
	my_alloc_vec(&result_2,p);

	/* zR6*kf4*y2etilde */
    mx_mult(p,q,p,zR6, kf4, result_1);
    mx_mult_diag1(p,p,result_1, y2etilde, result_2);
       
    bx1 = 0.0;
    for(i=0;i<p;i++)
        bx1 += y2etilde[i][i]-result_2[i];
    bx1 /= p;
	
	my_free(result_1,p); free(result_2); 
	
	return bx1;
}

/* Account for bx1_p*/
void p5M_bx12(double *bx1, double **zR6, double **kf4, double **y2etilde, int p, int q){ 
	
	int i;
	double **result_1, *result_2;
	my_alloc(&result_1,p,p);
	my_alloc_vec(&result_2,p);

	/* zR6*kf4*y2etilde */
	mx_mult(p,q,p,zR6, kf4, result_1);
    mx_mult_diag1(p,p,result_1, y2etilde, result_2);
        
    for(i=0;i<p;i++)
		bx1[i] = y2etilde[i][i]-result_2[i];
	
	my_free(result_1,p); free(result_2); 
}

/* CUC case */
double p5M_bx13(double **zR6, double **kf4, double **y2eg, double **po1, int p, int q){ 
	
	int i;
	double **temp, **result_1, *result_2, *result_3; 
	double bx1;
	my_alloc(&temp,q,p);
	my_alloc(&result_1,p,p);
	my_alloc_vec(&result_2,p);
	my_alloc_vec(&result_3,p);

	/* zR6*kf4*y2e */
    mx_mult(p,q,p,zR6, kf4, result_1);
    mx_mult_diag1(p,p,result_1, y2eg, result_2);
       
	/* zR6*po1*zR6' */
    mx_trans(p,q,zR6,temp);
    mx_mult(p,q,q,zR6,po1,result_1); 
    mx_mult_diag1(p,q,result_1,temp,result_3);
        
    bx1 = 0.0;
    for(i=0;i<p;i++)
        bx1 += y2eg[i][i]-2*result_2[i]+result_3[i];
    bx1 /= p;
	
	my_free(temp,q); my_free(result_1,p); free(result_2); 
	free(result_3);
	
	return bx1;
}

/* CUU case */
void p5M_bx1_cuu(double **bx1, double **zR6, double ***kf4, double ***y2eg, double 
	***po1, int p, int q, int G){ 
	
	int i, g;
	double **temp, **result_1, **result_2, **result_3;
	my_alloc(&temp,q,p);
	my_alloc(&result_1,p,p);
	my_alloc(&result_2,G,p);
	my_alloc(&result_3,G,p);

	/* zR6*kf4*y2e */
    for (g=0;g<G;g++){
		mx_mult(p,q,p,zR6, kf4[g], result_1);
    	mx_mult_diag1(p,p,result_1, y2eg[g], result_2[g]);
	}
       
	/* zR6*po1*zR6' */
    for(g=0;g<G;g++){
		mx_trans(p,q,zR6,temp);
        mx_mult(p,q,q,zR6,po1[g],result_1); 
        mx_mult_diag1(p,q,result_1,temp,result_3[g]);
	}
        
    for(g=0;g<G;g++)
	    for(i=0;i<p;i++)
             	bx1[g][i] = y2eg[g][i][i]-2*result_2[g][i]+result_3[g][i];
	
	my_free(temp,q); my_free(result_1,p); my_free(result_2,G); 
	my_free(result_3,G);
}

/* UCU*/
void p5M_bx1_ucu(double *bx1, double ***zR6, double ***kf4, double ***y2e, int p, int q, double *pi, int G){ 
	
	int i,g;
	double **result_1, **result_2;
	my_alloc(&result_1,p,p);
	my_alloc(&result_2,G,p);

	/* zR6*kf4*y2etilde */
    for (g=0; g<G; g++){
	    mx_mult(p,q,p,zR6[g], kf4[g], result_1);
        mx_mult_diag1(p,p,result_1, y2e[g], result_2[g]);
	}
	
    for(i=0;i<p;i++){
		bx1[i]=0.0;
		for(g=0;g<G;g++) bx1[i] += pi[g]*(y2e[g][i][i]-result_2[g][i]);
	}
	
	my_free(result_1,p); my_free(result_2,G); 
}

/* UCC*/
double p5M_bx1_ucc(double ***zR6, double ***kf4, double ***y2e, int p, int q, double *pi, int G){ 
	
	int i,g;
	double **result_1, **result_2, bx1;
	my_alloc(&result_1,p,p);
	my_alloc(&result_2,G,p);

	/* zR6*kf4*y2etilde */
	for (g=0; g<G; g++){
		mx_mult(p,q,p,zR6[g], kf4[g], result_1);
        mx_mult_diag1(p,p,result_1, y2e[g], result_2[g]);
	}
	
    bx1=0.0;
	for(g=0;g<G;g++)
	    for(i=0;i<p;i++)
			bx1 += pi[g]*(y2e[g][i][i]-result_2[g][i]);
	bx1 /= p;
	
	my_free(result_1,p); my_free(result_2,G);

	return bx1;
}

double p5M_det_sigma_NEW(double **zR6, double bx1, double ja2, int p, int q){

	int i, j;
	double **tmp, **tmp2, det[1];

	my_alloc(&tmp,p,p);
	my_alloc(&tmp2,p,p);

	p5M_kf41(tmp2, bx1, zR6, p, q);
	mx_mult(q,p,q,tmp2, zR6, tmp);
    for(i=0;i<q;i++){
		for(j=0;j<q;j++){
			tmp2[i][j] = tmp[i][j]*(-1.0);
			if(i==j) tmp2[i][i] += 1.0;
		}
	}
	GaussJordan(q,tmp2,tmp,det);

	my_free(tmp,p); my_free(tmp2,p);

	return (ja2 - log(det[0]));
}

double p5M_det_sigma_NEW2(double **zR6, double *bx1, double ja2, int p, int q){

	int i, j;
	double **tmp, **tmp2, det[1];

	my_alloc(&tmp,p,p);
	my_alloc(&tmp2,p,p);

	p5M_kf42(tmp2, bx1, zR6, p, q);
	mx_mult(q,p,q,tmp2, zR6, tmp);
    for(i=0;i<q;i++){
		for(j=0;j<q;j++){
			tmp2[i][j] = tmp[i][j]*(-1.0);
			if(i==j) tmp2[i][i] += 1.0;
		}
	}
	GaussJordan(q,tmp2,tmp,det);

	my_free(tmp,p); my_free(tmp2,p);

	return (ja2 - log(det[0]));
}

double p5M_det_sigma(double **zR6, double **sigma, double ja2, int p, int q){

	int i,j;
	double siginv, det[1];
	double **r_1, **r_2, **r_3, **tmp;
	
	my_alloc(&tmp,p,p);
	my_alloc(&r_1,p,p);
	my_alloc(&r_2,p,p);
	my_alloc(&r_3,p,p);
	
	/* |I-zR6'*sigma*zR6| */
    mx_trans(p, q, zR6, tmp);
    mx_mult(q, p, p, tmp, sigma, r_1);
    mx_mult(q,p,q,r_1,zR6,r_3);
    for(i=0;i<q;i++){
       	for(j=0;j<q;j++){
        	r_3[i][j] = 0-r_3[i][j];
            if(i==j) r_3[i][i] += 1.0;
        }
    }
    GaussJordan(q,r_3,r_2,det);
    siginv = ja2 - log(det[0]);
	
	my_free(tmp,p); my_free(r_1,p); my_free(r_2,p); my_free (r_3,p);
	/*if (siginv<0) printf("-");*/
	return siginv;
}

/***************** p5M Z values WITHOUT Sigma ************************/
int p5M_z(double **v, double **x, double **z, double **zR6, double bx1, double **mu, 
	      double *pi, double *max_v, double log_c, int N, int G, int p, int q){
    
    int i,g;
    double d,a,e,d_alt;        
   
    for(i=0; i<N; i++){
        d=0.0;
        for(g=0; g<G; g++){
			a = woodbury(x[i], zR6, bx1, mu[g], p, q);
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
double woodbury(double *x, double **zR6, double bx1, double *mu, int p, int q){
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
	lhs /= bx1; 

    /* Compute RHS */
	for (j=0;j<p;j++)
		lvec[j] = (x[j]-mu[j])/bx1;
	vec_mx_mult(p,q,lvec,zR6,tvec);
 	/* zR6'/bx1 */
    mx_trans(p, q, zR6, temp);
    for(i=0;i<q;i++)
		for(j=0;j<p;j++)
            temp[i][j]/=bx1;
 	/* zR6'/bx1* zR6 */
    mx_mult(q, p, q, temp, zR6, result);
    /* (I + zR6'/bx1*zR6)^{-1} */
    for(i=0; i<q;i++) result[i][i] += 1.0;
    GaussJordan(q, result, cp, det);
    mx_trans(p, q, zR6, temp);
    mx_mult(q, q, p, cp, temp, result);
	vec_mx_mult(q,p,tvec,result,cvec);
	rhs = 0.0;
	for (j=0;j<p;j++)
		rhs += cvec[j]*(x[j]-mu[j]);
	rhs /= bx1; 
	
	free(lvec); free(cvec); free(tvec); my_free(result,q); my_free(temp,q); my_free(cp,q);

	return (lhs-rhs);
}

/***************** p5M Z values WITHOUT Sigma ************************/
int p5M_z2(double **v, double **x, double **z, double **zR6, double *bx1, double **mu, 
	      double *pi, double *max_v, double log_c, int N, int G, int p, int q){
    
    int i,g;
    double d,a,e,d_alt;        
    
    for(i=0; i<N; i++){
        d=0.0;
        for(g=0; g<G; g++){
			a = woodbury2(x[i], zR6, bx1, mu[g], p, q);
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
double woodbury2(double *x, double **zR6, double *bx1, double *mu, int p, int q){
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
		lhs += (x[j]-mu[j])*(x[j]-mu[j])/bx1[j];

    /* Compute RHS */
	for (j=0;j<p;j++)
		lvec[j] = (x[j]-mu[j])/bx1[j];
	vec_mx_mult(p,q,lvec,zR6,tvec);
 	/* zR6'/bx1 */
    mx_trans(p, q, zR6, temp);
    for(i=0;i<q;i++)
		for(j=0;j<p;j++)
            temp[i][j]/=bx1[j];
 	/* zR6'/bx1* zR6 */
    mx_mult(q, p, q, temp, zR6, result);
    /* (I + zR6'/bx1*zR6)^{-1} */
    for(i=0; i<q;i++) result[i][i] += 1.0;
    GaussJordan(q, result, cp, det);
    mx_trans(p, q, zR6, temp);
    mx_mult(q, q, p, cp, temp, result);
	vec_mx_mult(q,p,tvec,result,cvec);
	rhs = 0.0;
	for (j=0;j<p;j++)
		rhs += cvec[j]*(x[j]-mu[j])/bx1[j];
	
	free(lvec); free(cvec); free(tvec); my_free(result,q); my_free(temp,q); my_free(cp,q);

	return (lhs-rhs);
}

/***************** p5M Z values WITHOUT Sigma ************************/
int p5M_z3(double **v, double **x, double **z, double **zR6, double *bx1, double **mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g;
    double d,a,e,d_alt;        
    
    for(i=0; i<N; i++){
        d=0.0;
        for(g=0; g<G; g++){
			a = woodbury(x[i], zR6, bx1[g], mu[g], p, q);
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

/***************** p5M Z values WITHOUT Sigma ************************/
int p5M_z4(double **v, double **x, double **z, double **zR6, double **bx1, double **mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g;
    double d,a,e,d_alt;        
    
    for(i=0; i<N; i++){
        d=0.0;
        for(g=0; g<G; g++){
			a = woodbury2(x[i], zR6, bx1[g], mu[g], p, q);
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

int p5M_z5(double **v, double **x, double **z, double ***zR6, double bx1, double **mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g;
    double d,a,e,d_alt;        
   
    for(i=0; i<N; i++){
        d=0.0;
        for(g=0; g<G; g++){
			a = woodbury(x[i], zR6[g], bx1, mu[g], p, q);
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

int p5M_z6(double **v, double **x, double **z, double ***zR6, double *bx1, double **mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g;
    double d,a,e,d_alt;        
    
    for(i=0; i<N; i++){
        d=0.0;
        for(g=0; g<G; g++){
			a = woodbury2(x[i], zR6[g], bx1, mu[g], p, q);
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

int p5M_z7(double **v, double **x, double **z, double ***zR6, double *bx1, double **mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g;
    double d,a,e,d_alt;        
    
    for(i=0; i<N; i++){
        d=0.0;
        for(g=0; g<G; g++){
			a = woodbury(x[i], zR6[g], bx1[g], mu[g], p, q);
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

int p5M_z8(double **v, double **x, double **z, double ***zR6, double **bx1, double **mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g;
    double d,a,e,d_alt;        
    
    for(i=0; i<N; i++){
        d=0.0;
        for(g=0; g<G; g++){
			a = woodbury2(x[i], zR6[g], bx1[g], mu[g], p, q);
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

int p5M_z9(double **v, double **x, double **z, double **zR6, double *dsw, double *g9p, double **mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g,j;
    double d,a,e,d_alt, *bx1;

	my_alloc_vec(&bx1,p);
    
    for(i=0; i<N; i++){
        d=0.0;
        for(g=0; g<G; g++){
			for(j=0;j<p;j++) bx1[j]=dsw[g]*g9p[j];
			a = woodbury2(x[i], zR6, bx1, mu[g], p, q);
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

int p5M_z10(double **v, double **x, double **z, double ***zR6, double *dsw, double *g9p, double **mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g,j;
    double d,a,e,d_alt, *bx1;

	my_alloc_vec(&bx1,p);
    
    for(i=0; i<N; i++){
        d=0.0;
        for(g=0; g<G; g++){
			for(j=0;j<p;j++) bx1[j]=dsw[g]*g9p[j];
			a = woodbury2(x[i], zR6[g], bx1, mu[g], p, q);
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

int p5M_z11(double **v, double **x, double **z, double **zR6, double dsw, double **g9p, double **mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g,j;
    double d,a,e,d_alt, *bx1;

	my_alloc_vec(&bx1,p);
    
    for(i=0; i<N; i++){
        d=0.0;
        for(g=0; g<G; g++){
			for(j=0;j<p;j++) bx1[j]=dsw*g9p[g][j];
			a = woodbury2(x[i], zR6, bx1, mu[g], p, q);
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

int p5M_z12(double **v, double **x, double **z, double ***zR6, double dsw, double **g9p, double **mu, 
	      double *pi, double *max_v, double *log_c, int N, int G, int p, int q){
    
    int i,g,j;
    double d,a,e,d_alt, *bx1;

	my_alloc_vec(&bx1,p);
    
    for(i=0; i<N; i++){
        d=0.0;
        for(g=0; g<G; g++){
			for(j=0;j<p;j++) bx1[j]=dsw*g9p[g][j];
			a = woodbury2(x[i], zR6[g], bx1, mu[g], p, q);
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

double p5M_dsw(double **zR6, double *g9p, double **kf4, double **y2eg, double **po1, int p, int q){ 
	
	int i;
	double **temp, **result_1, *result_2, *result_3; 
	double dsw;

	my_alloc(&temp,q,p);
	my_alloc(&result_1,p,p);
	my_alloc_vec(&result_2,p);
	my_alloc_vec(&result_3,p);

	/* zR6*kf4*y2e */
    mx_mult(p,q,p,zR6, kf4, result_1);
    mx_mult_diag1(p,p,result_1, y2eg, result_2);
       
	/* zR6*po1*zR6' */
    mx_trans(p,q,zR6,temp);
    mx_mult(p,q,q,zR6,po1,result_1); 
    mx_mult_diag1(p,q,result_1,temp,result_3);
        
    dsw = 0.0;
    for(i=0;i<p;i++)
        dsw += (y2eg[i][i]-2*result_2[i]+result_3[i])/g9p[i];
    dsw /= p;
	
	my_free(temp,q); my_free(result_1,p); free(result_2);free(result_3);
	
	return dsw;
}

double p5M_dsw2(double **zR6, double *g9p, double **kf4, double **y2eg, int p, int q){ 
	
	int i;
	double **result_1, *result_2; 
	double dsw;

	my_alloc(&result_1,p,p);
	my_alloc_vec(&result_2,p);

	/* zR6*kf4*y2e */
    mx_mult(p,q,p,zR6, kf4, result_1);
    mx_mult_diag1(p,p,result_1, y2eg, result_2);
       
    dsw = 0.0;
    for(i=0;i<p;i++)
        dsw += (y2eg[i][i]-result_2[i])/g9p[i];
    dsw /= p;
	
	my_free(result_1,p); free(result_2); 
	
	return dsw;
}


void p5M_g9p(double *g9p, double **zR6, double *dsw, double ***kf4, double ***y2e, double ***po1, 
		double *n, int p, int q, int N, int G){

	int i,g;
	double **temp, **result_1, **result_2, **result_3, *temp1;
	double lagrange = 0.0; /* Lagrange Multiplier */
	my_alloc(&temp,q,p);
	my_alloc(&result_1,p,p);
	my_alloc(&result_2,G,p);
	my_alloc(&result_3,G,p);
	my_alloc_vec(&temp1,p);

	/* zR6*kf4_g*y2e_g */
    for(g=0;g<G;g++){
		mx_mult(p,q,p,zR6, kf4[g], result_1);
       	mx_mult_diag1(p,p,result_1, y2e[g], result_2[g]);
    }
	
	/* zR6*po1_g*zR6' */
    for(g=0;g<G;g++){
		mx_trans(p,q,zR6,temp);
       	mx_mult(p,q,q,zR6,po1[g],result_1); 
       	mx_mult_diag1(p,q,result_1,temp,result_3[g]);
	}
	
    for(i=0;i<p;i++){
        temp1[i] = 0.0;
	    for(g=0;g<G;g++)
		    temp1[i] += (y2e[g][i][i]-2.0*result_2[g][i]+result_3[g][i])*n[g]/dsw[g];
	    lagrange += log(temp1[i]);
	}	
	
	/* Compute Lagrange Multiplier */
	lagrange /= (double)p;
	lagrange = exp(lagrange);
	lagrange -= (double)N; 
	lagrange /= (double)2; 

	for(i=0;i<p;i++) g9p[i] = temp1[i]/((double)N+(2.0*lagrange));
		
	my_free(temp,q); my_free(result_1,p); my_free(result_2,G); 
	my_free(result_3,G); free(temp1);
}

void p5M_g9p2(double *g9p, double ***zR6, double *dsw, double ***kf4, double ***y2e, double ***po1, 
		double *n, int p, int q, int N, int G){

	int i,g;
	double **temp, **result_1, **result_2, **result_3, *temp1, toprint, lagrange_alt=0.0;
	my_alloc(&temp,q,p);
	my_alloc(&result_1,p,p);
	my_alloc(&result_2,G,p);
	my_alloc(&result_3,G,p);
	my_alloc_vec(&temp1,p);

	/* zR6*kf4_g*y2e_g */
    for(g=0;g<G;g++){
		mx_mult(p,q,p,zR6[g], kf4[g], result_1);
       	mx_mult_diag1(p,p,result_1, y2e[g], result_2[g]);
    }
	
	/* zR6*po1_g*zR6' */
    for(g=0;g<G;g++){
		mx_trans(p,q,zR6[g],temp);
        mx_mult(p,q,q,zR6[g],po1[g],result_1); 
        mx_mult_diag1(p,q,result_1,temp,result_3[g]);
	}
	
    for(i=0;i<p;i++){
        temp1[i] = 0.0;
	    for(g=0;g<G;g++)
		temp1[i] += ((y2e[g][i][i]-2*result_2[g][i]+result_3[g][i])*n[g]/dsw[g]);
	    lagrange_alt += log(temp1[i]);
	}	
	
	lagrange_alt /= p;
	toprint = exp(lagrange_alt);
	lagrange_alt = toprint;
	lagrange_alt -= N; 
	lagrange_alt /= 2; 
	
	for(i=0;i<p;i++)
	    g9p[i] = temp1[i]/(N+(2*lagrange_alt)); 
	
	my_free(temp,q); my_free(result_1,p); my_free(result_2,G); 
	my_free(result_3,G); free(temp1);
}

/* dsw not a vector */
void p5M_g9p3(double *g9p, double **zR6, double dsw, double **kf4, double **y2e, double **po1, 
		double n, int p, int q){

	int i;
	double **temp, **result_1, *result_2, *result_3, *temp1;
	double lagrange = 0.0; /* Lagrange Multiplier THIS WAS CHANGED FROM 1.0 NOV 2008*/
	my_alloc(&temp,p,p);
	my_alloc(&result_1,p,p);
	my_alloc_vec(&result_2,p);
	my_alloc_vec(&result_3,p);
	my_alloc_vec(&temp1,p);

	/* zR6*kf4_g*y2e_g */
    mx_mult(p,q,p,zR6, kf4, result_1);
    mx_mult_diag1(p,p,result_1, y2e, result_2);
       		
	/* zR6*po1_g*zR6' */
    mx_trans(p,q,zR6,temp);
    mx_mult(p,q,q,zR6,po1,result_1); 
    mx_mult_diag1(p,q,result_1,temp,result_3);
		
    for(i=0;i<p;i++) temp1[i] = y2e[i][i]-2*result_2[i]+result_3[i];
	
	for(i=0;i<p;i++) lagrange += log(temp1[i]);
		
	/* Compute Lagrange Multiplier */
	lagrange /= p; 
	lagrange = exp(lagrange);
	lagrange /= dsw;
	lagrange -= 1; 
	lagrange *= (n/2); 

	for(i=0;i<p;i++) g9p[i] = temp1[i]/((1+(2*lagrange/n))*dsw);
		
	my_free(temp,p); my_free(result_1,p); free(result_2); 
	free(result_3); free(temp1);
}
