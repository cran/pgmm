#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "functions.h"
#include "time.h"
funcType funcs[13] = {NULL, aecm,aecm2,aecm3,aecm4,aecm5,aecm6,aecm7,aecm8,aecm9,aecm10,aecm11,aecm12};
funcType2 funcs2[13] = {NULL, claecm,claecm2,claecm3,claecm4,claecm5,claecm6,claecm7,claecm8,claecm9,claecm10,claecm11,claecm12};

void pgmm_c(double *x1, double *z1, double *bic, int *cls, int *q, int *p, int *G, int *N, int *model, int *clust, double *lambda, double *psi, double *tol){
    funcType func;
    funcType2 func2;
    int NN, pp, GG;
    NN = *N;
    pp = *p;
    GG = *G;
    double *x = malloc(sizeof(double)*NN*pp);
    double *z = malloc(sizeof(double)*NN*GG);
    get_data(x1,x,*N,*p);
    get_data(z1,z,*N,*G);
	
	if (*clust!=0){
		func=funcs[*model];
		*bic = func(z, x, cls, *q, *p, *G, *N, lambda, psi,*tol);
	}else{
		func2=funcs2[*model];
		*bic = func2(z, x, *q, *p, *G, *N, lambda, psi,*tol);
	}
	
    give_data(z1,z,*N,*G);
      free(x); free(z);    
}

int convergtest_NEW(double *l, double *at, double *v_max, double *v, int N, int it, int G, double TOL){
	int i,g, flag=0;
	double sum, l_inf;

	l[it]=0.0;
    for(i=0; i<N; i++){
        sum=0.0; 
        for(g=0; g<G; g++){
            sum += exp(v[i*G+g]-v_max[i]);
	    }
		l[it] += log(sum)+v_max[i];
	    if(isnan(l[it])||isinf(l[it])) return -1;
    }
	
    if(it > 0)
		if(l[it]<l[it-1])
			return -1;   
        
	if(it > 2){
        at[it-1]=(l[it]-l[it-1])/(l[it-1]-l[it-2]);
	    if(at[it-1]<1.0){
            l_inf = l[it-1]+(l[it]-l[it-1])/(1-at[it-1]); 
			if(fabs(l_inf - l[it])<TOL) flag=1; 
	    }     
	}
	return flag;
}

void get_data(double *x1, double *x, int n, int m){
   int i,j,k;
   i=0;j=0;k=0; 
   for (i=0;i<n;i++){
	  for (j=0;j<m;j++){
		 x[i*m+j]=x1[k];
		 k++;
	  }
   }
}

void get_data2(double *x1, double **x, int G, int n, int m){
   int g,i,k=0;
   for (g=0;g<G;g++){
	for (i=0;i<n*m;i++){
		 x[g][i]=x1[k];
		 k++;
	  }
	}
}

void give_data(double *x1, double *x, int n, int m){
	int i,j,k;
	i=0;j=0;k=0; 
	for (i=0;i<n;i++){
		for (j=0;j<m;j++){
			x1[k]=x[i*m+j];
			k++;
		}
	}
}		

void known_z(int *class, double *z, int N, int G){
	int i,g;
	for(i=0;i<N;i++){
		if(class[i]!=0){
			for(g=1;g<=G;g++){
				z[i*G+g-1]=0.0;
				if(g==class[i]) z[i*G+g-1]=1.0;
			}
		}
	}
}
