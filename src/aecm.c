#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<R.h>
#include "functions.h"

double aecm(double **z, double **x, int *cls, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol){
	
	double **zR6, **mu, bic, **kf4, **po1, **y2etilde, **v, *max_v, *l, *at, *pi, *n; 
	int it=0,stop=0,paras;
	double bx1, ja2, log_c=0.0, q2h;

   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc(&y2etilde,p,p);		
	my_alloc_vec(&max_v,N);
   	my_alloc(&v,N,G);		
   	my_alloc(&zR6,p,q);
   	my_alloc(&kf4,q,p);
   	my_alloc(&po1,q,q);
   	my_alloc(&mu,G,p);
	
	bx1 = *bx1_vec;	
	get_data(smf,zR6,p,q);
	
    while(stop==0){            
    
		p5M_n(n, z, G, N);
	
	    p5M_pi(pi, n, G, N);
	  
	    updau8m(mu, n, x, z, G, N, p);

		if(it>0){ 
			p5M_z(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
			known_z(cls,z,N,G);
		}

	    p5M_stilde(y2etilde,x, z, mu, G, N, p);

	    p5M_kf41(kf4, bx1, zR6, p, q);
        
   	    p5M_po1(po1, kf4, zR6, y2etilde, p, q); 
    
        p5M_zR6(zR6, kf4, y2etilde, po1, p, q);
 
        bx1 = p5M_bx1(zR6, kf4, y2etilde, p, q); 
   
 	    ja2=p*log(bx1);
	    
	    q2h = p5M_det_sigma_NEW(zR6, bx1, ja2, p, q);
	  
	    log_c = (p/2.0)*log(2.0*M_PI) + 0.5*q2h;
       
		p5M_z(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
		known_z(cls,z,N,G);
 
        stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++;
		/*printf("ll=%f\t",l[it-1]);*/
	}
	/*printf("ll=%f\n",l[it-1]);*/

	paras = G-1 + G*p + p*q - q*(q-1)/2 + 1;
    	
   	bic = 2.0*l[it-1] - paras*log(N);
	
	/* printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it-1); */
    
   	my_free(zR6,p); my_free(mu,G); free(n); my_free(kf4,q); 
	my_free(po1,q); my_free(y2etilde,p); free(l); free(at); free(pi);
       
	return bic;
} 

double aecm2(double **z, double **x, int *cls, int q, int p, int G, int N, double *smf, double *bx1, double tol){

   	double **zR6, **mu, bic, **w, **kf4, **po1, 
	       **y2etilde, **v, *max_v, *l, *at, *pi, *det, *n; 
	int i,it=0,stop=0;
	double a, ja2,q2h, log_c=0.0;

	/* printf("G=%d \t q=%d\n",G,q);*/

   	/* Allocate memory - matrices and vectors */ 
   	my_alloc_vec(&det,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc(&y2etilde,p,p);		
   	my_alloc(&zR6,p,q);			
   	my_alloc(&kf4,q,p);
   	my_alloc(&po1,q,q);
   	my_alloc(&mu,G,p);
   	my_alloc(&w,G,N);
	my_alloc_vec(&max_v,N);
   	my_alloc(&v,N,G);	
		
	get_data(smf,zR6,p,q);
		
   	while(stop==0){            
	    
	    p5M_n(n, z, G, N);
	    
	    p5M_pi(pi, n, G, N);
	    
	    updau8m(mu, n, x, z, G, N, p); 

        if (it>0){ 
			p5M_z2(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
			known_z(cls,z,N,G);
		}
        
		p5M_stilde(y2etilde, x, z, mu, G, N, p);        
        
	    p5M_kf42(kf4, bx1, zR6, p, q);

	    p5M_po1(po1, kf4, zR6, y2etilde, p, q);

	    p5M_zR6(zR6, kf4, y2etilde, po1, p, q);
			    
	    p5M_bx12(bx1, zR6, kf4, y2etilde, p, q);

		ja2 = 0.0; for(i=0;i<p;i++) ja2 += log(bx1[i]);

	    q2h = p5M_det_sigma_NEW2(zR6, bx1, ja2, p, q);
	    
		log_c = (p/2.0)*log(2.0*M_PI) + 0.5*q2h; 
   
		p5M_z2(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
		known_z(cls,z,N,G);
 
        stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	
	    it++; 
    }
	/*printf("ll=%f\n",l[it-1]);*/
    
	/* no of parameters {pq - q(q - 1)/2} + p */
    a = G-1 + G*p + p*q - q*(q-1)/2 + p;
    
    /* BIC = 2log_likelihood - mlog n; */
    bic = 2.0*l[it-1] - a*log(N);
    
    /* printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it-1); */   
    
    /* Deallocate memory */
    my_free(zR6,p); my_free(mu,G); my_free(w,G); free(n); free(det); 
    my_free(kf4,q); my_free(po1,q); my_free(y2etilde,p); free(l); 
    free(at); free(pi);   
	return bic;
}

double aecm3(double **z, double **x, int *cls, int q, int p, int G, int N, double *smf, double *bx1, double tol){

	double **zR6, **mu, bic, ***kf4, ***po1, ***y2e, 
	       *l, *at, *pi, *n, *ja2, *q2h, *log_det,
			**v, *max_v, a; 
	int g,it=0,stop=0;

	/* printf("G=%d \t q=%d\n",G,q);*/

    my_alloc_vec(&log_det,G);
    my_alloc_vec(&ja2,G);
    my_alloc_vec(&q2h,G);
   	my_alloc_vec(&pi,G);
    my_alloc_vec(&n,G);
    my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
    my_alloc_3d(&y2e,G,p,p);		
    my_alloc(&zR6,p,q);		
    my_alloc_3d(&kf4,G,q,p);		
    my_alloc_3d(&po1,G,q,q);
    my_alloc(&mu,G,p);
	my_alloc_vec(&max_v,N);
   	my_alloc(&v,N,G);		
    
	get_data(smf,zR6,p,q);
		
    while(stop==0){            
           
        p5M_n(n, z, G, N);
	
	    p5M_pi(pi, n, G, N);
	
	    updau8m(mu, n, x, z, G, N, p); 
                          
		if(it>0){ 
			p5M_z3(v, x, z, zR6, bx1, mu, pi, max_v, log_det, N, G, p, q);
			known_z(cls,z,N,G);
        }

	    p5M_sg(y2e, x, z, mu, n, p, G, N);
	    
	    for(g=0;g<G;g++) p5M_kf41(kf4[g], bx1[g], zR6, p, q);    
	
    	for (g=0;g<G;g++) p5M_po1(po1[g], kf4[g], zR6, y2e[g], p, q);
        
	    p5M_zR62(zR6, kf4, y2e, po1, n, bx1, p, q, G);
	    
        for(g=0;g<G;g++) bx1[g] = p5M_bx13(zR6, kf4[g], y2e[g], po1[g], p, q);
    
        for(g=0;g<G;g++){
	    	ja2[g] = p*log(bx1[g]);
			q2h[g] = p5M_det_sigma_NEW(zR6, bx1[g], ja2[g], p, q);
			log_det[g] = (p/2.0)*log(2*M_PI) + 0.5*q2h[g];
		} 
	    
		p5M_z3(v, x, z, zR6, bx1, mu, pi, max_v, log_det, N, G, p, q);
		known_z(cls,z,N,G);

	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
		it++;
    }
	
	/*printf("ll=%f\n",l[it-1]);*/
	
   	a = G-1 + G*p+ p*q - q*(q-1)/2 + G;
    
   	bic = 2.0*l[it-1] - a*log(N);
    
   	/* printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it); */    
    
   	/* Deallocate memory */
   	my_free(zR6,p); my_free(mu,G); my_free(v,G); free(n); free(log_det); free(max_v); 
	my_free_3d(kf4,G,q); my_free_3d(po1,G,q); my_free_3d(y2e,G,p); 
	free(l); free(at); free(pi); free(ja2);
    
	return bic;
}

double aecm4(double **z, double **x, int *cls, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol){

	double **zR6, **mu, bic, **w, ***kf4, ***po1, ***y2e, *l, *at, *pi, *n,
	       **bx1, *ja2, *q2h, *log_c, **v, *max_v; 
	int it=0,stop=0,paras, g, j;
	
	/* printf("G=%d \t q=%d\n",G,q);*/

   	/* Allocate memory - matrices and vectors */ 
   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&y2e,G,p,p);		
   	my_alloc(&zR6,p,q);
   	my_alloc_3d(&kf4,G,q,p);
   	my_alloc_3d(&po1,G,q,q);
   	my_alloc(&mu,G,p);
   	my_alloc(&w,G,N);
	my_alloc(&bx1,G,p);
	my_alloc_vec(&ja2,G);
	my_alloc_vec(&q2h,G);
	my_alloc_vec(&log_c,G);
	
	get_data(bx1_vec,bx1,G,p);
	get_data(smf,zR6,p,q);
        
   	while(stop==0){            
   
	    p5M_n(n, z, G, N);
	
	    p5M_pi(pi, n, G, N);
	
	    updau8m(mu, n, x, z, G, N, p); 

        if(it>0){ 
			p5M_z4(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
			known_z(cls,z,N,G);
		} 
 
	    p5M_sg(y2e,x, z, mu, n, p, G, N);
        
	    for(g=0;g<G;g++) p5M_kf42(kf4[g], bx1[g], zR6, p, q);
        
   	    for(g=0;g<G;g++) p5M_po1(po1[g], kf4[g], zR6, y2e[g], p, q); 
    
	    p5M_zR6_cuu(zR6, kf4, y2e, po1, n, bx1, p, q, G);
	 
        p5M_bx1_cuu(bx1, zR6, kf4, y2e, po1, p, q, G);

	    for (g=0;g<G;g++){
			ja2[g]=0.0;
			for(j=0;j<p;j++) ja2[g] += log(bx1[g][j]);
	    }
	
	    for(g=0;g<G;g++) q2h[g] = p5M_det_sigma_NEW2(zR6, bx1[g], ja2[g], p, q);

		for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*q2h[g];

        p5M_z4(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
		known_z(cls,z,N,G);
            
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++;
	}
	
	/* no of parameters {pq - q(q - 1)/2} + Gp */
	paras = G-1 + G*p + p*q - q*(q-1)/2 + G*p;
    
   	/* BIC = 2log_likelihood - mlog n; */
   	bic = 2.0*l[it-1] - paras*log(N);
    
	/*if((G==3)&&(q==4)) printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it);*/
    
    /* Deallocate memory */
    my_free(zR6,p); my_free(mu,G); my_free(w,G); free(n); my_free_3d(kf4,G,q); 
	my_free_3d(po1,G,q); my_free_3d(y2e,G,p); free(l); free(at); free(pi);
	free(q2h); free(log_c); free(ja2);
	return bic;
} 

double aecm5(double **z, double **x, int*cls, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol){

	double ***zR6, **mu, bic, **w, ***kf4, ***po1, **v, *max_v, ***y2e, *l, *at, *pi, 
	*n, *q2h, *log_c, ja2, bx1; 
	int it=0,stop=0,paras, g,i;

	/* printf("G=%d \t q=%d\n",G,q);*/

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
   	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&q2h,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&y2e,G,p,p);		
   	my_alloc_3d(&zR6,G,p,q);
   	my_alloc_3d(&kf4,G,q,p);
   	my_alloc_3d(&po1,G,q,q);
   	my_alloc(&mu,G,p);
   	my_alloc(&w,G,N);
	
	bx1=*bx1_vec;
	get_data2(smf,zR6,G,p,q);
       
   	while(stop==0){            
    
	    p5M_n(n, z, G, N);
	
	    p5M_pi(pi, n, G, N);
	
	    updau8m(mu, n, x, z, G, N, p); 

		if(it>0){ 
			p5M_z5(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
			known_z(cls,z,N,G);
        }
 
	    p5M_sg(y2e,x, z, mu, n, p, G, N);
        
	    for(g=0;g<G;g++) p5M_kf41(kf4[g], bx1, zR6[g], p, q);
        
   	    for(g=0;g<G;g++) p5M_po1(po1[g], kf4[g], zR6[g], y2e[g], p, q); 
    
        for (g=0;g<G;g++) p5M_zR6(zR6[g], kf4[g], y2e[g], po1[g], p, q);
	 
	    bx1 = p5M_bx1_ucc(zR6, kf4, y2e, p, q, pi, G);
   
		ja2 = 0.0;
		for(i=0;i<p;i++) ja2 += log(bx1);
	    
	    for(g=0;g<G;g++) q2h[g] = p5M_det_sigma_NEW(zR6[g], bx1, ja2, p, q);
	    
		for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*q2h[g];
       
		p5M_z5(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
		known_z(cls,z,N,G);
            
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++; 
	}
	/*printf("ll=%f\n",l[it-1]);*/

	paras = G-1 + G*p + G*(p*q - q*(q-1)/2) + 1;
    
   	bic = 2.0*l[it-1] - paras*log(N);
    
	/* printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it); */   
    
   	my_free_3d(zR6,G,p); my_free(mu,G); my_free(w,G); free(n); 
	my_free_3d(kf4,G,q); my_free_3d(po1,G,q); my_free_3d(y2e,G,p); 
	free(l); free(at); free(pi);free(q2h); free(log_c);
	return bic;
}

double aecm6(double **z, double **x, int *cls, int q, int p, int G, int N, double *smf, double *bx1, double tol){

	double ***zR6, **mu, bic, ***kf4, ***po1,**v,*max_v, ***y2e, *l, *at, *pi, 
	*n, *q2h,*log_c, ja2; 
	int it=0,stop=0,paras, g, j;

	/* printf("G=%d \t q=%d\n",G,q);*/

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
	my_alloc_vec(&q2h,G);
	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&y2e,G,p,p);		
   	my_alloc_3d(&zR6,G,p,q);
   	my_alloc_3d(&kf4,G,q,p);
   	my_alloc_3d(&po1,G,q,q);
   	my_alloc(&mu,G,p);
	
	get_data2(smf,zR6,G,p,q);
      
   	while(stop==0){            
    
	    p5M_n(n, z, G, N);
	
	    p5M_pi(pi, n, G, N);
	
	    updau8m(mu, n, x, z, G, N, p); 

        if(it>0){ 
			p5M_z6(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
			known_z(cls,z,N,G);
		} 

	    p5M_sg(y2e,x, z, mu, n, p, G, N);
        
	    for(g=0;g<G;g++) p5M_kf42(kf4[g], bx1, zR6[g], p, q);
        
   	    for(g=0;g<G;g++) p5M_po1(po1[g], kf4[g], zR6[g], y2e[g], p, q); 
    
        for (g=0;g<G;g++) p5M_zR6(zR6[g], kf4[g], y2e[g], po1[g], p, q);
	 
	    p5M_bx1_ucu(bx1, zR6, kf4, y2e, p, q, pi, G);
                    
		ja2 = 0.0;
		for(j=0;j<p;j++) ja2 += log(bx1[j]);
	    
	    for(g=0;g<G;g++) q2h[g] = p5M_det_sigma_NEW2(zR6[g], bx1, ja2, p, q);
	    
		for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*q2h[g];

        p5M_z6(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
		known_z(cls,z,N,G);
            
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++; 
	}
	/*printf("ll=%f\n",l[it-1]);*/
	
	paras = G-1 + G*p + G*(p*q - q*(q-1)/2) + p;
    
   	bic = 2.0*l[it-1] - paras*log(N);
    
	/* printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it); */   
    
   	my_free_3d(zR6,G,p); my_free(mu,G); my_free(v,N); free(n);free(max_v); 
	my_free_3d(kf4,G,q); my_free_3d(po1,G,q); my_free_3d(y2e,G,p); 
	free(l); free(at); free(pi); free(q2h); free(log_c);
	return bic;
}

double aecm7(double **z, double **x, int *cls, int q, int p, int G, int N, double *smf, double *bx1, double tol){

	double ***zR6, **mu, bic, ***kf4, ***po1, **v, *max_v, ***y2e, *l, *at, *pi, 
		*n, *ja2, *q2h, *log_c; 
	int it=0,stop=0,paras, g;

	/* printf("G=%d \t q=%d\n",G,q); */

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
	my_alloc_vec(&ja2,G);
	my_alloc_vec(&q2h,G);
	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&y2e,G,p,p);		
   	my_alloc_3d(&zR6,G,p,q);
   	my_alloc_3d(&kf4,G,q,p);
   	my_alloc_3d(&po1,G,q,q);
   	my_alloc(&mu,G,p);
	
	get_data2(smf,zR6,G,p,q);

   	while(stop==0){            
    
	    p5M_n(n, z, G, N);
	
	    p5M_pi(pi, n, G, N);
	
	    updau8m(mu, n, x, z, G, N, p); 

        if(it>0){ 
			p5M_z7(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
			known_z(cls,z,N,G);
		}

	    p5M_sg(y2e,x, z, mu, n, p, G, N);
        
	    for(g=0;g<G;g++) p5M_kf41(kf4[g], bx1[g], zR6[g], p, q);
        
   	    for(g=0;g<G;g++) p5M_po1(po1[g], kf4[g], zR6[g], y2e[g], p, q); 
    
        for (g=0;g<G;g++) p5M_zR6(zR6[g], kf4[g], y2e[g], po1[g], p, q);
	 
        for (g=0;g<G;g++) bx1[g]=p5M_bx1(zR6[g], kf4[g], y2e[g], p, q); 
           
	    for(g=0;g<G;g++) ja2[g] = p*log(bx1[g]);
	    
	    for(g=0;g<G;g++) q2h[g] = p5M_det_sigma_NEW(zR6[g], bx1[g], ja2[g], p, q);
	   
		for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*q2h[g];

		p5M_z7(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
		known_z(cls,z,N,G);
            
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++;
	}
	/*printf("ll=%f\n",l[it-1]);*/

	paras = G-1 + G*p + G*(p*q - q*(q-1)/2) + G;
    
   	bic = 2.0*l[it-1] - paras*log(N);
    
	/* printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it); */   
    
    my_free_3d(zR6,G,p); my_free(mu,G); my_free(v,N); free(n); 
	my_free_3d(kf4,G,q); free(ja2); my_free_3d(po1,G,q); 
	my_free_3d(y2e,G,p); free(l); free(at); free(pi);
	free(q2h); free(log_c);
	
	return bic;
}

double aecm8(double **z, double **x, int *cls, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol){

	double ***zR6, **mu, bic, ***kf4, ***po1, **v, *max_v, ***y2e, *l, *at, *pi, 
	*n, **bx1, *log_c,*q2h,*ja2; 
	int it=0,stop=0,paras, g, j;

	/* printf("G=%d \t q=%d\n",G,q); */

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&y2e,G,p,p);		
   	my_alloc_3d(&zR6,G,p,q);
   	my_alloc_3d(&kf4,G,q,p);
   	my_alloc_3d(&po1,G,q,q);
   	my_alloc(&mu,G,p);
	my_alloc(&bx1,G,p);
	my_alloc_vec(&ja2,G);
	my_alloc_vec(&q2h,G);
	my_alloc_vec(&log_c,G);
	
	get_data(bx1_vec,bx1,G,p);
	get_data2(smf,zR6,G,p,q);
	       
   	while(stop==0){            
    
	    p5M_n(n, z, G, N);
	
	    p5M_pi(pi, n, G, N);
	
	    updau8m(mu, n, x, z, G, N, p);

        if(it>0){ 
			p5M_z8(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
			known_z(cls,z,N,G);
		}	
       
	    p5M_sg(y2e,x, z, mu, n, p, G, N);
	    
	    for(g=0;g<G;g++) p5M_kf42(kf4[g], bx1[g], zR6[g], p, q);
        
   	    for(g=0;g<G;g++) p5M_po1(po1[g], kf4[g], zR6[g], y2e[g], p, q); 
    
        for (g=0;g<G;g++) p5M_zR6(zR6[g], kf4[g], y2e[g], po1[g], p, q);
	 
        for (g=0;g<G;g++) p5M_bx12(bx1[g], zR6[g], kf4[g], y2e[g], p, q); 
           
	    for(g=0;g<G;g++){
		    ja2[g]=0.0;
		    for(j=0;j<p;j++) ja2[g] += log(bx1[g][j]);
	    }
	    
		for(g=0;g<G;g++) q2h[g] = p5M_det_sigma_NEW2(zR6[g], bx1[g], ja2[g], p, q);
	    	    
		for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*q2h[g];

		p5M_z8(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
		known_z(cls,z,N,G);
            
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++; 
	}
	/*printf("ll=%f\n",l[it-1]);*/
	paras = G-1 + G*p + G*(p*q - q*(q-1)/2) + G*p;
    
   	bic = 2.0*l[it-1] - paras*log(N);
    
	/* printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it); */   
    
   	my_free_3d(zR6,G,p); my_free(mu,G); my_free(v,N); free(n); 
	my_free_3d(kf4,G,q); free(ja2); my_free_3d(po1,G,q); 
	my_free_3d(y2e,G,p); free(l); free(at); free(pi);
	free(q2h); free(log_c); my_free(bx1,G); free(max_v);
	return bic;
}

double aecm9(double **z, double **x, int *cls, int q, int p, int G, int N, double *smf, double *dsw, double tol){

	double **zR6, **mu, bic, ***kf4, ***po1, ***y2e, **v, *max_v, *l, *at, *pi, 
	*n, *g9p, *log_c, *ja2, *q2h, a, *bx1; 
	int g,i,it=0,stop=0;
	       
	/* printf("G=%d \t q=%d\n",G,q);*/

   	/* Allocate memory - matrices and vectors */ 
   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
	my_alloc_vec(&ja2,G);
	my_alloc_vec(&q2h,G);
	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&y2e,G,p,p);		
   	my_alloc(&zR6,p,q);		
   	my_alloc_3d(&kf4,G,q,p);		
   	my_alloc_3d(&po1,G,q,q);
   	my_alloc(&mu,G,p);
	my_alloc_vec(&g9p,p);
	my_alloc_vec(&bx1,p);
	
	get_data(smf,zR6,p,q);
	for(i=0;i<p;i++) g9p[i] = 1.0; 
	
    while(stop==0){            
           
        p5M_n(n, z, G, N);
	
	    p5M_pi(pi, n, G, N);
	
	    updau8m(mu, n, x, z, G, N, p); 
                          
        if(it>0){ 
			p5M_z9(v, x, z, zR6, dsw, g9p, mu, pi, max_v, log_c, N, G,p,q);
			known_z(cls,z,N,G);
		}	
		
	    p5M_sg(y2e, x, z, mu, n, p, G, N);

	    for(g=0;g<G;g++){ 
			for (i=0;i<p;i++) bx1[i] = dsw[g]*g9p[i];
			p5M_kf42(kf4[g], bx1, zR6, p, q);    
		}	

    	for (g=0;g<G;g++) p5M_po1(po1[g], kf4[g], zR6, y2e[g], p, q);
        
	    p5M_zR62(zR6, kf4, y2e, po1, n, dsw, p, q, G);
	    
        for(g=0;g<G;g++) dsw[g] = p5M_dsw(zR6, g9p, kf4[g], y2e[g], po1[g], p, q);

	    p5M_g9p(g9p, zR6, dsw, kf4, y2e, po1, n, p, q, N, G); 
      
        for(g=0;g<G;g++){
			for (i=0;i<p;i++) bx1[i] = dsw[g]*g9p[i];
	    	ja2[g] = p*log(dsw[g]);
			q2h[g] = p5M_det_sigma_NEW2(zR6, bx1, ja2[g], p, q);
			log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*q2h[g];
		}

	    p5M_z9(v, x, z, zR6, dsw, g9p, mu, pi, max_v, log_c, N, G,p,q);
        known_z(cls,z,N,G);
		
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
        it++;  
    }
	/*printf("ll=%f\n",l[it-1]);*/

    /* no of parameters */
    a = G-1 + G*p+ p*q - q*(q-1)/2 + G + (p-1);
    
    /* BIC = 2log_likelihood - mlog n; */
    bic = 2.0*l[it-1] - a*log(N);
    
    /* printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it); */
    
    /* Deallocate memory */
    my_free(zR6,p); my_free(mu,G); my_free(v,N); free(n); free(log_c); 
	my_free_3d(kf4,G,q); my_free_3d(po1,G,q); my_free_3d(y2e,G,p); 
	free(l); free(at); free(pi); free(ja2); free(g9p); free(q2h);
	
	return bic;
}

double aecm10(double **z, double **x, int *cls, int q, int p, int G, int N, double *smf, double *dsw, double tol){

	double ***zR6, **mu, bic, ***kf4, ***po1, ***y2e,**v, *max_v, *l, *at, *pi, 
	*n, *g9p, *q2h, *log_c, *ja2, a, *bx1; 
	int g,i,it=0,stop=0;
	       
	/* printf("G=%d \t q=%d\n",G,q);*/

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
	my_alloc_vec(&ja2,G);
	my_alloc_vec(&q2h,G);
	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&y2e,G,p,p);		
   	my_alloc_3d(&zR6,G,p,q);		
   	my_alloc_3d(&kf4,G,q,p);		
   	my_alloc_3d(&po1,G,q,q);
   	my_alloc(&mu,G,p);
	my_alloc_vec(&g9p,p);
	my_alloc_vec(&bx1,p);
	
	get_data2(smf,zR6,G,p,q);
	for(i=0;i<p;i++) g9p[i] = 1.0; 
	
   	while(stop==0){            
           
        p5M_n(n, z, G, N);
	
	    p5M_pi(pi, n, G, N);
	
	    updau8m(mu, n, x, z, G, N, p); 
                          
        if(it>0){ 
			p5M_z10(v, x, z, zR6, dsw, g9p, mu, pi, max_v, log_c, N, G,p,q);
            known_z(cls,z,N,G);
		}	
		
	    p5M_sg(y2e, x, z, mu, n, p, G, N);

	    for(g=0;g<G;g++){ 
			for (i=0;i<p;i++) bx1[i] = dsw[g]*g9p[i];
			p5M_kf42(kf4[g], bx1, zR6[g], p, q);    
		}
	
    	for (g=0;g<G;g++) p5M_po1(po1[g], kf4[g], zR6[g], y2e[g], p, q);
        
	    for(g=0;g<G;g++) p5M_zR6(zR6[g], kf4[g], y2e[g], po1[g], p, q);
	    
        for(g=0;g<G;g++) dsw[g] = p5M_dsw2(zR6[g], g9p, kf4[g], y2e[g], p, q);
	   	
	    p5M_g9p2(g9p, zR6, dsw, kf4, y2e, po1, n, p, q, N, G); 
	    
        for(g=0;g<G;g++){
			for (i=0;i<p;i++) bx1[i] = dsw[g]*g9p[i];
			ja2[g] = p*log(dsw[g]);
	    	q2h[g] = p5M_det_sigma_NEW2(zR6[g], bx1, ja2[g], p, q);
			log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*q2h[g];
		}

	    p5M_z10(v, x, z, zR6, dsw, g9p, mu, pi, max_v, log_c, N, G,p,q);
        known_z(cls,z,N,G);
		
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
        it++;
   	}
	/*printf("ll=%f\n",l[it-1]);*/
   
   	a = G-1 + G*p+ G*(p*q - q*(q-1)/2) + G + (p-1);
    
   	bic = 2.0*l[it-1] - a*log(N);
    
   	/* printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it); */
    
   	my_free_3d(zR6,G,p); my_free(mu,G); my_free(v,N); free(n); free(log_c); 
	my_free_3d(kf4,G,q); my_free_3d(po1,G,q); my_free_3d(y2e,G,p); 
	free(l); free(at); free(pi); free(ja2); free(g9p); free(q2h);

	return bic;
}

double aecm11(double **z, double **x, int *cls, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol){

	double **zR6, **mu, bic, ***kf4, ***po1, ***y2e, **v, *max_v, *l, *at, *pi, 
	*n, **g9p, *log_c, *q2h, ja2, a, dsw, *bx1; 
	int g,i,it=0,stop=0;

	/* printf("G=%d \t q=%d\n",G,q);*/

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&y2e,G,p,p);		
   	my_alloc(&zR6,p,q);		
   	my_alloc_3d(&kf4,G,q,p);		
	my_alloc_3d(&po1,G,q,q);
   	my_alloc(&mu,G,p);
	my_alloc(&g9p,G,p);
	my_alloc_vec(&q2h,G);
	my_alloc_vec(&log_c,G);
	my_alloc_vec(&bx1,p);
	
	dsw = *bx1_vec;
	get_data(smf,zR6,p,q);
	
	for(g=0;g<G;g++)
		for(i=0;i<p;i++) 
			g9p[g][i] = 1.0; 
	
   	while(stop==0){            
           
        p5M_n(n, z, G, N);
	
	    p5M_pi(pi, n, G, N);
	
	    updau8m(mu, n, x, z, G, N, p); 
                          
        if(it>0){ 
			p5M_z11(v, x, z, zR6, dsw, g9p, mu, pi, max_v, log_c, N, G,p,q);
            known_z(cls,z,N,G);
		}	
		
	    p5M_sg(y2e, x, z, mu, n, p, G, N);

	    for(g=0;g<G;g++){ 
			for (i=0;i<p;i++) bx1[i] = dsw*g9p[g][i];
			p5M_kf42(kf4[g], bx1, zR6, p, q);    
		}	
	
    	for (g=0;g<G;g++) p5M_po1(po1[g], kf4[g], zR6, y2e[g], p, q);
        
	    p5M_zR6_cuu(zR6, kf4, y2e, po1, n, g9p, p, q, G);
	    
        dsw =0.0;
	    for(g=0;g<G;g++) dsw += pi[g]*p5M_dsw(zR6, g9p[g], kf4[g], y2e[g], po1[g], p, q);
 	   	
	    for(g=0;g<G;g++) p5M_g9p3(g9p[g], zR6, dsw, kf4[g], y2e[g], po1[g], n[g], p, q); 
		
		ja2 = p*log(dsw);

        for(g=0;g<G;g++){
			for (i=0;i<p;i++) bx1[i] = dsw*g9p[g][i];
	    	q2h[g] = p5M_det_sigma_NEW2(zR6, bx1, ja2, p, q);
			log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*q2h[g];
		}

	    p5M_z11(v, x, z, zR6, dsw, g9p, mu, pi, max_v, log_c, N, G,p,q);
        known_z(cls,z,N,G);

	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
        it++;
    }
	/*printf("ll=%f\n",l[it-1]);*/

   	a = G-1 + G*p+ p*q - q*(q-1)/2 + 1 + G*(p-1);
    
   	bic = 2.0*l[it-1] - a*log(N);
    
   	/* printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it); */
    
   	my_free(zR6,p); my_free(mu,G); my_free(v,N); free(n); free(log_c); my_free_3d(kf4,G,q); my_free_3d(po1,G,q); 
	my_free_3d(y2e,G,p); free(l); free(at); free(pi); my_free(g9p,G); free(q2h); 
	
	return bic;
} 

double aecm12(double **z, double **x, int *cls, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol){

	double ***zR6, **mu, bic, ***kf4, ***po1, ***y2e, **v, *max_v, *l, *at, *pi, 
	*n, **g9p, *q2h, *log_c, ja2, a, dsw, *bx1; 
	int g,i,it=0,stop=0;
	       
	/* printf("G=%d \t q=%d\n",G,q);*/

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
	my_alloc_vec(&q2h,G);
	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&y2e,G,p,p);		
   	my_alloc_3d(&zR6,G,p,q);		
   	my_alloc_3d(&kf4,G,q,p);		
	my_alloc_3d(&po1,G,q,q);
   	my_alloc(&mu,G,p);
	my_alloc(&g9p,G,p);
	my_alloc_vec(&bx1,p);
	
	dsw = *bx1_vec;
	get_data2(smf,zR6,G,p,q);
	
	for(g=0;g<G;g++)
		for(i=0;i<p;i++) 
			g9p[g][i] = 1.0; 
	
   	while(stop==0){            
           
		p5M_n(n, z, G, N);
	
	    p5M_pi(pi, n, G, N);
	
	    updau8m(mu, n, x, z, G, N, p); 
                          
	    if(it>0){
			p5M_z12(v, x, z, zR6, dsw, g9p, mu, pi, max_v, log_c, N, G,p,q);
            known_z(cls,z,N,G);
		}	
		
	    p5M_sg(y2e, x, z, mu, n, p, G, N);
	    
	    for(g=0;g<G;g++){ 
			for (i=0;i<p;i++) bx1[i] = dsw*g9p[g][i];
			p5M_kf42(kf4[g], bx1, zR6[g], p, q);    
		}	
	
    	for (g=0;g<G;g++) p5M_po1(po1[g], kf4[g], zR6[g], y2e[g], p, q);
        
    	for (g=0;g<G;g++) p5M_zR6(zR6[g], kf4[g], y2e[g], po1[g], p, q);
		    
        dsw =0.0;
	    for(g=0;g<G;g++) dsw += pi[g]*p5M_dsw2(zR6[g], g9p[g], kf4[g], y2e[g], p, q);
	    
	    for(g=0;g<G;g++) p5M_g9p3(g9p[g], zR6[g], dsw, kf4[g], y2e[g], po1[g], n[g], p, q); 
	    
		ja2 = p*log(dsw);        

		for(g=0;g<G;g++){
			for (i=0;i<p;i++) bx1[i] = dsw*g9p[g][i];
			q2h[g] = p5M_det_sigma_NEW2(zR6[g], bx1, ja2, p, q);
			log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*q2h[g];
		}

	    p5M_z12(v, x, z, zR6, dsw, g9p, mu, pi, max_v, log_c, N, G,p,q);
        known_z(cls,z,N,G);

	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
        it++;
    }
	/*printf("ll=%f\n",l[it-1]);*/
   
   	a = G-1 + G*p+ G*(p*q - q*(q-1)/2) + 1+ G*(p-1);
    
   	bic = 2.0*l[it-1] - a*log(N);
    
   	/* printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it); */
    
   	my_free_3d(zR6,G,p); my_free(mu,G); my_free(v,N); free(n);  
	my_free_3d(kf4,G,q); my_free_3d(po1,G,q); my_free_3d(y2e,G,p); 
	free(l); free(at); free(pi); my_free(g9p,G); free(log_c);
	free(q2h);

	return bic;
}

double claecm(double **z, double **x, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol){
	
	double **zR6, **mu, bic, **kf4, **po1, **y2etilde, **v, *max_v, *l, *at, *pi, *n; 
	int it=0,stop=0,paras;
	double bx1, ja2, log_c=0.0, q2h;
	
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc(&y2etilde,p,p);		
	my_alloc_vec(&max_v,N);
   	my_alloc(&v,N,G);		
   	my_alloc(&zR6,p,q);
   	my_alloc(&kf4,q,p);
   	my_alloc(&po1,q,q);
   	my_alloc(&mu,G,p);
	
	bx1 = *bx1_vec;	
	get_data(smf,zR6,p,q);
	
    while(stop==0){            
    
		p5M_n(n, z, G, N);
	
	    p5M_pi(pi, n, G, N);
	  
	    updau8m(mu, n, x, z, G, N, p);

		if(it>0){ 
			p5M_z(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
		}

	    p5M_stilde(y2etilde,x, z, mu, G, N, p);

	    p5M_kf41(kf4, bx1, zR6, p, q);
        
   	    p5M_po1(po1, kf4, zR6, y2etilde, p, q); 
    
        p5M_zR6(zR6, kf4, y2etilde, po1, p, q);
 
        bx1 = p5M_bx1(zR6, kf4, y2etilde, p, q); 
   
 	    ja2=p*log(bx1);
	    
	    q2h = p5M_det_sigma_NEW(zR6, bx1, ja2, p, q);
	  
	    log_c = (p/2.0)*log(2.0*M_PI) + 0.5*q2h;
       
		p5M_z(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
		 
        stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++;
		/*printf("ll=%f\t",l[it-1]);*/
	}
	/*printf("ll=%f\n",l[it-1]);*/

	paras = G-1 + G*p + p*q - q*(q-1)/2 + 1;
    	
   	bic = 2.0*l[it-1] - paras*log(N);
	
	/* printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it-1); */
	    
   	my_free(zR6,p); my_free(mu,G); free(n); my_free(kf4,q); 
	my_free(po1,q); my_free(y2etilde,p); free(l); free(at); free(pi);
       
	return bic;
} 

double claecm2(double **z, double **x, int q, int p, int G, int N, double *smf, double *bx1, double tol){

   	double **zR6, **mu, bic, **w, **kf4, **po1, 
	       **y2etilde, **v, *max_v, *l, *at, *pi, *det, *n; 
	int i,it=0,stop=0;
	double a, ja2,q2h, log_c=0.0;

	/* printf("G=%d \t q=%d\n",G,q);*/

   	/* Allocate memory - matrices and vectors */ 
   	my_alloc_vec(&det,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc(&y2etilde,p,p);		
   	my_alloc(&zR6,p,q);			
   	my_alloc(&kf4,q,p);
   	my_alloc(&po1,q,q);
   	my_alloc(&mu,G,p);
   	my_alloc(&w,G,N);
	my_alloc_vec(&max_v,N);
   	my_alloc(&v,N,G);		

	get_data(smf,zR6,p,q);
	
   	while(stop==0){            
	    
	    p5M_n(n, z, G, N);
	    
	    p5M_pi(pi, n, G, N);
	    
	    updau8m(mu, n, x, z, G, N, p); 

        if (it>0){ 
			p5M_z2(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
		}
        
		p5M_stilde(y2etilde, x, z, mu, G, N, p);        
        
	    p5M_kf42(kf4, bx1, zR6, p, q);

	    p5M_po1(po1, kf4, zR6, y2etilde, p, q);

	    p5M_zR6(zR6, kf4, y2etilde, po1, p, q);
			    
	    p5M_bx12(bx1, zR6, kf4, y2etilde, p, q);

		ja2 = 0.0; for(i=0;i<p;i++) ja2 += log(bx1[i]);

	    q2h = p5M_det_sigma_NEW2(zR6, bx1, ja2, p, q);
	    
		log_c = (p/2.0)*log(2.0*M_PI) + 0.5*q2h; 
   
		p5M_z2(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
 
        stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	
	    it++; 
    }
	/*printf("ll=%f\n",l[it-1]);*/
    
	/* no of parameters {pq - q(q - 1)/2} + p */
    a = G-1 + G*p + p*q - q*(q-1)/2 + p;
    
    /* BIC = 2log_likelihood - mlog n; */
    bic = 2.0*l[it-1] - a*log(N);
    
    /* printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it-1); */   
    
    /* Deallocate memory */
    my_free(zR6,p); my_free(mu,G); my_free(w,G); free(n); free(det); 
    my_free(kf4,q); my_free(po1,q); my_free(y2etilde,p); free(l); 
    free(at); free(pi);   
	return bic;
}

double claecm3(double **z, double **x, int q, int p, int G, int N, double *smf, double *bx1, double tol){

	double **zR6, **mu, bic, ***kf4, ***po1, ***y2e, 
	       *l, *at, *pi, *n, *ja2, *q2h, *log_det,
			**v, *max_v, a; 
	int g,it=0,stop=0;

	/* printf("G=%d \t q=%d\n",G,q);*/

    my_alloc_vec(&log_det,G);
    my_alloc_vec(&ja2,G);
    my_alloc_vec(&q2h,G);
   	my_alloc_vec(&pi,G);
    my_alloc_vec(&n,G);
    my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
    my_alloc_3d(&y2e,G,p,p);		
    my_alloc(&zR6,p,q);		
    my_alloc_3d(&kf4,G,q,p);		
    my_alloc_3d(&po1,G,q,q);
    my_alloc(&mu,G,p);
	my_alloc_vec(&max_v,N);
   	my_alloc(&v,N,G);		
        
	get_data(smf,zR6,p,q);
	
    while(stop==0){            
           
        p5M_n(n, z, G, N);
	
	    p5M_pi(pi, n, G, N);
	
	    updau8m(mu, n, x, z, G, N, p); 
                          
		if(it>0){ 
			p5M_z3(v, x, z, zR6, bx1, mu, pi, max_v, log_det, N, G, p, q);
        }

	    p5M_sg(y2e, x, z, mu, n, p, G, N);
	    
	    for(g=0;g<G;g++) p5M_kf41(kf4[g], bx1[g], zR6, p, q);    
	
    	for (g=0;g<G;g++) p5M_po1(po1[g], kf4[g], zR6, y2e[g], p, q);
        
	    p5M_zR62(zR6, kf4, y2e, po1, n, bx1, p, q, G);
	    
        for(g=0;g<G;g++) bx1[g] = p5M_bx13(zR6, kf4[g], y2e[g], po1[g], p, q);
    
        for(g=0;g<G;g++){
	    	ja2[g] = p*log(bx1[g]);
			q2h[g] = p5M_det_sigma_NEW(zR6, bx1[g], ja2[g], p, q);
			log_det[g] = (p/2.0)*log(2*M_PI) + 0.5*q2h[g];
		} 
	    
		p5M_z3(v, x, z, zR6, bx1, mu, pi, max_v, log_det, N, G, p, q);

	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
		it++;
    }
	
	/*printf("ll=%f\n",l[it-1]);*/
	
   	a = G-1 + G*p+ p*q - q*(q-1)/2 + G;
    
   	bic = 2.0*l[it-1] - a*log(N);
    
   	/* printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it); */    
    
   	/* Deallocate memory */
   	my_free(zR6,p); my_free(mu,G); my_free(v,G); free(n); free(log_det); free(max_v); 
	my_free_3d(kf4,G,q); my_free_3d(po1,G,q); my_free_3d(y2e,G,p); 
	free(l); free(at); free(pi); free(ja2);
    
	return bic;
}

double claecm4(double **z, double **x, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol){

	double **zR6, **mu, bic, **w, ***kf4, ***po1, ***y2e, *l, *at, *pi, *n,
	       **bx1, *ja2, *q2h, *log_c, **v, *max_v; 
	int it=0,stop=0,paras, g, j;
	
   	/* Allocate memory - matrices and vectors */ 
   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&y2e,G,p,p);		
   	my_alloc(&zR6,p,q);
   	my_alloc_3d(&kf4,G,q,p);
   	my_alloc_3d(&po1,G,q,q);
   	my_alloc(&mu,G,p);
   	my_alloc(&w,G,N);
	my_alloc(&bx1,G,p);
	my_alloc_vec(&ja2,G);
	my_alloc_vec(&q2h,G);
	my_alloc_vec(&log_c,G);
	
	get_data(bx1_vec,bx1,G,p);
	get_data(smf,zR6,p,q);
        
   	while(stop==0){            
   
	    p5M_n(n, z, G, N);
	
	    p5M_pi(pi, n, G, N);
	
	    updau8m(mu, n, x, z, G, N, p); 

        if(it>0){ 
			p5M_z4(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
		} 
 
	    p5M_sg(y2e,x, z, mu, n, p, G, N);
        
	    for(g=0;g<G;g++) p5M_kf42(kf4[g], bx1[g], zR6, p, q);
        
   	    for(g=0;g<G;g++) p5M_po1(po1[g], kf4[g], zR6, y2e[g], p, q); 
    
	    p5M_zR6_cuu(zR6, kf4, y2e, po1, n, bx1, p, q, G);
	 
        p5M_bx1_cuu(bx1, zR6, kf4, y2e, po1, p, q, G);

	    for (g=0;g<G;g++){
			ja2[g]=0.0;
			for(j=0;j<p;j++) ja2[g] += log(bx1[g][j]);
	    }
	
	    for(g=0;g<G;g++) q2h[g] = p5M_det_sigma_NEW2(zR6, bx1[g], ja2[g], p, q);

		for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*q2h[g];

        p5M_z4(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
            
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++;
	}
	
	/* no of parameters {pq - q(q - 1)/2} + Gp */
	paras = G-1 + G*p + p*q - q*(q-1)/2 + G*p;
    
   	/* BIC = 2log_likelihood - mlog n; */
   	bic = 2.0*l[it-1] - paras*log(N);
    
	/*if((G==3)&&(q==4)) printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it);*/
    
    /* Deallocate memory */
    my_free(zR6,p); my_free(mu,G); my_free(w,G); free(n); my_free_3d(kf4,G,q); 
	my_free_3d(po1,G,q); my_free_3d(y2e,G,p); free(l); free(at); free(pi);
	free(q2h); free(log_c); free(ja2);
	return bic;
} 

double claecm5(double **z, double **x, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol){

	double ***zR6, **mu, bic, **w, ***kf4, ***po1, **v, *max_v, ***y2e, *l, *at, *pi, 
	*n, *q2h, *log_c, ja2, bx1; 
	int it=0,stop=0,paras, g,i;

	/*Rprintf("G=%d \t q=%d\n",G,q);*/

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
   	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&q2h,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&y2e,G,p,p);		
   	my_alloc_3d(&zR6,G,p,q);
   	my_alloc_3d(&kf4,G,q,p);
   	my_alloc_3d(&po1,G,q,q);
   	my_alloc(&mu,G,p);
   	my_alloc(&w,G,N);
	
	bx1=*bx1_vec;
	get_data2(smf,zR6,G,p,q);
       
   	while(stop==0){            
    
	    p5M_n(n, z, G, N);
	
	    p5M_pi(pi, n, G, N);
	
	    updau8m(mu, n, x, z, G, N, p); 

		if(it>0){ 
			p5M_z5(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
        }
 
	    p5M_sg(y2e,x, z, mu, n, p, G, N);
        
	    for(g=0;g<G;g++) p5M_kf41(kf4[g], bx1, zR6[g], p, q);
        
   	    for(g=0;g<G;g++) p5M_po1(po1[g], kf4[g], zR6[g], y2e[g], p, q); 
    
        for (g=0;g<G;g++) p5M_zR6(zR6[g], kf4[g], y2e[g], po1[g], p, q);
	 
	    bx1 = p5M_bx1_ucc(zR6, kf4, y2e, p, q, pi, G);
   
		ja2 = 0.0;
		for(i=0;i<p;i++) ja2 += log(bx1);
	    
	    for(g=0;g<G;g++) q2h[g] = p5M_det_sigma_NEW(zR6[g], bx1, ja2, p, q);
	    
		for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*q2h[g];
       
		p5M_z5(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
            
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++; 
	}
	/*printf("ll=%f\n",l[it-1]);*/

	paras = G-1 + G*p + G*(p*q - q*(q-1)/2) + 1;
    
   	bic = 2.0*l[it-1] - paras*log(N);
    
	/* printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it); */   
    
   	my_free_3d(zR6,G,p); my_free(mu,G); my_free(w,G); free(n); 
	my_free_3d(kf4,G,q); my_free_3d(po1,G,q); my_free_3d(y2e,G,p); 
	free(l); free(at); free(pi);free(q2h); free(log_c);
	return bic;
}

double claecm6(double **z, double **x, int q, int p, int G, int N, double *smf, double *bx1, double tol){

	double ***zR6, **mu, bic, ***kf4, ***po1,**v,*max_v, ***y2e, *l, *at, *pi, 
	*n, *q2h,*log_c, ja2; 
	int it=0,stop=0,paras, g, j;

	/* printf("G=%d \t q=%d\n",G,q);*/

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
	my_alloc_vec(&q2h,G);
	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&y2e,G,p,p);		
   	my_alloc_3d(&zR6,G,p,q);
   	my_alloc_3d(&kf4,G,q,p);
   	my_alloc_3d(&po1,G,q,q);
   	my_alloc(&mu,G,p);
	
	get_data2(smf,zR6,G,p,q);
  
   	while(stop==0){            
    
	    p5M_n(n, z, G, N);
	
	    p5M_pi(pi, n, G, N);
	
	    updau8m(mu, n, x, z, G, N, p); 

        if(it>0){ 
			p5M_z6(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
		} 

	    p5M_sg(y2e,x, z, mu, n, p, G, N);
        
	    for(g=0;g<G;g++) p5M_kf42(kf4[g], bx1, zR6[g], p, q);
        
   	    for(g=0;g<G;g++) p5M_po1(po1[g], kf4[g], zR6[g], y2e[g], p, q); 
    
        for (g=0;g<G;g++) p5M_zR6(zR6[g], kf4[g], y2e[g], po1[g], p, q);
	 
	    p5M_bx1_ucu(bx1, zR6, kf4, y2e, p, q, pi, G);
                    
		ja2 = 0.0;
		for(j=0;j<p;j++) ja2 += log(bx1[j]);
	    
	    for(g=0;g<G;g++) q2h[g] = p5M_det_sigma_NEW2(zR6[g], bx1, ja2, p, q);
	    
		for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*q2h[g];

        p5M_z6(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
            
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++; 
	}
	/*printf("ll=%f\n",l[it-1]);*/
	
	paras = G-1 + G*p + G*(p*q - q*(q-1)/2) + p;
    
   	bic = 2.0*l[it-1] - paras*log(N);
    
	/* printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it); */   
    
   	my_free_3d(zR6,G,p); my_free(mu,G); my_free(v,N); free(n);free(max_v); 
	my_free_3d(kf4,G,q); my_free_3d(po1,G,q); my_free_3d(y2e,G,p); 
	free(l); free(at); free(pi); free(q2h); free(log_c);
	return bic;
}

double claecm7(double **z, double **x, int q, int p, int G, int N, double *smf, double *bx1, double tol){

	double ***zR6, **mu, bic, ***kf4, ***po1, **v, *max_v, ***y2e, *l, *at, *pi, 
		*n, *ja2, *q2h, *log_c; 
	int it=0,stop=0,paras, g;

	/* printf("G=%d \t q=%d\n",G,q); */

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
	my_alloc_vec(&ja2,G);
	my_alloc_vec(&q2h,G);
	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&y2e,G,p,p);		
   	my_alloc_3d(&zR6,G,p,q);
   	my_alloc_3d(&kf4,G,q,p);
   	my_alloc_3d(&po1,G,q,q);
   	my_alloc(&mu,G,p);
	
	get_data2(smf,zR6,G,p,q);

   	while(stop==0){            
    
	    p5M_n(n, z, G, N);
	
	    p5M_pi(pi, n, G, N);
	
	    updau8m(mu, n, x, z, G, N, p); 

        if(it>0){ 
			p5M_z7(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
		}

	    p5M_sg(y2e,x, z, mu, n, p, G, N);
        
	    for(g=0;g<G;g++) p5M_kf41(kf4[g], bx1[g], zR6[g], p, q);
        
   	    for(g=0;g<G;g++) p5M_po1(po1[g], kf4[g], zR6[g], y2e[g], p, q); 
    
        for (g=0;g<G;g++) p5M_zR6(zR6[g], kf4[g], y2e[g], po1[g], p, q);
	 
        for (g=0;g<G;g++) bx1[g]=p5M_bx1(zR6[g], kf4[g], y2e[g], p, q); 
           
	    for(g=0;g<G;g++) ja2[g] = p*log(bx1[g]);
	    
	    for(g=0;g<G;g++) q2h[g] = p5M_det_sigma_NEW(zR6[g], bx1[g], ja2[g], p, q);
	   
		for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*q2h[g];

		p5M_z7(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
            
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++;
	}
	/*printf("ll=%f\n",l[it-1]);*/

	paras = G-1 + G*p + G*(p*q - q*(q-1)/2) + G;
    
   	bic = 2.0*l[it-1] - paras*log(N);
    
	/* printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it); */   
    
    my_free_3d(zR6,G,p); my_free(mu,G); my_free(v,N); free(n); 
	my_free_3d(kf4,G,q); free(ja2); my_free_3d(po1,G,q); 
	my_free_3d(y2e,G,p); free(l); free(at); free(pi);
	free(q2h); free(log_c);
	
	return bic;
}

double claecm8(double **z, double **x, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol){

	double ***zR6, **mu, bic, ***kf4, ***po1, **v, *max_v, ***y2e, *l, *at, *pi, 
	*n, **bx1, *log_c,*q2h,*ja2; 
	int it=0,stop=0,paras, g, j;

	/* printf("G=%d \t q=%d\n",G,q); */

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&y2e,G,p,p);		
   	my_alloc_3d(&zR6,G,p,q);
   	my_alloc_3d(&kf4,G,q,p);
   	my_alloc_3d(&po1,G,q,q);
   	my_alloc(&mu,G,p);
	my_alloc(&bx1,G,p);
	my_alloc_vec(&ja2,G);
	my_alloc_vec(&q2h,G);
	my_alloc_vec(&log_c,G);
	
	get_data(bx1_vec,bx1,G,p);
	get_data2(smf,zR6,G,p,q);
	
   	while(stop==0){            
    
	    p5M_n(n, z, G, N);
	
	    p5M_pi(pi, n, G, N);
	
	    updau8m(mu, n, x, z, G, N, p);

        if(it>0){ 
			p5M_z8(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
		}	
       
	    p5M_sg(y2e,x, z, mu, n, p, G, N);
	    
	    for(g=0;g<G;g++) p5M_kf42(kf4[g], bx1[g], zR6[g], p, q);
        
   	    for(g=0;g<G;g++) p5M_po1(po1[g], kf4[g], zR6[g], y2e[g], p, q); 
    
        for (g=0;g<G;g++) p5M_zR6(zR6[g], kf4[g], y2e[g], po1[g], p, q);
	 
        for (g=0;g<G;g++) p5M_bx12(bx1[g], zR6[g], kf4[g], y2e[g], p, q); 
           
	    for(g=0;g<G;g++){
		    ja2[g]=0.0;
		    for(j=0;j<p;j++) ja2[g] += log(bx1[g][j]);
	    }
	    
		for(g=0;g<G;g++) q2h[g] = p5M_det_sigma_NEW2(zR6[g], bx1[g], ja2[g], p, q);
	    	    
		for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*q2h[g];

		p5M_z8(v, x, z, zR6, bx1, mu, pi, max_v, log_c, N, G, p, q);
            
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
	    it++; 
	}
	/*printf("ll=%f\n",l[it-1]);*/
	paras = G-1 + G*p + G*(p*q - q*(q-1)/2) + G*p;
    
   	bic = 2.0*l[it-1] - paras*log(N);
    
	/* printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it); */   
    
   	my_free_3d(zR6,G,p); my_free(mu,G); my_free(v,N); free(n); 
	my_free_3d(kf4,G,q); free(ja2); my_free_3d(po1,G,q); 
	my_free_3d(y2e,G,p); free(l); free(at); free(pi);
	free(q2h); free(log_c); my_free(bx1,G); free(max_v);
	return bic;
}

double claecm9(double **z, double **x, int q, int p, int G, int N, double *smf, double *dsw, double tol){

	double **zR6, **mu, bic, ***kf4, ***po1, ***y2e, **v, *max_v, *l, *at, *pi, 
	*n, *g9p, *log_c, *ja2, *q2h, a, *bx1; 
	int g,i,it=0,stop=0;
	       
	/* printf("G=%d \t q=%d\n",G,q);*/

   	/* Allocate memory - matrices and vectors */ 
   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
	my_alloc_vec(&ja2,G);
	my_alloc_vec(&q2h,G);
	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&y2e,G,p,p);		
   	my_alloc(&zR6,p,q);		
   	my_alloc_3d(&kf4,G,q,p);		
   	my_alloc_3d(&po1,G,q,q);
   	my_alloc(&mu,G,p);
	my_alloc_vec(&g9p,p);
	my_alloc_vec(&bx1,p);
	
	get_data(smf,zR6,p,q);
	for(i=0;i<p;i++) g9p[i] = 1.0; 
	
    while(stop==0){            
           
        p5M_n(n, z, G, N);
	
	    p5M_pi(pi, n, G, N);
	
	    updau8m(mu, n, x, z, G, N, p); 
                          
        if(it>0){ 
			p5M_z9(v, x, z, zR6, dsw, g9p, mu, pi, max_v, log_c, N, G,p,q);
		}	
		
	    p5M_sg(y2e, x, z, mu, n, p, G, N);

	    for(g=0;g<G;g++){ 
			for (i=0;i<p;i++) bx1[i] = dsw[g]*g9p[i];
			p5M_kf42(kf4[g], bx1, zR6, p, q);    
		}	

    	for (g=0;g<G;g++) p5M_po1(po1[g], kf4[g], zR6, y2e[g], p, q);
        
	    p5M_zR62(zR6, kf4, y2e, po1, n, dsw, p, q, G);
	    
        for(g=0;g<G;g++) dsw[g] = p5M_dsw(zR6, g9p, kf4[g], y2e[g], po1[g], p, q);

	    p5M_g9p(g9p, zR6, dsw, kf4, y2e, po1, n, p, q, N, G); 
      
        for(g=0;g<G;g++){
			for (i=0;i<p;i++) bx1[i] = dsw[g]*g9p[i];
	    	ja2[g] = p*log(dsw[g]);
			q2h[g] = p5M_det_sigma_NEW2(zR6, bx1, ja2[g], p, q);
			log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*q2h[g];
		}

	    p5M_z9(v, x, z, zR6, dsw, g9p, mu, pi, max_v, log_c, N, G,p,q);
		
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
        it++;  
    }
	/*printf("ll=%f\n",l[it-1]);*/

    /* no of parameters */
    a = G-1 + G*p+ p*q - q*(q-1)/2 + G + (p-1);
    
    /* BIC = 2log_likelihood - mlog n; */
    bic = 2.0*l[it-1] - a*log(N);
    
    /* printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it); */
    
    /* Deallocate memory */
    my_free(zR6,p); my_free(mu,G); my_free(v,N); free(n); free(log_c); 
	my_free_3d(kf4,G,q); my_free_3d(po1,G,q); my_free_3d(y2e,G,p); 
	free(l); free(at); free(pi); free(ja2); free(g9p); free(q2h);
	
	return bic;
}

double claecm10(double **z, double **x, int q, int p, int G, int N, double *smf, double *dsw, double tol){

	double ***zR6, **mu, bic, ***kf4, ***po1, ***y2e,**v, *max_v, *l, *at, *pi, 
	*n, *g9p, *q2h, *log_c, *ja2, a, *bx1; 
	int g,i,it=0,stop=0;
	       
	/* printf("G=%d \t q=%d\n",G,q);*/

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
	my_alloc_vec(&ja2,G);
	my_alloc_vec(&q2h,G);
	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&y2e,G,p,p);		
   	my_alloc_3d(&zR6,G,p,q);		
   	my_alloc_3d(&kf4,G,q,p);		
   	my_alloc_3d(&po1,G,q,q);
   	my_alloc(&mu,G,p);
	my_alloc_vec(&g9p,p);
	my_alloc_vec(&bx1,p);
	
	get_data2(smf,zR6,G,p,q);
	for(i=0;i<p;i++) g9p[i] = 1.0; 
	
   	while(stop==0){            
           
        p5M_n(n, z, G, N);
	
	    p5M_pi(pi, n, G, N);
	
	    updau8m(mu, n, x, z, G, N, p); 
                          
        if(it>0){ 
			p5M_z10(v, x, z, zR6, dsw, g9p, mu, pi, max_v, log_c, N, G,p,q);
		}	
		
	    p5M_sg(y2e, x, z, mu, n, p, G, N);

	    for(g=0;g<G;g++){ 
			for (i=0;i<p;i++) bx1[i] = dsw[g]*g9p[i];
			p5M_kf42(kf4[g], bx1, zR6[g], p, q);    
		}
	
    	for (g=0;g<G;g++) p5M_po1(po1[g], kf4[g], zR6[g], y2e[g], p, q);
        
	    for(g=0;g<G;g++) p5M_zR6(zR6[g], kf4[g], y2e[g], po1[g], p, q);
	    
        for(g=0;g<G;g++) dsw[g] = p5M_dsw2(zR6[g], g9p, kf4[g], y2e[g], p, q);
	   	
	    p5M_g9p2(g9p, zR6, dsw, kf4, y2e, po1, n, p, q, N, G); 
	    
        for(g=0;g<G;g++){
			for (i=0;i<p;i++) bx1[i] = dsw[g]*g9p[i];
			ja2[g] = p*log(dsw[g]);
	    	q2h[g] = p5M_det_sigma_NEW2(zR6[g], bx1, ja2[g], p, q);
			log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*q2h[g];
		}

	    p5M_z10(v, x, z, zR6, dsw, g9p, mu, pi, max_v, log_c, N, G,p,q);
		
	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
        it++;
   	}
	/*printf("ll=%f\n",l[it-1]);*/
   
   	a = G-1 + G*p+ G*(p*q - q*(q-1)/2) + G + (p-1);
    
   	bic = 2.0*l[it-1] - a*log(N);
    
   	/* printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it); */
    
   	my_free_3d(zR6,G,p); my_free(mu,G); my_free(v,N); free(n); free(log_c); 
	my_free_3d(kf4,G,q); my_free_3d(po1,G,q); my_free_3d(y2e,G,p); 
	free(l); free(at); free(pi); free(ja2); free(g9p); free(q2h);

	return bic;
}

double claecm11(double **z, double **x, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol){

	double **zR6, **mu, bic, ***kf4, ***po1, ***y2e, **v, *max_v, *l, *at, *pi, 
	*n, **g9p, *log_c, *q2h, ja2, a, dsw, *bx1; 
	int g,i,it=0,stop=0;

	/* printf("G=%d \t q=%d\n",G,q);*/

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&y2e,G,p,p);		
   	my_alloc(&zR6,p,q);		
   	my_alloc_3d(&kf4,G,q,p);		
	my_alloc_3d(&po1,G,q,q);
   	my_alloc(&mu,G,p);
	my_alloc(&g9p,G,p);
	my_alloc_vec(&q2h,G);
	my_alloc_vec(&log_c,G);
	my_alloc_vec(&bx1,p);
	
	dsw = *bx1_vec;
	get_data(smf,zR6,p,q);
	
	for(g=0;g<G;g++)
		for(i=0;i<p;i++) 
			g9p[g][i] = 1.0; 
	
   	while(stop==0){            
           
        p5M_n(n, z, G, N);
	
	    p5M_pi(pi, n, G, N);
	
	    updau8m(mu, n, x, z, G, N, p); 
                          
        if(it>0){ 
			p5M_z11(v, x, z, zR6, dsw, g9p, mu, pi, max_v, log_c, N, G,p,q);
		}	
		
	    p5M_sg(y2e, x, z, mu, n, p, G, N);

	    for(g=0;g<G;g++){ 
			for (i=0;i<p;i++) bx1[i] = dsw*g9p[g][i];
			p5M_kf42(kf4[g], bx1, zR6, p, q);    
		}	
	
    	for (g=0;g<G;g++) p5M_po1(po1[g], kf4[g], zR6, y2e[g], p, q);
        
	    p5M_zR6_cuu(zR6, kf4, y2e, po1, n, g9p, p, q, G);
	    
        dsw =0.0;
	    for(g=0;g<G;g++) dsw += pi[g]*p5M_dsw(zR6, g9p[g], kf4[g], y2e[g], po1[g], p, q);
 	   	
	    for(g=0;g<G;g++) p5M_g9p3(g9p[g], zR6, dsw, kf4[g], y2e[g], po1[g], n[g], p, q); 
		
		ja2 = p*log(dsw);

        for(g=0;g<G;g++){
			for (i=0;i<p;i++) bx1[i] = dsw*g9p[g][i];
	    	q2h[g] = p5M_det_sigma_NEW2(zR6, bx1, ja2, p, q);
			log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*q2h[g];
		}

	    p5M_z11(v, x, z, zR6, dsw, g9p, mu, pi, max_v, log_c, N, G,p,q);

	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
        it++;
    }
	/*printf("ll=%f\n",l[it-1]);*/

   	a = G-1 + G*p+ p*q - q*(q-1)/2 + 1 + G*(p-1);
    
   	bic = 2.0*l[it-1] - a*log(N);
    
   	/* printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it); */
    
   	my_free(zR6,p); my_free(mu,G); my_free(v,N); free(n); free(log_c); my_free_3d(kf4,G,q); my_free_3d(po1,G,q); 
	my_free_3d(y2e,G,p); free(l); free(at); free(pi); my_free(g9p,G); free(q2h); 
	
	return bic;
} 

double claecm12(double **z, double **x, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol){

	double ***zR6, **mu, bic, ***kf4, ***po1, ***y2e, **v, *max_v, *l, *at, *pi, 
	*n, **g9p, *q2h, *log_c, ja2, a, dsw, *bx1; 
	int g,i,it=0,stop=0;
	       
	/* printf("G=%d \t q=%d\n",G,q);*/

   	my_alloc_vec(&max_v,N);
	my_alloc(&v, N, G);
	my_alloc_vec(&q2h,G);
	my_alloc_vec(&log_c,G);
   	my_alloc_vec(&pi,G);
   	my_alloc_vec(&n,G);
   	my_alloc_vec(&at,150000);
	my_alloc_vec(&l,150000);
   	my_alloc_3d(&y2e,G,p,p);		
   	my_alloc_3d(&zR6,G,p,q);		
   	my_alloc_3d(&kf4,G,q,p);		
	my_alloc_3d(&po1,G,q,q);
   	my_alloc(&mu,G,p);
	my_alloc(&g9p,G,p);
	my_alloc_vec(&bx1,p);
	
	dsw = *bx1_vec;
	get_data2(smf,zR6,G,p,q);
		
	for(g=0;g<G;g++)
		for(i=0;i<p;i++) 
			g9p[g][i] = 1.0; 
	
   	while(stop==0){            
           
		p5M_n(n, z, G, N);
	
	    p5M_pi(pi, n, G, N);
	
	    updau8m(mu, n, x, z, G, N, p); 
                          
	    if(it>0){
			p5M_z12(v, x, z, zR6, dsw, g9p, mu, pi, max_v, log_c, N, G,p,q);
		}	
		
	    p5M_sg(y2e, x, z, mu, n, p, G, N);
	    
	    for(g=0;g<G;g++){ 
			for (i=0;i<p;i++) bx1[i] = dsw*g9p[g][i];
			p5M_kf42(kf4[g], bx1, zR6[g], p, q);    
		}	
	
    	for (g=0;g<G;g++) p5M_po1(po1[g], kf4[g], zR6[g], y2e[g], p, q);
        
    	for (g=0;g<G;g++) p5M_zR6(zR6[g], kf4[g], y2e[g], po1[g], p, q);
		    
        dsw =0.0;
	    for(g=0;g<G;g++) dsw += pi[g]*p5M_dsw2(zR6[g], g9p[g], kf4[g], y2e[g], p, q);
	    
	    for(g=0;g<G;g++) p5M_g9p3(g9p[g], zR6[g], dsw, kf4[g], y2e[g], po1[g], n[g], p, q); 
	    
		ja2 = p*log(dsw);        

		for(g=0;g<G;g++){
			for (i=0;i<p;i++) bx1[i] = dsw*g9p[g][i];
			q2h[g] = p5M_det_sigma_NEW2(zR6[g], bx1, ja2, p, q);
			log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*q2h[g];
		}

	    p5M_z12(v, x, z, zR6, dsw, g9p, mu, pi, max_v, log_c, N, G,p,q);

	    stop = convergtest_NEW(l, at, max_v, v, N, it, G, tol);
        it++;
    }
	/*printf("ll=%f\n",l[it-1]);*/
   
   	a = G-1 + G*p+ G*(p*q - q*(q-1)/2) + 1+ G*(p-1);
    
   	bic = 2.0*l[it-1] - a*log(N);
    
   	/* printf("BIC[G=%d][q=%d]=%f\t iteration %d\n",G,q,bic,it); */
    
   	my_free_3d(zR6,G,p); my_free(mu,G); my_free(v,N); free(n);  
	my_free_3d(kf4,G,q); my_free_3d(po1,G,q); my_free_3d(y2e,G,p); 
	free(l); free(at); free(pi); my_free(g9p,G); free(log_c);
	free(q2h);

	return bic;
}
