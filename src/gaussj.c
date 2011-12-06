#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include<R.h>
#include "functions.h"

void generate_identity(int N, double **c21){
    int i, j;
    for (i = 0; i < N; i ++)
        for (j = 0; j < N; j ++){
            c21[i][j] = 0;
            if(i==j) c21[i][i]=1;
        }
}

void GaussJordan(int N, double **c21, double **a8e, double *det){
    int c, r, r_max, j, sign;
    double temp, v_v, v_max, q7j = 0.0;
    generate_identity(N, a8e); det[0]=1; sign=0;
    
    for (c = 0; c < N; c++){
        r_max = c;
        v_max = fabs(c21[c][c]);
        for (r = c + 1; r < N; r++){
            v_v = fabs(c21[r][c]);
            if (v_v > v_max){
                r_max = r;
                v_max = v_v;
            }
        } 
        if (r_max != c){
            for (j = c; j < N; j++){
                temp = c21[c][j];
                c21[c][j] = c21[r_max][j];
                c21[r_max][j] = temp;
            } 
            for (j = 0; j < N; j++){
                temp = a8e[c][j];
                a8e[c][j] = a8e[r_max][j];
                a8e[r_max][j] = temp;
            }
            sign++;
        }       
        q7j = c21[c][c];
        
        det[0]*= q7j;
        for (j = c; j < N; j++){
            c21[c][j] /= q7j;
        } 
        for (j = 0; j < N; j++){
            a8e[c][j] /= q7j;
        } 
        for (r = c + 1; r < N; r++){
            q7j = c21[r][c];
            for (j = c; j < N; j++)
                c21[r][j] -= q7j * c21[c][j];
            for (j = 0; j < N; j++)
                a8e[r][j] -= q7j * a8e[c][j];
        }
    }
        
    if(sign%2 != 0) det[0] *= -1; 
    
    for (c = 1; c < N; c++){
        for (j = 0; j < N-c; j++){
            q7j = c21[j][N-c];
            for (r = 0; r < N; r++){
                c21[j][r] =c21[j][r] - q7j * c21[N-c][r];
                a8e[j][r] =a8e[j][r] - q7j * a8e[N-c][r];
            }
        }
    }
}

/* Function to find mx R = AB */
void mx_mult(int m, int n, int q, double **a, double **b, double **r){
    int i,j,k;
    
    if(((m>40)&&(q>40))||(n>40)){
        /*Rprintf("Computing Strassen's\n");*/
        str_mx_mult(m, n, q, a, 0, 0, b, 0, 0, r);
    } else {
        for(i=0; i<m; i++)
            for(j=0; j<q; j++) {
                r[i][j]=0.0;
                for(k=0; k<n; k++)
                    r[i][j] += a[i][k]*b[k][j];
            }
    }
}

/* NEW Function to find the DIAGONAL of a mx R = AB */
void mx_mult_diag1(int m, int n, double **a, double **b, double *r){
    int i,k;
    for(i=0; i<m; i++){
			r[i]=0.0;
            for(k=0; k<n; k++)
                r[i] += a[i][k]*b[k][i];
	}
}

/* Function to find the DIAGONAL of a mx R = AB */
void mx_mult_diag(int m, int n, double **a, double **b, double **r){
    int i,k;
    for(i=0; i<m; i++){
			r[i][i]=0.0;
            for(k=0; k<n; k++)
                r[i][i] += a[i][k]*b[k][i];
	}
}

/* Function to find r = bA */
void vec_mx_mult(int n, int q, double *a, double **b, double *r){
	int j,k;
        for(j=0; j<q; j++){
            r[j]=0.0;
            	for(k=0; k<n; k++)
                    r[j] += a[k]*b[k][j];
        }
}
/* Function to find r = Ab */
void mx_vec_mult(int n, int q, double *a, double **b, double *r){
	int j,k;
        for(j=0; j<q; j++){
            r[j]=0.0;
            	for(k=0; k<n; k++)
                    r[j] += a[k]*b[j][k];
        }
}

/* Function to find the transpose of an m*n c21 A; R=A' */
void mx_trans(int m, int n, double **a, double **r){
    int i,j;
    for(i=0; i<n; i++)
        for(j=0; j<m; j++)
            r[i][j]=a[j][i];
}

int maximum(double *z){
    int j,k; k=0;
    for(j=1;j<5;j++)        
        if(z[j]>z[k]) k=j;   
    return k;
}

int maxi(double *z, int G){
    int g,k; k=0;
    
    if(G==1) return 0;
    
    for(g=1;g<G;g++)        
        if(z[g]>z[k]) k=g;   
    return k;
}

/* Function to allocated memory to a row x col array */
void mx_alloc(int row, int col, double **a){
    int i;
    a = (double **)calloc(row, sizeof(double*));      
    for (i=0; i<row; i++) a[i] = (double *)calloc(col, sizeof(double));    
}

/* Function to allocated memory to a row x col array */
void vec_alloc(int row, double *a){
    a = (double *)calloc(row, sizeof(double*));      
}

double maximum_array(double *array, int k){
    int i;
    double max = array[0];

    for(i=1;i<k;i++)
	if(array[i]> max) max=array[i];

    return max;
}

double** init_mx(int nr, int nc) {
    int i,j;
    double **A = (double**)malloc(sizeof(double*)*nr);
    
    for(i=0; i < nr; i++) { 
        A[i] = (double*)malloc(sizeof(double)*nc);
        for(j=0; j < nc; j++)
            A[i][j] = 0;
    }
    return A;
}

void release_mx(int nr, double **A) {
    int i;
    for(i=0; i < nr; i++)
        free(A[i]);
    free(A);
}

void std_mx_mult(int m, int n, int q, double **A, int ax, int ay, double **B, int bx, int by, double **R) {
    int i,j,k;
    for(i=0; i<m; i++)
        for(j=0; j<q; j++){
            R[i][j]=0.0;
            for(k=0; k<n; k++)
                R[i][j] += A[ax+i][ay+k]*B[bx+k][by+j];
        }
}

void str_mx_sum(int diff, int nr, int nc, double **A, int ax, int ay, double **B, int bx, int by, double **R) {
    int i, j;
    int fac = 1;
    
    if(diff) 
        fac = -1;
    for(i=0; i < nr; i++)
        for(j=0; j < nc; j++)
            R[i][j] = A[ax+i][ay+j] + fac*B[bx+i][by+j];
}

void str_mx_mult(int m, int n, int q, double **A, int ax, int ay, double **B, int bx, int by, double **R) {
    int i, j, k;
    double sum;
    double **S1, **S2, **S3, **S4;
    double **T1, **T2, **T3, **T4;
    double **P1, **P2, **P3, **P4, **P5, **P6, **P7;
    double u1, u2, u3, u4, u5, u6, u7;
    
    if(((m>40)&&(q>40))||(n>40)) {
        std_mx_mult(m, n, q, A, ax, ay, B, bx, by, R);
        return;
    }
    
    S1 = init_mx(m/2, n/2);
    S2 = init_mx(m/2, n/2);
    S3 = init_mx(m/2, n/2);
    S4 = init_mx(m/2, n/2);
    
    T1 = init_mx(n/2, q/2);
    T2 = init_mx(n/2, q/2);
    T3 = init_mx(n/2, q/2);
    T4 = init_mx(n/2, q/2);
    
    P1 = init_mx(m/2, q/2);
    P2 = init_mx(m/2, q/2);
    P3 = init_mx(m/2, q/2);
    P4 = init_mx(m/2, q/2);
    P5 = init_mx(m/2, q/2);
    P6 = init_mx(m/2, q/2);
    P7 = init_mx(m/2, q/2);
    
    str_mx_sum(0, m/2, n/2,  A, ax+m/2,   ay+0,  A, ax+m/2, ay+n/2, S1);
    str_mx_sum(1, m/2, n/2, S1,      0,      0,  A,   ax+0,   ay+0, S2);
    str_mx_sum(1, m/2, n/2,  A,   ax+0,   ay+0,  A, ax+m/2,   ay+0, S3);
    str_mx_sum(1, m/2, n/2,  A,   ax+0, ay+n/2, S2,      0,      0, S4);
    
    str_mx_sum(1, n/2, q/2, B, bx+0 ,  by+q/2, B,  bx+0,   by+0, T1);
    str_mx_sum(1, n/2, q/2, B, bx+n/2, by+q/2, T1,    0,      0, T2);
    str_mx_sum(1, n/2, q/2, B, bx+n/2, by+q/2, B,  bx+0, by+q/2, T3);
    str_mx_sum(1, n/2, q/2, B, bx+n/2,   by+0, T2,    0,      0, T4);
    
    str_mx_mult(m/2, n/2, q/2,  A,     ax,     ay,  B,     bx,      by, P1);
    str_mx_mult(m/2, n/2, q/2,  A,     ax, ay+n/2,  B, bx+n/2,      by, P2);
    str_mx_mult(m/2, n/2, q/2, S1,      0,      0, T1,      0,       0, P3);
    str_mx_mult(m/2, n/2, q/2, S2,      0,      0, T2,      0,       0, P4);
    str_mx_mult(m/2, n/2, q/2, S3,      0,      0, T3,      0,       0, P5);
    str_mx_mult(m/2, n/2, q/2, S4,      0,      0,  B, bx+n/2,  by+q/2, P6);
    str_mx_mult(m/2, n/2, q/2,  A, ax+m/2, ay+n/2, T4,      0,       0, P7);
    
    for(i=0; i < m/2; i++) {
        for(j=0; j < q/2; j++) {
            u1 = P1[i][j]+P2[i][j];  
            u2 = P1[i][j]+P4[i][j];
            u3 = u2 + P5[i][j];
            u4 = u3 + P7[i][j];
            u5 = u3 + P3[i][j];
            u6 = u2 + P3[i][j]; 
            u7 = u6 + P6[i][j];
            
            R[i][j] = u1;
            R[i][j+q/2]= u7;
            R[i+m/2][j]= u4;
            R[i+m/2][j+q/2]= u5;
        }
    }
    
    if(n%2 == 1) {
        for(i=0; i < (m/2)*2; i++)
            for(k=0; k < (q/2)*2; k++)
                R[i][k] += A[ax+i][ay+n-1]*B[bx+n-1][by+k];
    }
    
    if(m%2 == 1) {
        for(k=0; k < q; k++) {
            sum=0;
            for(j=0; j < n; j++)
                sum += A[ax+m-1][ay+j]*B[bx+j][by+k];
            R[m-1][k] = sum;
        }
    }
    
    if(q%2 == 1) {
        for(i=0; i < m; i++) {
            sum=0;
            for(j=0; j < n; j++)
                sum += A[ax+i][ay+j]*B[bx+j][by+q-1];
            R[i][q-1] = sum;
        }
    }
    
    release_mx(m/2, S1);
    release_mx(m/2, S2);
    release_mx(m/2, S3);
    release_mx(m/2, S4);
    
    release_mx(n/2, T1);
    release_mx(n/2, T2);
    release_mx(n/2, T3);
    release_mx(n/2, T4);
    
    release_mx(m/2, P1);
    release_mx(m/2, P2);
    release_mx(m/2, P3);
    release_mx(m/2, P4);
    release_mx(m/2, P5);
    release_mx(m/2, P6);
    release_mx(m/2, P7);
}
