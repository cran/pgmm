#include <stdlib.h>
#include <stdio.h>
#include "functions.h"

/* Code for int matrices */
void my_alloc_int(int ***vec, int n, int m){
	int i;
	int **vec1;
	
	vec1= (int **)calloc(n, sizeof(int*));
	/*if(vec1==NULL){ printf("\t\tYou have failed!\n"); exit(4);}*/
	
	for (i=0; i<n; i++)
		vec1[i]= (int *)calloc(m, sizeof(int));
	
	*vec = vec1;
}

/* Code for matrices */
void my_alloc(double ***vec, int n, int m){
	int i;
	double **vec1;
	
	vec1= (double **)calloc(n, sizeof(double*));
	/*if(vec1==NULL){ printf("\t\tYou have failed!\n"); exit(4);}*/
	
	for (i=0; i<n; i++)
		vec1[i]= (double *)calloc(m, sizeof(double));
	
	*vec = vec1;
}

/* Code for vectors */
void my_alloc_vec(double **vec, int n){
	
	double *vec1;
	vec1= (double *)calloc(n, sizeof(double));
	*vec = vec1;
}

/* Code for 3-d arrays */
void my_alloc_3d(double ****vec, int n, int m, int s){
	int i,j;
	double ***vec1;
	vec1= (double ***)calloc(n, sizeof(double**));
	for (i=0; i<n; i++)
		vec1[i]= (double **)calloc(m, sizeof(double*));
	for (i=0; i<n; i++)
		for (j=0; j<m; j++)
			vec1[i][j]= (double *)calloc(s, sizeof(double));
	*vec = vec1;
}

/* Code for matrices */ 
void my_free(double **vec, int n){
	int i;
	
	for (i=0; i<n; i++)
		free(vec[i]);

	free(vec);
}
/* Code for 3-d arrays */
void my_free_3d(double ***vec, int n, int m){
	int i;
	
	for (i=0; i<n; i++)
		my_free(vec[i],m);

	free(vec);
}


