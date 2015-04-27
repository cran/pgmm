double (*testfunc)(double *, double *, int *,int,int,int,int,double *,double *,double);
typedef typeof(testfunc) funcType;

double (*testfunc2)(double *, double *,int,int,int,int,double *,double *,double);
typedef typeof(testfunc2) funcType2;

/* Functions that dynamically allocate memory. */
void my_alloc(double ***vec, int n, int m);
void my_alloc_3d(double ****vec, int n, int m, int s);
void my_alloc_vec(double **vec, int n);
void my_free(double **vec, int n);
void my_free_3d(double ***vec, int n, int m);
void mx_alloc(int row, int col, double **a);
void release_mx(int nr, double **A);

/* AECM Algorithms*/ 
double aecm(double *z, double *x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);
double aecm2(double *z, double *x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);
double aecm3(double *z, double *x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);
double aecm4(double *z, double *x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);
double aecm5(double *z, double *x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);
double aecm6(double *z, double *x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);
double aecm7(double *z, double *x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);
double aecm8(double *z, double *x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);
double aecm9(double *z, double *x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);
double aecm10(double *z, double *x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);
double aecm11(double *z, double *x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);
double aecm12(double *z, double *x, int *cls, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);
double claecm(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);
double claecm2(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);
double claecm3(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);
double claecm4(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);
double claecm5(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);
double claecm6(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);
double claecm7(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);
double claecm8(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);
double claecm9(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);
double claecm10(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);
double claecm11(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);
double claecm12(double *z, double *x, int q, int p, int G, int N, double *lam_vec, double *psi_vec, double tol);

/* Other utility functions. */
void lambda_store(double *lam_vec, double *lambda, int p, int q);
void lambda_storeG(double *lam_vec, double **lambda, int G, int p, int q);
void get_data(double *x1, double *x, int n, int m);
void get_data2(double *x1, double **x, int G, int n, int m);
void give_data(double *x1, double *x, int n, int m);
void generate_identity(int n, double *matrix);
void GaussJordan(int n, double *MATRIX, double *INVERSE, double *det);
void mx_mult(int m, int n, int q, double *a, double *b, double *r);
void mx_mult_diag(int m, int n, double *a, double *b, double *r);
void mx_mult_diag1(int m, int n, double *a, double *b, double *r);
void vec_mx_mult(int n, int q, double *a, double *b, double *r);
void mx_trans(int m, int n, double *a, double *r);
int maximum(double *z);
void vec_alloc(int row, double *a);
double** init_mx(int nr, int nc);
void std_mx_mult(int m, int n, int q, double *A, int ax, int ay, double *B, int bx, int by, double *R);
void str_mx_sum(int diff, int nr, int nc, double *A, int ax, int ay, double *B, int bx, int by, double *R);
void str_mx_mult(int m, int n, int q, double *A, int ax, int ay, double *B, int bx, int by, double *R);

/*********** Update Z values ************************/
int update_z(double *v, double *x, double *z, double *lambda, double psi, double *mu, double *pi, double *max_v, double c, int N, int G, int p, int q);
int update_z2(double *v, double *x, double *z, double *lambda, double *psi, double *mu, double *pi, double *max_v, double c, int N, int G, int p, int q);
int update_z3(double *v, double *x, double *z, double *lambda, double *psi, double *mu, double *pi, double *max_v, double *c, int N, int G, int p, int q);
int update_z4(double *v, double *x, double *z, double *lambda, double *psi, double *mu, double *pi, double *max_v, double *c, int N, int G, int p, int q);
int update_z5(double *v, double *x, double *z, double **lambda, double psi, double *mu, double *pi, double *max_v, double *c, int N, int G, int p, int q);
int update_z6(double *v, double *x, double *z, double **lambda, double *psi, double *mu, double *pi, double *max_v, double *c, int N, int G, int p, int q);
int update_z7(double *v, double *x, double *z, double **lambda, double *psi, double *mu, double *pi, double *max_v, double *c, int N, int G, int p, int q);
int update_z8(double *v, double *x, double *z, double **lambda, double *psi, double *mu, double *pi, double *max_v, double *c, int N, int G, int p, int q);
int update_z9(double *v, double *x, double *z, double *lambda, double *omega, double *delta,double *mu, double *pi, double *max_v, double *c, int N, int G, int p, int q);
int update_z10(double *v, double *x, double *z, double **lambda, double *omega, double *delta,double *mu, double *pi, double *max_v, double *c, int N, int G, int p, int q);
int update_z11(double *v, double *x, double *z, double *lambda, double omega, double *delta,double *mu, double *pi, double *max_v, double *c, int N, int G, int p, int q);
int update_z12(double *v, double *x, double *z, double **lambda, double omega, double *delta,double *mu, double *pi, double *max_v, double *c, int N, int G, int p, int q);

/************* Woodbury Trick ******************/
double woodbury2(double *x, double *lambda, double *psi, double *mu, int p, int q);
double woodbury(double *x, double *lambda, double psi, double *mu, int p, int q);

/************* Updating Functions ***************/
void update_n(double *n, double *z, int G, int N);
void update_pi(double *pi, double *n, int G, int N);
void update_mu(double *mu, double *n, double *x, double *z, int G, int N, int p);
void update_stilde(double *sampcovtilde, double *x, double *z, double *mu, int G, int N, int p);
void update_sg(double **sg, double *x, double *z, double *mu, double *n, int p, int G, int N);
void update_beta1(double *beta, double psi, double *lambda, int p, int q);
void update_beta2(double *beta, double *PSI, double *lambda, int p, int q);
void update_theta(double *theta, double *beta, double *lambda, double *sampcovtilde, int p, int q);
void update_lambda(double *lambda, double *beta, double *s, double *theta, int p, int q);
void update_lambda2(double *lambda, double **beta, double **s, double **theta, double *n, double *Psi, int p, int q, int G);
void update_lambda_cuu(double *lambda, double **beta, double **s, double **theta, double *n, double *Psi, int p, int q, int G);
double update_psi(double *lambda, double *beta, double *sampcovtilde, int p, int q);
void update_psi2(double *psi, double *lambda, double *beta, double *sampcovtilde, int p, int q);
double update_psi3(double *lambda, double *beta, double *sampcovg, double *theta, int p, int q);
void update_psi_ucu(double *psi, double **lambda, double **beta, double **sampcov, int p, int q, double *pi, int G);
double update_psi_ucc(double **lambda, double **beta, double **sampcov, int p, int q, double *pi, int G);
void update_psi_cuu(double *psi, double *lambda, double **beta, double **sampcovg, double **theta, int p, int q, int G);
double update_det_sigma(double *lambda, double *sigma, double psi, int p, int q);
double update_det_sigma_NEW(double *lambda, double psi, double log_detpsi, int p, int q);
double update_det_sigma_NEW2(double *lambda, double *psi, double log_detpsi, int p, int q);
double update_omega(double *lambda, double *delta, double *beta, double *sampcovg, double *theta, int p, int q);
double update_omega2(double *lambda, double *delta, double *beta, double *sampcovg, int p, int q);
void update_delta(double *delta, double *lambda, double *omega, double **beta, double **sampcov, double **theta, double *n, int p, int q, int N, int G);
void update_delta2(double *delta, double **lambda, double *omega, double **beta, double **sampcov, double **theta, double *n, int p, int q, int N, int G);
void update_delta3(double *delta, double *lambda, double omega, double *beta, double *sampcov, double *theta, double n, int p, int q);

/*********** Functions to test for convergence ***********/
int convergtest_NEW(double *l, double *at, double *v_max, double *v, int N, int it, int G, double tol);

/************* Others ***************/
int maxi(double *z, int G);
extern funcType funcs[13];
extern funcType2 funcs2[13];
double maximum_array(double *array, int k);
void known_z(int *class, double *z, int N, int G);
void printmx(double *A, int r, int c);
