double (*testfunc)(double **, double **, int *,int,int,int,int,double *,double *,double);
typedef typeof(testfunc) funcType;

double (*testfunc2)(double **, double **,int,int,int,int,double *,double *,double);
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
double aecm(double **z, double **x, int *cls, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);
double aecm2(double **z, double **x, int *cls, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);
double aecm3(double **z, double **x, int *cls, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);
double aecm4(double **z, double **x, int *cls, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);
double aecm5(double **z, double **x, int *cls, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);
double aecm6(double **z, double **x, int *cls, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);
double aecm7(double **z, double **x, int *cls, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);
double aecm8(double **z, double **x, int *cls, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);
double aecm9(double **z, double **x, int *cls, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);
double aecm10(double **z, double **x, int *cls, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);
double aecm11(double **z, double **x, int *cls, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);
double aecm12(double **z, double **x, int *cls, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);
double claecm(double **z, double **x, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);
double claecm2(double **z, double **x, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);
double claecm3(double **z, double **x, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);
double claecm4(double **z, double **x, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);
double claecm5(double **z, double **x, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);
double claecm6(double **z, double **x, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);
double claecm7(double **z, double **x, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);
double claecm8(double **z, double **x, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);
double claecm9(double **z, double **x, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);
double claecm10(double **z, double **x, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);
double claecm11(double **z, double **x, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);
double claecm12(double **z, double **x, int q, int p, int G, int N, double *smf, double *bx1_vec, double tol);

/* Other utility functions. */
void get_data(double *x1, double **x, int n, int m);
void get_data2(double *x1, double ***x, int G, int n, int m);
void give_data(double *x1, double **x, int n, int m);
void generate_identity(int n, double **matrix);
void GaussJordan(int n, double **MATRIX, double **INVERSE, double *det);
void mx_mult(int m, int n, int q, double **a, double **b, double **r);
void mx_mult_diag(int m, int n, double **a, double **b, double **r);
void mx_mult_diag1(int m, int n, double **a, double **b, double *r);
void vec_mx_mult(int n, int q, double *a, double **b, double *r);
void mx_trans(int m, int n, double **a, double **r);
int maximum(double *z);
void vec_alloc(int row, double *a);
double** init_mx(int nr, int nc);
void std_mx_mult(int m, int n, int q, double **A, int ax, int ay, double **B, int bx, int by, double **R);
void str_mx_sum(int diff, int nr, int nc, double **A, int ax, int ay, double **B, int bx, int by, double **R);
void str_mx_mult(int m, int n, int q, double **A, int ax, int ay, double **B, int bx, int by, double **R);

/*********** p5M Z values ************************/
int p5M_z(double **v, double **x, double **z, double **zR6, double bx1, double **mu, double *pi, double *max_v, double c, int N, int G, int p, int q);
int p5M_z2(double **v, double **x, double **z, double **zR6, double *bx1, double **mu, double *pi, double *max_v, double c, int N, int G, int p, int q);
int p5M_z3(double **v, double **x, double **z, double **zR6, double *bx1, double **mu, double *pi, double *max_v, double *c, int N, int G, int p, int q);
int p5M_z4(double **v, double **x, double **z, double **zR6, double **bx1, double **mu, double *pi, double *max_v, double *c, int N, int G, int p, int q);
int p5M_z5(double **v, double **x, double **z, double ***zR6, double bx1, double **mu, double *pi, double *max_v, double *c, int N, int G, int p, int q);
int p5M_z6(double **v, double **x, double **z, double ***zR6, double *bx1, double **mu, double *pi, double *max_v, double *c, int N, int G, int p, int q);
int p5M_z7(double **v, double **x, double **z, double ***zR6, double *bx1, double **mu, double *pi, double *max_v, double *c, int N, int G, int p, int q);
int p5M_z8(double **v, double **x, double **z, double ***zR6, double **bx1, double **mu, double *pi, double *max_v, double *c, int N, int G, int p, int q);
int p5M_z9(double **v, double **x, double **z, double **zR6, double *dsw, double *g9p,double **mu, double *pi, double *max_v, double *c, int N, int G, int p, int q);
int p5M_z10(double **v, double **x, double **z, double ***zR6, double *dsw, double *g9p,double **mu, double *pi, double *max_v, double *c, int N, int G, int p, int q);
int p5M_z11(double **v, double **x, double **z, double **zR6, double dsw, double **g9p,double **mu, double *pi, double *max_v, double *c, int N, int G, int p, int q);
int p5M_z12(double **v, double **x, double **z, double ***zR6, double dsw, double **g9p,double **mu, double *pi, double *max_v, double *c, int N, int G, int p, int q);

/************* Woodbury Trick ******************/
double woodbury2(double *x, double **zR6, double *bx1, double *mu, int p, int q);
double woodbury(double *x, double **zR6, double bx1, double *mu, int p, int q);

/************* Updating Functions ***************/
void p5M_n(double *n, double **z, int G, int N);
void p5M_pi(double *pi, double *n, int G, int N);
void updau8m(double **mu, double *n, double **x, double **z, int G, int N, int p);
void p5M_stilde(double **y2etilde, double **x, double **z, double **mu, int G, int N, int p);
void p5M_sg(double ***sg, double **x, double **z, double **mu, double *n, int p, int G, int N);
void p5M_kf41(double **kf4, double bx1, double **zR6, int p, int q);
void p5M_kf42(double **kf4, double *bx1, double **zR6, int p, int q);
void p5M_po1(double **po1, double **kf4, double **zR6, double **y2etilde, int p, int q);
void p5M_zR6(double **zR6, double **kf4, double **s, double **po1, int p, int q);
void p5M_zR62(double **zR6, double ***kf4, double ***s, double ***po1, double *n, double *bx1, int p, int q, int G);
void p5M_zR6_cuu(double **zR6, double ***kf4, double ***s, double ***po1, double *n, double **bx1, int p, int q, int G);
double p5M_bx1(double **zR6, double **kf4, double **y2etilde, int p, int q);
void p5M_bx12(double *bx1, double **zR6, double **kf4, double **y2etilde, int p, int q);
double p5M_bx13(double **zR6, double **kf4, double **y2eg, double **po1, int p, int q);
void p5M_bx1_ucu(double *bx1, double ***zR6, double ***kf4, double ***y2e, int p, int q, double *pi, int G);
double p5M_bx1_ucc(double ***zR6, double ***kf4, double ***y2e, int p, int q, double *pi, int G);
void p5M_bx1_cuu(double **bx1, double **zR6, double ***kf4, double ***y2eg, double ***po1, int p, int q, int G);
double p5M_det_sigma(double **zR6, double **sigma, double bx1, int p, int q);
double p5M_det_sigma_NEW(double **zR6, double bx1, double ja2, int p, int q);
double p5M_det_sigma_NEW2(double **zR6, double *bx1, double ja2, int p, int q);
double p5M_dsw(double **zR6, double *g9p, double **kf4, double **y2eg, double **po1, int p, int q);
double p5M_dsw2(double **zR6, double *g9p, double **kf4, double **y2eg, int p, int q);
void p5M_g9p(double *g9p, double **zR6, double *dsw, double ***kf4, double ***y2e, double ***po1, double *n, int p, int q, int N, int G);
void p5M_g9p2(double *g9p, double ***zR6, double *dsw, double ***kf4, double ***y2e, double ***po1, double *n, int p, int q, int N, int G);
void p5M_g9p3(double *g9p, double **zR6, double dsw, double **kf4, double **y2e, double **po1, double n, int p, int q);

/*********** Functions to test for convergence ***********/
int convergtest_NEW(double *l, double *at, double *v_max, double **v, int N, int it, int G, double tol);

/************* Others ***************/
int maxi(double *z, int G);
extern funcType funcs[13];
extern funcType2 funcs2[13];
double maximum_array(double *array, int k);
void known_z(int *class, double **z, int N, int G);
