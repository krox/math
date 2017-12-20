module bindings.openblas;
import std.complex;
extern(C):

/*Set the number of threads on runtime.*/
void openblas_set_num_threads(int num_threads);

/*Get the number of threads on runtime.*/
int openblas_get_num_threads();

/*Get the number of physical processors (cores).*/
int openblas_get_num_procs();

/*Get the build configure on runtime.*/
const(char)* openblas_get_config();

/*Get the CPU corename on runtime.*/
const(char)* openblas_get_corename();

/* Get the parallelization type which is used by OpenBLAS */
int openblas_get_parallel();
/* OpenBLAS is compiled for sequential use */
enum OPENBLAS_SEQUENTIAL = 0;
/* OpenBLAS is compiled using normal threading model */
enum OPENBLAS_THREAD = 1;
/* OpenBLAS is compiled using OpenMP threading model */
enum OPENBLAS_OPENMP = 2;

enum
{
    CblasRowMajor=101,
    CblasColMajor=102,

    CblasNoTrans=111,
    CblasTrans=112,
    CblasConjTrans=113,
    CblasConjNoTrans=114,

    CblasUpper=121,
    CblasLower=122,

    CblasNonUnit=131,
    CblasUnit=132,

    CblasLeft=141,
    CblasRight=142,
}

float cblas_sdsdot(int n, float alpha, const(float)* x, int incx, const(float)* y, int incy);
double cblas_dsdot (int n, const(float)* x, int incx, const(float)* y, int incy);
float cblas_sdot(int n, const(float)* x, int incx, const(float)* y, int incy);
double cblas_ddot(int n, const(double)* x, int incx, const(double)* y, int incy);

Complex!float cblas_cdotu(int n, const(float)* x, int incx, const(float)* y, int incy);
Complex!float cblas_cdotc(int n, const(float)* x, int incx, const(float)* y, int incy);
Complex!double cblas_zdotu(int n, const(double)* x, int incx, const(double)* y, int incy);
Complex!double cblas_zdotc(int n, const(double)* x, int incx, const(double)* y, int incy);

void cblas_cdotu_sub(int n, const(float)* x, int incx, const(float)* y, int incy, Complex!float *ret);
void cblas_cdotc_sub(int n, const(float)* x, int incx, const(float)* y, int incy, Complex!float *ret);
void cblas_zdotu_sub(int n, const(double)* x, int incx, const(double)* y, int incy, Complex!double *ret);
void cblas_zdotc_sub(int n, const(double)* x, int incx, const(double)* y, int incy, Complex!double *ret);

float cblas_sasum (int n, const(float)* x, int incx);
double cblas_dasum (int n, const(double)* x, int incx);
float cblas_scasum(int n, const(float)* x, int incx);
double cblas_dzasum(int n, const(double)* x, int incx);

float cblas_snrm2 (int N, const(float)* X, int incX);
double cblas_dnrm2 (int N, const(double)* X, int incX);
float cblas_scnrm2(int N, const(float)* X, int incX);
double cblas_dznrm2(int N, const(double)* X, int incX);

size_t cblas_isamax(int n, const(float)* x, int incx);
size_t cblas_idamax(int n, const(double)* x, int incx);
size_t cblas_icamax(int n, const(float)* x, int incx);
size_t cblas_izamax(int n, const(double)* x, int incx);

void cblas_saxpy(int n, float alpha, const(float)* x, int incx, float *y, int incy);
void cblas_daxpy(int n, double alpha, const(double)* x, int incx, double *y, int incy);
void cblas_caxpy(int n, const(float)* alpha, const(float)* x, int incx, float *y, int incy);
void cblas_zaxpy(int n, const(double)* alpha, const(double)* x, int incx, double *y, int incy);

void cblas_scopy(int n, const(float)* x, int incx, float *y, int incy);
void cblas_dcopy(int n, const(double)* x, int incx, double *y, int incy);
void cblas_ccopy(int n, const(float)* x, int incx, float *y, int incy);
void cblas_zcopy(int n, const(double)* x, int incx, double *y, int incy);

void cblas_sswap(int n, float *x, int incx, float *y, int incy);
void cblas_dswap(int n, double *x, int incx, double *y, int incy);
void cblas_cswap(int n, float *x, int incx, float *y, int incy);
void cblas_zswap(int n, double *x, int incx, double *y, int incy);

void cblas_srot(int N, float *X, int incX, float *Y, int incY, float c, float s);
void cblas_drot(int N, double *X, int incX, double *Y, int incY, double c, double s);

void cblas_srotg(float *a, float *b, float *c, float *s);
void cblas_drotg(double *a, double *b, double *c, double *s);

void cblas_srotm(int N, float *X, int incX, float *Y, int incY, const(float)* P);
void cblas_drotm(int N, double *X, int incX, double *Y, int incY, const(double)* P);

void cblas_srotmg(float *d1, float *d2, float *b1, float b2, float *P);
void cblas_drotmg(double *d1, double *d2, double *b1, double b2, double *P);

void cblas_sscal(int N, float alpha, float *X, int incX);
void cblas_dscal(int N, double alpha, double *X, int incX);
void cblas_cscal(int N, const(float)* alpha, float *X, int incX);
void cblas_zscal(int N, const(double)* alpha, double *X, int incX);
void cblas_csscal(int N, float alpha, float *X, int incX);
void cblas_zdscal(int N, double alpha, double *X, int incX);

void cblas_sgemv(int order, int trans, int m, int n, float alpha, const(float)* a, int lda, const(float)* x, int incx, float beta, float *y, int incy);
void cblas_dgemv(int order, int trans, int m, int n, double alpha, const(double)* a, int lda, const(double)* x, int incx, double beta, double *y, int incy);
void cblas_cgemv(int order, int trans, int m, int n, const(float)* alpha, const(float)* a, int lda, const(float)* x, int incx, const(float)* beta, float *y, int incy);
void cblas_zgemv(int order, int trans, int m, int n, const(double)* alpha, const(double)* a, int lda, const(double)* x, int incx, const(double)* beta, double *y, int incy);

void cblas_sger (int order, int M, int N, float alpha, const(float)* X, int incX, const(float)* Y, int incY, float *A, int lda);
void cblas_dger (int order, int M, int N, double alpha, const(double)* X, int incX, const(double)* Y, int incY, double *A, int lda);
void cblas_cgeru(int order, int M, int N, const(float)* alpha, const(float)* X, int incX, const(float)* Y, int incY, float *A, int lda);
void cblas_cgerc(int order, int M, int N, const(float)* alpha, const(float)* X, int incX, const(float)* Y, int incY, float *A, int lda);
void cblas_zgeru(int order, int M, int N, const(double)* alpha, const(double)* X, int incX, const(double)* Y, int incY, double *A, int lda);
void cblas_zgerc(int order, int M, int N, const(double)* alpha, const(double)* X, int incX, const(double)* Y, int incY, double *A, int lda);

void cblas_strsv(int order, int Uplo, int TransA, int Diag, int N, const(float)* A, int lda, float *X, int incX);
void cblas_dtrsv(int order, int Uplo, int TransA, int Diag, int N, const(double)* A, int lda, double *X, int incX);
void cblas_ctrsv(int order, int Uplo, int TransA, int Diag, int N, const(float)* A, int lda, float *X, int incX);
void cblas_ztrsv(int order, int Uplo, int TransA, int Diag, int N, const(double)* A, int lda, double *X, int incX);

void cblas_strmv(int order, int Uplo, int TransA, int Diag, int N, const(float)* A, int lda, float *X, int incX);
void cblas_dtrmv(int order, int Uplo, int TransA, int Diag, int N, const(double)* A, int lda, double *X, int incX);
void cblas_ctrmv(int order, int Uplo, int TransA, int Diag, int N, const(float)* A, int lda, float *X, int incX);
void cblas_ztrmv(int order, int Uplo, int TransA, int Diag, int N, const(double)* A, int lda, double *X, int incX);

void cblas_ssyr(int order, int Uplo, int N, float alpha, const(float)* X, int incX, float *A, int lda);
void cblas_dsyr(int order, int Uplo, int N, double alpha, const(double)* X, int incX, double *A, int lda);
void cblas_cher(int order, int Uplo, int N, float alpha, const(float)* X, int incX, float *A, int lda);
void cblas_zher(int order, int Uplo, int N, double alpha, const(double)* X, int incX, double *A, int lda);

void cblas_ssyr2(int order, int Uplo,int N, float alpha, const(float)* X, int incX, const(float)* Y, int incY, float *A, int lda);
void cblas_dsyr2(int order, int Uplo, int N, double alpha, const(double)* X, int incX, const(double)* Y, int incY, double *A, int lda);
void cblas_cher2(int order, int Uplo, int N, const(float)* alpha, const(float)* X, int incX, const(float)* Y, int incY, float *A, int lda);
void cblas_zher2(int order, int Uplo, int N, const(double)* alpha, const(double)* X, int incX, const(double)* Y, int incY, double *A, int lda);

void cblas_sgbmv(int order, int TransA, int M, int N, int KL, int KU, float alpha, const(float)* A, int lda, const(float)* X, int incX, float beta, float *Y, int incY);
void cblas_dgbmv(int order, int TransA, int M, int N, int KL, int KU, double alpha, const(double)* A, int lda, const(double)* X, int incX, double beta, double *Y, int incY);
void cblas_cgbmv(int order, int TransA, int M, int N, int KL, int KU, const(float)* alpha, const(float)* A, int lda, const(float)* X, int incX, const(float)* beta, float *Y, int incY);
void cblas_zgbmv(int order, int TransA, int M, int N, int KL, int KU, const(double)* alpha, const(double)* A, int lda, const(double)* X, int incX, const(double)* beta, double *Y, int incY);

void cblas_ssbmv(int order, int Uplo, int N, int K, float alpha, const(float)* A, int lda, const(float)* X, int incX, float beta, float *Y, int incY);
void cblas_dsbmv(int order, int Uplo, int N, int K, double alpha, const(double)* A, int lda, const(double)* X, int incX, double beta, double *Y, int incY);

void cblas_stbmv(int order, int Uplo, int TransA, int Diag, int N, int K, const(float)* A, int lda, float *X, int incX);
void cblas_dtbmv(int order, int Uplo, int TransA, int Diag, int N, int K, const(double)* A, int lda, double *X, int incX);
void cblas_ctbmv(int order, int Uplo, int TransA, int Diag, int N, int K, const(float)* A, int lda, float *X, int incX);
void cblas_ztbmv(int order, int Uplo, int TransA, int Diag, int N, int K, const(double)* A, int lda, double *X, int incX);

void cblas_stbsv(int order, int Uplo, int TransA, int Diag, int N, int K, const(float)* A, int lda, float *X, int incX);
void cblas_dtbsv(int order, int Uplo, int TransA, int Diag, int N, int K, const(double)* A, int lda, double *X, int incX);
void cblas_ctbsv(int order, int Uplo, int TransA, int Diag, int N, int K, const(float)* A, int lda, float *X, int incX);
void cblas_ztbsv(int order, int Uplo, int TransA, int Diag, int N, int K, const(double)* A, int lda, double *X, int incX);

void cblas_stpmv(int order, int Uplo, int TransA, int Diag, int N, const(float)* Ap, float *X, int incX);
void cblas_dtpmv(int order, int Uplo, int TransA, int Diag, int N, const(double)* Ap, double *X, int incX);
void cblas_ctpmv(int order, int Uplo, int TransA, int Diag, int N, const(float)* Ap, float *X, int incX);
void cblas_ztpmv(int order, int Uplo, int TransA, int Diag, int N, const(double)* Ap, double *X, int incX);

void cblas_stpsv(int order, int Uplo, int TransA, int Diag, int N, const(float)* Ap, float *X, int incX);
void cblas_dtpsv(int order, int Uplo, int TransA, int Diag, int N, const(double)* Ap, double *X, int incX);
void cblas_ctpsv(int order, int Uplo, int TransA, int Diag, int N, const(float)* Ap, float *X, int incX);
void cblas_ztpsv(int order, int Uplo, int TransA, int Diag, int N, const(double)* Ap, double *X, int incX);

void cblas_ssymv(int order, int Uplo, int N, float alpha, const(float)* A, int lda, const(float)* X, int incX, float beta, float *Y, int incY);
void cblas_dsymv(int order, int Uplo, int N, double alpha, const(double)* A, int lda, const(double)* X, int incX, double beta, double *Y, int incY);
void cblas_chemv(int order, int Uplo, int N, const(float)* alpha, const(float)* A, int lda, const(float)* X, int incX, const(float)* beta, float *Y, int incY);
void cblas_zhemv(int order, int Uplo, int N, const(double)* alpha, const(double)* A, int lda, const(double)* X, int incX, const(double)* beta, double *Y, int incY);


void cblas_sspmv(int order, int Uplo, int N, float alpha, const(float)* Ap, const(float)* X, int incX, float beta, float *Y, int incY);
void cblas_dspmv(int order, int Uplo, int N, double alpha, const(double)* Ap, const(double)* X, int incX, double beta, double *Y, int incY);

void cblas_sspr(int order, int Uplo, int N, float alpha, const(float)* X, int incX, float *Ap);
void cblas_dspr(int order, int Uplo, int N, double alpha, const(double)* X, int incX, double *Ap);

void cblas_chpr(int order, int Uplo, int N, float alpha, const(float)* X, int incX, float *A);
void cblas_zhpr(int order, int Uplo, int N, double alpha, const(double)* X,int incX, double *A);

void cblas_sspr2(int order, int Uplo, int N, float alpha, const(float)* X, int incX, const(float)* Y, int incY, float *A);
void cblas_dspr2(int order, int Uplo, int N, double alpha, const(double)* X, int incX, const(double)* Y, int incY, double *A);
void cblas_chpr2(int order, int Uplo, int N, const(float)* alpha, const(float)* X, int incX, const(float)* Y, int incY, float *Ap);
void cblas_zhpr2(int order, int Uplo, int N, const(double)* alpha, const(double)* X, int incX, const(double)* Y, int incY, double *Ap);

void cblas_chbmv(int order, int Uplo, int N, int K, const(float)* alpha, const(float)* A, int lda, const(float)* X, int incX, const(float)* beta, float *Y, int incY);
void cblas_zhbmv(int order, int Uplo, int N, int K, const(double)* alpha, const(double)* A, int lda, const(double)* X, int incX, const(double)* beta, double *Y, int incY);

void cblas_chpmv(int order, int Uplo, int N, const(float)* alpha, const(float)* Ap, const(float)* X, int incX, const(float)* beta, float *Y, int incY);
void cblas_zhpmv(int order, int Uplo, int N, const(double)* alpha, const(double)* Ap, const(double)* X, int incX, const(double)* beta, double *Y, int incY);

void cblas_sgemm(int Order, int TransA, int TransB, int M, int N, int K, float alpha, const(float)* A, int lda, const(float)* B, int ldb, float beta, float *C, int ldc);
void cblas_dgemm(int Order, int TransA, int TransB, int M, int N, int K, double alpha, const(double)* A, int lda, const(double)* B, int ldb, double beta, double *C, int ldc);
void cblas_cgemm(int Order, int TransA, int TransB, int M, int N, int K, const(float)* alpha, const(float)* A, int lda, const(float)* B, int ldb, const(float)* beta, float *C, int ldc);
void cblas_cgemm3m(int Order, int TransA, int TransB, int M, int N, int K, const(float)* alpha, const(float)* A, int lda, const(float)* B, int ldb, const(float)* beta, float *C, int ldc);
void cblas_zgemm(int Order, int TransA, int TransB, int M, int N, int K, const(double)* alpha, const(double)* A, int lda, const(double)* B, int ldb, const(double)* beta, double *C, int ldc);
void cblas_zgemm3m(int Order, int TransA, int TransB, int M, int N, int K, const(double)* alpha, const(double)* A, int lda, const(double)* B, int ldb, const(double)* beta, double *C, int ldc);


void cblas_ssymm(int Order, int Side, int Uplo, int M, int N, float alpha, const(float)* A, int lda, const(float)* B, int ldb, float beta, float *C, int ldc);
void cblas_dsymm(int Order, int Side, int Uplo, int M, int N, double alpha, const(double)* A, int lda, const(double)* B, int ldb, double beta, double *C, int ldc);
void cblas_csymm(int Order, int Side, int Uplo, int M, int N, const(float)* alpha, const(float)* A, int lda, const(float)* B, int ldb, const(float)* beta, float *C, int ldc);
void cblas_zsymm(int Order, int Side, int Uplo, int M, int N, const(double)* alpha, const(double)* A, int lda, const(double)* B, int ldb, const(double)* beta, double *C, int ldc);

void cblas_ssyrk(int Order, int Uplo, int Trans, int N, int K, float alpha, const(float)* A, int lda, float beta, float *C, int ldc);
void cblas_dsyrk(int Order, int Uplo, int Trans, int N, int K, double alpha, const(double)* A, int lda, double beta, double *C, int ldc);
void cblas_csyrk(int Order, int Uplo, int Trans, int N, int K, const(float)* alpha, const(float)* A, int lda, const(float)* beta, float *C, int ldc);
void cblas_zsyrk(int Order, int Uplo, int Trans, int N, int K, const(double)* alpha, const(double)* A, int lda, const(double)* beta, double *C, int ldc);

void cblas_ssyr2k(int Order, int Uplo, int Trans, int N, int K, float alpha, const(float)* A, int lda, const(float)* B, int ldb, float beta, float *C, int ldc);
void cblas_dsyr2k(int Order, int Uplo, int Trans, int N, int K, double alpha, const(double)* A, int lda, const(double)* B, int ldb, double beta, double *C, int ldc);
void cblas_csyr2k(int Order, int Uplo, int Trans, int N, int K, const(float)* alpha, const(float)* A, int lda, const(float)* B, int ldb, const(float)* beta, float *C, int ldc);
void cblas_zsyr2k(int Order, int Uplo, int Trans, int N, int K, const(double)* alpha, const(double)* A, int lda, const(double)* B, int ldb, const(double)* beta, double *C, int ldc);

void cblas_strmm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, float alpha, const(float)* A, int lda, float *B, int ldb);
void cblas_dtrmm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, double alpha, const(double)* A, int lda, double *B, int ldb);
void cblas_ctrmm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, const(float)* alpha, const(float)* A, int lda, float *B, int ldb);
void cblas_ztrmm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, const(double)* alpha, const(double)* A, int lda, double *B, int ldb);

void cblas_strsm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, float alpha, const(float)* A, int lda, float *B, int ldb);
void cblas_dtrsm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, double alpha, const(double)* A, int lda, double *B, int ldb);
void cblas_ctrsm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, const(float)* alpha, const(float)* A, int lda, float *B, int ldb);
void cblas_ztrsm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, const(double)* alpha, const(double)* A, int lda, double *B, int ldb);

void cblas_chemm(int Order, int Side, int Uplo, int M, int N, const(float)* alpha, const(float)* A, int lda, const(float)* B, int ldb, const(float)* beta, float *C, int ldc);
void cblas_zhemm(int Order, int Side, int Uplo, int M, int N, const(double)* alpha, const(double)* A, int lda, const(double)* B, int ldb, const(double)* beta, double *C, int ldc);

void cblas_cherk(int Order, int Uplo, int Trans, int N, int K, float alpha, const(float)* A, int lda, float beta, float *C, int ldc);
void cblas_zherk(int Order, int Uplo, int Trans, int N, int K, double alpha, const(double)* A, int lda, double beta, double *C, int ldc);

void cblas_cher2k(int Order, int Uplo, int Trans, int N, int K, const(float)* alpha, const(float)* A, int lda, const(float)* B, int ldb, float beta, float *C, int ldc);
void cblas_zher2k(int Order, int Uplo, int Trans, int N, int K, const(double)* alpha, const(double)* A, int lda, const(double)* B, int ldb, double beta, double *C, int ldc);

//void cblas_xerbla(int p, char *rout, char *form, ...);

/*** BLAS extensions ***/

void cblas_saxpby(int n, float alpha, const(float)* x, int incx,float beta, float *y, int incy);
void cblas_daxpby(int n, double alpha, const(double)* x, int incx,double beta, double *y, int incy);
void cblas_caxpby(int n, const(float)* alpha, const(float)* x, int incx,const(float)* beta, float *y, int incy);
void cblas_zaxpby(int n, const(double)* alpha, const(double)* x, int incx,const(double)* beta, double *y, int incy);

void cblas_somatcopy(int CORDER, int CTRANS, int crows, int ccols, float calpha, const(float)* a, int clda, float *b, int cldb);
void cblas_domatcopy(int CORDER, int CTRANS, int crows, int ccols, double calpha, const(double)* a, int clda, double *b, int cldb);
void cblas_comatcopy(int CORDER, int CTRANS, int crows, int ccols, float* calpha, float* a, int clda, float*b, int cldb);
void cblas_zomatcopy(int CORDER, int CTRANS, int crows, int ccols, double* calpha, double* a, int clda, double *b, int cldb);

void cblas_simatcopy(int CORDER, int CTRANS, int crows, int ccols, float calpha, float *a, int clda, int cldb);
void cblas_dimatcopy(int CORDER, int CTRANS, int crows, int ccols, double calpha, double *a, int clda, int cldb);
void cblas_cimatcopy(int CORDER, int CTRANS, int crows, int ccols, float* calpha, float* a, int clda, int cldb);
void cblas_zimatcopy(int CORDER, int CTRANS, int crows, int ccols, double* calpha, double* a, int clda, int cldb);

void cblas_sgeadd(int CORDER,int crows, int ccols, float calpha, float *a, int clda, float cbeta, float *c, int cldc);
void cblas_dgeadd(int CORDER,int crows, int ccols, double calpha, double *a, int clda, double cbeta, double *c, int cldc);
void cblas_cgeadd(int CORDER,int crows, int ccols, const(float)* calpha, float *a, int clda, const(float)* cbeta, float *c, int cldc);
void cblas_zgeadd(int CORDER,int crows, int ccols, const(double)* calpha, double *a, int clda, const(double)* cbeta, double *c, int cldc);
