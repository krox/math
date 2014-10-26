module math.lapacke;

import std.complex;

extern(C):

enum LAPACK_ROW_MAJOR = 101;
enum LAPACK_COL_MAJOR = 102;

/// solve A*x = b
int LAPACKE_sgesv( int matrix_order, int n, int nrhs, float* a, int lda, int* ipiv, float* b, int ldb );
int LAPACKE_dgesv( int matrix_order, int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb );
int LAPACKE_cgesv( int matrix_order, int n, int nrhs, Complex!float* a, int lda, int* ipiv, Complex!float* b, int ldb );
int LAPACKE_zgesv( int matrix_order, int n, int nrhs, Complex!double* a, int lda, int* ipiv, Complex!double* b, int ldb );

/// eigenvalues of a general matrix A
int LAPACKE_sgeev( int matrix_order, char jobvl, char jobvr, int n, float* a, int lda, float* wr, float* wi, float* vl, int ldvl, float* vr, int ldvr );
int LAPACKE_dgeev( int matrix_order, char jobvl, char jobvr, int n, double* a, int lda, double* wr, double* wi, double* vl, int ldvl, double* vr, int ldvr );
int LAPACKE_cgeev( int matrix_order, char jobvl, char jobvr, int n, Complex!float* a, int lda, Complex!float* w, Complex!float* vl, int ldvl, Complex!float* vr, int ldvr );
int LAPACKE_zgeev( int matrix_order, char jobvl, char jobvr, int n, Complex!double* a, int lda, Complex!double* w, Complex!double* vl, int ldvl, Complex!double* vr, int ldvr );
