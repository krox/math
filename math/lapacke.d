module math.lapacke;

import std.complex;
import jive.array;

extern(C):

enum LAPACK_ROW_MAJOR = 101;
enum LAPACK_COL_MAJOR = 102;

template ComplexType(T)
{
	static if(is(Complex!T))
		alias ComplexType = Complex!T;
	else
		alias ComplexType = T;
}

alias RealType(T:Complex!T) = T;
alias RealType(T) = T;

/// //////////////////////////////////////////////////////////////////////////
/// solve linear equations A*x = b
/// //////////////////////////////////////////////////////////////////////////

/// dense matrices
int LAPACKE_sgesv( int matrix_order, int n, int nrhs, float* a, int lda, int* ipiv, float* b, int ldb );
int LAPACKE_dgesv( int matrix_order, int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb );
int LAPACKE_cgesv( int matrix_order, int n, int nrhs, Complex!float* a, int lda, int* ipiv, Complex!float* b, int ldb );
int LAPACKE_zgesv( int matrix_order, int n, int nrhs, Complex!double* a, int lda, int* ipiv, Complex!double* b, int ldb );

// band matrix
int LAPACKE_sgbsv( int matrix_order, int n, int kl, int ku, int nrhs, float* ab, int ldab, int* ipiv, float* b, int ldb );
int LAPACKE_dgbsv( int matrix_order, int n, int kl, int ku, int nrhs, double* ab, int ldab, int* ipiv, double* b, int ldb );
int LAPACKE_cgbsv( int matrix_order, int n, int kl, int ku, int nrhs, Complex!float* ab, int ldab, int* ipiv, Complex!float* b, int ldb );
int LAPACKE_zgbsv( int matrix_order, int n, int kl, int ku, int nrhs, Complex!double* ab, int ldab, int* ipiv, Complex!double* b, int ldb );

/// //////////////////////////////////////////////////////////////////////////
/// compute (complex) eigenvalues of a general matrix A
/// //////////////////////////////////////////////////////////////////////////

// dense matrix
int LAPACKE_sgeev( int matrix_order, char jobvl, char jobvr, int n, float* a, int lda, float* wr, float* wi, float* vl, int ldvl, float* vr, int ldvr );
int LAPACKE_dgeev( int matrix_order, char jobvl, char jobvr, int n, double* a, int lda, double* wr, double* wi, double* vl, int ldvl, double* vr, int ldvr );
int LAPACKE_cgeev( int matrix_order, char jobvl, char jobvr, int n, Complex!float* a, int lda, Complex!float* w, Complex!float* vl, int ldvl, Complex!float* vr, int ldvr );
int LAPACKE_zgeev( int matrix_order, char jobvl, char jobvr, int n, Complex!double* a, int lda, Complex!double* w, Complex!double* vl, int ldvl, Complex!double* vr, int ldvr );

int my_LAPACKE_sgeev( int matrix_order, char jobvl, char jobvr, int n, float* a, int lda, Complex!float* w, float* vl, int ldvl, float* vr, int ldvr )
{
	auto wr = Array!float(n); // TODO: it is possible to not use these temporary arrays (O(n log n) time should be doable, O(n) time is really tough)
	auto wi = Array!float(n);
	int r = LAPACKE_sgeev(matrix_order, jobvl, jobvr, n, a, lda, wr.ptr, wi.ptr, vl, ldvl, vr, ldvr);
	for(int i = 0; i < n; ++i)
		w[i] = Complex!float(wr[i], wi[i]);
	return r;
}

int my_LAPACKE_dgeev( int matrix_order, char jobvl, char jobvr, int n, double* a, int lda, Complex!double* w, double* vl, int ldvl, double* vr, int ldvr )
{
	auto wr = Array!double(n); // TODO: it is possible to not use these temporary arrays (O(n log n) time should be doable, O(n) time is really tough)
	auto wi = Array!double(n);
	int r = LAPACKE_dgeev(matrix_order, jobvl, jobvr, n, a, lda, wr.ptr, wi.ptr, vl, ldvl, vr, ldvr);
	for(int i = 0; i < n; ++i)
		w[i] = Complex!double(wr[i], wi[i]);

	return r;
}

// band matrix
/* these routine dont seem to exist in lapack for some reason */

/// //////////////////////////////////////////////////////////////////////////
/// compute (real) eigenvalues of a hermitian matrix A
/// //////////////////////////////////////////////////////////////////////////

// dense matrix
int LAPACKE_ssyev( int matrix_order, char jobz, char uplo, int n, float* a, int lda, float* w );
int LAPACKE_dsyev( int matrix_order, char jobz, char uplo, int n, double* a, int lda, double* w );
int LAPACKE_cheev( int matrix_order, char jobz, char uplo, int n, Complex!float* a, int lda, float* w );
int LAPACKE_zheev( int matrix_order, char jobz, char uplo, int n, Complex!double* a, int lda, double* w );

// band matrix
int LAPACKE_ssbev( int matrix_order, char jobz, char uplo, int n, int kd, float* ab, int ldab, float* w, float* z, int ldz );
int LAPACKE_dsbev( int matrix_order, char jobz, char uplo, int n, int kd, double* ab, int ldab, double* w, double* z, int ldz );
int LAPACKE_chbev( int matrix_order, char jobz, char uplo, int n, int kd, Complex!float* ab, int ldab, float* w, Complex!float* z, int ldz );
int LAPACKE_zhbev( int matrix_order, char jobz, char uplo, int n, int kd, Complex!double* ab, int ldab, double* w, Complex!double* z, int ldz );
