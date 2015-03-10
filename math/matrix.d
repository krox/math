module math.matrix;

private import std.traits;
private import std.conv : to;
private import std.exception : assumeUnique;
private import std.typecons;
private import std.complex;
private import std.functional : binaryFun;
private import std.algorithm : min, max;

import jive.array;
private import math.lapacke;

class Matrix(T)
{
	abstract size_t height() const @property;

	abstract size_t width() const @property;

	abstract T opIndex(size_t i, size_t j) const;

	final T opIndex(size_t i) const
	{
		if(width == 1)
			return this[i,0];
		if(height == 1)
			return this[0,i];

		throw new Exception("not a vector");
	}

	final override string toString() const @property
	{
		string s;
		auto strings = Array2!string(height, width);
		auto pitch = Array!size_t(width, 0);

		for(size_t i = 0; i < height; ++i)
			for(size_t j = 0; j < width; ++j)
			{
				strings[i,j] = to!string(this[i,j]);
				pitch[j] = max(pitch[j], strings[i,j].length);
			}

		for(size_t i = 0; i < height; ++i)
		{
			if(i == 0)
				s ~= "⎛ ";
			else if(i == height-1)
				s ~= "⎝ ";
			else
				s ~= "⎜ ";


			for(size_t j = 0; j < width; ++j)
			{
				s ~= strings[i,j];
				for(int k = 0; k < pitch[j]+1-strings[i,j].length; ++k)
					s ~= " ";
			}

			if(i == 0)
				s ~= "⎞\n";
			else if(i == height-1)
				s ~= "⎠";
			else
				s ~= "⎟\n";
		}
		return s;
	}

	final Matrix opBinary(string op)(Matrix b) const
		if(op == "+" || op == "-")
	{
		if(width != b.width || height != b.height)
			throw new Exception("matrix dimension mismatch");

		auto a = Slice2!T(height, b.width);
		for(size_t i = 0; i < height; ++i)
			for(size_t j = 0; j < b.width; ++j)
			{
				a[i,j] = mixin("this[i,j] "~op~" b[i,j]");
			}
		return Matrix(a.assumeUnique);
	}

	final Matrix opBinary(string op)(Matrix b) const
		if(op == "*")
	{
		if(width != b.height)
			throw new Exception("matrix dimension mismatch");

		auto a = Slice2!T(height, b.width);
		for(size_t i = 0; i < height; ++i)
			for(size_t j = 0; j < b.width; ++j)
			{
				a[i,j] = this[i,0]*b[0,j];
				for(size_t k = 1; k < width; ++k)
					a[i,j] = a[i,j] + this[i,k]*b[k,j];
			}
		return Matrix(a.assumeUnique);
	}

	final Matrix pow(long exp)
	{
		assert(width > 0);
		if(width != height)
			throw new Exception("can not take a power of a non-square matrix");
		if(exp <= 0)
			throw new Exception("non-positive powers not implemented (yet?)");

		Matrix r;
		Matrix base = this;

		while(exp)
		{
			if(exp & 1)
			{
				if(r.width == 0)
					r = base;
				else
					r = r * base;
			}
			exp >>= 1;
			base = base * base;
		}
		return r;
	}

	/** solve the linear equations this * x = b */
	Matrix!T solve(Matrix!T b)
	{
		auto n = cast(int)this.width;
		if(this.height != n || b.height != n)
			throw new Exception("invalid matrix dimensions");

		auto rhs = b.dup;
		matSolve!T(this.dup, rhs);
		return Matrix!T(rhs.assumeUnique);
	}

	/** compute (a base of) the null-space of this */
	Matrix!T nullSpace()
	{
		return Matrix!T(matNullSpace!T(this.dup).assumeUnique);
	}

	/** compute (a base of) the eigen-space of an eigenvalue lambda */
	Matrix!T eigenSpace(T lambda)
	{
		return Matrix!T(matEigenSpace!T(this.dup, lambda).assumeUnique);
	}

	/** compute (complex) eigenvalues of this */
	Array!(ComplexType!T) eigenvalues()
	{
		auto n = cast(int)this.width;
		if(this.height != n)
			throw new Exception("invalid matrix dimensions");

		auto a = this.dup;
		auto w = Array!(ComplexType!T)(n);
		int r;

		     static if(is(T == float))       r = my_LAPACKE_sgeev(LAPACK_COL_MAJOR, 'N', 'N', n, a.ptr, n, w.ptr, null, n, null, n);
		else static if(is(T == double))	     r = my_LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'N', n, a.ptr, n, w.ptr, null, n, null, n);
		else static if(is(T == Complex!float))  r = LAPACKE_cgeev(LAPACK_COL_MAJOR, 'N', 'N', n, a.ptr, n, w.ptr, null, n, null, n);
		else static if(is(T == Complex!double)) r = LAPACKE_zgeev(LAPACK_COL_MAJOR, 'N', 'N', n, a.ptr, n, w.ptr, null, n, null, n);
		else throw new Exception("unsupported matrix type");

		if(r != 0)
			throw new Exception("lapack error");

		return w;
	}

	/** compute (real) eigenvalues of this assuming matrix is hermitian */
	Array!(RealType!T) hermitianEigenvalues()
	{
		auto n = cast(int)this.width;
		if(this.height != n)
			throw new Exception("invalid matrix dimensions");

		auto a = this.dup;
		auto w = Array!(RealType!T)(n);
		int r;

		     static if(is(T == float))          r = LAPACKE_ssyev(LAPACK_COL_MAJOR, 'N', 'U', n, a.ptr, n, w.ptr);
		else static if(is(T == double))         r = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'N', 'U', n, a.ptr, n, w.ptr);
		else static if(is(T == Complex!float))  r = LAPACKE_cheev(LAPACK_COL_MAJOR, 'N', 'U', n, a.ptr, n, w.ptr);
		else static if(is(T == Complex!double)) r = LAPACKE_zheev(LAPACK_COL_MAJOR, 'N', 'U', n, a.ptr, n, w.ptr);
		else throw new Exception("unsupported matrix type");

		if(r != 0)
			throw new Exception("lapack error");

		return w;
	}

	static auto opCall(Slice!(immutable(T), 2) data)
	{
		return new DenseMatrix!T(data);
	}

	static auto opCall(size_t height, size_t width, immutable(T)[] data)
	{
		return new DenseMatrix!T(Slice!(immutable(T), 2)(height, width, data));
	}

	static auto build(alias fun)(size_t h, size_t w)
	{
		auto data = Slice2!T(h, w);
		for(size_t j = 0; j < w; ++j)
			for(size_t i = 0; i < h; ++i)
				data[i,j] = binaryFun!(fun,"i","j")(i, j);
		return new DenseMatrix!T(data.assumeUnique);
	}

	/** only explicitly genrates upper/right half */
	static auto buildSymmetric(alias fun)(size_t n)
	{
		auto data = Slice2!T(n, n);
		for(size_t j = 0; j < n; ++j)
			for(size_t i = 0; i <= j; ++i)
			{
				data[i,j] = binaryFun!(fun,"i","j")(i, j);
				if(i != j)
					data[j,i] = data[i,j];
			}
		return new DenseMatrix!T(data.assumeUnique);
	}

	static auto buildBand(alias fun)(size_t n, int kl, int ku)
	{
		assert(kl >= 0 && ku >= 0);
		auto data = Slice2!T(kl+ku+1, n);
		for(size_t j = 0; j < n; ++j)
			for(size_t i = max(ku,j)-ku; i <= min(n-1,j+kl); ++i)
				data[i+ku-j,j] = binaryFun!(fun,"i","j")(i, j);
		return new BandMatrix!T(kl, ku, data.assumeUnique);
	}

	/** only explicitly genrates upper/right half */
	static auto buildSymmetricBand(alias fun)(size_t n, int kd)
	{
		assert(kd >= 0);
		auto data = Slice2!T(2*kd+1, n);
		for(size_t j = 0; j < n; ++j)
			for(size_t i = max(kd,j)-kd; i <= j; ++i)
			{
				data[i+kd-j,j] = binaryFun!(fun,"i","j")(i, j);
				if(i != j)
					data[j+kd-i, i] = data[i+kd-j, j];
			}
		return new BandMatrix!T(kd, kd, data.assumeUnique);
	}

	Slice2!T dup() const @property
	{
		auto r = Slice2!T(height, width);
		foreach(i,j, ref x; r)
			x = this[i,j];
		return r;
	}
}

/**
 * simple storage of a (possibly non-square) matrix.
 * column major order in a contigous array.
 */
final class DenseMatrix(T) : Matrix!T
{
	// TODO: decide if this should be Array2 instead of Slice2
	private Slice2!(immutable(T)) data;

	this(Slice2!(immutable(T)) data)
	{
		this.data = data;
	}

	override size_t height() const @property
	{
		return data.size[0];
	}

	override size_t width() const @property
	{
		return data.size[1];
	}

	override T opIndex(size_t i, size_t j) const
	{
		return data[i,j];
	}
}

/** efficient storage of a square band matrix */
final class BandMatrix(T) : Matrix!T
{
	const int kl, ku; // number of sub/super-diagonals

	private Slice2!(immutable(T)) data;

	this(int kl, int ku, Slice2!(immutable(T)) data)
	{
		this.data = data;
		this.kl = kl;
		this.ku = ku;

		if(kl+ku+1 != data.size[0] || kl >= height || ku >= height)
			throw new Exception("invalid dimensions of band matrix");
	}

	override size_t height() const @property
	{
		return data.size[1]; // width == height
	}

	override size_t width() const @property
	{
		return data.size[1];
	}

	override T opIndex(size_t i, size_t j) const
	{
		if(i+ku < j || i > j+kl)
			return T(0);

		return data[i-j+ku,j];
	}

	/** solve the linear equations this * x = b */
	override Matrix!T solve(Matrix!T b)
	{
		auto n = cast(int)this.width;
		auto nrhs = cast(int)b.width;
		if(this.height != n || b.height != n)
			throw new Exception("invalid matrix dimensions");

		auto a = this.data.data.dup;
		auto rhs = b.dup;
		auto p = new int[n];
		int r;

		     static if(is(T == float))          r = LAPACKE_sgbsv(LAPACK_COL_MAJOR, n, kl, ku, nrhs, a.ptr, 1+kl+ku, p.ptr, rhs.ptr, n);
		else static if(is(T == double))         r = LAPACKE_dgbsv(LAPACK_COL_MAJOR, n, kl, ku, nrhs, a.ptr, 1+kl+ku, p.ptr, rhs.ptr, n);
		else static if(is(T == Complex!float))  r = LAPACKE_cgbsv(LAPACK_COL_MAJOR, n, kl, ku, nrhs, a.ptr, 1+kl+ku, p.ptr, rhs.ptr, n);
		else static if(is(T == Complex!double)) r = LAPACKE_zgbsv(LAPACK_COL_MAJOR, n, kl, ku, nrhs, a.ptr, 1+kl+ku, p.ptr, rhs.ptr, n);
		else throw new Exception("unsupported matrix type");

		if(r != 0)
			throw new Exception("lapack error");

		return Matrix!T(rhs.assumeUnique);
	}

	/** compute (complex) eigenvalues of this */
	override Array!(ComplexType!T) eigenvalues()
	{
		throw new Exception("non-symmetric eigenvalues on band matrices not supported");
	}

	/** compute (real) eigenvalues of this, assuming matrix is hermitian */
	override Array!(RealType!T) hermitianEigenvalues()
	{
		auto n = cast(int)this.width;
		if(this.height != n)
			throw new Exception("invalid matrix dimensions");

		auto a = this.data.data.dup;
		auto w = Array!(RealType!T)(n);
		int r;

	         static if(is(T == float))          r = LAPACKE_ssbev(LAPACK_COL_MAJOR, 'N', 'U', n, ku, a.ptr, 1+kl+ku, w.ptr, null, n);
		else static if(is(T == double))         r = LAPACKE_dsbev(LAPACK_COL_MAJOR, 'N', 'U', n, ku, a.ptr, 1+kl+ku, w.ptr, null, n);
		else static if(is(T == Complex!float))  r = LAPACKE_chbev(LAPACK_COL_MAJOR, 'N', 'U', n, ku, a.ptr, 1+kl+ku, w.ptr, null, n);
		else static if(is(T == Complex!double)) r = LAPACKE_zhbev(LAPACK_COL_MAJOR, 'N', 'U', n, ku, a.ptr, 1+kl+ku, w.ptr, null, n);
		else throw new Exception("unsupported matrix type");

		if(r != 0)
			throw new Exception("lapack error");

		return w;
	}
}

/**
 * gauss a (possibly non-square) matrix.
 * afterwards column i has pivot in row r[i] (or r[i] = -1)
 */
Array!int matGauss(T)(Slice2!T m, Slice2!T b)
{
	assert(m.size[0] == b.size[0]);

	auto p = Array!int(m.size[0], -1); // "column i has pivot in row p[i]"

	for(int row = 0; row < m.size[0]; ++row)
	{
		// find a non-zero column to use as pivot
		int col = 0;
		while(col < m.size[1] && m[row,col] == 0)
			++col;
		if(col >= m.size[1]) // row is zero
			continue;
		p[col] = row;

		// normalize the row
		auto inv = 1/m[row,col];
		for(int i = 0; i < m.size[1]; ++i)
			m[row,i] = m[row,i] * inv;
		for(int i = 0; i < b.size[1]; ++i)
			b[row,i] = b[row,i] * inv;

		// substract it from every other row
		for(int i = 0; i < m.size[0]; ++i)
			if(i != row)
			{
				auto f = m[i,col];
				if(f == 0)
					continue;
				for(int j = 0; j < m.size[1]; ++j)
					m[i,j] = m[i,j] - f*m[row,j];
				for(int j = 0; j < b.size[1]; ++j)
					b[i,j] = b[i,j] - f*b[row,j];
			}
	}

	return p;
}

/**
 * backend for solving linear equations m*x = b.
 * m gets destroyed, result is stored into b.
 * currently only implemented for m square and invertible
 * uses lapack for float/double, and simple gauss for everything else
 * (which is only appropriate for exact types)
 */
void matSolve(T)(Slice2!T m, Slice2!T b)
	if(isFloatingPoint!(RealType!T))
{
	auto p = new int[n];
	int r;

	static if(is(T == float))
		r = LAPACKE_sgesv(LAPACK_COL_MAJOR, n, nrhs, a.ptr, n, p.ptr, b.ptr, n);
	else static if(is(T == double))
		r = LAPACKE_dgesv(LAPACK_COL_MAJOR, n, nrhs, a.ptr, n, p.ptr, b.ptr, n);
	else static if(is(T == Complex!float))
		r = LAPACKE_cgesv(LAPACK_COL_MAJOR, n, nrhs, a.ptr, n, p.ptr, b.ptr, n);
	else static if(is(T == Complex!double))
		r = LAPACKE_zgesv(LAPACK_COL_MAJOR, n, nrhs, a.ptr, n, p.ptr, b.ptr, n);
	else throw new Exception("unsupported matrix type");

	if(r != 0)
		throw new Exception("lapack error");
}

/** ditto */
void matSolve(T)(Slice2!T m, Slice2!T b)
	if(!isFloatingPoint!(RealType!T))
{
	int n = cast(int)m.size[0];
	int nrhs = cast(int)b.size[1];
	assert(m.size[1] == n && b.size[0] == n);

	// gauss the matrix and check it is not singular
	auto p = matGauss!T(m, b);
	foreach(x; p)
		if(x == -1)
			throw new Exception("matrix not invertible");

	// permute rhs according to b[i] = b[p2[i]]
	for(int start = 0; start < n; ++start)
	{
		// nothing to do or already done
		if(p[start] == start || p[start] == -1)
			continue;

		// permute the cycle
		for(int j = 0; j < nrhs; ++j)
		{
			T tmp = b[start, j];
			int i;
			for(i = start; p[i] != start; i = p[i])
				b[i,j] = b[p[i],j];
			b[i,j] = tmp;
		}

		// mark the cycle as done
		for(int i = start; i != -1; )
		{
			int t = p[i];
			p[i] = -1;
			i = t;
		}
	}
}

/**
 * backend for computing the nullspace of a matrix.
 * m gets destroyed, result is newly allocated.
 * uses simple gauss, so it is only suitable for exact types
 */
Slice2!T matNullSpace(T)(Slice2!T m)
{
	// gauss which does the actual work
	auto p = matGauss!T(m, Slice2!T(m.size[0],0,null));

	// count zero-rows, which is equal to the dimension of the nullspace
	int d = 0;
	for(int col = 0; col < m.size[1]; ++col)
		if(p[col] == -1)
			++d;

	// extract nullspace-base from columns that do not contain a pivot
	auto ns = Slice2!T(m.size[1], d, T(0));
	int k = 0;
	for(int col = 0; col < m.size[1]; ++col)
		if(p[col] == -1)
		{
			ns[col,k] = T(1);
			for(int c = 0; c < m.size[1]; ++c)
				if(p[c] != -1)
					ns[c,k] = -m[p[c],col];
			++k;
		}

	return ns;
}

/**
 * backend for computing an eigenspace of a matrix.
 * m gets destroyed, result is newly allocated.
 * uses simple gauss, so it is only suitable for exact types
 */
Slice2!T matEigenSpace(T)(Slice2!T m, T lambda)
{
	if(m.size[0] != m.size[1])
		throw new Exception("there is no eigenspace of a non-square matrix");

	for(int i = 0; i < m.size[0]; ++i)
		m[i,i] = m[i,i] - lambda;
	return matNullSpace!T(m);
}

struct Mat(T, size_t N, size_t M)
{
	T[N][M] c;

	enum Mat zero = constant(0);
	enum Mat identity = _identity();

	static Mat constant(T v) pure
	{
		Mat r;
		foreach(n; 0..N)
			foreach(m; 0..M)
				r[n,m] = v;
		return r;
	}

	static Mat _identity() pure
	{
		Mat r;
		foreach(n; 0..N)
			foreach(m; 0..M)
				r[n,m] = n==m ? 1 : 0;
		return r;
	}

	ref T opIndex(size_t n, size_t m)
	{
		return c[m][n];
	}

	ref const(T) opIndex(size_t n, size_t m) const
	{
		return c[m][n];
	}

	static if(M == 1)
	{
		ref inout(T) opIndex(size_t n) inout
		{
			return c[0][n];
		}
	}

	static if(N == 1)
	{
		ref inout(T) opIndex(size_t m) inout
		{
			return c[m][0];
		}
	}

	void opAddAssign(const Mat b)
	{
		foreach(n; 0..N)
			foreach(m; 0..M)
				 this[n,m] += b[n,m];
	}

	void opSubAssign(const Mat b)
	{
		foreach(n; 0..N)
			foreach(m; 0..M)
				 this[n,m] -= b[n,m];
	}

	void opMulAssign(T s)
	{
		foreach(n; 0..N)
			foreach(m; 0..M)
				 this[n,m] *= s;
	}

	Mat opMul(T s) const
	{
		Mat r = this;
		r *= s;
		return r;
	}

	static if(isFloatingPoint!T)
	{
		void opDivAssign(T s)
		{
			opMulAssign(1.0/s);
		}

		Mat opDiv(T s) const
		{
			return opMul(1.0/s);
		}
	}

	Mat opAdd(const Mat b) const
	{
		Mat r = this;
		r += b;
		return r;
	}

	Mat!(T,N,K) opMul(size_t K)(ref const Mat!(T,M,K) rhs) const
	{
		auto r = Mat!(T,N,K).zero;
		foreach(n; 0..N)
			foreach(k; 0..K)
				foreach(m; 0..M)
					r[n,k] += this[n,m]*rhs[m,k];
		return r;
	}

	static if(N==M) Mat pow(long exp) const
	{
		assert(exp >= 0);
		Mat r = identity;
		Mat base = this;

		while(exp)
		{
			if(exp & 1)
				r = r * base;
			exp >>= 1;
			base = base * base;
		}
		return r;
	}

	string toString() const @property
	{
		return to!string(c);	// TODO: transpose this
	}
}

alias Vec(T, size_t N)  = Mat!(T, N, 1);

T cross(T)(Vec!(T,2) a, Vec!(T,2) b)
{
    return a[0] * b[1] - a[1] * b[0];
}

bool ccw(T)(Vec!(T,2) a, Vec!(T,2) b)
{
	return cross!T(a, b) > 0;
}

bool ccw(T)(Vec!(T,2) a, Vec!(T,2) b, Vec!(T,2) c)
{
	return cross(b-a, c-a) > 0;
}
