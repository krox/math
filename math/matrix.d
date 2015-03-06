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
		for(size_t i = 0; i < height; ++i)
		{
			if(i == 0)
				s ~= "⎛ ";
			else if(i == height-1)
				s ~= "⎝ ";
			else
				s ~= "⎜ ";


			for(size_t j = 0; j < width; ++j)
				s ~= to!string(this[i,j]) ~ " ";

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
		return Matrix(height, width, a.assumeUnique);
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
		auto nrhs = cast(int)b.width;
		if(this.height != n || b.height != n)
			throw new Exception("invalid matrix dimensions");

		auto a = this.dup;
		auto rhs = b.dup;
		auto p = new int[n];

		int r;

		     static if(is(T == float))          r = LAPACKE_sgesv(LAPACK_COL_MAJOR, n, nrhs, a.ptr, n, p.ptr, rhs.ptr, n);
		else static if(is(T == double))         r = LAPACKE_dgesv(LAPACK_COL_MAJOR, n, nrhs, a.ptr, n, p.ptr, rhs.ptr, n);
		else static if(is(T == Complex!float))  r = LAPACKE_cgesv(LAPACK_COL_MAJOR, n, nrhs, a.ptr, n, p.ptr, rhs.ptr, n);
		else static if(is(T == Complex!double)) r = LAPACKE_zgesv(LAPACK_COL_MAJOR, n, nrhs, a.ptr, n, p.ptr, rhs.ptr, n);
		else throw new Exception("unsupported matrix type");

		if(r != 0)
			throw new Exception("lapack error");

		return Matrix!T(rhs.assumeUnique);
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

final class DenseMatrix(T) : Matrix!T
{
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
